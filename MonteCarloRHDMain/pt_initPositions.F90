!!****if* source/Particles/localAPI/pt_initPositions
!!
!! NAME
!!    pt_initPositions
!!
!! SYNOPSIS
!!
!!    call pt_initPositions(integer(in)  :: blockID,
!!                          logical(out) :: success)
!!
!! DESCRIPTION
!!
!!  Initializes particle locations for one block in the grid.
!!
!!  The initialization of locations may fail, if adding the
!!  number of particles requested (by the choice of initialization
!!  method and runtime parameters) would exceed the maximum allowed
!!  number of particles allowed in any MPI task given by runtime
!!  parameter pr_maxPerProc. In that case, an impementation of this
!!  interface may:
!!    o abort the simulation (by calling Driver_abortFlash); or
!!    o return with success set to FALSE.
!!  If returning with success=FALSE, pt_numLocal may have been
!!  reset to 0.
!!
!! ARGUMENTS
!!
!!  blockID:        local block ID containing particles to create
!!
!!  success:        returns .TRUE. if positions for all particles
!!                  that should be assigned to this block have been
!!                  successfully initialized.
!!
!! SIDE EFFECTS
!!
!!  Updates the variable pt_numLocal (in Particles_data), but adding
!!  to the previous value (if successful) the number of particles that
!!  were placed in the block given by blockID. If unsuccessful, pt_numLocal
!!  may be reset to 0.
!!
!!  The values of position attributes of the newly initialized particles
!!  will be updated.
!!
!!  The values of blk and proc attributes of newly initialized particles
!!  will also be updated appropriately.
!!
!!  The values of velocity attributes of particles MAY also be updated
!!  by some implementations of this interface.
!!
!! NOTES
!!
!!  The values of the tag attribute of particles will not be initialize
!!  here, that is left for later and typically done later in the same
!!  invocation of particles_initpositions that calls this interface.
!!
!! SEE ALSO
!!
!!  Particles_data
!!  Particles_initPositions
!!***


subroutine pt_initPositions (blockID,success)
  use Particles_data, only : pt_numLocal, pt_initradfield_num, pt_is_veldp,&
                             pt_maxPerProc, particles, pt_meshMe,&
                             pt_indexList, pt_indexCount, pt_samp_Tgas,&
                             pt_interp_vel_LT
  use pt_interface, only :  pt_updateTypeDS
  use Grid_data, only : gr_geometry
  use Grid_interface, only : Grid_getBlkIndexLimits, Grid_getBlkBoundBox,&
                             Grid_getDeltas, Grid_getSingleCellVol,&
                             Grid_getBlkPtr, Grid_releaseBlkPtr,&
                             Grid_outsideBoundBox,&
                             Grid_sortParticles, Grid_moveParticles
  use Simulation_data, only : a_rad
  use new_mcp, only : sample_cell_position, sample_iso_velocity,&
                      sample_time, sample_energy
  use relativity, only : compute_gamma, transform_comoving_to_lab
  use transport, only : interp_gas_vel_at_mcp
  use spherical, only : get_cartesian_velocity
  use Driver_interface, only : Driver_abortFlash

  implicit none
#include "constants.h"
#include "Flash.h"

  integer, INTENT(in) :: blockID
  logical,intent(OUT) :: success

  ! aux variables
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  real, dimension(2,MDIM) :: bndBox
  real, dimension(MDIM) :: deltaCell
  integer, dimension(MDIM) :: cellID
  real, pointer :: solnVec(:,:,:,:)
  integer, dimension(MAXBLOCKS,NPART_TYPES) :: particlesPerBlk
  real :: Tgas, dV, totalE, weight, mcpweight
  integer :: ii, p, i, j, k
  real, dimension(MDIM) :: newxyz, newvel
  real, dimension(MDIM) :: newpos
  logical :: isoutside
  integer, dimension(MDIM) :: neghdir
  integer :: num_LT_off
  real :: newenergy, dshift
  real :: gamm_fac, dV_0

  ! cell coordinates for vel interpolation
  integer, parameter :: ycoordsize = NYB + 2*NGUARD
  real, dimension(ycoordsize) :: ycoords
  integer, parameter :: xcoordsize = NXB + 2*NGUARD
  real, dimension(xcoordsize) :: xcoords
  integer, parameter :: zcoordsize = NZB + 2*NGUARD
  real, dimension(zcoordsize) :: zcoords
  real, dimension(MDIM) :: v_gas, v_gas_cart
  real, dimension(MDIM) :: cell_coords

  ! Commenting out stub code
  !success = .true. ! DEV: returns true because this stub creates no particles,
                   ! therefore all of those zero particles were created successfully
  !return

  ! Implementation for MCRHD's initial thermal radiation field
  success = .false. ! default
  p = pt_numLocal

  if (pt_initradfield_num <= 0) then
    success = .true.    ! successfully create no particles
    return
  end if

  ! Non-zero initial radiation field MCP number starts

  ! Gather cell indicies in current block
  call Grid_getBlkIndexLimits(blockID, blkLimits, blkLimitsGC)

  ! Get the grid geometry of this block
  call Grid_getBlkBoundBox(blockID,bndBox)
  call Grid_getDeltas(blockID, deltaCell)
  call Grid_getCellCoords(IAXIS, blockID, CENTER, .TRUE.,&
                          xcoords, xcoordsize)
  ycoords = 1.0
  zcoords = 1.0
  if (NDIM > 1) then
    call Grid_getCellCoords(JAXIS, blockID, CENTER, .TRUE.,&
                            ycoords, ycoordsize)

    if (NDIM > 2) then
      call Grid_getCellCoords(KAXIS, blockID, CENTER, .TRUE.,&
                              zcoords, zcoordsize)
    end if
  end if

  call Grid_getBlkPtr(blockID, solnVec, CENTER)

  num_LT_off = 0
  ! Do the loop over cells
  do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
    cellID(3) = k

    do j = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
      cellID(2) = j

      do i = blkLimits(LOW,IAXIS), blkLimits(HIGH, IAXIS)
        cellID(1) = i

        ! get gas velocity in current cell for gamma calc.
        v_gas = solnVec(VELX_VAR:VELZ_VAR, cellID(IAXIS),&
                         cellID(JAXIS), cellID(KAXIS))

        ! convert spherical velocity to Cartesian velocity
        if (gr_geometry == SPHERICAL) then
          cell_coords = (/ xcoords(cellID(IAXIS)),&
                           ycoords(cellID(JAXIS)),&
                           zcoords(cellID(KAXIS)) /)
          call get_cartesian_velocity(cell_coords, v_gas,&
                                      v_gas_cart)
          v_gas = v_gas_cart
        end if

        call Grid_getSingleCellVol(blockID, EXTERIOR, cellID, dV)
        dV_0 = dV  ! default the same
        ! convert lab-frame volume to CMF volume
        if (pt_is_veldp) then
          call compute_gamma(v_gas, gamm_fac)
          dV_0 = gamm_fac * dV
        end if

        Tgas = solnVec(TEMP_VAR, i, j, k)

        ! radiation energy in CMF
        totalE = a_rad * (Tgas**4) * dV_0
        weight = totalE / pt_initradfield_num

        ! Start loop to sample MCPs
        do ii = 1, pt_initradfield_num
          p = p + 1

          call sample_cell_position(bndBox, deltaCell, cellID, newxyz)

          call sample_iso_velocity(newvel)

          call sample_energy(solnVec, cellID, pt_samp_Tgas, -1.0, newenergy)

          ! initialize the particle
          ! Directly edit particles array
          particles(BLK_PART_PROP,p) = real(blockID)
          particles(PROC_PART_PROP,p) = real(pt_meshMe)

          ! TREM should be -1.0 because it's an initial field
          particles(TREM_PART_PROP,p) = -1.0
          particles(ISNW_PART_PROP,p) = 1.0

          particles(POSX_PART_PROP:POSZ_PART_PROP,p) = newxyz
          particles(VELX_PART_PROP:VELZ_PART_PROP,p) = newvel

          particles(ENER_PART_PROP,p) = newenergy

          mcpweight = weight / newenergy

          particles(NPIN_PART_PROP,p) = mcpweight
          particles(NUMP_PART_PROP,p) = mcpweight
          particles(NUM0_PART_PROP,p) = mcpweight

          ! Both TYPE and TAG are taken care of in Particles_initPositions

          ! transform back to lab frame if needed
          if (pt_is_veldp) then
            if (pt_interp_vel_LT) then
              call interp_gas_vel_at_mcp(bndBox, deltaCell, newxyz,&
                                         solnVec, v_gas)
              cell_coords = newxyz  ! put cell location at MCP
            end if

            ! convert spherical velocity to Cartesian velocity
            if (gr_geometry == SPHERICAL) then
              call get_cartesian_velocity(cell_coords, v_gas,&
                                          v_gas_cart)
              v_gas = v_gas_cart
            end if

            ! call some conversion function to convert
            call transform_comoving_to_lab(v_gas,&
                          -1.0, particles(:,p), dshift)  ! dtNew = -1 as dummy

            ! check whether the mcp is still contained in current block
            newpos = particles(POSX_PART_PROP:POSZ_PART_PROP, p)
            call Grid_outsideBoundBox(newpos, bndBox, isoutside, neghdir)
            if (isoutside) num_LT_off = num_LT_off + 1
          end if

        end do ! done MCP sampling

      end do
    end do
  end do

  call Grid_releaseBlkPtr(blockID, solnVec)

  ! Update final particle count
  pt_numLocal = p

  ! Do a moveParticles if isveldp is active.
  ! MCPs could be Lorentz transformed off the parent block
  if (pt_is_veldp) then
    print *, "Rank", pt_meshMe, "has", num_LT_off, "Lorentz-boosted off block."
    call Grid_moveParticles(particles, NPART_PROPS, pt_maxPerProc,&
                            pt_numLocal, pt_indexList, pt_indexCount, .false.)
#ifdef TYPE_PART_PROP
    call Grid_sortParticles(particles, NPART_PROPS, pt_numLocal, NPART_TYPES,&
                            pt_maxPerProc, particlesPerBlk, BLK_PART_PROP,&
                            TYPE_PART_PROP)
#else
    call Grid_sortParticles(particles, NPART_PROPS, pt_numLocal, NPART_TYPES,&
                            pt_maxPerProc, particlesPerBlk, BLK_PART_PROP)
#endif
    call pt_updateTypeDS(particlesPerBlk)
  end if ! only do this when LT is on

  ! success check
  if (pt_numLocal > pt_maxPerProc) then
    call Driver_abortFlash("pt_initPositions: initial radiation field exceeds&
                            maximum allowed number pt_maxPerProc.")
  else
    success = .true.
  end if

!----------------------------------------------------------------------
  
end subroutine pt_initPositions


