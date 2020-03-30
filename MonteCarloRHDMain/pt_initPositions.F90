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
                             pt_maxPerProc, particles, pt_meshMe
  use Grid_interface, only : Grid_getBlkIndexLimits, Grid_getBlkBoundBox,&
                             Grid_getDeltas, Grid_getSingleCellVol,&
                             Grid_getBlkPtr, Grid_releaseBlkPtr
  use Simulation_data, only : a_rad
  use new_mcp, only : sample_cell_position, sample_iso_velocity,&
                      sample_time, sample_energy
  use relativity, only : transform_comoving_to_lab
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
  real :: Tgas, dV, totalE, weight, mcpweight
  integer :: ii, p, i, j, k
  real, dimension(MDIM) :: newxyz, newvel
  real :: newenergy, dshift

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

  call Grid_getBlkPtr(blockID, solnVec, CENTER)

  ! Do the loop over cells
  do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
    cellID(3) = k

    do j = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
      cellID(2) = j

      do i = blkLimits(LOW,IAXIS), blkLimits(HIGH, IAXIS)
        cellID(1) = i

        call Grid_getSingleCellVol(blockID, EXTERIOR, cellID, dV)

        Tgas = solnVec(TEMP_VAR, i, j, k)
        totalE = a_rad * (Tgas**4) * dV
        weight = totalE / pt_initradfield_num

        ! Start loop to sample MCPs
        do ii = 1, pt_initradfield_num
          p = p + 1

          call sample_cell_position(bndBox, deltaCell, cellID, newxyz)

          call sample_iso_velocity(newvel)

          call sample_energy(solnVec, cellID, newenergy)

          ! initialize the particle
          ! Directly edit particles array
          particles(BLK_PART_PROP,p) = real(blockID)
          particles(PROC_PART_PROP,p) = real(pt_meshMe)

          ! TREM should be -1.0 because it's an initial field
          particles(TREM_PART_PROP,p) = -1.0

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
            ! call some conversion function to convert
            call transform_comoving_to_lab(cellID, solnVec,&
                          particles(:,p), dshift)
          end if

        end do ! done MCP sampling

      end do
    end do
  end do

  call Grid_releaseBlkPtr(blockID, solnVec)

  ! Update final particle count
  pt_numLocal = p

  ! success check
  if (pt_numLocal > pt_maxPerProc) then
    call Driver_abortFlash("pt_initPositions: initial radiation field exceeds&
                            maximum allowed number pt_maxPerProc.")
  else
    success = .true.
  end if

!----------------------------------------------------------------------
  
end subroutine pt_initPositions


