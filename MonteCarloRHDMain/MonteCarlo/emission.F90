module emission

  implicit none
  contains


! Thermal and face emission
subroutine emit_mcps(particles, p_count, dtNew, ind)
  use Particles_data, only : pt_maxPerProc, pt_maxnewnum, pt_numLocal,&
                             pt_meshMe, pt_typeInfo
  use Particles_interface, only : Particles_addNew 
  use Timers_interface, only : Timers_start, Timers_stop
  use Grid_interface, only : Grid_getBlkPtr, Grid_getListOfBlocks,&
                             Grid_releaseBlkPtr, Grid_sortParticles
  use Driver_interface, only : Driver_abortFlash
  use ionization, only : calc_recomb_emissivity
  implicit none
#include "Flash.h"
#include "constants.h"
#include "Particles.h"

  ! Input/Output 
  real,dimension(NPART_PROPS,pt_maxPerProc),intent(inout) :: particles
  integer, intent(in) :: p_count
  real, intent(in)  :: dtNew
  integer, intent(in) :: ind

  ! Other parameters
  integer, dimension(MAXBLOCKS) :: blkList
  integer :: numofblks
  real, pointer :: solnVec(:,:,:,:)

  ! aux parameters
  integer :: b
  logical :: success 
  integer :: num_thermal, num_face, old_pt_numLocal

  integer :: now_num
  real, dimension(MDIM,pt_maxnewnum) :: now_pos, now_vel
  real, dimension(pt_maxnewnum) :: now_time, now_energy, now_weight

  integer :: new_num
  real, dimension(MDIM,pt_maxnewnum) :: new_pos, new_vel
  real, dimension(pt_maxnewnum) :: new_time, new_energy, new_weight

  integer :: i_old_begin, i_old_end
  integer, dimension(MAXBLOCKS,NPART_TYPES) :: particlesPerBlk

  call Timers_start("MCP Emission")

  call Grid_getListOfBlocks(LEAF, blkList, numofblks)

  ! Initializing new mcp arrays
  new_pos = 0.0d0
  new_vel = 0.0d0
  new_time = 0.0d0
  new_energy = 0.0d0
  new_weight = 0.0d0
  new_num = 0

  ! Loop through leaf blocks
  do b = 1, numofblks
    call Grid_getBlkPtr(blkList(b), solnVec, CENTER)

    ! Reset Monte Carlo tally variables
    call reset_deposition_var(solnVec)

    ! Thermal radiation
    call thermal_emission(blkList(b), solnVec, dtNew,&
           now_num, now_pos, now_time, &
           now_energy, now_vel, now_weight)
    if (now_num .GT. 0) then
      new_pos(1:MDIM,new_num+1:new_num+now_num) = now_pos(1:MDIM,1:now_num)
      new_vel(1:MDIM,new_num+1:new_num+now_num) = now_vel(1:MDIM,1:now_num)
      new_time(new_num+1:new_num+now_num) = now_time(1:now_num)
      new_energy(new_num+1:new_num+now_num) = now_energy(1:now_num)
      new_weight(new_num+1:new_num+now_num) = now_weight(1:now_num)

      new_num = new_num + now_num
      num_thermal = new_num
    end if

    ! Face emission
    call face_emission(blkList(b), solnVec, dtNew,&
           now_num, now_pos, now_time, &
           now_energy, now_vel, now_weight)
    if (now_num .GT. 0) then
      new_pos(1:MDIM,new_num+1:new_num+now_num) = now_pos(1:MDIM,1:now_num)
      new_vel(1:MDIM,new_num+1:new_num+now_num) = now_vel(1:MDIM,1:now_num)
      new_time(new_num+1:new_num+now_num) = now_time(1:now_num)
      new_energy(new_num+1:new_num+now_num) = now_energy(1:now_num)
      new_weight(new_num+1:new_num+now_num) = now_weight(1:now_num)

      new_num = new_num + now_num
      num_face = new_num
    end if

    ! Recombination emission
    call calc_recomb_emissivity(blkList(b), solnVec, dtNew)
    ! Case B approximation, no MCP generation required

    call Grid_releaseBlkPtr(blkList(b), solnVec)
  end do

  ! Appending new MCPs to particles array, pt_numLocal updated
  old_pt_numLocal = pt_numLocal
  print *, "Rank", pt_meshMe, "receives", new_num, "new MCPs."

  call Particles_addNew(new_num, new_pos(:,1:new_num), success)

  ! Putting the energy into the particles DS
  if (new_num .GT. 0) then
    particles(TREM_PART_PROP,old_pt_numLocal+1:pt_numLocal) = new_time(1:new_num)
    particles(ENER_PART_PROP,old_pt_numLocal+1:pt_numLocal) = new_energy(1:new_num)
    particles(VELX_PART_PROP:VELZ_PART_PROP,old_pt_numLocal+1:pt_numLocal) = &
                                                              new_vel(:,1:new_num)
    particles(NPIN_PART_PROP, old_pt_numLocal+1:pt_numLocal) = new_weight(1:new_num)
    particles(NUMP_PART_PROP, old_pt_numLocal+1:pt_numLocal) = new_weight(1:new_num)
    particles(NUM0_PART_PROP, old_pt_numLocal+1:pt_numLocal) = new_weight(1:new_num)

    ! Other bookkeeping attributes
    particles(ISNW_PART_PROP, old_pt_numLocal+1:pt_numLocal) = 1.0d0
    particles(ISAB_PART_PROP, old_pt_numLocal+1:pt_numLocal) = 0.0d0
    particles(ISCP_PART_PROP, old_pt_numLocal+1:pt_numLocal) = 0.0d0

#ifdef TYPE_PART_PROP
    ! Photons come first if sink particles are present
    particles(TYPE_PART_PROP, old_pt_numLocal+1:pt_numLocal) = pt_typeInfo(PART_TYPE, ind)
#endif

  end if

  ! Filling in key attributes for census MCPs.
  ! Here PART_LOCAL is the old number count of particles
  i_old_begin = pt_typeInfo(PART_TYPE_BEGIN,ind)
  i_old_end   = i_old_begin + pt_typeInfo(PART_LOCAL,ind) - 1
  particles(TREM_PART_PROP, i_old_begin:i_old_end) = dtNew
  particles(ISNW_PART_PROP, i_old_begin:i_old_end) = 0.0d0

  ! Sorting the old and new MCPs
  call Timers_start("Sort before RT")
#ifdef TYPE_PART_PROP
  call Grid_sortParticles(particles,NPART_PROPS,&
                          pt_numLocal,NPART_TYPES, pt_maxPerProc,&
                          particlesPerBlk,BLK_PART_PROP, TYPE_PART_PROP)
#else
  call Grid_sortParticles(particles,NPART_PROPS,&
                          pt_numLocal,NPART_TYPES, pt_maxPerProc,&
                          particlesPerBlk,BLK_PART_PROP)
#endif

  !! Update the pt_typeInfo
  call pt_updateTypeDS(particlesPerBlk)
  call Timers_stop("Sort before RT")

  call Timers_stop("MCP Emission")

end subroutine emit_mcps


subroutine thermal_emission(blockID, solnVec, dtNew,&
           now_num, now_pos, now_time, &
           now_energy, now_vel, now_weight)
  use Particles_data, only : pt_maxnewnum, pt_ThermalEmission,&
              pt_is_grey, pt_is_eff_scattering, pt_num_tmcps_tstep,&
              pt_is_veldp
  use Grid_interface, only : Grid_getBlkIndexLimits,&
              Grid_getBlkBoundBox, Grid_getDeltas
  use Simulation_data, only : R, sigma, a_rad
  use Eos_data, only : eos_singleSpeciesA
  use opacity, only : calc_abs_opac
  use Driver_interface, only : Driver_abortFlash
  use rhd, only : cellAddVar
  use new_mcp, only : sample_cell_position, sample_iso_velocity,&
                      sample_time, sample_energy
  use relativity, only : transform_comoving_to_lab

  implicit none
#include "Flash.h"
#include "constants.h"

  ! Input/Output
  integer, intent(in) :: blockID
  real, pointer :: solnVec(:,:,:,:)
  real, intent(in) :: dtNew
  integer, intent(out) :: now_num
  real, dimension(MDIM,pt_maxnewnum), intent(out) :: now_pos, now_vel
  real, dimension(pt_maxnewnum), intent(out) :: now_time, now_energy, now_weight

  ! aux parameters
  integer :: i, j, k, ii
  integer, dimension(MDIM) :: cellID
  integer, dimension(2, MDIM) :: blkLimits, blkLimitsGC
  real, dimension(LOW:HIGH, MDIM) :: bndBox
  real, dimension(MDIM) :: deltaCell
  real :: dV, temp, rho, gamc, ka, kp, dE, dEMIE, c_V, beta, fn
  real, parameter :: eps_dummy = 1.0d-10

  real, dimension(MDIM) :: newxyz, newvel
  real :: weight_per_mcp, dE_per_mcp, dtNow, newenergy
  real, dimension(NPART_PROPS) :: newparticle
  real :: dshift


  ! Initialization
  now_num = 0
  now_pos = 0.0d0
  now_vel = 0.0d0
  now_time = 0.0d0
  now_energy = 0.0d0
  now_weight = 0.0d0

  if (.not. pt_ThermalEmission) return

  ! Obtain block info
  call Grid_getBlkIndexLimits(blockID, blkLimits, blkLimitsGC)
  call Grid_getBlkBoundBox(blockID, bndBox)
  call Grid_getDeltas(blockID, deltaCell) 

  do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
    do j = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
      do i = blkLimits(LOW,IAXIS), blkLimits(HIGH, IAXIS)

        cellID = (/ i, j, k /)

        CALL Grid_getSingleCellVol(blockID, EXTERIOR, cellID, dV)

        ! Obtaining grid info of current cell
        temp = solnVec(TEMP_VAR, i, j, k)
        rho  = solnVec(DENS_VAR, i, j, k)
        gamc = solnVec(GAMC_VAR, i, j, k)

        dE = 0.0d0
        if (pt_is_grey) then
          call calc_abs_opac(cellID, solnVec, eps_dummy, ka)
          ! ka here is already in co-moving frame
          dE = 4.0 * ka * sigma * (temp**4.0) * dV * dtNew
          kp = ka
        else
          call Driver_abortFlash("thermal_emission: non-grey thermal emission&
                                  not yet implemented.")
        end if

        fn = 1.0d0
        if ((pt_is_eff_scattering) .and. (dE .ne. 0.0d0)) then
          c_V = R / ((gamc - 1.0d0)*eos_singleSpeciesA)

          ! IMC beta
          beta = (4.0d0 * a_rad * (temp**3)) / (rho * c_V)

          call comp_fleck_factor(kp, beta, dtNew, fn)
        end if
        solnVec(FLEC_VAR,i,j,k) = fn

        dE = fn * dE

        dEMIE = -dE / (rho * dV)

        ! Actual emission
        if (dE .gt. 0.0d0) then
          dE_per_mcp = dE / pt_num_tmcps_tstep

          call cellAddVar(solnVec, cellID, EMIE_VAR, dEMIE)

          ! Loop to create new MCPs
          do ii = 1, pt_num_tmcps_tstep

            call sample_cell_position(bndBox, deltaCell, cellID, newxyz)

            call sample_time(dtNew, dtNow)

            call sample_iso_velocity(newvel)

            call sample_energy(solnVec, cellID, newenergy)
            weight_per_mcp = dE_per_mcp / newenergy

            ! Initialization
            newparticle = 0.0
            ! Put attributes into a new particle array for
            ! easier Lorentz transformation
            newparticle(POSX_PART_PROP:POSZ_PART_PROP) = newxyz
            newparticle(VELX_PART_PROP:VELZ_PART_PROP) = newvel
            newparticle(ENER_PART_PROP) = newenergy
            newparticle(TREM_PART_PROP) = dtNow
            newparticle(NUMP_PART_PROP) = weight_per_mcp

            ! Convert to lab frame if velocity dependent is on
            if (pt_is_veldp) then
              ! call some conversion function to convert
              call transform_comoving_to_lab(cellID, solnVec,&
                            newparticle, dshift) 
            end if


            ! Record the MCP attributes
            now_time(now_num + ii)   = newparticle(TREM_PART_PROP)
            now_pos(:, now_num + ii) = newparticle(POSX_PART_PROP:&
                                                   POSZ_PART_PROP)
            now_energy(now_num + ii) = newparticle(ENER_PART_PROP)
            now_vel(:, now_num + ii) = newparticle(VELX_PART_PROP:&
                                                   VELZ_PART_PROP)
            now_weight(now_num + ii) = newparticle(NUMP_PART_PROP)
          end do

          now_num = now_num + pt_num_tmcps_tstep
          
        end if

      end do
    end do
  end do

end subroutine thermal_emission

subroutine point_emission(blockID, solnVec, dtNew,&
                          now_num, now_pos, now_time,&
                          now_energy, now_vel, now_weight)
  use Particles_data, only : pt_maxnewnum, pt_PointEmission, pt_PointPulse,&
                             pt_PointPulseErad,&
                             pt_PointLuminosity, pt_PointSrcPosOffset,&
                             pt_num_pmcps_tstep, pt_meshMe
  use Grid_interface, only : Grid_getDeltas, Grid_getBlkIDFromPos
  use Paramesh_comm_data, only : amr_mpi_meshComm
  use new_mcp, only : sample_iso_velocity,&
                      sample_time, sample_energy
  use transport, only : get_cellID
  implicit none
#include "constants.h"

  ! Input/Output
  integer, intent(in) :: blockID
  real, pointer :: solnVec(:,:,:,:)
  real, intent(in) :: dtNew
  integer, intent(out) :: now_num
  real, dimension(MDIM,pt_maxnewnum), intent(out) :: now_pos, now_vel
  real, dimension(pt_maxnewnum), intent(out) :: now_time, now_energy, now_weight

  ! saved attributes
  logical, save :: first_call = .true.

  ! Aux variables
  integer :: ii
  real, dimension(LOW:HIGH, MDIM) :: bndBox
  real, dimension(MDIM) :: deltaCell, newxyz, newvel
  integer :: originblkID, originprocID
  integer, dimension(MDIM) :: cellID
  real, dimension(MDIM) :: origin = (/ 0.0d0, 0.0d0, 0.0d0 /)
  real :: dE, weight_per_mcp, dE_per_mcp, dtNow, newenergy

  ! Exit when not using point emission
  if (.not. pt_PointEmission) then
    return
  end if 

  ! Exit after first call if in pulse mode
  if ((.not. first_call) .and. (pt_PointPulse)) then
    return
  end if

  ! Get local cell size
  call Grid_getBlkBoundBox(blockID, bndBox)
  call Grid_getDeltas(blockID, deltaCell)

  ! Offset the source location to avoid on-grid complications,
  ! default is 1.0e-3*dx
  do ii = 1, MDIM
    origin(ii) = origin(ii) + pt_PointSrcPosOffset*deltaCell(ii)
  end do

  ! Find out which block the origin resides on
  call Grid_getBlkIDFromPos(origin, originblkID, originprocID, amr_mpi_meshComm)

  ! Output initialization
  now_num = 0
  now_pos = 0.0d0
  now_vel = 0.0d0
  now_time = 0.0d0
  now_energy = 0.0d0
  now_weight = 0.0d0

  ! Only emit if (blockID == block where the origin is)
  if ((pt_meshMe == originprocID) .and. (blockID == originblkID)) then
    if (pt_PointPulse) then
      dE = pt_PointPulseErad
      first_call = .false.
    else ! continuous source
      dE = pt_PointLuminosity * dtNew
    end if

    dE_per_mcp = dE / pt_num_pmcps_tstep

    do ii = 1, pt_num_pmcps_tstep
      ! The full time step
      dtNow = dtNew 

      ! newxyz is just the origin
      newxyz = origin

      ! Isotropic radial velocity for MCP
      call sample_iso_velocity(newvel)

      ! Use cell information for sampling MCP frequency,
      ! constant eps in the gray case
      call get_cellID(bndBox, deltaCell, newxyz, cellID)
      call sample_energy(solnVec, cellID, newenergy)
      weight_per_mcp = dE_per_mcp / newenergy

      ! Record the MCP attributes
      now_time(now_num + ii) = dtNow
      now_pos(:, now_num + ii) = newxyz
      now_energy(now_num + ii) = newenergy
      now_vel(:, now_num + ii) = newvel
      now_weight(now_num + ii) = weight_per_mcp
    end do

    ! Done sampling the point source MCPs
    now_num = now_num + pt_num_pmcps_tstep
  end if

end subroutine point_emission

subroutine face_emission(blockID, solnVec, dtNew,&
           now_num, now_pos, now_time,&
           now_energy, now_vel, now_weight)
  use Particles_data, only : pt_maxnewnum, pt_FaceEmission, pt_is_FacePlanck,&
                             pt_is_grey, pt_num_fmcps_tstep,&
                             pt_FaceEmissionSide, pt_FaceEmissionAxis,&
                             pt_FacePlanckTemp, pt_constFaceFlux,&
                             pt_is_radial_face_vel, pt_is_iso_face_vel,&
                             pt_is_therm_face_vel, pt_smlpush
  use Grid_interface, only : Grid_getBlkBoundBox, Grid_getDeltas,&
                             Grid_getBlkBC, Grid_getBlkPhysicalSize
  use Simulation_data, only : sigma, clight
  use Driver_interface, only : Driver_abortFlash
  use new_mcp, only : sample_blk_position, sample_iso_velocity,&
                      sample_therm_face_velocity, &
                      sample_time, sample_energy
  use transport, only : get_cellID
  use spherical, only : get_cartesian_position
  ! debug
  use Driver_data, only : dr_meshme

  implicit none
#include "Flash.h"
#include "constants.h"

  ! Input/Output
  integer, intent(in) :: blockID
  real, pointer :: solnVec(:,:,:,:)
  real, intent(in) :: dtNew
  integer, intent(out) :: now_num
  real, dimension(MDIM,pt_maxnewnum), intent(out) :: now_pos, now_vel
  real, dimension(pt_maxnewnum), intent(out) :: now_time, now_energy, now_weight

  ! aux parameters
  integer, dimension(2,MDIM) :: faces, onBoundary
  logical, dimension(2,MDIM) :: isEmits
  real, dimension(LOW:HIGH, MDIM) :: bndBox
  real, dimension(MDIM) :: deltaCell

  real :: FaceArea
  real :: r_in, theta_in, theta_out, phi_in, phi_out
  integer :: i, ii, jj
  real :: flux, dE
  real, parameter :: eps_dummy = 1.0d-11

  real, dimension(MDIM) :: newxyz, newvel, cart_pos, r_hat
  integer, dimension(MDIM) :: cellID
  real :: weight_per_mcp, dE_per_mcp, dtNow, newenergy


  ! Initialization
  now_num = 0
  now_pos = 0.0d0
  now_vel = 0.0d0
  now_time = 0.0d0
  now_energy = 0.0d0
  now_weight = 0.0d0

  ! Quit if not including face emission
  if (.NOT. pt_FaceEmission) return

  if (pt_is_FacePlanck) then
    flux = sigma * (pt_FacePlanckTemp**4)
  else
    flux = pt_constFaceFlux
  end if

  ! Obtain block info
  call Grid_getBlkBoundBox(blockID, bndBox)
  call Grid_getDeltas(blockID, deltaCell)

  call Grid_getBlkBC(blockID, faces, onBoundary)

  ! Activate the chosen face
  isEmits = .false.
  if (onBoundary(pt_FaceEmissionSide, pt_FaceEmissionAxis) /= NOT_BOUNDARY) then
    isEmits(pt_FaceEmissionSide, pt_FaceEmissionAxis) = .true.

    ! Abort if not emitting from inner radial face
    if ((pt_FaceEmissionSide /= IAXIS) .or. (pt_FaceEmissionAxis /= 1)) then
      call Driver_abortFlash("face_emission: non-radial face emission&
                              not yet implemented.&
                              Check pt_FaceEmissionSide's value.")
    end if
  end if

  do ii = LOW, HIGH
    do jj = IAXIS, KAXIS

      ! Only select the active face
      if (isEmits(ii, jj)) then
        call get_face_area(ii, jj, bndBox, FaceArea)

        dE = flux * FaceArea * dtNew

        dE_per_mcp = dE / pt_num_fmcps_tstep

        do i = 1, pt_num_fmcps_tstep
          call sample_time(dtNew, dtNow)

          ! newxyz is in spherical cooridnates here
          call sample_blk_position(bndBox, newxyz)
          ! Impose the face value along FaceAxis
          newxyz(jj) = bndBox(ii, jj)
          if (ii .EQ. HIGH) newxyz(jj) = (1.0d0 - pt_smlpush)*newxyz(jj)
          if (ii .EQ. LOW) newxyz(jj) = (1.0d0 + pt_smlpush)*newxyz(jj)

          if (pt_is_radial_face_vel) then
            call get_cartesian_position(newxyz, cart_pos) 
            r_hat = cart_pos / sqrt(dot_product(cart_pos, cart_pos))
            newvel = clight * r_hat 
          else if (pt_is_iso_face_vel) then
            call sample_iso_velocity(newvel)
          else if (pt_is_therm_face_vel) then
            call get_cartesian_position(newxyz, cart_pos) 
            r_hat = cart_pos / sqrt(dot_product(cart_pos, cart_pos))
            call sample_therm_face_velocity(r_hat, newxyz(JAXIS), newvel)
          end if

          call get_cellID(bndBox, deltaCell, newxyz, cellID)

          call sample_energy(solnVec, cellID, newenergy)
          weight_per_mcp = dE_per_mcp / newenergy

          ! Record the MCP attributes
          now_time(now_num + i) = dtNow
          now_pos(:, now_num + i) = newxyz
          now_energy(now_num + i) = newenergy
          now_vel(:, now_num + i) = newvel
          now_weight(now_num + i) = weight_per_mcp

          ! End of sampling the face MCPs
        end do

        now_num = now_num + pt_num_fmcps_tstep
      end if
    end do
  end do

end subroutine face_emission


subroutine get_face_area(dir, axis, bndBox, FaceArea)
  use Driver_interface, only : Driver_abortFlash
  use Grid_data, only: gr_geometry
  implicit none
#include "constants.h"
#include "Flash.h"

  ! Input/Output
  integer, intent(in) :: dir, axis
  real, dimension(LOW:HIGH, MDIM) :: bndBox
  real, intent(out)   :: FaceArea

  ! aux variables
  integer :: i
  real, dimension(MDIM) :: blockSize
  real :: r_in, theta_in, theta_out, phi_in, phi_out

  ! Initialization
  FaceArea = 1.0d0
  if (gr_geometry == CARTESIAN) then
    blockSize(:) = bndBox(HIGH,:) - bndBox(LOW,:)
    ! FaceArea = 1.0 if (NDIM == 1)

    if (NDIM .eq. 2) then
      blockSize(KAXIS) = 1.0 ! set blockSize(KAXIS) to one for 2D
      do i = IAXIS, JAXIS
        if (i .NE. axis) FaceArea = FaceArea * blockSize(i)
      end do
    end if

    if (NDIM .eq. 3) then
      do i = IAXIS, KAXIS
        if (i .NE. axis) FaceArea = FaceArea * blockSize(i)
      end do
    end if

  else if (gr_geometry == SPHERICAL) then
    r_in = bndBox(LOW, IAXIS)

    theta_in  = bndBox(LOW, JAXIS)
    theta_out = bndBox(HIGH, JAXIS)

    phi_in  = bndBox(LOW, KAXIS)
    phi_out = bndBox(HIGH, KAXIS)


    if (NDIM == 3) then
      FaceArea = FaceArea * r_in*r_in&
                          * (phi_out - phi_in)&
                          * (cos(theta_in) - cos(theta_out))
    else
      call Driver_abortFlash("face_emission: not yet implemented&
                              for 1D/2D simulations.")
    end if
  else
    call Driver_abortFlash("get_face_area: unknown geometry requested.")
  end if

end subroutine get_face_area


subroutine comp_fleck_factor(kp, beta, dt, fn)
  use Particles_data, only : pt_es_alpha
  use Simulation_data, only : clight
  implicit none

  real, intent(in) :: kp, beta, dt
  real, intent(out) :: fn

  fn = 1.0d0 / (1.0d0 + (pt_es_alpha * dt * beta * clight * kp))

end subroutine comp_fleck_factor


! Zeroing the deposition variables on the grid
subroutine reset_deposition_var(solnVec)

  implicit none

#include "constants.h"
#include "Flash.h"

  real, pointer :: solnVec(:,:,:,:)

  solnVec(URAD_VAR,:,:,:) = 0.0
  solnVec(ABSE_VAR,:,:,:) = 0.0
  solnVec(EMIE_VAR,:,:,:) = 0.0

  solnVec(ABMX_VAR,:,:,:) = 0.0
  solnVec(ABMY_VAR,:,:,:) = 0.0
  solnVec(ABMZ_VAR,:,:,:) = 0.0

  solnVec(SCMX_VAR,:,:,:) = 0.0
  solnVec(SCMY_VAR,:,:,:) = 0.0
  solnVec(SCMZ_VAR,:,:,:) = 0.0

  solnVec(FLEC_VAR,:,:,:) = 1.0
#ifdef FLEP_VAR
  solnVec(FLEP_VAR,:,:,:) = 1.0
#endif
#ifdef IONR_VAR
  solnVec(IONR_VAR,:,:,:) = 0.0
#endif
#ifdef RECR_VAR
  solnVec(RECR_VAR,:,:,:) = 0.0
#endif
#ifdef HEAT_VAR
  solnVec(HEAT_VAR,:,:,:) = 0.0
#endif
#ifdef COOL_VAR
  solnVec(COOL_VAR,:,:,:) = 0.0
#endif

end subroutine reset_deposition_var

end module emission
