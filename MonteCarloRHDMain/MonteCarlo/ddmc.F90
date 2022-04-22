module ddmc

  implicit none
  contains

! subroutine to populate the SCRATCH face variables
! for leakage opacity (KLEAX, KLEAY, KLEAZ in Config)
subroutine calc_block_leakage_opac(blockID, scratchVec, solnVec)
  use Particles_data, only : pt_is_ddmc, pt_ddmc_tau_thres, pt_use_fromPos
  use opacity, only : calc_abs_opac, calc_sca_opac
  use Grid_interface, only : Grid_getBlkIndexLimits, Grid_getDeltas,&
                             Grid_getBlkBoundBox, Grid_getBlkBC,&
                             Grid_getCellCoords, Grid_getBlkIDFromPos
  use Paramesh_comm_data, only : amr_mpi_meshComm
  use gr_interface, only : gr_findNeghID
  implicit none
#include "constants.h"
#include "Flash.h"

  ! Input/Output
  integer, intent(in) :: blockID
  real, pointer :: scratchVec(:,:,:,:)
  real, pointer :: solnVec(:,:,:,:)

  ! aux variables
  integer, dimension(MDIM) :: cellID, cellID_negh
  integer, dimension(2, MDIM) :: blkLimits, blkLimitsGC
  real, dimension(LOW:HIGH, MDIM) :: bndBox
  real, dimension(MDIM) :: deltaCell
  integer :: i, j, k
  real, dimension(LOW:HIGH, MDIM) :: k_leak
  integer, dimension(LOW:HIGH,MDIM) :: faces
  integer, dimension(LOW:HIGH,MDIM) :: onBoundary
  integer :: axis

  if (.not. pt_is_ddmc) return

#ifdef FLASH_GREY_DDMC
  ! for grey DDMC
  ! DDMC is on, gather grid size info
  call Grid_getBlkIndexLimits(blockID, blkLimits, blkLimitsGC)
  call Grid_getBlkBoundBox(blockID, bndBox)
  call Grid_getDeltas(blockID, deltaCell)

  ! for correcting reflective BC
  call Grid_getBlkBC(blockID, faces, onBoundary)

  ! zero out all the leakage opacities
  scratchVec(KLEX_SCRATCH_GRID_VAR,:,:,:) = 0.0
  scratchVec(KLEY_SCRATCH_GRID_VAR,:,:,:) = 0.0
  scratchVec(KLEZ_SCRATCH_GRID_VAR,:,:,:) = 0.0

  ! Begin loop over interior cells
  ! The +/- 1 limits are to cover the faces at domain boundary
  ! so guard cells also have TRAP populated
  do k = blkLimits(LOW,KAXIS)-1, blkLimits(HIGH,KAXIS)+1
    do j = blkLimits(LOW, JAXIS)-1, blkLimits(HIGH, JAXIS)+1
      do i = blkLimits(LOW,IAXIS)-1, blkLimits(HIGH, IAXIS)+1
        cellID = (/ i, j, k /)

        ! dummy energy and dshift for grey case
        call calc_cell_leakage_opacities(blockID, cellID, solnVec,&
                                         bndBox, deltaCell,&
                                         1.0, 1.0, k_leak)

        ! assign k_leak to scratch grid variables
        scratchVec(KLEX_SCRATCH_GRID_VAR, i, j, k)   = k_leak(LOW,  IAXIS)
        scratchVec(KLEX_SCRATCH_GRID_VAR, i+1, j, k) = k_leak(HIGH, IAXIS)
        scratchVec(KLEY_SCRATCH_GRID_VAR, i, j, k)   = k_leak(LOW,  JAXIS)
        scratchVec(KLEY_SCRATCH_GRID_VAR, i, j+1, k) = k_leak(HIGH, JAXIS)
        scratchVec(KLEZ_SCRATCH_GRID_VAR, i, j, k)   = k_leak(LOW,  KAXIS)
        scratchVec(KLEZ_SCRATCH_GRID_VAR, i, j, k+1) = k_leak(HIGH, KAXIS)
      end do
    end do
  end do

  ! zero out leakage opacities towards reflective boundaries
  do axis = 1, NDIM
    if (onboundary(LOW, axis) == REFLECTING) then
      ! set scratch(5,:,:) = 0.0
      if (axis == IAXIS) scratchVec(KLEX_SCRATCH_GRID_VAR, GRID_ILO, :, :) = 0.0
      if (axis == JAXIS) scratchVec(KLEY_SCRATCH_GRID_VAR, :, GRID_JLO, :) = 0.0
      if (axis == KAXIS) scratchVec(KLEZ_SCRATCH_GRID_VAR, :, :, GRID_KLO) = 0.0
    else if (onboundary(HIGH, axis) == REFLECTING) then
      ! set scratch(12,:,:) = 0.0
      if (axis == IAXIS) scratchVec(KLEX_SCRATCH_GRID_VAR, GRID_IHI+1, :, :) = 0.0
      if (axis == JAXIS) scratchVec(KLEY_SCRATCH_GRID_VAR, :, GRID_JHI+1, :) = 0.0
      if (axis == KAXIS) scratchVec(KLEZ_SCRATCH_GRID_VAR, :, :, GRID_KHI+1) = 0.0
    end if
  end do
#endif

end subroutine calc_block_leakage_opac

! subroutine to get local cell widths of a cell
subroutine calc_cell_widths(cellID, bndBox, deltaCell, cell_widths)

  use Grid_data, only : gr_geometry
  implicit none

#include "constants.h"
#include "Flash.h"

  ! Input/outputs
  integer, dimension(MDIM), intent(in) :: cellID
  real, dimension(LOW:HIGH, MDIM), intent(in) :: bndBox
  real, dimension(MDIM), intent(in) :: deltaCell
  real, dimension(MDIM), intent(out) :: cell_widths

  ! aux variables
  real :: dx, r0 ! for spherical coord.
  real, dimension(LOW:HIGH) :: r_bnd
  real :: theta, sintheta
  integer :: axis

  select case (gr_geometry)

    case (CARTESIAN)
      cell_widths = deltaCell

    case (SPHERICAL)
      r_bnd(LOW)  = bndBox(LOW, IAXIS) +&
                      (cellID(IAXIS) - NGUARD - 1)*deltaCell(IAXIS)
      r_bnd(HIGH) = r_bnd(LOW) + deltaCell(IAXIS)

      r0 = 0.5*sum(r_bnd)  ! center radial center

      theta = bndBox(LOW, JAXIS) +&
                      (cellID(JAXIS) - NGUARD - 1)*deltaCell(JAXIS)+&
                      0.5*deltaCell(JAXIS)
      sintheta = sin(theta)

      do axis = IAXIS, NDIM ! over valid axes
        dx = deltaCell(axis)
        if (axis == JAXIS) then ! r*dtheta
          dx = dx*r0
        else if (axis == KAXIS) then ! r*sin(theta)*dphi
          dx = r0*sintheta*dx
        end if

        cell_widths(axis) = dx
      end do
  end select

end subroutine calc_cell_widths

! subroutine to compute leakage opacities for one cell
! the energy and dshift inputs are for non-grey, MCP-based case.
subroutine calc_cell_leakage_opacities(blockID, cellID, solnVec,&
                                       bndBox, deltaCell,&
                                       energy, dshift, k_leak)
  use Particles_data, only : pt_is_ddmc, pt_ddmc_tau_thres,&
                             pt_ddmc_lambda, pt_meshMe,&
                             pt_use_fromPos
  use Grid_data, only : gr_geometry
  use opacity, only : calc_abs_opac, calc_sca_opac
  use Grid_interface, only : Grid_getCellCoords, Grid_getBlkIDFromPos
  use Paramesh_comm_data, only : amr_mpi_meshComm
  use gr_interface, only : gr_findNeghID
  implicit none

#include "constants.h"
#include "Flash.h"

  ! Input/Output
  integer, intent(in) :: blockID
  integer, dimension(MDIM), intent(in) :: cellID
  real, pointer :: solnVec(:,:,:,:)
  real, dimension(LOW:HIGH, MDIM), intent(in) :: bndBox
  real, dimension(MDIM), intent(in) :: deltaCell
  real, intent(in) :: energy, dshift
  real, dimension(LOW:HIGH, MDIM), intent(out) :: k_leak

  ! aux variables
  integer, dimension(LOW:HIGH), parameter :: dir2dx = (/ -1, 1 /)
  integer :: axis, dir
  integer, dimension(MDIM) :: cellID_negh
  real :: ka, ka_negh, ks, ks_negh, ktot, ktot_negh
  real :: dx, tau, dx_negh, tau_negh
  real :: Xi, r0, geo_factor ! for spherical coord.
  real, dimension(LOW:HIGH) :: r_bnd
  real :: theta, sintheta
  real, dimension(MDIM) :: tau_3d
  logical, dimension(MDIM) :: tau_mask
  integer :: num_ddmc_axes

  ! bnd determination
  integer :: local_trap_index, negh_trap_index
  real :: local_dx_corr, negh_dx_corr

  k_leak = 0.0
  geo_factor = 1.0
  if (.not. pt_is_ddmc) return

#ifdef FLASH_GREY_DDMC
  ! x direction, current cell
  call calc_abs_opac(cellID, solnVec, energy, dshift, ka)
  call calc_sca_opac(cellID, solnVec, energy, dshift, ks)
  ktot = ka + ks

  ! Compute geometric factor for spherical coord.
  r_bnd(LOW)  = bndBox(LOW, IAXIS) + (cellID(IAXIS) - NGUARD - 1)*deltaCell(IAXIS)
  r_bnd(HIGH) = r_bnd(LOW) + deltaCell(IAXIS)
  r0 = 0.5*sum(r_bnd)  ! center center
  dx = deltaCell(IAXIS) ! get dr for geometric factor
  Xi = 1.0 + dx*dx/(12.0*r0*r0)

  if (gr_geometry == SPHERICAL) then
    theta = bndBox(LOW, JAXIS) + (cellID(JAXIS) - NGUARD - 1)*deltaCell(JAXIS) +&
                   0.5*deltaCell(JAXIS)
    sintheta = sin(theta)
  end if

  tau_3d = 0.0 ! not active by default
  do axis = IAXIS, NDIM ! over valid axes

    ! current cell width
    dx = deltaCell(axis)
    if (gr_geometry == SPHERICAL) then
      if (axis == JAXIS) then ! r*dtheta
        dx = dx*r0
      else if (axis == KAXIS) then ! r*sin(theta)*dphi
        dx = r0*sintheta*dx
      end if
    end if

    ! also query trap index in case local cell is a guard cell
    call get_negh_trap_index(blockID, bndBox, cellID, local_trap_index)
    local_dx_corr = 1.0
    ! if index = +/-1, correct the dx_negh to reflect
    if (local_trap_index == 1) then
      local_dx_corr = 0.5
    else if (local_trap_index == -1) then
      local_dx_corr = 2.0
!    else if (local_trap_index < -1) then
!      local_dx_corr = 0.0 ! mark as non-DDMC
    end if

    tau = ktot * dx * local_dx_corr
    tau_3d(axis) = tau
  
    !if (tau < pt_ddmc_tau_thres) cycle ! no leakage from this axis

    do dir = LOW, HIGH  ! loop over two neghbors
      cellID_negh = cellID
      cellID_negh(axis) = cellID_negh(axis) + dir2dx(dir)

      call calc_abs_opac(cellID_negh, solnVec, energy, dshift, ka_negh)
      call calc_sca_opac(cellID_negh, solnVec, energy, dshift, ks_negh)

      ktot_negh = ka_negh + ks_negh
      dx_negh = deltaCell(axis)
      if (gr_geometry == SPHERICAL) then
        if (axis == JAXIS) then ! r*dtheta
          dx_negh = dx_negh*r0
        else if (axis == KAXIS) then ! r*sin(theta)*dphi
          dx_negh = r0*sintheta*dx
        end if
      end if

      ! insert boundary query here
      call get_negh_trap_index(blockID, bndBox, cellID_negh, negh_trap_index)
      negh_dx_corr = 1.0
      ! if index = +/-1, correct the dx_negh to reflect
      if (negh_trap_index == 1) then
        negh_dx_corr = 0.5
      else if (negh_trap_index == -1) then
        negh_dx_corr = 2.0
!      else if (negh_trap_index < -1) then
!        negh_dx_corr = 0.0
      end if

      tau_negh = ktot_negh * dx_negh * negh_dx_corr

      if (tau_negh < pt_ddmc_tau_thres) then
        tau_negh =  2.0*pt_ddmc_lambda ! override with 
      end if

      ! geometric factor for radial direction in spherical coord.
      if ((gr_geometry == SPHERICAL) .and. (axis == IAXIS)) then
        geo_factor = (r_bnd(dir)/r0)**2/Xi
      end if

      ! Compute leakage opacity
      k_leak(dir, axis) = (2.0/3.0/dx/local_dx_corr)*(1.0/(tau + tau_negh))*geo_factor

      ! debugging density jump case
      !if ((axis == IAXIS) .and. (dir == HIGH) .and. (cellID(1) == 12) .and.&
      !    (tau > pt_ddmc_tau_thres) .and. (tau_negh < pt_ddmc_tau_thres)) then
      !  print *, "cellID", cellID
      !  print *, "cellID_negh", cellID_negh
      !  print *, "dens", solnVec(DENS_VAR, cellID(1):cellID(1)+1, cellID(2), cellID(3))
      !  print *, "kleakCALC", k_leak(:,IAXIS)
      !end if
    end do ! left and right neighbor
  end do ! up to NDIM

  ! DDMC-active if number of axis exceeding tau_thres == NDIM
  tau_mask = (tau_3d >= pt_ddmc_tau_thres)
  num_ddmc_axes = count(tau_mask, 1)
  !if (num_ddmc_axes == NDIM) then ! leakage possible from this axis
  if ((num_ddmc_axes == NDIM) .and. (abs(negh_trap_index) <= 1)) then ! leakage possible from this axis
    solnVec(TRAP_VAR, cellID(IAXIS), cellID(JAXIS), cellID(KAXIS)) = 1.0
    solnVec(TAUX_VAR, cellID(IAXIS), cellID(JAXIS), cellID(KAXIS)) = tau_3d(1)
    solnVec(TAUY_VAR, cellID(IAXIS), cellID(JAXIS), cellID(KAXIS)) = tau_3d(2)
    solnVec(TAUZ_VAR, cellID(IAXIS), cellID(JAXIS), cellID(KAXIS)) = tau_3d(3)
  end if
#endif

end subroutine calc_cell_leakage_opacities


! subroutine to get neighboring cell's trap/refinement condition
subroutine get_negh_trap_index(blockID, bndBox, negh_cellID, negh_trap_index)
  use Particles_data, only : pt_meshMe, pt_use_fromPos
  use Grid_interface, only : Grid_getCellCoords, Grid_getBlkIDFromPos,&
                             Grid_getBlkRefineLevel, Grid_getBlkNeighLevels
  use Paramesh_comm_data, only : amr_mpi_meshComm
  use gr_interface, only : gr_findNeghID
  use Driver_interface,  ONLY : Driver_abortFlash
  use tree, only : lrefine_max
  implicit none

#include "constants.h"
#include "Flash.h"

  ! Input/Output
  integer, intent(in) :: blockID
  real, dimension(LOW:HIGH, MDIM), intent(in) :: bndBox
  integer, dimension(MDIM), intent(in) :: negh_cellID
  integer, intent(out) :: negh_trap_index

  ! bnd determination
  integer, parameter :: ycoordsize = NYB + 2*NGUARD
  real, dimension(ycoordsize) :: ycoords
  integer, parameter :: xcoordsize = NXB + 2*NGUARD
  real, dimension(xcoordsize) :: xcoords
  integer, parameter :: zcoordsize = NZB + 2*NGUARD
  real, dimension(zcoordsize) :: zcoords
  real, dimension(MDIM) :: negh_coords
  logical :: isoutside
  integer, dimension(MDIM) :: neghdir
  integer,dimension(BLKNO:PROCNO) :: neghID
  integer :: old_rflvl, new_rflvl
  integer :: levels(-1:1, -K2D:K2D , -K3D:K3D)

  ! get cell info for neighbor case
  ! cell locations for checking neghbor's refinement level
  ! xyz here can be Cartesian or Spherical
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

  ! get neghbor cell's coordinates
  ! gather cell position
  negh_coords = (/ xcoords(negh_cellID(IAXIS)),&
                   ycoords(negh_cellID(JAXIS)),&
                   zcoords(negh_cellID(KAXIS)) /)

  call Grid_outsideBoundBox(negh_coords, bndBox,&
                            isoutside, neghdir)

  ! if outside (i.e. in a guard cell), check proc,
  ! if same proc check refine_level
  ! only allow TRAP if same proc and same refinelevel
  negh_trap_index = 0
!  if (isoutside) then
!    if (pt_use_fromPos) then
!      call Grid_getBlkIDFromPos(negh_coords,&
!                                neghID(BLKNO), neghID(PROCNO),&
!                                amr_mpi_meshComm)
!    else
!      call gr_findNeghID(blockID, negh_coords, neghdir, neghID)
!    end if

!    if (neghID(PROCNO) == pt_meshMe) then
!      call Grid_getBlkRefineLevel(blockID, old_rflvl)
!      call Grid_getBlkRefineLevel(neghID(BLKNO), new_rflvl)

!      negh_trap_index = new_rflvl - old_rflvl
!    else
!      negh_trap_index = -99 ! if not 
!    end if
!  end if

  ! new implementation to use neghdir and Grid_getBlkNeighLevels
  if (isoutside) then
    call Grid_getBlkNeighLevels(blockID, levels)

    old_rflvl = levels(0,0,0)  ! local block refinement
    new_rflvl = levels(neghdir(IAXIS)-2, neghdir(JAXIS)-2, neghdir(KAXIS)-2)
    if ((new_rflvl < 1) .or. (new_rflvl > lrefine_max)) then
      call Driver_abortFlash("get_negh_trap_index: unknown new_rflvl.")
    end if

    negh_trap_index = new_rflvl - old_rflvl
  end if

end subroutine get_negh_trap_index


! subroutine to retrieve leakage opacity from scratch face variable
subroutine get_grey_leakage_opac(cellID, scratchVec, k_leak)

  implicit none
#include "constants.h"
#include "Flash.h"

  ! Input/Output
  integer, dimension(MDIM) :: cellID
  real, pointer :: scratchVec(:,:,:,:)
  real, dimension(LOW:HIGH, MDIM) :: k_leak

#ifdef FLASH_GREY_DDMC
  ! directly copy from scratch grid for leakage opacity
  k_leak(LOW:HIGH, IAXIS) = scratchVec(KLEX_SCRATCH_GRID_VAR,&
                                       cellID(IAXIS):cellID(IAXIS)+1,&
                                       cellID(JAXIS),&
                                       cellID(KAXIS))
  k_leak(LOW:HIGH, JAXIS) = scratchVec(KLEY_SCRATCH_GRID_VAR,&
                                       cellID(IAXIS),&
                                       cellID(JAXIS):cellID(JAXIS)+1,&
                                       cellID(KAXIS))
  k_leak(LOW:HIGH, KAXIS) = scratchVec(KLEZ_SCRATCH_GRID_VAR,&
                                       cellID(IAXIS),&
                                       cellID(JAXIS),&
                                       cellID(KAXIS):cellID(KAXIS)+1)
#endif

end subroutine get_grey_leakage_opac

! subroutine to leak mcp to a neighboring cell
subroutine leak_mcp_to_neighbor(bndBox, deltaCell, target_cellID, now_cellID,&
                                is_trap_negh, leak_axis, leak_dir, particle)
  use Particles_data, only : pt_smlpush
  use new_mcp, only : sample_cell_position, sample_iso_velocity
  implicit none
#include "constants.h"
#include "Flash.h"

  ! Input/Output
  real, dimension(LOW:HIGH, MDIM), intent(in) :: bndBox
  real, dimension(MDIM), intent(in) :: deltaCell
  integer, dimension(MDIM), intent(in) :: target_cellID, now_cellID
  logical, intent(in) :: is_trap_negh
  integer, intent(out) :: leak_axis, leak_dir
  real, dimension(NPART_PROPS), intent(inout) :: particle

  ! aux variables
  real, dimension(MDIM) :: newxyz, newvel
  real :: leak_bnd

  call sample_cell_position(bndBox, deltaCell, target_cellID, newxyz)
  call sample_iso_velocity(newvel)
  call get_leak_info(now_cellID, target_cellID, leak_axis, leak_dir)

  ! neighbor is not DDMC-active, fix position and velocity
  if (.not. is_trap_negh) then
    leak_bnd = now_cellID(leak_axis) - NGUARD + minval((/ 0, leak_dir /), 1)&
                + sign(pt_smlpush, real(leak_dir))
    leak_bnd = bndBox(LOW, leak_axis) + leak_bnd * deltaCell(leak_axis)
    newxyz(leak_axis) = leak_bnd

    !call sample_velocity_ddmc2imc(leak_dir, leak_axis, newvel)
    call sample_velocity_ddmc2imc_general(now_cellID, leak_axis, leak_dir,&
                                            bndBox, deltaCell,&
                                            newvel)

    ! debugging
    !print *, "leaking to MC", leak_axis, leak_dir
    !print *, "newxyz", newxyz
    !print *, "newvel", newvel
  end if

  particle(POSX_PART_PROP:POSZ_PART_PROP) = newxyz
  particle(VELX_PART_PROP:VELZ_PART_PROP) = newvel

end subroutine leak_mcp_to_neighbor


! subroutine to construct leak_axis and leak_dir from cellID and xcellID
subroutine get_leak_info(cellID, xcellID, leak_axis, leak_dir)
  use Driver_interface,  ONLY : Driver_abortFlash
  implicit none
#include "constants.h"
#include "Flash.h"

  ! Input/Output
  integer, dimension(MDIM), intent(in) :: cellID, xcellID
  integer, intent(out) :: leak_axis, leak_dir

  ! aux variables
  integer, dimension(MDIM) :: delta_cellID
  logical, dimension(MDIM) :: leak_mask
  integer :: i, num_leak_axes

  delta_cellID = 0
  delta_cellID = xcellID - cellID

  leak_mask = (delta_cellID /= 0)
  num_leak_axes = count(leak_mask, 1)

  ! assert one leakage direction
  if (num_leak_axes /= 1) then
    call Driver_abortFlash("get_leak_info: leakage in more than one directions!")
  end if

  do i = 1, MDIM
    if (delta_cellID(i) /= 0) then
      leak_axis = i
      leak_dir  = sign(1, delta_cellID(i))
    end if
  end do 

end subroutine get_leak_info


!!! Subroutine to sample the asymptotically correct angular distribution
!!! for IMC-DDMC interfaces.
!!! Based on Davidson et al. 2006
subroutine sample_velocity_ddmc2imc(dir, axis, velvec)
  use Simulation_data, only : clight
  use Particles_data, only : pt_is_ddmc
  use random, only : rand
  use Driver_interface, only : Driver_abortFlash
  implicit none

#include "Flash.h"
#include "constants.h"

  ! Input/Output
  integer, intent(in) :: dir, axis
  real, dimension(MDIM), intent(out) :: velvec

  ! Temporary variables
  real :: mu, xi1, xi2
  real :: theta, costheta, sintheta
  real :: phi,   cosphi,   sinphi

  real :: v_forward, v2, v3

  ! Initialization of output
  velvec = 0.0d0

  if (pt_is_ddmc) then
    ! First, sample the component normal to the face
    xi1 = rand()
    xi2 = rand()

    if (xi1 .lt. 0.5d0) then
      mu = sqrt(xi2)
    else
      mu = xi2 ** (1.0/3.0)
    endif
    ! mu must be positive here

    if (dir .eq. -1) mu = -mu

    theta = acos(mu)
    costheta = cos(theta)
    sintheta = sin(theta)

    phi    = 2.0d0 * PI * rand()
    cosphi = cos(phi)
    sinphi = sin(phi)

    ! Constructing the output velocity
    v_forward = clight * mu
    v2 = clight * sintheta * cosphi
    v3 = clight * sintheta * sinphi

    velvec(axis) = v_forward

    ! The order of the non-forward direction does not matter
    if (axis .EQ. IAXIS) then
      velvec(JAXIS) = v2
      velvec(KAXIS) = v3
    else if (axis .EQ. JAXIS) then
      velvec(IAXIS) = v2
      velvec(KAXIS) = v3
    else if (axis .EQ. KAXIS) then
      velvec(IAXIS) = v2
      velvec(JAXIS) = v3
    else
      call Driver_abortFlash("Error: Unknown axis for DDMC2IMC vel. sampling")
    end if

  else
    call Driver_abortFlash("Error: DDMC2IMC angular sampling, but DDMC is off.")
  end if

  !if (mu .gt. 0.0d0)  print *, "newvell", dir, -mu
end subroutine sample_velocity_ddmc2imc

! subroutine basically the same as sample_velocity_ddmc2imc
subroutine sample_velocity_ddmc2imc_general(cellID, axis, dir,&
                                            bndBox, deltaCell,&
                                            velvec)
  use Simulation_data, only : clight
  use Grid_data, only : gr_geometry
  use Particles_data, only : pt_is_ddmc
  use random, only : rand
  use Driver_interface, only : Driver_abortFlash
  implicit none

#include "Flash.h"
#include "constants.h"

  ! Input/Output
  integer, dimension(MDIM), intent(in) :: cellID
  integer, intent(in) :: axis, dir
  real, dimension(LOW:HIGH, MDIM), intent(in) :: bndBox
  real, dimension(MDIM), intent(in) :: deltaCell
  real, dimension(MDIM), intent(out) :: velvec

  ! Temporary variables
  real :: mu, xi1, xi2
  real :: radius
  real :: theta, costheta, sintheta
  real :: phi,   cosphi,   sinphi
  real, dimension(MDIM) :: r_hat, theta_hat, phi_hat

  real :: v_forward, v2, v3

  ! Initialization of output
  velvec = 0.0d0

  ! First, sample the component normal to the face
  xi1 = rand()
  xi2 = rand()

  if (xi1 .lt. 0.5d0) then
    mu = sqrt(xi2)
  else
    mu = xi2 ** (1.0/3.0)
  endif
  ! mu must be positive here

  !if (dir .eq. -1) mu = -mu
  ! forward must be positive

  theta = acos(mu)
  costheta = cos(theta)
  sintheta = sin(theta)

  phi    = 2.0d0 * PI * rand()
  cosphi = cos(phi)
  sinphi = sin(phi)

  ! Constructing the output velocity
  v_forward = clight * mu
  v2 = clight * sintheta * cosphi
  v3 = clight * sintheta * sinphi

  ! infer surface normal's unit vector
  select case(gr_geometry)

    case (CARTESIAN)
      ! negative if leaking left
      velvec(axis) = real(dir)*v_forward

      ! The order of the non-forward direction does not matter
      if (axis .EQ. IAXIS) then
        velvec(JAXIS) = v2
        velvec(KAXIS) = v3
      else if (axis .EQ. JAXIS) then
        velvec(IAXIS) = v2
        velvec(KAXIS) = v3
      else if (axis .EQ. KAXIS) then
        velvec(IAXIS) = v2
        velvec(JAXIS) = v3
      else
        call Driver_abortFlash("Error: Unknown axis for DDMC2IMC vel. sampling")
      end if

    case (SPHERICAL)

      ! compute cell center coordinate
      radius = bndBox(LOW, IAXIS) +&
                 (cellID(IAXIS) - NGUARD - 1)*deltaCell(IAXIS)+&
                   0.5*deltaCell(IAXIS)
      theta  = bndBox(LOW, JAXIS) +&
                 (cellID(JAXIS) - NGUARD - 1)*deltaCell(JAXIS)+&
                   0.5*deltaCell(JAXIS)
      phi    = bndBox(LOW, KAXIS) +&
                 (cellID(KAXIS) - NGUARD - 1)*deltaCell(KAXIS)+&
                   0.5*deltaCell(KAXIS)

      sintheta = sin(theta)
      costheta = cos(theta)
      sinphi   = sin(phi)
      cosphi   = cos(phi)

      r_hat = (/ sintheta*cosphi,&
                 sintheta*sinphi,&
                 costheta /)

      theta_hat = (/ costheta*cosphi,&
                     costheta*sinphi,&
                     -sintheta /)

      phi_hat   = (/ -sinphi,&
                     cosphi,&
                     0.0d0 /)

      if (axis == IAXIS) then

        velvec = real(dir)*v_forward*r_hat +&
                 v2*theta_hat +&
                 v3*phi_hat
      else if (axis == JAXIS) then
        velvec = v2*r_hat +&
                 real(dir)*v_forward*theta_hat +&
                 v3*phi_hat
      else if (axis == KAXIS) then
        velvec = v2*r_hat +&
                 v3*theta_hat +&
                 real(dir)*v_forward*phi_hat
      end if

  end select
  !call Driver_abortFlash("Error: DDMC2IMC angular sampling, but DDMC is off.")

end subroutine sample_velocity_ddmc2imc_general



! subroutine to perform energy/urad deposition in case of DDMC leakage
subroutine deposit_energy_ddmc(solnVec, cellID, v_gas, particle, mcp_fate,&
                                dt, dtNew, fleck, k_a, k_s, dvol)
  use Simulation_data, only : clight
  use Particles_data, only : pt_is_corrdl, pt_is_deposit_urad,&
                             pt_is_deposit_energy, pt_ABS_ID,&
                             pt_is_photoionization
  use rhd, only : cellAddVar
  implicit none
#include "Flash.h"
#include "constants.h"

  ! Input/Output
  real, pointer :: solnVec(:,:,:,:)
  integer, dimension(MDIM), intent(in) :: cellID
  real, dimension(MDIM), intent(in) :: v_gas
  real, dimension(NPART_PROPS), intent(inout) :: particle
  integer, intent(in) :: mcp_fate
  real, intent(in) :: dt, dtNew, fleck, k_a, k_s, dvol

  ! aux variables
  real :: rho
  real :: k_ea, k_es, dl_corr
  real :: old_weight, new_weight, avg_weight, mcp_eps, mcp_energy
  real :: dl, dtau, delta_e, u_rad

  rho = solnVec(DENS_VAR, cellID(IAXIS), cellID(JAXIS), cellID(KAXIS))

  old_weight  = particle(NUMP_PART_PROP) 
  mcp_eps     = particle(ENER_PART_PROP)

  ! Inputs of k_a and k_s are assumed to be in the comoving frame
  k_ea = fleck * k_a
  k_es = (1.0d0 - fleck) * k_a

  dl = clight * dt !* dshift ! convert dt to comoving frame
  dtau = k_ea * dl

  ! Correction for absorption during the integrated path 
  dl_corr = 1.0d0
  if (pt_is_corrdl .and. (k_ea > 0.0)) then
    dl_corr = (1.0d0 - exp(-dtau)) / dtau
  end if

  dl = dl
  dtau = dtau

  if (pt_is_deposit_energy) then
    new_weight = old_weight * exp(-dtau)

    if ((mcp_fate == pt_ABS_ID) .and. (.not. pt_is_photoionization)) then
      new_weight = 0.0d0
    end if

    ! no dshift needed, mcp_eps already in comoving frame
    ! dvol combined with dtnew (in lab, in rhd.F90) = dvol_0 * dt_0
    delta_e = (mcp_eps * (old_weight - new_weight)) / (rho * dvol)
    call cellAddVar(solnVec, cellID, ABSE_VAR, delta_e)

    ! Update MCP weights
    particle(NUMP_PART_PROP) = new_weight
    particle(NUM0_PART_PROP) = new_weight
  end if

  if (pt_is_deposit_urad) then
    avg_weight = old_weight*dl_corr
    mcp_energy = mcp_eps * avg_weight
    u_rad = (mcp_energy / dvol) * (dl / (clight * dtNew))
    call cellAddVar(solnVec, cellID, URAD_VAR, u_rad)
  end if

end subroutine deposit_energy_ddmc

! subroutine to deposit momentum to faces in ddmc
subroutine deposit_face_momentum_ddmc(blockID, cellID, dt, dtNew,&
                                      fleck, k_a, k_s,&
                                      particle,&
                                      is_remote, leak_axis, leak_dir)
  use Simulation_data, only : clight
  use Grid_interface, only : Grid_getBlkPtr, Grid_releaseBlkPtr
  use Particles_data, only : pt_is_corrdl

  implicit none
#include "Flash.h"
#include "constants.h"

  ! Input/Output
  integer, intent(in) :: blockID
  integer, dimension(MDIM), intent(in) :: cellID ! current cellID of MCP
  real, intent(in) :: dt, dtNew, fleck, k_a, k_s
  real, dimension(NPART_PROPS), intent(in) :: particle
  logical, intent(in) :: is_remote
  integer, intent(in) :: leak_axis, leak_dir

  ! aux variables
  real :: mcp_energy, mcp_eps, mcp_nump, mcp_avg_energy
  real, pointer :: faceVec(:,:,:,:)
#ifdef FLASH_GREY_DDMC
  integer, dimension(MDIM), parameter :: facegrid = (/ DDFX_FACE_VAR,&
                                                       DDFY_FACE_VAR,&
                                                       DDFZ_FACE_VAR /)
#endif
  integer, dimension(MDIM), parameter :: dir = (/ FACEX, FACEY, FACEZ /)
  integer, dimension(MDIM) :: faceID
  integer :: face_shift
  real :: k_total, dl, dtau, dl_corr
  real :: kEovercdt

  face_shift = 0

  ! find out face location
  faceID = cellID
  ! local right leakage
  if ((.not. is_remote) .and. (leak_dir > 0)) then
    face_shift = 1
  else if ((is_remote) .and. (leak_dir < 0)) then
    ! remote left leakage
    face_shift = 1
  end if
  faceID(leak_axis) = faceID(leak_axis) + face_shift

  k_total = k_a + k_s

  dl = clight * dt !* dshift ! convert dt to comoving frame
  dtau = k_total * dl

  ! Correction for absorption during the integrated path 
  dl_corr = 1.0d0 ! remote deposition should have dtau = 0
  if (pt_is_corrdl .and. (dtau > 0.0)) then
    dl_corr = (1.0d0 - exp(-dtau)) / dtau
  end if

  ! MCP energy
  mcp_eps  = particle(ENER_PART_PROP)
  mcp_nump = particle(NUMP_PART_PROP)
  mcp_energy = mcp_eps * mcp_nump
  mcp_avg_energy = mcp_energy * dl_corr

  ! accumulate k_tot*E_mcp/c/dt. Divide by face area and rho later
  kEovercdt = sign(k_total * mcp_avg_energy / clight / dtNew, real(leak_dir))

  ! debugging
!  print *, "depoDDMC", is_remote, leak_axis, leak_dir
!  print *, "cellID", cellID
!  print *, "faceID", faceID
!  print *, "k_tot", k_total
!  print *, "dl_corr", dl_corr

  ! retrieve and add to face variable
  call Grid_getBlkPtr(blockID, faceVec, dir(leak_axis))

#ifdef FLASH_GREY_DDMC
  faceVec(facegrid(leak_axis), faceID(IAXIS), faceID(JAXIS), faceID(KAXIS)) =&
    faceVec(facegrid(leak_axis), faceID(IAXIS), faceID(JAXIS), faceID(KAXIS)) + &
    kEovercdt
#endif

  call Grid_releaseBlkPtr(blockID, faceVec)
  ! end deposition of face grid variable

end subroutine deposit_face_momentum_ddmc

! subroutine to divide the DDF[X,Y,Z] face variables by the local face area
subroutine divide_by_face_area(blockID)
  use Grid_interface, only : Grid_getBlkBoundBox, Grid_getDeltas,&
                             Grid_getBlkIndexLimits,&
                             Grid_getBlkPtr, Grid_releaseBlkPtr
  use Grid_data, only : gr_geometry
  use spherical, only : get_spherical_cell_face_area
  implicit none
#include "Flash.h"
#include "constants.h"

  ! Input/Output
  integer, intent(in) :: blockID

  ! aux variables
  integer, dimension(2, MDIM) :: blkLimits, blkLimitsGC
  real, dimension(LOW:HIGH, MDIM) :: bndBox
  real, dimension(MDIM) :: deltaCell
#ifdef FLASH_GREY_DDMC
  integer, dimension(MDIM), parameter :: facegrid = (/ DDFX_FACE_VAR,&
                                                       DDFY_FACE_VAR,&
                                                       DDFZ_FACE_VAR /)
#endif
  integer, dimension(MDIM), parameter :: dir = (/ FACEX, FACEY, FACEZ /)
  real, pointer :: faceVec(:,:,:,:)
  integer :: ff, i, j, k, now_face
  real, dimension(MDIM) :: cart_face_areas
  integer, dimension(MDIM) :: faceID
  real :: face_area


  ! gather block info
  call Grid_getBlkBoundBox(blockID, bndBox)
  call Grid_getDeltas(blockID, deltaCell)
  cart_face_areas(IAXIS) = deltaCell(JAXIS)*deltaCell(KAXIS)
  cart_face_areas(JAXIS) = deltaCell(IAXIS)*deltaCell(KAXIS)
  cart_face_areas(KAXIS) = deltaCell(IAXIS)*deltaCell(JAXIS)

  ! indexing should be the same in FACEX, FACEY, and FACEZ

  do ff = IAXIS, KAXIS
    now_face = dir(ff)

    ! get face data structure limits
    call Grid_getBlkIndexLimits(blockID, blkLimits, blkLimitsGC, now_face)

    ! get face grid pointer
    call Grid_getBlkPtr(blockID, faceVec, dir(ff))

    ! go over interior cells' faces
    do k = blkLimitsGC(LOW,KAXIS), blkLimitsGC(HIGH,KAXIS)
      do j = blkLimitsGC(LOW, JAXIS), blkLimitsGC(HIGH, JAXIS)
        do i = blkLimitsGC(LOW,IAXIS), blkLimitsGC(HIGH, IAXIS)

          ! get current face area along direction ff
          if (gr_geometry == CARTESIAN) then
            face_area = cart_face_areas(ff)
          else if (gr_geometry == SPHERICAL) then
            faceID = (/ i, j, k /)
            call get_spherical_cell_face_area(bndBox, deltaCell,&
                                              ff, faceID, face_area)
          end if

          ! divide by face area
#ifdef FLASH_GREY_DDMC
          faceVec(facegrid(ff), i, j, k) = faceVec(facegrid(ff), i, j, k)&
                                                  / face_area

#endif

        end do
      end do
    end do
  
    ! return face grid pointer
    call Grid_releaseBlkPtr(blockID, faceVec)

  end do

end subroutine divide_by_face_area

! subroutine to store out leakage info for remote momentum depositio
subroutine store_ddmc_remote_momentum_deposition(leak_axis, leak_dir,&
                                                 is_trap_negh, particle)

  implicit none
#include "Flash.h"
#include "constants.h"

  ! Input/Output
  integer, intent(in) :: leak_axis, leak_dir
  logical, intent(in) :: is_trap_negh
  real, dimension(NPART_PROPS), intent(inout) :: particle

  ! aux variable
  integer :: particle_fate

  ! no need to continue if neighbor is not DDMC-active
  if (.not. is_trap_negh) return

  particle_fate = sign(leak_axis, leak_dir)
  particle(FATE_PART_PROP) = real(particle_fate)

end subroutine store_ddmc_remote_momentum_deposition

! subroutine to clean up remote deposition particle attribute
subroutine reset_particle_remote_depo(particle)
  implicit none
#include "Flash.h"
#include "constants.h"

  ! Input/Output
  real, dimension(NPART_PROPS), intent(inout) :: particle

  particle(FATE_PART_PROP) = 0.0

end subroutine reset_particle_remote_depo


! subroutine to average the face momentum deposition and populate them
! to DDMX, DDMY, DDMZ
subroutine avg_face_flux_and_populate_arad(blockID)
  use Particles_data, only : face_negh
  use Grid_interface, only : Grid_getBlkBoundBox, Grid_getDeltas,&
                             Grid_getBlkIndexLimits,&
                             Grid_getBlkPtr, Grid_releaseBlkPtr
  implicit none
#include "Flash.h"
#include "constants.h"

  ! Input/Output
  integer, intent(in) :: blockID

  ! aux variables
  integer, dimension(2, MDIM) :: blkLimits, blkLimitsGC
  real, dimension(LOW:HIGH, MDIM) :: bndBox
  real, dimension(MDIM) :: deltaCell
  integer, dimension(MDIM), parameter :: dir = (/ FACEX, FACEY, FACEZ /)
#ifdef FLASH_GREY_DDMC
  integer, dimension(MDIM), parameter :: facegrid = (/ DDFX_FACE_VAR,&
                                                       DDFY_FACE_VAR,&
                                                       DDFZ_FACE_VAR /)
  integer, dimension(MDIM), parameter :: ddm_var = (/ DDMX_VAR, DDMY_VAR, DDMZ_VAR /)
#endif
  real, pointer :: solnVec(:,:,:,:)
  real, pointer :: faceVec(:,:,:,:)
  integer :: i, j, k, ff
  real :: rho, avg_face_val
  integer, dimension(MDIM) :: adj_face

  call Grid_getBlkPtr(blockID, solnVec, CENTER)

  do ff = IAXIS, KAXIS
    ! prepare for indexing offset for the next face
    adj_face = face_negh(:, ff)

    ! get face grid pointer
    call Grid_getBlkPtr(blockID, faceVec, dir(ff))

    ! Get cell-center grid value's index limits
    call Grid_getBlkIndexLimits(blockID, blkLimits, blkLimitsGC)

    ! average L and R face values and populate DDM*
    do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
      do j = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
        do i = blkLimits(LOW,IAXIS), blkLimits(HIGH, IAXIS)

          rho = solnVec(DENS_VAR, i, j, k)
#ifdef FLASH_GREY_DDMC
          ! average L and R face values
          avg_face_val = 0.5*(faceVec(facegrid(ff),i,j,k) + &
                              faceVec(facegrid(ff),&
                                         i+adj_face(IAXIS),&
                                         j+adj_face(JAXIS),&
                                         k+adj_face(KAXIS)))

          ! divide it by cell-centered density
          ! At this point, the unit is cm/s/s
          avg_face_val = avg_face_val / rho

          ! populate DDM*_VAR with the average face value
          solnVec(ddm_var(ff), i, j, k) = avg_face_val
#endif
        end do
      end do
    end do

    ! return face grid pointer
    call Grid_releaseBlkPtr(blockID, faceVec)

  end do

  call Grid_releaseBlkPtr(blockID, solnVec)

end subroutine avg_face_flux_and_populate_arad

! subroutine to compute mu for conversion probability calc
! It depends on grid geometry.
subroutine compute_mcp_mu_toward_cell_bnd(particle, axis, dir, mu)
  use Grid_data, only : gr_geometry
  use Simulation_data, only : clight

  implicit none
#include "Flash.h"
#include "constants.h"

  ! Input/Output
  real, dimension(NPART_PROPS), intent(in) :: particle
  integer, intent(in) :: axis, dir
  real, intent(out) :: mu

  ! aux_variables
  real, dimension(MDIM) :: partvel
  real, dimension(MDIM, MDIM) :: vel_hat
  real :: r, theta, phi
  real :: sin_theta, cos_theta, sin_phi, cos_phi

  partvel = particle(VELX_PART_PROP:VELZ_PART_PROP)

  select case (gr_geometry)

    case (CARTESIAN)
      ! simple case: just v_i / c
      mu = abs(partvel(axis)) / clight

    case (SPHERICAL)
      r     = particle(POSX_PART_PROP)
      theta = particle(POSY_PART_PROP)
      phi   = particle(POSZ_PART_PROP)

      sin_theta = sin(theta)
      cos_theta = cos(theta)
      sin_phi   = sin(phi)
      cos_phi   = cos(phi)

      ! r_hat
      vel_hat(IAXIS, :) = (/ sin_theta*cos_phi,&
                             sin_theta*sin_phi,&
                             cos_theta /)
      ! theta_hat
      vel_hat(JAXIS, :) = (/ cos_theta*cos_phi,&
                             cos_theta*sin_phi,&
                             -sin_theta /)
      ! phi_hat
      vel_hat(KAXIS, :) = (/ -sin_phi,&
                              cos_phi,&
                              0.0d0 /)

      mu = abs(dot_product(partvel, vel_hat(axis, :)) / clight)

  end select

end subroutine compute_mcp_mu_toward_cell_bnd


! subroutine to compute conversion probability for MC2DDMC transition
subroutine compute_mc2ddmc_convert_prob(blockID, bndBox, deltaCell, solnVec,&
                                        cellID_negh, particle, axis, dir,&
                                        p_convert, is_converted)
  use Particles_data, only : pt_ddmc_lambda
  use opacity, only : calc_abs_opac, calc_sca_opac
  use random, only : rand

  implicit none
#include "Flash.h"
#include "constants.h"

  ! Input/Output
  integer, intent(in) :: blockID
  real, dimension(LOW:HIGH, MDIM), intent(in) :: bndBox
  real, dimension(MDIM), intent(in) :: deltaCell
  real, pointer :: solnVec(:,:,:,:)
  integer, dimension(MDIM), intent(in) :: cellID_negh
  real, dimension(NPART_PROPS), intent(in) :: particle
  integer, intent(in) :: axis, dir
  real, intent(out) :: p_convert
  logical, intent(out) :: is_converted

  ! aux variables
  real :: mu, ka_negh, ks_negh, ktot_negh
  real, dimension(MDIM) :: cell_widths 
  real :: tau_negh
  real :: xi
  integer :: negh_trap_index
  real :: negh_dx_corr

  ! use dummy 1.0 for energy and dshift
  ! gather neighbor opacities
  call calc_abs_opac(cellID_negh, solnVec, 1.0, 1.0, ka_negh)
  call calc_sca_opac(cellID_negh, solnVec, 1.0, 1.0, ks_negh)
  ktot_negh = ka_negh + ks_negh

  ! get cell widths
  call calc_cell_widths(cellID_negh, bndBox, deltaCell, cell_widths)

  call compute_mcp_mu_toward_cell_bnd(particle, axis, dir, mu)

  call get_negh_trap_index(blockID, bndBox, cellID_negh, negh_trap_index)

  negh_dx_corr = 1.0
  ! if index = +/-1, correct the dx_negh to reflect
  if (negh_trap_index == 1) then
    negh_dx_corr = 0.5
  else if (negh_trap_index == -1) then
    negh_dx_corr = 2.0
  end if

  tau_negh = ktot_negh * cell_widths(axis) * negh_dx_corr ! dir not used

  p_convert = 4.0/(3.0*tau_negh + 6.0*pt_ddmc_lambda)
  p_convert = p_convert*(1.0 + 1.5*mu)

  ! sampling imc2ddmc conversion
  xi = rand()
  is_converted = (xi < p_convert)

end subroutine compute_mc2ddmc_convert_prob

end module ddmc
