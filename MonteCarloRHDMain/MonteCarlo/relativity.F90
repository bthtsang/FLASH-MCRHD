module relativity

  implicit none
  contains

! Higher level subroutine to apply Lorentz transformation
subroutine transform_lab_to_comoving(cellID, solnVec, particle, dshift)

  implicit none
#include "constants.h"
#include "Flash.h"

  integer, dimension(MDIM), intent(in) :: cellID
  real, pointer :: solnVec(:,:,:,:)
  real, dimension(NPART_PROPS), intent(inout) :: particle
  real, intent(out) :: dshift

  call lorentz_transformation(cellID, solnVec, .false., particle, dshift)

end subroutine transform_lab_to_comoving


! Higher level subroutine to apply Lorentz transformation
subroutine transform_comoving_to_lab(cellID, solnVec, particle, dshift)

  implicit none
#include "constants.h"
#include "Flash.h"
  integer, dimension(MDIM), intent(in) :: cellID
  real, pointer :: solnVec(:,:,:,:)
  real, dimension(NPART_PROPS), intent(inout) :: particle
  real, intent(out) :: dshift

  call lorentz_transformation(cellID, solnVec, .true., particle, dshift)

end subroutine transform_comoving_to_lab


! Lorentz transformation for a MCP
! It transforms the MCP's
! propagation direction (VEL*_PART_PROP) and
! frequency (ENER_PART_PROP).
! The weight (number of photons) should conserve under
! the transformation (NUMP_PART_PROP).
subroutine lorentz_transformation(cellID, solnVec, tolab, particle, dshift)
  use Driver_interface, only : Driver_abortFlash
  use Simulation_data, only : clight
  use spherical, only : get_cartesian_position, get_spherical_position
  use Grid_data, only: gr_geometry
  implicit none
#include "constants.h"
#include "Flash.h"

  ! Input/output
  integer, dimension(MDIM), intent(in) :: cellID
  real, pointer :: solnVec(:,:,:,:)
  logical, intent(in) :: tolab
  real, dimension(NPART_PROPS), intent(inout) :: particle
  real, intent(out) :: dshift

  ! aux variables
  real, dimension(MDIM) :: v_mcp, v_gas, v_mcp_new
  real, dimension(MDIM) :: x_mcp, x_mcp_cart, x_mcp_sph, x_mcp_new
  real :: beta2, vdd, gamm
  real :: v_norm
  real :: vdx, v2
  real :: t_mcp, t_mcp_new

  ! Extract MCP and gas (cell) velocity
  v_mcp = particle(VELX_PART_PROP:VELZ_PART_PROP)
  v_gas = solnVec(cellID(IAXIS), cellID(JAXIS), cellID(KAXIS), VELX_VAR:VELZ_VAR)

  ! Extract MCP's position and time
  t_mcp = particle(TREM_PART_PROP)
  x_mcp = particle(POSX_PART_PROP:POSZ_PART_PROP)
  if (gr_geometry == SPHERICAL) then
    call get_cartesian_position(x_mcp, x_mcp_cart)
    x_mcp = x_mcp_cart
  end if

  if (NDIM < 3) then
    call Driver_abortFlash("relativity::lorentz_transformation")
  end if

  if (tolab) then
    v_gas = -v_gas
  end if

  ! Computing relativistic quantities
  v2 = dot_product(v_gas, v_gas) 
  beta2 = v2/clight/clight
  vdd  = dot_product(v_gas, v_mcp)
  gamm = 1.0/sqrt(1.0 - beta2)

  ! dshift is the epsilon_0/epsilon factor in the LT,
  ! which is also useful for converting kappa's between frames
  dshift = gamm*(1.0 - vdd/clight)

  ! Compute new position
  vdx = dot_product(v_gas, x_mcp)
  x_mcp_new(1) = x_mcp(1) + ((gamm - 1.0)*vdx/v2 - gamm*t_mcp)*v_gas(1)
  x_mcp_new(2) = x_mcp(2) + ((gamm - 1.0)*vdx/v2 - gamm*t_mcp)*v_gas(2)
  x_mcp_new(3) = x_mcp(3) + ((gamm - 1.0)*vdx/v2 - gamm*t_mcp)*v_gas(3)
  if (gr_geometry == SPHERICAL) then
    call get_spherical_position(x_mcp_new, x_mcp_sph)
    x_mcp_new = x_mcp_sph
  end if
  ! Compute new time
  t_mcp_new = gamm*(t_mcp - vdx/clight/clight)

  ! Update photon propagation direction to new frame
  v_mcp_new(1) = (v_mcp(1) - gamm*v_gas(1)/clight*(1.0 - gamm*vdd/clight/(gamm+1.0)))/dshift
  v_mcp_new(2) = (v_mcp(2) - gamm*v_gas(2)/clight*(1.0 - gamm*vdd/clight/(gamm+1.0)))/dshift
  v_mcp_new(3) = (v_mcp(3) - gamm*v_gas(3)/clight*(1.0 - gamm*vdd/clight/(gamm+1.0)))/dshift
  ! Re-normalize velocity vector just in case
  v_norm = sqrt(dot_product(v_mcp_new, v_mcp_new))
  v_mcp_new(1) = v_mcp_new(1) / v_norm * clight
  v_mcp_new(2) = v_mcp_new(2) / v_norm * clight
  v_mcp_new(3) = v_mcp_new(3) / v_norm * clight

  ! Update new position and time
  particle(POSX_PART_PROP:POSZ_PART_PROP) = x_mcp_new
  particle(TREM_PART_PROP) = t_mcp_new

  ! Update photon energy to new frame 
  particle(ENER_PART_PROP) = dshift*particle(ENER_PART_PROP)
  ! Update new propagation direction
  particle(VELX_PART_PROP:VELZ_PART_PROP) = v_mcp_new

end subroutine lorentz_transformation

! Subroutine to compute the epsilon/epsilon_0 value
! without modifying the particle's attributes. 
! This is mostly meant for opacity calculation
subroutine compute_dshift(cellID, solnVec, tolab, particle, dshift)
  use Simulation_data, only : clight

  implicit none
#include "Flash.h"

  ! Input/output
  integer, dimension(MDIM), intent(in) :: cellID
  real, pointer :: solnVec(:,:,:,:)
  logical, intent(in) :: tolab
  real, dimension(NPART_PROPS), intent(in) :: particle
  real, intent(out) :: dshift

  ! aux variables
  real, dimension(MDIM) :: v_mcp, v_gas
  real :: v2, beta2, vdd, gamm

  v_mcp = particle(VELX_PART_PROP:VELZ_PART_PROP)
  v_gas = solnVec(cellID(IAXIS), cellID(JAXIS), cellID(KAXIS), VELX_VAR:VELZ_VAR)

  if (tolab) then
    v_gas = -v_gas
  end if

  ! Computing relativistic quantities
  v2 = dot_product(v_gas, v_gas) 
  beta2 = v2/clight/clight
  vdd  = dot_product(v_gas, v_mcp)
  gamm = 1.0/sqrt(1.0 - beta2)
  dshift = gamm*(1.0 - vdd/clight)

end subroutine 

end module relativity
