module scattering 

  implicit none
  contains

subroutine scatter_mcp(solnVec, cellID, dt, particle,&
                       is_elastic, is_isotropic)
  use Driver_interface, only : Driver_abortFlash
  use Particles_data, only : pt_is_veldp, pt_samp_Tgas
  use new_mcp, only : sample_iso_velocity, sample_energy
  use relativity, only : transform_lab_to_comoving,&
                         transform_comoving_to_lab
  implicit none

#include "constants.h"
#include "Flash.h"

  ! Input/output
  real, pointer :: solnVec(:,:,:,:)
  integer, dimension(MDIM), intent(in) :: cellID
  real, intent(in) :: dt
  real, dimension(NPART_PROPS), intent(inout) :: particle
  logical, intent(in) :: is_elastic, is_isotropic

  ! aux variables
  real :: dshift, temp 
  real, dimension(MDIM) :: old_vel, new_vel
  real :: mcp_eps_bs, mcp_w_now_bs, mcp_w_ini_bs, mcp_energy_bs
  real :: mcp_eps_as, mcp_w_now_as, mcp_w_ini_as, mcp_energy_as

  temp = solnVec(TEMP_VAR, cellID(IAXIS), cellID(JAXIS), cellID(KAXIS))

  if (pt_is_veldp) then
    call transform_lab_to_comoving(cellID, solnVec, dt, particle, dshift)
  end if

  ! Getting MCP information
  old_vel = particle(VELX_PART_PROP:VELZ_PART_PROP)
  mcp_eps_bs    = particle(ENER_PART_PROP)
  mcp_w_now_bs  = particle(NUMP_PART_PROP)
  mcp_w_ini_bs  = particle(NPIN_PART_PROP)
  mcp_energy_bs = mcp_eps_bs * mcp_w_now_bs

  ! Sample photon energy
  if (is_elastic) then
    mcp_eps_as = mcp_eps_bs
  else 
    ! Sample new photon energy thermally
    call sample_energy(solnVec, cellID, pt_samp_Tgas, temp, mcp_eps_as)
  end if
  ! Conserve photon number if photon energy shifted
  mcp_energy_as = mcp_energy_bs ! MCP energy should conserve
  mcp_w_now_as = mcp_energy_as / mcp_eps_as
  mcp_w_ini_as = (mcp_w_ini_bs / mcp_w_now_bs) * mcp_w_now_as

  ! Sample velocity, remember this is in the comoving frame
  if (is_isotropic) then
    call sample_iso_velocity(new_vel)
  else 
    call Driver_abortFlash("scatter_mcp: anisotropic scattering&
                            has not been implemented.")
  end if

  ! Update MCP attributes  
  particle(ENER_PART_PROP) = mcp_eps_as
  particle(VELX_PART_PROP:VELZ_PART_PROP) = new_vel

  particle(NUMP_PART_PROP) = mcp_w_now_as
  particle(NUM0_PART_PROP) = mcp_w_now_as
  particle(NPIN_PART_PROP) = mcp_w_ini_as

  if (pt_is_veldp) then
    call transform_comoving_to_lab(cellID, solnVec, dt, particle, dshift)
  end if

end subroutine scatter_mcp

end module scattering
