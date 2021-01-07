module opacity 

  implicit none
  contains

subroutine calc_abs_opac(cellID, solnVec, energy, ka)
  use Particles_data, only : pt_is_grey, pt_grey_abs_opac, pt_dens_threshold,&
                             pt_is_kt_opac
  use Driver_interface,  ONLY : Driver_abortFlash
  implicit none

#include "constants.h"
#include "Flash.h"

  ! Input/Output
  integer, dimension(MDIM), intent(in) :: cellID
  real, pointer :: solnVec(:,:,:,:)
  real, intent(in) :: energy
  real, intent(out) :: ka

  ! aux variables
  real :: rho, temp, temp_opac

  rho = solnVec(DENS_VAR, cellID(IAXIS), cellID(JAXIS), cellID(KAXIS))
  temp = solnVec(TEMP_VAR, cellID(IAXIS), cellID(JAXIS), cellID(KAXIS))

  if (pt_is_grey) then
    ka = pt_grey_abs_opac

    if (pt_is_kt_opac) then
      temp_opac = temp
      if (temp >= 150.0) temp_opac = 150.0
      ka = 0.10*(temp_opac/10.0)**2
    end if 

    ka = ka * rho

    if (rho .lt. pt_dens_threshold) ka = 0.0d0
  else
    call Driver_abortFlash("calc_abs_opac:&
                           non-grey emission not implemented yet!")
  end if

end subroutine calc_abs_opac


subroutine calc_sca_opac(cellID, solnVec, energy, ks)
  use Particles_data, only : pt_is_grey, pt_grey_sca_opac, pt_dens_threshold
  use Driver_interface,  ONLY : Driver_abortFlash
  implicit none

#include "constants.h"
#include "Flash.h"

  ! Input/Output
  integer, dimension(MDIM), intent(in) :: cellID
  real, pointer :: solnVec(:,:,:,:)
  real, intent(in) :: energy
  real, intent(out) :: ks

  ! aux variables
  real :: rho, temp

  rho = solnVec(DENS_VAR, cellID(IAXIS), cellID(JAXIS), cellID(KAXIS))
  temp = solnVec(TEMP_VAR, cellID(IAXIS), cellID(JAXIS), cellID(KAXIS))

  if (pt_is_grey) then
    ks = pt_grey_sca_opac
    ks = ks * rho

    if (rho .lt. pt_dens_threshold) ks = 0.0d0
  else
    call Driver_abortFlash("calc_sca_opac:&
                            non-grey emission not implemented yet!")
  end if

end subroutine calc_sca_opac


subroutine get_fleck(cellID, solnVec, fleck)
  use Particles_data, only : pt_is_eff_scattering
  implicit none
#include "Flash.h"

  integer, dimension(MDIM), intent(in) :: cellID
  real, pointer :: solnVec(:,:,:,:)
  real, intent(out) :: fleck

  ! The grid should be populated with fleck during 
  ! the thermal emission call
  fleck = solnVec(FLEC_VAR, cellID(IAXIS), cellID(JAXIS), cellID(KAXIS))

  if (.not. pt_is_eff_scattering) fleck = 1.0d0

end subroutine get_fleck

end module opacity
