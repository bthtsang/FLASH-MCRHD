module rhd

  implicit none
  contains

!!! Subroutine to add a value to a cell's solnVec with specified field
subroutine cellAddVar(dataPtr, ijkindex, var, value)

  implicit none
#include "constants.h"
  real, pointer :: dataPtr(:,:,:,:)
  integer, dimension(MDIM), intent(in) :: ijkindex
  integer, intent(in) :: var
  real, intent(in) :: value

  dataPtr(var, ijkindex(1), ijkindex(2), ijkindex(3)) = &
         dataPtr(var, ijkindex(1), ijkindex(2), ijkindex(3)) + &
         value

end subroutine cellAddVar

!!! Subroutine to set gas temperature according to ionization fraction
subroutine set_temp_with_ion_fraction(now_h1, now_T)
  use Driver_interface, only : Driver_abortFlash
  use Grid_data, only : gr_geometry, gr_eosMode
  use Particles_data, only : pt_T_n, pt_T_i
  implicit none

  real, intent(in) :: now_h1
  real, intent(out) :: now_T

  ! Double check you are using dens_temp for EOS
  ! Turn off this error for spherical due to PPM error
  !if ((gr_geometry /= SPHERICAL) .and. (gr_eosMode /= MODE_DENS_TEMP)) then
  if ((gr_eosMode /= MODE_DENS_TEMP)) then
    call Driver_abortFlash("set_temp_with_ion_fraction: wrong EOS mode!")
  end if

  ! Perform linear interpolation
  ! Function in Henney et al. 2005
  now_T = pt_T_i - (pt_T_i - pt_T_n)*now_h1/(2.0 - now_h1)

end subroutine set_temp_with_ion_fraction


!!! Subroutine to explicitly apply the radiation source terms
!!! to the grid variables
subroutine apply_rad_source_terms(dt)
  use Particles_data, only : pt_is_thermally_coupled,&
                             pt_is_dynamically_coupled,&
                             pt_temp_floor, pt_is_apply_recombination,&
                             pt_is_pi_heating, pt_marshak_eos
  use Grid_interface, only : Grid_getBlkIndexLimits, Grid_getListOfBlocks,&
                             Grid_getBlkPtr, Grid_releaseBlkPtr
  use Grid_data, only : gr_geometry, gr_eosMode
  use Eos_interface, only : Eos_wrapped
  use Simulation_data, only : a_rad
  use Multispecies_interface, only : Multispecies_getProperty
  use Driver_interface, only : Driver_abortFlash
  use Timers_interface, only : Timers_start, Timers_stop

  implicit none
#include "Flash.h"
#include "constants.h"
#include "Multispecies.h"
  real, intent(in) :: dt

  ! aux variables
  integer, dimension(MAXBLOCKS) :: blkList
  integer :: numofblks, blockID
  integer :: b, i, j, k
  integer, dimension(2, MDIM) :: blkLimits, blkLimitsGC
  real, pointer :: solnVec(:,:,:,:)
  real :: Ae
  real :: ekin_old, ekin_new
  real :: dvxdt, dvydt, dvzdt
  integer, dimension(2, MDIM) :: EoscellID
  real :: temp, rho, eint, ekin, ener
  real :: x_n, T_interp
  ! Debug
  real :: oldtemp, oldeint, netheating, oldele

  call Timers_start("MCP RHD coupling")

  ! Get the list of leaf blocks
  call Grid_getListOfBlocks(LEAF, blkList, numofblks)

  ! Loop through leaf blocks on the current rank
  do b = 1, numofblks
    blockID = blkList(b)

    ! Get cell index limits
    call Grid_getBlkIndexLimits(blockID, blkLimits, blkLimitsGC)
    call Grid_getBlkPtr(blockID, solnVec, CENTER)

    do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
      do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
        do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)

          ! For marshak EOS, construct energies differently
          if (pt_marshak_eos) then
            rho = solnVec(DENS_VAR,i,j,k)
            temp = solnVec(TEMP_VAR,i,j,k)

            eint = a_rad*temp**4 / rho

            ! override the variable EINT_VAR here [erg/g]
            solnVec(EINT_VAR,i,j,k) = eint

            ekin = 0.5d0*dot_product(solnVec(VELX_VAR:VELZ_VAR,i,j,k),&
                                      solnVec(VELX_VAR:VELZ_VAR,i,j,k))
            ener = eint + ekin
            ! override the variable ENER_VAR here [erg/g]
            solnVec(ENER_VAR,i,j,k) = ener
          end if

          ! Apply radiation energy source
          if (pt_is_thermally_coupled) then
            oldeint = solnVec(EINT_VAR,i,j,k)
            netheating = solnVec(ABSE_VAR,i,j,k) + solnVec(EMIE_VAR,i,j,k)

            ! Thermal heating and cooling
            solnVec(EINT_VAR,i,j,k) = solnVec(EINT_VAR,i,j,k)&
                                       + solnVec(ABSE_VAR,i,j,k)&
                                       + solnVec(EMIE_VAR,i,j,k)
            solnVec(ENER_VAR,i,j,k) = solnVec(ENER_VAR,i,j,k)&
                                       + solnVec(ABSE_VAR,i,j,k)&
                                       + solnVec(EMIE_VAR,i,j,k)

            !if (netheating > oldeint) then
            !  print *, "huge_heating", oldeint
            !  print *, "net heating", netheating, solnVec(ABSE_VAR,i,j,k), solnVec(EMIE_VAR,i,j,k)
            !  print *, "now eint", solnVec(EINT_VAR,i,j,k)
            !end if

            ! Photoionization heating and cooling
#ifdef HEAT_VAR
#ifdef COOL_VAR

            ! debugging
            oldtemp = solnVec(TEMP_VAR,i,j,k)
            oldeint = solnVec(EINT_VAR,i,j,k)
            netheating = solnVec(HEAT_VAR,i,j,k) + solnVec(COOL_VAR,i,j,k)
            oldele  = solnVec(ELE_SPEC, i, j,k)


            if (pt_is_apply_recombination) then
              solnVec(EINT_VAR,i,j,k) = solnVec(EINT_VAR,i,j,k)&
                                         + solnVec(HEAT_VAR,i,j,k)&
                                         + solnVec(COOL_VAR,i,j,k)
              solnVec(ENER_VAR,i,j,k) = solnVec(ENER_VAR,i,j,k)&
                                         + solnVec(HEAT_VAR,i,j,k)&
                                         + solnVec(COOL_VAR,i,j,k)
            else 
              solnVec(EINT_VAR,i,j,k) = solnVec(EINT_VAR,i,j,k)&
                                         + solnVec(HEAT_VAR,i,j,k)
              solnVec(ENER_VAR,i,j,k) = solnVec(ENER_VAR,i,j,k)&
                                         + solnVec(HEAT_VAR,i,j,k)
            end if

            EoscellID(:, IAXIS) = i
            EoscellID(:, JAXIS) = j
            EoscellID(:, KAXIS) = k
            call Eos_wrapped(gr_eosMode, EoscellID, blockID)

            temp = solnVec(TEMP_VAR,i,j,k)
            eint = solnVec(EINT_VAR,i,j,k)
            if (temp > 1.0e6) then
              print *, "Overheating..", oldtemp, temp
              print *, "Eint..", oldeint, eint
              print *, "heat/cool:", netheating, solnVec(HEAT_VAR,i,j,k), solnVec(COOL_VAR,i,j,k)
              print *, "abse/emie:", solnVec(ABSE_VAR,i,j,k), solnVec(EMIE_VAR,i,j,k)
              print *, "ion state:", solnVec(IONR_VAR,i,j,k), solnVec(RECR_VAR,i,j,k)
              print *, "species:", solnVec(H1_SPEC,i,j,k), solnVec(HP_SPEC,i,j,k), solnVec(ELE_SPEC,i,j,k)
              print *, "vel", solnVec(VELX_VAR:VELZ_VAR,i,j,k)
              print *, "dens", solnVec(DENS_VAR,i,j,k)

              call Driver_abortFlash("rhd: overheating in cell, please follow up.")
            end if
            if (oldele > 1.0) then
              print*, "ELEwrong", oldele
              print*, "N12", solnVec(H1_SPEC,i,j,k), solnVec(HP_SPEC,i,j,k) 
              print*, "Temps", oldtemp, temp
              print*, "eints", oldeint, eint
              print*, "heats", netheating, netheating, solnVec(HEAT_VAR,i,j,k), solnVec(COOL_VAR,i,j,k)
              print*, "ions", solnVec(IONR_VAR,i,j,k), solnVec(RECR_VAR,i,j,k)
            end if
#endif 
#endif
          end if

          ! Apply radiation momentum source
          if (pt_is_dynamically_coupled) then
            ekin_old = 0.5d0*dot_product(solnVec(VELX_VAR:VELZ_VAR,i,j,k),&
                                         solnVec(VELX_VAR:VELZ_VAR,i,j,k))

            dvxdt = (solnVec(ABMX_VAR,i,j,k) + solnVec(SCMX_VAR,i,j,k))
            solnVec(VELX_VAR,i,j,k) = solnVec(VELX_VAR,i,j,k) + dvxdt * dt

            if (NDIM > 1) then
              dvydt = (solnVec(ABMY_VAR,i,j,k) + solnVec(SCMY_VAR,i,j,k))
              solnVec(VELY_VAR,i,j,k) = solnVec(VELY_VAR,i,j,k) + dvydt * dt

              if (NDIM == 3) then
                dvzdt = (solnVec(ABMZ_VAR,i,j,k) + solnVec(SCMZ_VAR,i,j,k))
                solnVec(VELZ_VAR,i,j,k) = solnVec(VELZ_VAR,i,j,k) + dvzdt * dt
              end if
            end if

            ! KE after radiation source terms
            ekin_new = 0.5d0*dot_product(solnVec(VELX_VAR:VELZ_VAR,i,j,k),&
                                         solnVec(VELX_VAR:VELZ_VAR,i,j,k))

            ! Energy source term due to momentum sources
            solnVec(ENER_VAR,i,j,k) = solnVec(ENER_VAR,i,j,k)&
                                       - ekin_old + ekin_new
          end if

          ! Simple case of not treating photoionization heating
          ! self-consistently
#ifdef FLASH_MCRHD_IONIZATION
          if (.not. pt_is_pi_heating) then
            x_n = solnVec(H1_SPEC, i, j, k)

            call set_temp_with_ion_fraction(x_n, T_interp)
            solnVec(TEMP_VAR, i, j, k) = T_interp

            ! Override with the appropriate EINT if in spherical coord.
            ! The DENS_IE mode then will predict the correct TEMP
            if (gr_geometry == SPHERICAL) then
              call Eos_wrapped(MODE_DENS_TEMP, blkLimits, blockID)
            end if
          end if
#endif

        end do
      end do
    end do

    ! Update cell temperature after applying the source terms
    call Eos_wrapped(gr_eosMode, blkLimits, blockID)

    ! Impose temperature floor
    if (pt_temp_floor > 0.0) then
      do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
          do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)    
            if (pt_temp_floor > 0.0) then
              temp = solnVec(TEMP_VAR,i,j,k)

#ifdef FLASH_MCRHD_IONIZATION
              if (temp > 1.0e6) then
                print *, "Overheating..", temp
                print *, "heat/cool:", solnVec(HEAT_VAR,i,j,k), solnVec(COOL_VAR,i,j,k)
                print *, "ion state:", solnVec(IONR_VAR,i,j,k), solnVec(RECR_VAR,i,j,k)
                print *, "species:", solnVec(H1_SPEC,i,j,k), solnVec(HP_SPEC,i,j,k), solnVec(ELE_SPEC,i,j,k)
              end if 
#endif

              if (temp < pt_temp_floor) then
                EoscellID(:, IAXIS) = i
                EoscellID(:, JAXIS) = j
                EoscellID(:, KAXIS) = k
                solnVec(TEMP_VAR,i,j,k) = pt_temp_floor

                ! Override with the new temperature
                call Eos_wrapped(MODE_DENS_TEMP, EoscellID, blockID)

                if (pt_marshak_eos) then
                  call Driver_abortFlash("rhd: eint modified by temp. floor but&
                                          marshak EOS is on. Please follow up.")
                end if
              end if
            end if
          end do
        end do
      end do
    end if

    ! Override temperature update with Marshak test EOS if needed
    if (pt_marshak_eos) then
      do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
          do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)

            eint = solnVec(EINT_VAR,i,j,k)
            rho  = solnVec(DENS_VAR,i,j,k)

            temp = (rho*eint/a_rad)**(0.25)

            ! Override temperature computed from Eos_wrapped
            solnVec(TEMP_VAR,i,j,k) = temp
            ! EINT and ENER should not be changed
          end do
        end do
      end do
    end if

    call Grid_releaseBlkPtr(blockID, solnVec)
  end do

  call Timers_stop("MCP RHD coupling")

end subroutine apply_rad_source_terms

end module rhd
