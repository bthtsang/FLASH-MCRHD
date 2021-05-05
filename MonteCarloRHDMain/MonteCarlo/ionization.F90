module ionization 

  implicit none
  contains

! Photoionization absorption coefficient
subroutine calc_pi_opac(cellID, solnVec, energy, kpi, nH, nH1)
  use Particles_data, only : pt_is_photoionization, pt_dens_threshold,&
                             pt_nH1_threshold
  use Driver_interface, only : Driver_abortFlash
  use Simulation_data, only : mH, ev2erg, sigma_pi
  use Multispecies_interface, only : Multispecies_getProperty
  implicit none

#include "constants.h"
#include "Flash.h"
#include "Multispecies.h"

  ! Input/Output
  integer, dimension(MDIM), intent(in) :: cellID
  real, pointer :: solnVec(:,:,:,:)
  real, intent(in) :: energy
  real, intent(out) :: kpi, nH, nH1

  ! aux variables
  real :: rho, temp
  real :: h1, Ah1, energy_eV
  real :: hp, Ahp

  kpi = 0.0d0
  nH1 = 1.0d0
  if (.not. pt_is_photoionization) return 

#ifdef FLASH_MCRHD_IONIZATION
#ifdef FLASH_MULTISPECIES
  rho = solnVec(DENS_VAR, cellID(IAXIS), cellID(JAXIS), cellID(KAXIS))
  temp = solnVec(TEMP_VAR, cellID(IAXIS), cellID(JAXIS), cellID(KAXIS))

  !if ((NSPECIES /= 0) .and. (rho >= pt_dens_threshold)) then
  if (NSPECIES /= 0) then
    h1 = solnVec(H1_SPEC, cellID(IAXIS), cellID(JAXIS), cellID(KAXIS))
    hp = solnVec(HP_SPEC, cellID(IAXIS), cellID(JAXIS), cellID(KAXIS))
    call Multispecies_getProperty(H1_SPEC, A, Ah1)
    nH1 = rho*h1/(Ah1*mH)
    nH =  rho*(h1+hp)/(Ah1*mH)
    energy_eV = energy / ev2erg

    !if ((energy_eV >= 13.6) .and. (nH1 > pt_nH1_threshold)) then
    if ((energy_eV >= 13.6) .and. (h1 > pt_nH1_threshold)) then
      kpi = nH1*sigma_pi ! assume temperature-independent
    end if
  else ! NSPECIES=0, not having separate species for HI and HII
    call Driver_abortFlash("calc_pi_opac: multispecies is not defined &
                            but photoionization opacity is requested! Aborting.")
  end if 
#endif
#endif

end subroutine calc_pi_opac


! Collisional ionization coefficient
subroutine calc_k_coll(cellID, solnVec, k_coll)
  use Particles_data, only : pt_is_coll_ionization, pt_dens_threshold
  use Simulation_data, only : kB, ev2erg
  implicit none

#include "constants.h"
#include "Flash.h"

  ! Input/Output
  integer, dimension(MDIM), intent(in) :: cellID
  real, pointer :: solnVec(:,:,:,:)
  real, intent(out) :: k_coll

  ! aux variables
  real :: rho, temp

  k_coll = 0.0d0
  if (.not. pt_is_coll_ionization) return

  rho = solnVec(DENS_VAR, cellID(IAXIS), cellID(JAXIS), cellID(KAXIS))
  temp = solnVec(TEMP_VAR, cellID(IAXIS), cellID(JAXIS), cellID(KAXIS))

  if (rho >= pt_dens_threshold) then
    k_coll = 5.84d-11*sqrt(temp)*exp(-13.6*ev2erg/(kB*temp))
  end if

end subroutine calc_k_coll


! subroutine to compute recombination coefficient alpha_B
! for now, constant at 2.59d-13
! Collisional ionization coefficient
subroutine calc_recomb_coeff(cellID, solnVec, is_caseB, alpha_X)
  use Particles_data, only : pt_is_photoionization, pt_dens_threshold
  implicit none

#include "constants.h"
#include "Flash.h"

  ! Input/Output
  integer, dimension(MDIM), intent(in) :: cellID
  real, pointer :: solnVec(:,:,:,:)
  logical, intent(in) :: is_caseB
  real, intent(out) :: alpha_X

  ! aux variables
  real :: rho, temp

  alpha_X = 0.0d0
  if (.not. pt_is_photoionization) return

  rho = solnVec(DENS_VAR, cellID(IAXIS), cellID(JAXIS), cellID(KAXIS))
  temp = solnVec(TEMP_VAR, cellID(IAXIS), cellID(JAXIS), cellID(KAXIS))

  ! adopt a constant value for now
  if (is_caseB) then
    alpha_X = 2.59d-13!*(temp/1.0d4)**(-0.7)
  else
    alpha_X = 4.18d-13
  end if

end subroutine calc_recomb_coeff


! subroutine to compute the phototionization fleck factor
subroutine calc_ionization_fleck(nH, nH1, k_coll, alpha_B, dtnow,&
                                 gammap, fleckp)
  use Particles_data, only : pt_is_photoionization,&
                             pt_is_es_photoionization
  implicit none
#include "constants.h"
#include "Flash.h"

  ! Input/Output
  real, intent(in) :: nH, nH1, k_coll, alpha_B, dtnow
  real, intent(out) :: gammap, fleckp 

  ! aux variables
  real :: term1, term2
  
  gammap = 1.0d0
  fleckp = 1.0d0
  if ((.not. pt_is_photoionization) .or.&
      (.not. pt_is_es_photoionization)) return

  term1 = (nH - 2.0*nH1)*k_coll
  term2 = 2.0*(nH - nH1)*alpha_B

  gammap = 1.0 / (1.0 + (term1 + term2)*dtnow)
  fleckp = (1.0 + term1*dtnow) * gammap

  if (fleckp < 1.0e-10) then
    print *, "calc_ion_terms", term1, term2, nH, nH1
    print *, "coeffs", k_coll, alpha_B
    print *, "calc_ion_fleck", gammap, fleckp, dtnow
  end if

end subroutine calc_ionization_fleck

! Subroutine to populate recombination and the associated
! cooling rate. 
subroutine calc_recomb_emissivity(blkID, solnVec, dtNew,&
                                  now_num, now_pos, now_time, &
                                  now_energy, now_vel, now_weight)
  use Particles_data, only : pt_is_photoionization, pt_is_coll_ionization,&
                             pt_is_es_photoionization, pt_maxnewnum,&
                             pt_is_apply_recombination,&
                             pt_nH1_threshold, pt_is_caseB, pt_is_caseA_radeqm,&
                             pt_num_r1mcps_tstep, pt_is_veldp
  use Grid_interface, only : Grid_getBlkIndexLimits, Grid_getBlkBoundBox,&
                             Grid_getDeltas, Grid_getSingleCellVol
  use new_mcp, only : sample_blk_position, sample_iso_velocity,&
                      sample_cell_position,&
                      sample_therm_face_velocity,&
                      sample_cart_therm_face_velocity,&
                      sample_time, sample_energy
  use Multispecies_interface, only : Multispecies_getProperty
  use Simulation_data, only : mH, ev2erg, kB
  use relativity, only : transform_comoving_to_lab
  implicit none

#include "constants.h"
#include "Flash.h"

  ! Input/Output
  integer, intent(in) :: blkID
  real, pointer :: solnVec(:,:,:,:)
  real, intent(in) :: dtNew
  integer, intent(out) :: now_num
  real, dimension(MDIM,pt_maxnewnum), intent(out) :: now_pos, now_vel
  real, dimension(pt_maxnewnum), intent(out) :: now_time, now_energy, now_weight

  ! aux variables
  integer, dimension(2, MDIM) :: blkLimits, blkLimitsGC
  real, dimension(LOW:HIGH, MDIM) :: bndBox
  real, dimension(MDIM) :: deltaCell
  integer, dimension(MDIM) :: cellID
  integer :: i, j, k
  real :: rho, temp
  real :: h1, hp, ele, Ah1, Ae
  real :: nH, nH1, nHp, ne
  real :: k_coll, alpha_X, gamma_n, fleck_n
  real :: alpha_A, alpha_B
  real :: recomb_rate, cooling_rate
  integer :: ii
  real, dimension(NPART_PROPS) :: newparticle
  real :: dV, dE, dE_per_mcp, dtNow
  real :: weight_per_mcp, dshift
  real, dimension(MDIM) :: newxyz, newvel
  real :: newenergy

  real :: Lambda_rec
  real :: avg_recomb_energy

  ! Initialization
  now_num = 0
  now_pos = 0.0d0
  now_vel = 0.0d0
  now_time = 0.0d0
  now_energy = 0.0d0
  now_weight = 0.0d0

  ! exit if photoionization is off
  if (.not. pt_is_photoionization) return

  avg_recomb_energy = 2.4*ev2erg

#ifdef FLASH_MCRHD_IONIZATION
#ifdef H1_SPEC
  call Grid_getBlkIndexLimits(blkID, blkLimits, blkLimitsGC)
  call Grid_getBlkBoundBox(blkID, bndBox)
  call Grid_getDeltas(blkID, deltaCell)
  do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
    do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
      do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
        cellID = (/ i, j, k /)
        call Grid_getSingleCellVol(blkID, EXTERIOR, cellID, dV)

        ! Thermodynamical quantities
        rho  = solnVec(DENS_VAR, i, j, k)
        temp = solnVec(TEMP_VAR, i, j, k)

        ! ionization states
        h1   = solnVec(H1_SPEC, i, j, k)
        hp   = solnVec(HP_SPEC, i, j, k)
        ele  = solnVec(ELE_SPEC, i, j, k)
        call Multispecies_getProperty(H1_SPEC, A, Ah1)
        call Multispecies_getProperty(ELE_SPEC, A, Ae)
        nH  = rho*(h1+hp)/(Ah1*mH)
        nH1 = rho*h1/(Ah1*mH)
        nHp = nH - nH1
        ne  = rho*ele/(Ae*mH)

        if (temp < 10.0) then
          print *, "FailT", temp
          print *, "Failrho", rho
          print *, "FailV", solnVec(VELX_VAR:VELZ_VAR, i, j, k)
          print *, "Failedijk", i, j, k
          print *, "Failedns", h1, hp, ele
        end if

        call calc_recomb_coeff(cellID, solnVec, pt_is_caseB, alpha_X)        
        call calc_k_coll(cellID, solnVec, k_coll)

        call calc_ionization_fleck(nH, nH1, k_coll, alpha_X, dtNew,&
                                     gamma_n, fleck_n)
        if (pt_is_caseA_radeqm) then
          call calc_recomb_coeff(cellID, solnVec, .false., alpha_A)
          call calc_recomb_coeff(cellID, solnVec, .true.,  alpha_B)
          gamma_n = alpha_B / alpha_A
          fleck_n = alpha_B / alpha_A
        end if

        if (fleck_n < 1.0d-26) then
          print *, "tinyfleckp", nH, nH1, fleck_n
          print *, "cellID", i, j, k
          print *, "rho", rho
          print *, "temp", temp 
          print *, "h1,hp", h1, hp
          print *, "consts", Ah1, mH
        end if

        solnVec(GAMP_VAR, i, j, k) = gamma_n
        solnVec(FLEP_VAR, i, j, k) = fleck_n

        recomb_rate = 0.0
        cooling_rate = 0.0
        !if (nHp > pt_nH1_threshold) then
        if (hp > pt_nH1_threshold) then
          if (pt_is_es_photoionization) then
            !recomb_rate = (gamma_n - 1.0)*nH1& 
            !            + gamma_n*(nH*nH - nH1*nH1)*alpha_X*dtNew&
            !            - gamma_n*nH1*nH1*k_coll*dtNew
            recomb_rate = gamma_n*(nH - nH1)**2*alpha_X*dtNew
          else ! non implicit approach
            recomb_rate = alpha_X*nHp*nHp ! number/cm^3/s
            !recomb_rate = (nH*nH - nH1*nH1)*alpha_B&
            !            - nH1*nH1*k_coll
            ! Convert to change in mass fraction (dimensionless)
            recomb_rate = recomb_rate*dtNew ! number/cm^3
          end if
          recomb_rate = recomb_rate/nH ! dimensionless
        end if

        if (recomb_rate > hp) then
          recomb_rate = hp
        end if

        ! Setting the associated recombination cooling rate
        cooling_rate = avg_recomb_energy*recomb_rate/mH

        ! Sample recombination radiation field if not in case B approx.
        ! and if the cell has non-zero recombination rate
        ! Additionally, requires rad eqm for case A to be OFF.
        if ((.not. pt_is_caseB) .and. (recomb_rate > 0.0) .and.&
            (.not. pt_is_caseA_radeqm)) then
          dE = cooling_rate * rho * dV
          ! Separate out the ionizing fraction of the recombination
          dE = dE * (alpha_A - alpha_B) / alpha_A

          dE_per_mcp = dE / pt_num_r1mcps_tstep

          ! Loop to create new MCPs
          do ii = 1, pt_num_r1mcps_tstep

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

          now_num = now_num + pt_num_r1mcps_tstep

        end if

        !if (hp > pt_nH1_threshold) then
        !  Lambda_rec = 6.1e-10*kB*temp*temp**(-0.89)
        !  cooling_rate = Lambda_rec*nHp*nHp*dtnow/rho
        !end if

        ! Debug
      !  if (abs(cooling_rate) > solnVec(EINT_VAR, i, j, k)) then
      !    print *, "Overcooling:", solnVec(EINT_VAR, i, j, k), solnVec(TEMP_VAR,i, j, k)
      !    print *, "Info:", h1, hp, recomb_rate, cooling_rate
      !    print *, "nH info:", nH, nH1, nHp
      !    print *, "spec info:", solnVec(H1_SPEC, i, j, k), solnVec(HP_SPEC, i, j, k)
      !    print *, "ion state:", solnVec(IONR_VAR, i, j, k), solnVec(RECR_VAR, i, j, k)
      !  end if


        if (pt_is_coll_ionization) then
          recomb_rate = recomb_rate * (1.0 + k_coll*dtnow*nH)
        end if
        solnVec(RECR_VAR, i, j, k) = recomb_rate
        solnVec(COOL_VAR, i, j, k) = -cooling_rate ! -sign for convention

        if (pt_is_apply_recombination) then
          solnVec(H1_SPEC, i, j, k) = solnVec(H1_SPEC, i, j, k)&
                                    + recomb_rate
          solnVec(HP_SPEC, i, j, k) = solnVec(HP_SPEC, i, j, k)&
                                    - recomb_rate
          solnVec(ELE_SPEC, i, j, k) = solnVec(HP_SPEC, i, j, k)*Ae

          ! Make sure mass fraction \in [0., 1.]
          call sanitize_mass_fraction(cellID, solnVec)
        end if
      end do
    end do
  end do
#endif
#endif

end subroutine calc_recomb_emissivity


! Subroutine to update the ionization state due to absorption of 
! ionizing photons.
subroutine deposit_ionizing_radiation(solnVec, cellID, particle,&
                                        mcp_fate, dt, dtNew,&
                                        fleckp, k_ion, N_H1, dvol,&
                                        is_empty_cell_event)
  use Simulation_data, only : clight, mH, ev2erg
  use Particles_data, only : pt_is_photoionization, pt_ABS_ID,&
                             pt_is_corrdl, pt_is_cont_photo
  use Multispecies_interface, only : Multispecies_getProperty
  use rhd, only : cellAddVar
  implicit none
#include "constants.h"
#include "Flash.h"
#include "Multispecies.h"

  ! Input/output
  real, pointer :: solnVec(:,:,:,:)
  integer, dimension(MDIM), intent(in) :: cellID
  real, dimension(NPART_PROPS), intent(inout) :: particle
  integer, intent(in) :: mcp_fate
  real, intent(in) :: dt, dtNew, fleckp, k_ion, N_H1, dvol
  logical, intent(in) :: is_empty_cell_event

  ! aux variables
  integer :: i, j, k
  real :: rho, Ae, Ah1, nH1
  real :: k_ea_ion, nump_init, nump_old, nump_new
  real :: energy, d_trav, old_h1, old_hp, old_ele
  real :: delta_n_ion, delta_e_ion
  real :: dtau, dl_corr

  if (.not. pt_is_photoionization) return

#ifdef FLASH_MCRHD_IONIZATION
#ifdef FLASH_MULTISPECIES
  call Multispecies_getProperty(ELE_SPEC, A, Ae)
  i = cellID(IAXIS)
  j = cellID(JAXIS)
  k = cellID(KAXIS)

  ! Tally the absorption of ionizing radiation
  k_ea_ion = fleckp * k_ion
  d_trav = clight * dt
  nump_old  = particle(NUMP_PART_PROP)
  nump_init = particle(NPIN_PART_PROP)
  energy = particle(ENER_PART_PROP)
  rho = solnVec(DENS_VAR, i, j, k)
  old_h1 = solnVec(H1_SPEC, i, j, k)
  old_hp = solnVec(HP_SPEC, i, j, k)
  old_ele = solnVec(ELE_SPEC, i, j, k)

  dtau = k_ea_ion*d_trav
  ! New version
  if (pt_is_cont_photo) then
    ! nump_new = nump_old - (number of neutral H destroyed)
    !          = nump_old - (N_H1_initial - N_H1_final)
    ! With new scheme, N_H1_final = N_H1_initial / (k_ion*dl + 1)
    nump_new = nump_old - N_H1*dtau/ (dtau + 1.0)
   else
    ! Old version
    nump_new = nump_old*exp(-dtau) 
   end if

  ! Correction for absorption during the integrated path 
!  dl_corr = 1.0d0
!  dtau = k_ea_ion*d_trav  ! original tau without correction
!  if (pt_is_corrdl .and. (dtau > 2.0)) then
!    dl_corr = (1.0d0 - exp(-dtau)) / dtau
!    nump_new = nump_old*dl_corr
!  end if

  if ((mcp_fate == pt_ABS_ID) .and. (.not. is_empty_cell_event)) then
    nump_new = 0.0d0
  end if

  ! Tally the absorbed radiation
  delta_n_ion = -(nump_old - nump_new) / dvol
  delta_n_ion = delta_n_ion*mH/rho
  delta_e_ion = (energy - 13.6*ev2erg)*(nump_old - nump_new) / (rho*dvol)

  !if (abs(delta_n_ion) .gt. old_h1) then
  if (.not. is_empty_cell_event) then
    ! Update ionization/mass fraction
    solnVec(H1_SPEC, i, j, k)  = solnVec(H1_SPEC, i, j, k) + delta_n_ion
    solnVec(HP_SPEC, i, j, k)  = solnVec(HP_SPEC, i, j, k) - delta_n_ion
    solnVec(ELE_SPEC, i, j, k) = solnVec(HP_SPEC, i, j, k) * Ae
  else
    delta_n_ion = -old_h1 ! remove all H1 in cell
    delta_e_ion = abs(delta_n_ion)*(energy-13.6*ev2erg) / mH ! rho*dV canceled
    nump_new = nump_old - abs(delta_n_ion)*dvol*rho/mH

    ! Impose the fully ionized case
    solnVec(H1_SPEC, i, j, k)  = 0.0
    solnVec(HP_SPEC, i, j, k)  = 1.0/(1.0 + Ae)
    solnVec(ELE_SPEC, i, j, k) = Ae/(1.0+Ae)
  end if

  ! Check for negative mass fraction?
  if ((solnVec(H1_SPEC, i, j, k) < 0.0) .or.&
      (solnVec(HP_SPEC, i, j, k) < 0.0) .or.&
      (solnVec(ELE_SPEC, i, j, k) < 0.0)) then
    print *, "Abundances:"
    print *, "mcp_fate", mcp_fate
    print *, "H1", old_h1, solnVec(H1_SPEC, i, j, k)
    print *, "HP", old_hp, solnVec(HP_SPEC, i, j, k)
    print *, "ELE", old_ele, solnVec(ELE_SPEC, i, j, k)
    print *, "delta_n_ion", delta_n_ion
    print *, "k_ea_ion", k_ea_ion, fleckp
    print *, "d_trav", d_trav
    print *, "numps", nump_init, nump_old, nump_new
    print *, "is_empty_cell?", is_empty_cell_event
    call Multispecies_getProperty(H1_SPEC, A, Ah1)
    nH1 = rho*old_h1/(Ah1*mH)
    print *, "N_H1 (init)", nH1*dvol
    print *, "dN", delta_n_ion*rho/mH*dvol

    call Driver_abortFlash("deposit_ionizing_rad: negative mass fraction!")
  end if

  if (solnVec(ELE_SPEC,i,j,k) >1.0) then
    print *, "wrongELE-depo", solnVec(ELE_SPEC,i,j,k)
  end if

  ! Accumulate heating/ionization terms
  call cellAddVar(solnVec, cellID, IONR_VAR, delta_n_ion)
  call cellAddVar(solnVec, cellID, HEAT_VAR, delta_e_ion)

  ! Update particle's weight info
  particle(NUMP_PART_PROP) = nump_new

  if (nump_new < 0.0d0) then
    print *, "negative mcp weight"
    print *, "old_h1", old_h1
    print *, "delta_n_ion/delta_e_ion", delta_n_ion, delta_e_ion
    print *, "k_ea_ion*d_trav", k_ea_ion, d_trav
  end if
#endif
#endif

end subroutine deposit_ionizing_radiation


subroutine eff_scatter_ionizing_mcp(solnVec, cellID, particle, is_caseB)
  use Simulation_data, only : ev2erg
  use Driver_interface, only : Driver_abortFlash
  use new_mcp, only : sample_iso_velocity, sample_energy
  use Particles_data, only : pt_is_rm_mcps_caseB
  implicit none

#include "constants.h"
#include "Flash.h"

  ! Input/output
  real, pointer :: solnVec(:,:,:,:)
  integer, dimension(MDIM), intent(in) :: cellID
  real, dimension(NPART_PROPS), intent(inout) :: particle
  logical, intent(in) :: is_caseB

  ! aux variables
  real :: temp
  real, dimension(MDIM) :: old_vel, new_vel
  real :: mcp_eps_bs, mcp_w_now_bs, mcp_w_ini_bs, mcp_energy_bs
  real :: mcp_eps_as, mcp_w_now_as, mcp_w_ini_as, mcp_energy_as

  temp = solnVec(TEMP_VAR, cellID(IAXIS), cellID(JAXIS), cellID(KAXIS))

  ! Getting MCP information
  old_vel = particle(VELX_PART_PROP:VELZ_PART_PROP)
  mcp_eps_bs    = particle(ENER_PART_PROP) !  > 13.6 eV
  mcp_w_now_bs  = particle(NUMP_PART_PROP)
  mcp_w_ini_bs  = particle(NPIN_PART_PROP)
  mcp_energy_bs = mcp_eps_bs * mcp_w_now_bs

  ! Sample photon energy
  if (is_caseB) then
    mcp_eps_as = 10.0*ev2erg  ! reprocessed to non-ionizing
  else
    !call Driver_abortFlash("eff_scatter_ionizing_mcp: non case-B &
    !                        approx. not implemented yet.")
    mcp_eps_as = mcp_eps_bs  ! case A: stay ionizing
  end if

  ! Conserve photon number if photon energy shifted
  mcp_energy_as = mcp_energy_bs ! MCP energy should conserve
  mcp_w_now_as = mcp_energy_as / mcp_eps_as
  mcp_w_ini_as = (mcp_w_ini_bs / mcp_w_now_bs) * mcp_w_now_as

  ! Sample velocity
  call sample_iso_velocity(new_vel)

  ! Update MCP attributes  
  particle(ENER_PART_PROP) = mcp_eps_as
  particle(VELX_PART_PROP:VELZ_PART_PROP) = new_vel

  particle(NUMP_PART_PROP) = mcp_w_now_as
  particle(NUM0_PART_PROP) = mcp_w_now_as
  particle(NPIN_PART_PROP) = mcp_w_ini_as

  if (pt_is_rm_mcps_caseB) then
    particle(ISAB_PART_PROP) = 1.0
    particle(TREM_PART_PROP) = -1.0
  end if

end subroutine eff_scatter_ionizing_mcp


! Make sure after empty cell event H1 is actually 0.0
subroutine sanitize_fully_ionized_cell(cellID, solnVec)
  use Particles_data, only : pt_is_photoionization
  use Driver_interface, only : Driver_abortFlash
  implicit none

#include "constants.h"
#include "Flash.h"
#include "Multispecies.h"

  ! Input/Output
  integer, dimension(MDIM), intent(in) :: cellID
  real, pointer :: solnVec(:,:,:,:)

  ! aux variables
  real :: h1

  if (.not. pt_is_photoionization) return

#ifdef FLASH_MCRHD_IONIZATION
#ifdef FLASH_MULTISPECIES
  if (NSPECIES /= 0) then
    h1 = solnVec(H1_SPEC, cellID(IAXIS), cellID(JAXIS), cellID(KAXIS))

    if (h1 > 0.0) then
      call Driver_abortFlash("sanitize_fully_ionized_cell:&
                              target cell is not fully ionized, aborting...")
    end if
  else ! NSPECIES=0, not having separate species for HI and HII
    call Driver_abortFlash("sanitize_fully_ionized_cell:&
                            multispecies is not defined but ionization state&
                            sanitziation is requested! Aborting...")
  end if
#endif
#endif

end subroutine sanitize_fully_ionized_cell


! Make sure mass fraction is between 0.0 and 1.0
subroutine sanitize_mass_fraction(cellID, solnVec)
  use Multispecies_interface, only : Multispecies_getProperty
  use Driver_interface, only : Driver_abortFlash
  implicit none
#include "constants.h"
#include "Flash.h"
#include "Multispecies.h"

  ! Input/Output
  integer, dimension(MDIM), intent(in) :: cellID
  real, pointer :: solnVec(:,:,:,:)

  ! Aux variable
  integer :: i, j, k
  real:: h1_frac, hp_frac, Ae

  i = cellID(IAXIS)
  j = cellID(JAXIS)
  k = cellID(KAXIS)

#ifdef FLASH_MCRHD_IONIZATION
#ifdef FLASH_MULTISPECIES
  h1_frac = solnVec(H1_SPEC, i, j, k)
  hp_frac = solnVec(HP_SPEC, i, j, k)

  call Multispecies_getProperty(ELE_SPEC, A, Ae)

  ! Check neutral fraction
  if (h1_frac > 1.0) then
    h1_frac = 1.0
  else if (h1_frac < 0.0) then
    h1_frac = 0.0
  end if

  ! Check ionized fraction:w
  if (hp_frac > 1.0) then
    hp_frac = 1.0
  else if (hp_frac < 0.0) then
    hp_frac = 0.0
  end if
  solnVec(H1_SPEC, i, j, k) = h1_frac
  solnVec(HP_SPEC, i, j, k) = hp_frac
  solnVec(ELE_SPEC, i, j, k) = hp_frac*Ae
#endif
#endif
  
end subroutine sanitize_mass_fraction

end module ionization
