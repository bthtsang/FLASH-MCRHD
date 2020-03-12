!!****if* source/Particles/localAPI/pt_advanceCustom
!!
!! NAME
!!
!!  pt_advanceCharged
!!
!! SYNOPSIS
!!
!!  call pt_advanceCharged(real(in)    :: dtold,
!!                        real(in)    :: dtnew,
!!                        real(inout) :: particlesunused(NPART_PROPS,p_countUnused),
!!                        integer(in) :: p_countunused)
!!
!! DESCRIPTION
!!
!!   Advances particles in time
!!
!! ARGUMENTS
!!
!!   dtOld : previous time interval 
!!   dtNew : current time interval
!!   particlesunused -- particles on which to operate
!!   p_countunused - the number of particles in the list to advance
!!
!!
!!
!!***

subroutine pt_advanceCustom(dtOld, dtNew, particles, p_count, ind)

  use Particles_data, only : useParticles, pt_maxPerProc   
  use Timers_interface, only : Timers_start, Timers_stop  
  use emission, only : emit_mcps
  use transport, only : transport_mcps
  use rhd, only : apply_rad_source_terms
  implicit none


#include "Flash.h"
#include "constants.h"
 
  ! Input/Output 
  real, intent(in)  :: dtOld, dtNew
  integer, intent(in) :: p_count, ind
  real,dimension(NPART_PROPS,pt_maxPerProc),intent(inout) :: particles

  ! Other parameters
  integer, dimension(MAXBLOCKS) :: blkList
  integer :: numofblks
  real, pointer :: solnVec(:,:,:,:)

  ! Quit if particle module is not in use
  if (.not. useParticles) return


  ! Start performing RT
  call Timers_start("Radiation Transport")

  ! Emission of radiation, including both thermal, point, and face
  call emit_mcps(particles, p_count, dtNew, ind)

  ! Radiation transport
  call transport_mcps(dtOld, dtNew, particles, p_count, pt_maxPerProc, ind)

  ! Deposit radiation source terms
  call apply_rad_source_terms(dtNew)

  ! End of current RT step
  call Timers_stop("Radiation Transport")
  
end subroutine pt_advanceCustom
