!!****if* source/Particles/ParticlesMain/Particles_finalize
!!
!! NAME
!!    Particles_finalize
!!
!! SYNOPSIS
!!    Particles_finalize( )
!!
!! DESCRIPTION
!!
!!    Finalize routine for the particle unit.  Removes memory usage
!!      set in Particles_init.
!!
!! ARGUMENTS
!!
!! PARAMETERS
!!
!!
!!***
  
subroutine Particles_finalize()
  use Particles_data, ONLY : particles,useParticles,&
                             pt_is_grey,&
                             pt_energy_grid, pt_energy_centers,&
                             pt_delta_energy
  use pt_interface, ONLY : pt_initFinalize
  implicit none
  if(useParticles) deallocate(particles)
  if(useParticles .and. (.not. pt_is_grey)) then
    deallocate(pt_energy_grid)
    deallocate(pt_energy_centers)
    deallocate(pt_delta_energy)
  end if
  call pt_initFinalize()
end subroutine Particles_finalize
