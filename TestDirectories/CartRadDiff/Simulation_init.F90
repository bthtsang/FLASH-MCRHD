!!****if* source/Simulation/SimulationMain/Sedov/Simulation_init
!!
!! NAME
!!
!!  Simulation_init
!!
!!
!! SYNOPSIS
!!
!!  Simulation_init()
!!
!! DESCRIPTION
!!
!!  Initializes all the data specified in Simulation_data.
!!  It calls RuntimeParameters_get rotuine for initialization.
!!
!! ARGUMENTS
!!
!!   
!!
!! PARAMETERS
!!
!!***

subroutine Simulation_init()
  use Particles_data, ONLY : pt_meshMe
  use Simulation_data 
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Driver_data, ONLY : dr_meshMe
  use reader, ONLY : read_array_from_file
  USE Grid_interface, ONLY : Grid_getDomainBoundBox 
  USE PhysicalConstants_interface, ONLY : PhysicalConstants_get

  implicit none
#include "constants.h"
#include "Flash.h"

  call RuntimeParameters_get('sim_pAmbient', sim_pAmbient)
  call RuntimeParameters_get('sim_rhoAmbient', sim_rhoAmbient)
  call RuntimeParameters_get('sim_tempAmbient', sim_tempAmbient)
  call RuntimeParameters_get('gamma', sim_gamma)
  call RuntimeParameters_get('xmin',sim_xMin)
  call RuntimeParameters_get('ymin',sim_yMin)
  call RuntimeParameters_get('zmin',sim_zMin)
  call RuntimeParameters_get('xmax',sim_xMax)
  call RuntimeParameters_get('ymax',sim_yMax)
  call RuntimeParameters_get('zmax',sim_zMax)
  
  call PhysicalConstants_get("proton mass", mH)
  call PhysicalConstants_get("Boltzmann", kB)
  !CALL PhysicalConstants_get("Planck", h_planck)
  call PhysicalConstants_get("speed of light", clight)
  call PhysicalConstants_get("Stefan-Boltzmann", sigma)
  call PhysicalConstants_get("ideal gas constant", R)
  !CALL PhysicalConstants_get("Newton", G_grav)

  a_rad = (4.0d0 * sigma) / clight

end subroutine Simulation_init
