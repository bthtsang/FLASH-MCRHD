!!****f* source/Simulation/Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!!
!! SYNOPSIS
!!
!!  Simulation_initBlock(integer(IN) :: blockId)
!!
!!
!!
!! DESCRIPTION
!!  This routine applies initial conditions of a specific simulation
!!  to the specified block.
!!
!! 
!! ARGUMENTS
!!
!!  blockId -         the number of the block to update
!!
!! PARAMETERS
!!
!!  eosModeInit -     after this routine sets up initial conditions,
!!                    the grid package calls Eos to insure the
!!                    values are thermodynamically consistent.  This
!!                    parameter controls the mode of application of
!!                    the Eos.  Its default is "dens_ie", and it can
!!                    be one of "dens_ie", "dens_pres", "dens_temp".
!!                    Setting this value to dens_temp, for instance,
!!                    would make it possible to leave this routine
!!                    having just initialized density and temperature,
!!                    leaving Eos to calculate the rest.
!!
!! SEE ALSO
!!
!!  Eos_wrapped
!!***

subroutine Simulation_initBlock(blockId)
  ! the inclusion of pAmbient and rhoAmbient for hydro/eos usage
  use Simulation_data, ONLY: sim_pAmbient, sim_rhoAmbient, sim_tempAmbient, sim_gamma,&
                             R, sim_bulkvelaxis, sim_bulkvelval
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, Grid_putPointData,&
                             Grid_getCellCoords, Grid_getMinCellSize
  USE Eos_data, ONLY : eos_singleSpeciesA

  implicit none
#include "constants.h"
#include "Flash.h"
 
  integer, intent(in) :: blockId
  integer  ::  i, j, k
  integer,dimension(2,MDIM) :: blkLimits,blkLimitsGC
  integer,dimension(MDIM) :: axis
  ! for hydro/eos
  REAL :: vx, vy, vz, p, ek, e, rho, gamma
  REAL :: temp
  ! for cross-cell advancement checking
  REAL, DIMENSION(LOW:HIGH, MDIM) :: bndBox
  REAL, DIMENSION(MDIM) :: deltaCell
  INTEGER, PARAMETER :: ycoordsize = NYB + 2*NGUARD
  REAL, DIMENSION(ycoordsize) :: ycoords
  INTEGER, PARAMETER :: xcoordsize = NXB + 2*NGUARD
  REAL, DIMENSION(xcoordsize) :: xcoords
  INTEGER, PARAMETER :: zcoordsize = NZB + 2*NGUARD
  REAL, DIMENSION(zcoordsize) :: zcoords

  CALL Grid_getCellCoords(IAXIS, blockId, CENTER, .TRUE., xcoords, xcoordsize)
  if (NDIM > 1) then
    CALL Grid_getCellCoords(JAXIS, blockId, CENTER, .TRUE., ycoords, ycoordsize)

    if (NDIM > 2) then
      CALL Grid_getCellCoords(KAXIS, blockId, CENTER, .TRUE., zcoords, zcoordsize)
    end if
  end if

  call Grid_getBlkIndexLimits(blockId,blkLimits,blkLimitsGC)

  do k = blkLimitsGC(LOW,KAXIS), blkLimitsGC(HIGH,KAXIS)
    do j = blkLimitsGC(LOW, JAXIS), blkLimitsGC(HIGH, JAXIS)
      do i = blkLimitsGC(LOW,IAXIS), blkLimitsGC(HIGH, IAXIS)

        axis(IAXIS)=i
        axis(JAXIS)=j    
        axis(KAXIS)=k

        vx = 0.0d0
        vy = 0.0d0
        vz = 0.0d0

        ! Initialize velocity in the selected direction
        if (sim_bulkvelaxis == IAXIS) then
         vx = sim_bulkvelval
        else if (sim_bulkvelaxis == IAXIS) then
         vy = sim_bulkvelval
        else if (sim_bulkvelaxis == IAXIS) then
         vz = sim_bulkvelval
        end if

        ! gas kinetic energy
        ek = 0.5 * (vx*vx + vy*vy + vz*vz)

        rho = sim_rhoAmbient  ! fixed nH to be 100.
        temp = sim_tempAmbient
        p = (R/eos_singleSpeciesA)*rho*temp ! just dummy value, init with dens_temp
        e = p / (sim_gamma - 1.0)
        e = e / rho
        e = e + ek

        call Grid_putPointData(blockId, CENTER, DENS_VAR, EXTERIOR, axis, rho)
        call Grid_putPointData(blockId, CENTER, TEMP_VAR, EXTERIOR, axis, temp)
        call Grid_putPointData(blockId, CENTER, PRES_VAR, EXTERIOR, axis, p)
        call Grid_putPointData(blockId, CENTER, EINT_VAR, EXTERIOR, axis, e - ek)
        call Grid_putPointData(blockId, CENTER, ENER_VAR, EXTERIOR, axis, e)
        call Grid_putPointData(blockId, CENTER, GAME_VAR, EXTERIOR, axis, sim_gamma)
        call Grid_putPointData(blockId, CENTER, GAMC_VAR, EXTERIOR, axis, sim_gamma)
        call Grid_putPointData(blockId, CENTER, VELX_VAR, EXTERIOR, axis, vx)
        call Grid_putPointData(blockId, CENTER, VELY_VAR, EXTERIOR, axis, vy)
        call Grid_putPointData(blockId, CENTER, VELZ_VAR, EXTERIOR, axis, vz)
      enddo
    enddo
  enddo
 
  return
end subroutine Simulation_initBlock
