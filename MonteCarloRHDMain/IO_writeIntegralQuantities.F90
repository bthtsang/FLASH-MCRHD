!!****if* source/IO/IOMain/IO_writeIntegralQuantities
!!
!!
!!  NAME
!!    IO_writeIntegralQuantities
!!
!!  SYNOPSIS
!!    call IO_writeIntegralQuantities(integer(in) :: isFirst,
!!                                    real(in)    :: simTime)
!!
!!  DESCRIPTION
!!
!!   Compute the values of integral quantities (eg. total energy)
!!   and write them to an ASCII file.  If this is the initial step,
!!   create the file and write a header to it before writing the data.
!!
!!   Presently, this supports 1, 2, and 3-d Cartesian geometry and 2-d
!!   cylindrical geometry (r,z).  More geometries can be added by
!!   modifying the volume of each zone (dvol).
!!
!!   Users should modify this routine if they want to store any
!!   quantities other than default values in the flash.dat.  Make sure
!!   to modify the nGlobalSum parameter to match the number of
!!   quantities written.  Also make sure to modify the header to match
!!   the names of quantities with those calculated in the lsum and
!!   gsum arrays.
!!  
!!  ARGUMENTS
!!    
!!   isFirst - if 1 then write header info plus data, otherwise just write data
!!   simTime - simulation time
!!
!!
!!***

!!REORDER(4):solnData

subroutine IO_writeIntegralQuantities ( isFirst, simTime)

  use IO_data, ONLY : io_restart, io_statsFileName, io_globalComm
  use Grid_interface, ONLY : Grid_getListOfBlocks, &
    Grid_getBlkIndexLimits, Grid_getBlkPtr, Grid_getSingleCellVol, &
    Grid_releaseBlkPtr

   use IO_data, ONLY : io_globalMe, io_writeMscalarIntegrals

  use Particles_data, only : pt_typeInfo, particles, pt_meshMe
  implicit none

#include "Flash_mpi.h"
#include "constants.h"
#include "Flash.h"
#include "Particles.h"
  
  real, intent(in) :: simTime

  integer, intent(in) :: isFirst

  integer :: lb, count
  
  integer :: funit = 99
  integer :: error
  integer :: nGlobalSumUsed, iSum
  
  character (len=MAX_STRING_LENGTH), save :: fname 
  
  integer :: blockList(MAXBLOCKS)

  integer :: blkLimits(HIGH, MDIM), blkLimitsGC(HIGH, MDIM)

#ifdef MAGP_VAR
  integer, parameter ::  nGlobalSumProp = 8              ! Number of globally-summed regular quantities
#else
  integer, parameter ::  nGlobalSumProp = 7 + 4          ! Number of globally-summed regular quantities
#endif
  integer, parameter ::  nGlobalSum = nGlobalSumProp + NMASS_SCALARS ! Number of globally-summed quantities
  real :: gsum(nGlobalSum) !Global summed quantities
  real :: lsum(nGlobalSum) !Global summed quantities

  integer :: ivar
  integer :: i, j, k
  real :: dvol             !, del(MDIM)
  real, DIMENSION(:,:,:,:), POINTER :: solnData

  integer :: point(MDIM)
  integer :: ioStat

! MCRHD additions
  integer :: i_begin, i_count, i_end, p
  integer :: p_blk
  real    :: p_nump, p_ener

  if (io_writeMscalarIntegrals) then
     nGlobalSumUsed = nGlobalSum
  else
     nGlobalSumUsed = nGlobalSumProp
  end if

  ! Sum quantities over all locally held leaf-node blocks.
  gsum(1:nGlobalSumUsed) = 0.
  lsum(1:nGlobalSumUsed) = 0.
  
  call Grid_getListOfBlocks(LEAF, blockList, count)

  do lb = 1, count
     !get the index limits of the block
     call Grid_getBlkIndexLimits(blockList(lb), blkLimits, blkLimitsGC)

     ! get a pointer to the current block of data
     call Grid_getBlkPtr(blockList(lb), solnData)

     ! Sum contributions from the indicated blkLimits of cells.
     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
              
              point(IAXIS) = i
              point(JAXIS) = j
              point(KAXIS) = k

!! Get the cell volume for a single cell
              call Grid_getSingleCellVol(blockList(lb), EXTERIOR, point, dvol)
     
              ! mass   
#ifdef DENS_VAR
              lsum(1) = lsum(1) + solnData(DENS_VAR,i,j,k)*dvol 
#endif           


#ifdef DENS_VAR
#ifdef VELX_VAR      
              ! momentum
              lsum(2) = lsum(2) + solnData(DENS_VAR,i,j,k) * & 
                   &                                solnData(VELX_VAR,i,j,k)*dvol 
           
#endif
#ifdef VELY_VAR      

              lsum(3) = lsum(3) + solnData(DENS_VAR,i,j,k) * & 
                   &                                solnData(VELY_VAR,i,j,k)*dvol
           
#endif
#ifdef VELZ_VAR      
              lsum(4) = lsum(4) + solnData(DENS_VAR,i,j,k) * & 
                   &                                solnData(VELZ_VAR,i,j,k)*dvol
#endif

              ! total energy
#ifdef ENER_VAR
              lsum(5) = lsum(5) + solnData(ENER_VAR,i,j,k) * & 
                   &                                solnData(DENS_VAR,i,j,k)*dvol
#ifdef MAGP_VAR
              ! total plasma energy
!!$              lsum(5) = lsum(5) + (solnData(ENER_VAR,i,j,k) * & 
!!$                   &    solnData(DENS_VAR,i,j,k) + solnData(MAGP_VAR,i,j,k))*dvol

              lsum(5) = lsum(5) + solnData(MAGP_VAR,i,j,k)*dvol
#endif
#endif

           
#ifdef VELX_VAR      
#ifdef VELY_VAR      
#ifdef VELZ_VAR      
              ! kinetic energy
              lsum(6) = lsum(6) + 0.5*solnData(DENS_VAR,i,j,k) * & 
                   &                             (solnData(VELX_VAR,i,j,k)**2+ & 
                   &                              solnData(VELY_VAR,i,j,k)**2+ & 
                   &                              solnData(VELZ_VAR,i,j,k)**2)*dvol           

#endif
#endif
#endif


#ifdef EINT_VAR
              ! internal energy
              lsum(7) = lsum(7) + solnData(DENS_VAR,i,j,k) * & 
                   &                                solnData(EINT_VAR,i,j,k)*dvol
#endif
#endif ! ifdef DENS_VAR

#ifdef URAD_VAR
              ! radiation energy
              lsum(8) = lsum(8) + solnData(URAD_VAR,i,j,k)*dvol
#endif

#ifdef ABSE_VAR
              ! thermal radiation energy absorbed
              lsum(9) = lsum(9) + solnData(ABSE_VAR,i,j,k)*&
                                  solnData(DENS_VAR,i,j,k)*dvol
#endif
#ifdef HEAT_VAR
              ! thermal radiation energy absorbed
              lsum(9) = lsum(9) + solnData(HEAT_VAR,i,j,k)*&
                                  solnData(DENS_VAR,i,j,k)*dvol
#endif

#ifdef EMIE_VAR
              ! thermal radiation energy emitted
              lsum(10) = lsum(10) + solnData(EMIE_VAR,i,j,k)*&
                                    solnData(DENS_VAR,i,j,k)*dvol
#endif
#ifdef COOL_VAR
              ! thermal radiation energy absorbed
              lsum(10) = lsum(10) + solnData(COOL_VAR,i,j,k)*&
                                    solnData(DENS_VAR,i,j,k)*dvol
#endif

#ifdef MAGP_VAR
              ! magnetic energy
              lsum(12) = lsum(12) + solnData(MAGP_VAR,i,j,k)*dvol
#endif

#ifdef DENS_VAR
              if (io_writeMscalarIntegrals) then
                 iSum = nGlobalSumProp
!!$                 do ivar=MASS_SCALARS_BEGIN,MASS_SCALARS_END
                    lsum(iSum+1:iSum+NMASS_SCALARS) = &
                         lsum(iSum+1:iSum+NMASS_SCALARS) + &
                           solnData(DENS_VAR,i,j,k) * & 
                           solnData(MASS_SCALARS_BEGIN: &
                                    MASS_SCALARS_END,i,j,k)*dvol
!!$                 end do
              end if
#endif
           enddo
        enddo
     enddo
     call Grid_releaseBlkPtr(blockList(lb), solnData)

#ifdef PHOTON_PART_TYPE
     ! Added by Benny to track also total MCP radiation energy
     ! The pt_typeInfo data structure tracks the PART_LOCAL number, 
     ! which is the number of type-TYPE MCPs that are residing in
     ! the local leaf blocks. 'particles' array contains local MCPs,
     ! so no PROC_PART_PROP filtering is needed.
     i_begin = pt_typeInfo(PART_TYPE_BEGIN, PHOTON_PART_TYPE)
     i_count = pt_typeInfo(PART_LOCAL, PHOTON_PART_TYPE)
     i_end   = i_begin + i_count - 1

     do p = i_begin, i_end
       p_blk = int(particles(BLK_PART_PROP, p))

       ! count for the current blkID
       if (p_blk == blockList(lb)) then
         p_nump = particles(NUMP_PART_PROP, p)
         p_ener = particles(ENER_PART_PROP, p)
         lsum(11) = lsum(11) + p_nump*p_ener
       end if
     end do
     ! End E_MCPs summation
#endif

  enddo
  

  
  ! Now the MASTER_PE sums the local contributions from all of
  ! the processors and writes the total to a file.
  
  call MPI_Reduce (lsum, gsum, nGlobalSumUsed, FLASH_REAL, MPI_SUM, & 
       &                MASTER_PE, io_globalComm, error)
  

  if (io_globalMe  == MASTER_PE) then
     
     ! create the file from scratch if it is a not a restart simulation, 
     ! otherwise append to the end of the file
     
     !No matter what, we are opening the file. Check to see if already there
     ioStat = 0
     open(funit, file=trim(io_statsFileName), position='APPEND', status='OLD', iostat=ioStat)
     if(ioStat .NE. 0) then
        !print *, 'FILE FOUND'
        open(funit, file=trim(io_statsFileName), position='APPEND')
     endif
     
     if (isFirst .EQ. 1 .AND. (.NOT. io_restart .or. ioStat .NE. 0)) then
        
#ifndef MAGP_VAR
        write (funit, 10)               &
             '#time                     ', &
             'mass                      ', &
             'x-momentum                ', &
             'y-momentum                ', & 
             'z-momentum                ', &
             'E_total                   ', &
             'E_kinetic                 ', &
             'E_internal                ', &
             'E_radiation               ', &
             'E_absorbed                ', &
             'E_emitted                 ', &
             'E_MCPs                    ', &
             (msName(ivar),ivar=MASS_SCALARS_BEGIN,&
              min(MASS_SCALARS_END,&
                  MASS_SCALARS_BEGIN+nGlobalSumUsed-nGlobalSumProp-1))

#else
        
        write (funit, 10)               &
             '#time                     ', &
             'mass                      ', &
             'x-momentum                ', &
             'y-momentum                ', & 
             'z-momentum                ', &
             'E_total                   ', &
             'E_kinetic                 ', &
             'E_internal                ', &
             'MagEnergy                 ', &
             (msName(ivar),ivar=MASS_SCALARS_BEGIN,&
              min(MASS_SCALARS_END,&
                  MASS_SCALARS_BEGIN+nGlobalSumUsed-nGlobalSumProp-1))
#endif
        
10         format (2x,50(a25, :, 1X))

     else if(isFirst .EQ. 1) then
        write (funit, 11) 
11      format('# simulation restarted')
     endif
     
     ! Write the global sums to the file.
     write (funit, 12) simtime, gsum(1:nGlobalSumUsed)

12   format (1x, 50(es25.18, :, 1x))
 
     close (funit)          ! Close the file.
     
  endif
  
  call MPI_Barrier (io_globalComm, error)
  
  !=============================================================================
  
  return

  contains
    character(len=25) function msName(ivar)
      integer,intent(in) :: ivar
      character(len=25) :: str
      call Simulation_mapIntToStr(ivar,str,MAPBLOCK_UNK)
      msName = str
    end function msName
end subroutine IO_writeIntegralQuantities



