!!****if* source/Grid/GridParticles/GridParticlesMove/paramesh/gr_ptMoveOffBlk
!!
!! NAME
!!
!!  gr_ptMoveOffBlk
!!
!! SYNOPSIS
!!
!!  gr_ptMoveOffBlk(real(INOUT)   :: dataBuf(propCount,bufferDim2),
!!                 integer(in)     :: propCount,
!!                 integer(INOUT) :: localCount,
!!                 integer(IN)    :: bufferDim2,
!!                 real(INOUT)    :: destBuf(propCount,bufferDim2),
!!                 integer(INOUT) :: numDest)
!!
!! DESCRIPTION
!!     
!!    This routine is used in moving the nonstationary data elements 
!!    associated with structures like particles and rays when a data element
!!    moves off a block without re gridding. Here every element currently 
!!    on the processor is examined to see if it still belongs to the same block.
!!    If it does not, it is further examimned to see if it has moved out of the physical boundary.
!!    If is out of physical boundary, it may either leave the domain, stay on the same block,
!!    or be moved to destBuf (which holds elements to be passed to the next processor), depending
!!    on the boundary conditions. If it is still in the physical domain, it may have
!!    moved to another block on the same processor, in which case only its BLK
!!    needs to change, otherwise it is moved to destBuf.
!!
!! ARGUMENTS
!!
!!     dataBuf -           A 2 dimensional array of elements and their property
!!                           values
!!     propCount - number of element attributes
!!     localCount -   The number of valid elements in the element array that
!!                           are stored on this processor.
!!     bufferDim2 -   The second dimension of the dataBuf and destBuf buffer arrays.
!!     destBuf   -           temporary storage for elements that need to move off processor
!!     numDest   -           number of elements in destBuf
!! 
!!
!!***

#ifdef DEBUG_ALL
#define DEBUG_PARTICLES
#endif
subroutine gr_ptMoveOffBlk(dataBuf,propCount,localCount,bufferDim2,destBuf,numDest)

#include "constants.h"
#include "Flash.h"
  use gr_ptData, ONLY : gr_ptBlkList, gr_ptBlk, gr_ptProc,&
       gr_ptPosx,gr_ptPosy,gr_ptPosz, gr_ptKeepLostParticles
  use Grid_data, ONLY : gr_meshMe
  use Grid_interface, ONLY : Grid_getListOfBlocks, &
    Grid_getBlkBoundBox, Grid_outsideBoundBox
  use gr_ptInterface, ONLY : gr_ptApplyBC
  use gr_interface, ONLY:  gr_findNeghID, gr_findBlock

  implicit none
  integer, intent(IN) :: propCount,bufferDim2
  integer,intent(INOUT)::localCount,numDest
  real,dimension(propCount,bufferDim2),intent(INOUT)::dataBuf,destBuf
  integer :: blockID,lostElements,blkCount,newBlockID
  real,dimension(LOW:HIGH,MDIM) :: bndBox
  real,dimension(MDIM) :: pos
  integer,dimension(MDIM) :: Negh, pxyz
  integer,dimension(BLKNO:PROCNO) :: neghID
  logical :: moved
  integer i,k,j

  numDest = 0
  lostElements=0
  pxyz = (/gr_ptPosx,gr_ptPosy,gr_ptPosz/)
  pos = 0
  
  call Grid_getListOfBlocks(LEAF,gr_ptBlkList,blkCount)
  i=1
  j=localCount
  k=0
  do while(k < localCount)
     blockID=int(dataBuf(gr_ptBlk,i))  !! Find the ID of the block the element
                       
     !! belonged to before time advance
     if((blockID/=NONEXISTENT).and.(blockID/=LOST)) then
!
!     x,y,z positions in databuf are not necessarily contiguous and so we copy each
!     array element separately (no array slices)
!
        pos(1:NDIM) = dataBuf(pxyz(1:NDIM),i)

        k=k+1
        call Grid_getBlkBoundBox(blockID,bndBox)    !! get its bound box, and find if 
        !! the element has moved out
        call Grid_outsideBoundBox(pos,bndBox,moved,Negh)

        if(.not.moved) then
           i=i+1
        else

           !We know "Negh" so find the relevant "neghID".  
           !We may later ignore "neghID" after considering boundary conditions. 
           call gr_findNeghID(blockID,pos,Negh,neghID)

           !! if elements moved out of its block, find if it moved out 
           !! of external boundary
           !! Update the logical variable "moved" in gr_ptApplyBCs as it is dependant 
           !! upon the external boundary condition, i.e. set to .false. in the
           !! following situations:
           !! 1. Periodic boundary conditions and one block.  
           !! 2. Reflecting boundary conditions.
           call gr_ptApplyBC(dataBuf(:,i),propCount, blockID,lostElements,moved,Negh,bndBox)

           if(dataBuf(gr_ptBlk,i)==NONEXISTENT) then
              dataBuf(:,i)=dataBuf(:,j)
              j=j-1
           else
              if(.not.moved)then
                 i=i+1
              else
                 !if we got here the particle has moved to
                 !another block which exists, and its identity can
                 !be found. 


#ifdef BITTREE
                 !The element position may have been updated in gr_ptApplyBC 
                 !so get fresh values for "pos". This is needed to be done
                 !only when bittree is in operation because in dealing with
                 !periodic condition, the neghID found before updating
                 !the particle position is valid if surrblks are used
                 !but using the bittree it is not valid. Bittree based
                 !gr_findNeghID only finds the block which has the pos
                 !so the block before boundary conditions are applied
                 !would be wrong with periodic conditions.
                 pos(1:NDIM) = dataBuf(pxyz(1:NDIM),i)
                 call gr_findNeghID(blockID,pos,Negh,neghID)
#endif
                 if(neghID(BLKNO)==NONEXISTENT)then
                    dataBuf(gr_ptBlk,i)=UNKNOWN
                    !NOTE: This code is expected to get invoked when we build
                    !FLASH with PARAMESH2.  This is because PARAMESH2 does not 
                    !store complete block neighbor metadata, meaning the 
                    !gr_findNeghID call may return neghID(BLKNO)=NONEXISTENT.
#ifndef FLASH_GRID_PARAMESH2
                    print*,'moveOffBlk:',gr_meshMe,k,blockID,negh(1:NDIM),pos(1:NDIM)
#endif

                    !The element position may have been updated in gr_ptApplyBC 
                    !so get fresh values for "pos".
                    pos(1:NDIM) = dataBuf(pxyz(1:NDIM),i)
                    neghID(BLKNO)=UNKNOWN
                    call gr_findBlock(gr_ptBlkList,blkCount,pos,neghID(BLKNO))
                    if(neghID(BLKNO)/=NONEXISTENT)neghID(PROCNO)=gr_meshMe

#ifdef DEBUG_PARTICLES
                 else if(neghID(BLKNO) .LE. 0 .OR. neghID(PROCNO) < 0)then
9991                format(I5,': gr_findNeghID(',I4,',pos,..)-> Negh=',&
                         I1,I2,I2,', neghID=',I4,I5)
                    print 9991,gr_meshMe,blockID,Negh,neghID
9992                FORMAT(I5,': SUSPICIOUS dataBuf(gr_ptBlk,',I6,')<-',I9)
                    print 9992,gr_meshMe,i,neghID(PROCNO),neghID(BLKNO)
                    dataBuf(gr_ptProc,i)=neghID(PROCNO)
                    dataBuf(gr_ptBlk,i) = neghID(BLKNO)
#endif
                 else
                    dataBuf(gr_ptProc,i)=neghID(PROCNO)
                    dataBuf(gr_ptBlk,i) = neghID(BLKNO)
                 end if
                 if(neghID(PROCNO)==gr_meshMe) then
                    dataBuf(gr_ptBlk,i)=neghID(BLKNO)
                    dataBuf(gr_ptProc,i)=neghID(PROCNO)
                    i=i+1
                 else
                    numDest=numDest+1                    !! put it in the destBuf
                    destBuf(:,numDest)=dataBuf(:,i)
                    dataBuf(:,i)=dataBuf(:,j)
                    j=j-1
                 end if
              end if
           end if
           
        endif
     ! Benny added for saving LOST MCPs
     else if (blockID==LOST) then
       k=k+1
       i=i+1
     end if
  end do
  localCount=localCount-numDest
  if(.not.gr_ptKeepLostParticles)localCount=localCount-lostElements
end subroutine gr_ptMoveOffBlk
