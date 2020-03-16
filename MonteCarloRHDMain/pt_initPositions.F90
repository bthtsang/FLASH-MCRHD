!!****if* source/Particles/localAPI/pt_initPositions
!!
!! NAME
!!    pt_initPositions
!!
!! SYNOPSIS
!!
!!    call pt_initPositions(integer(in)  :: blockID,
!!                          logical(out) :: success)
!!
!! DESCRIPTION
!!
!!  Initializes particle locations for one block in the grid.
!!
!!  The initialization of locations may fail, if adding the
!!  number of particles requested (by the choice of initialization
!!  method and runtime parameters) would exceed the maximum allowed
!!  number of particles allowed in any MPI task given by runtime
!!  parameter pr_maxPerProc. In that case, an impementation of this
!!  interface may:
!!    o abort the simulation (by calling Driver_abortFlash); or
!!    o return with success set to FALSE.
!!  If returning with success=FALSE, pt_numLocal may have been
!!  reset to 0.
!!
!! ARGUMENTS
!!
!!  blockID:        local block ID containing particles to create
!!
!!  success:        returns .TRUE. if positions for all particles
!!                  that should be assigned to this block have been
!!                  successfully initialized.
!!
!! SIDE EFFECTS
!!
!!  Updates the variable pt_numLocal (in Particles_data), but adding
!!  to the previous value (if successful) the number of particles that
!!  were placed in the block given by blockID. If unsuccessful, pt_numLocal
!!  may be reset to 0.
!!
!!  The values of position attributes of the newly initialized particles
!!  will be updated.
!!
!!  The values of blk and proc attributes of newly initialized particles
!!  will also be updated appropriately.
!!
!!  The values of velocity attributes of particles MAY also be updated
!!  by some implementations of this interface.
!!
!! NOTES
!!
!!  The values of the tag attribute of particles will not be initialize
!!  here, that is left for later and typically done later in the same
!!  invocation of particles_initpositions that calls this interface.
!!
!! SEE ALSO
!!
!!  Particles_data
!!  Particles_initPositions
!!***


subroutine pt_initPositions (blockID,success)


  implicit none

  integer, INTENT(in) :: blockID
  logical,intent(OUT) :: success

  success = .true. ! DEV: returns true because this stub creates no particles,
                   ! therefore all of those zero particles were created successfully
  return

!----------------------------------------------------------------------
  
end subroutine pt_initPositions


