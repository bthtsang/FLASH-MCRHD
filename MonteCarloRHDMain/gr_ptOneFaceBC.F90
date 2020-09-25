!!****if* source/Grid/GridParticles/gr_ptOneFaceBC
!!
!! NAME
!!
!!  gr_ptOneFaceBC
!!
!! SYNOPSIS
!!
!!  gr_ptOneFaceBC(real(INOUT)     :: particle(propCount),
!!                     integer(IN) :: propCount,
!!                     integer(IN) :: axis,
!!                     integer(IN) :: face,
!!                     integer(INOUT)  :: lostParticles)
!!
!! DESCRIPTION
!!   This routine applies the correct boundary conditions to a particle
!!   that is known to have hit the physical boundary on a given axis
!!   and face
!!
!!   On return, the particle may be changed in one of the following ways:
!!    o  The particle is marked for dropping by setting its BLK_PART_PROP to NONEXISTENT.
!!       When this happens, lostParticles is also incremented by 1.
!!    o  The POS{X,Y,Z}_PART_PROP and possibly other properties like
!!       POSPRED{X,Y,Z}_PART_PROP, VEL{X,Y,Z}_PART_PROP, and/or VELPRED{X,Y,Z}_PART_PROP
!!       are changed so as to implement the action of a "periodic" or "reflecting"
!!       or similar boundary.
!!
!!   It is left to the caller to take appropriate action based on these changes, such as
!!   freeing storage space for a dropped particle or assigning a particle to a different
!!   block (possibly on a different processor).
!!
!! ARGUMENTS
!!
!!   particle  : Data structure holding the properties of one particle.
!!   propCount : The count of fields in the particles data structure
!!   axis : The axis on which to apply BC
!!   face      : indicates which (lower or upper) face boundary to apply
!!
!!  lostParticles : counter for particles that go permanently missing.
!!
!! NOTES
!!   If a setup is using user defined boundary conditions for the fluid, it is very likely
!!   to need a customized version of this routine in the setup directory.
!!
!!   The current implementation does not deal correctly with boundaries at inner boundary
!!   blocks defined with Simulation_defineDomain. Use of this GridParticles implementation
!!   together with Simulation_defineDomain may thus result in undefined behavior when a
!!   particles moves across a boundary block boundary.
!!
!!   In this implementation the following boundary conditions are specifically recognized:
!!
!!       Boundary Condition         Action
!!   ================================================================================
!!       OUTFLOW                    Drop particle
!!       DIODE
!!   --------------------------------------------------------------------------------
!!       REFLECTING                 Reflect particle (normally into the same block)
!!       HYDROSTATIC_NVREFL
!!   --------------------------------------------------------------------------------
!!       PERIODIC                   Move particle (to the opposite side of the domain)
!!   --------------------------------------------------------------------------------
!!
!!   All unrecognized boundary conditions result in dropping, as for OUTFLOW.
!!
!!***

subroutine gr_ptOneFaceBC(particle,propCount, axis, face, blockID, lostParticles)
!! , moved)

#include "constants.h"
#include "Flash.h"

  use Driver_interface, ONLY : Driver_abortFlash
#ifdef FLASH_GRID_PARAMESH
  use tree, ONLY : lrefine
#else
  use Grid_data,ONLY : gr_axisNumProcs
#endif


  use Grid_data,ONLY : gr_imin,gr_imax,gr_jmin,gr_jmax,&
       gr_kmin,gr_kmax,gr_domainBC,&
       gr_geometry ! Added by Benny

  use gr_ptData, ONLY : gr_ptBlk, gr_ptProc, gr_ptTag, &
       gr_ptPosx, gr_ptPosy, gr_ptPosz,&
       gr_ptPos2x, gr_ptPos2y, gr_ptPos2z,&
       gr_ptVelx, gr_ptVely, gr_ptVelz,&
       gr_ptVel2x, gr_ptVel2y, gr_ptVel2z,&
       gr_ptPosTmp, gr_ptVel, gr_ptVelTmp, gr_ptKeepLostParticles

  implicit none

 
  integer, intent(IN) :: propCount
  real,dimension(propCount),intent(INOUT)::particle
  integer, intent(IN) :: axis, face, blockID
  integer, intent(INOUT) :: lostParticles
!!  logical, intent(INOUT) :: moved
  real :: dist
  integer :: pos,vel,posPred,velPred
  real,dimension(LOW:HIGH) :: corner
  logical :: predictor, singlePeriodicBlock

  ! Added by Benny for spherical coordinate system
  real, dimension(MDIM) :: vel_reflected

  ! Benny debugging
  real, dimension(MDIM) :: old_mcp_pos, old_mcp_vel


  predictor=gr_ptPosTmp

  if (axis==IAXIS) then
     corner(LOW)=gr_imin
     corner(HIGH)=gr_imax
     pos=gr_ptPosx
     vel=gr_ptVelx
     posPred=gr_ptPos2x
     velPred=gr_ptVel2x
  elseif(axis==JAXIS) then
     if(NDIM<2)then
        call Driver_abortFlash("gr_ptOneFaceBC, NDIM<2, axis is JAXIS")
     else
        corner(LOW)=gr_jmin
        corner(HIGH)=gr_jmax
        pos=gr_ptPosy
        vel=gr_ptVely
        posPred=gr_ptPos2y
        velPred=gr_ptVel2y
     end if
  elseif(axis==KAXIS) then
     if(NDIM<3)then
        call Driver_abortFlash("gr_ptOneFaceBC, NDIM<3, axis is KAXIS")
     else
        corner(LOW)=gr_kmin
        corner(HIGH)=gr_kmax
        pos=gr_ptPosz
        vel=gr_ptVelz
        posPred=gr_ptPos2z
        velPred=gr_ptVel2z
     end if
  end if
!!  moved=.false.

  !if (axis == JAXIS) then
  !  print *, "crossing BC", axis, particle(pos)
  !  call Driver_abortFlash("gr_ptOneFaceBC, crossing phi boundary")
  !end if
  
  if(gr_domainBC(face,axis)==OUTFLOW) then
     if(gr_ptKeepLostParticles) then
        particle(gr_ptBlk)=LOST
     else
        particle(gr_ptBlk)=NONEXISTENT
        lostParticles=lostParticles+1
     end if
  elseif(gr_domainBC(face,axis)==REFLECTING .OR. &
       gr_domainBC(face,axis)==HYDROSTATIC_NVREFL) then
     !! with reflecting conditions, the particle stays within
     !! the same block, the position changes
     old_mcp_pos = particle(gr_ptPosx:gr_ptPosz)
     old_mcp_vel = particle(gr_ptVelx:gr_ptVelz)
     particle(pos)= 2.0*corner(face)-particle(pos)
     !if (particle(TAG_PART_PROP)==82.) write(*,*)'particle(pos) a=',particle(pos)



     if (gr_geometry == CARTESIAN) then
       particle(vel)= -particle(vel)
     else if (gr_geometry == SPHERICAL) then
       ! flip the radial velocity using the subroutine below
       call reflect_velocity(particle(gr_ptPosx:gr_ptPosz),&
                             particle(gr_ptVelx:gr_ptVelz),&
                             axis, face, vel_reflected)
       particle(gr_ptVelx:gr_ptVelz) = vel_reflected

       ! Benny debugging
!       if (axis == 2) then
!         print *, "theta Bcrossing", face
!         print *, "old pos", old_mcp_pos
!         print *, "new pos", particle(gr_ptPosx:gr_ptPosz) 
!         print *, "old vel", old_mcp_vel
!         print *, "new vel", particle(gr_ptVelx:gr_ptVelz) 
!       end if
     end if

     if(predictor)then
        particle(posPred)= 2.0*corner(face)-particle(posPred)
        particle(velPred)= -particle(velPred)
     end if
  elseif(gr_domainBC(face,axis)==PERIODIC) then
     dist=(-1)**(face-1)*(corner(HIGH)-corner(LOW))
     particle(pos)=dist+particle(pos)
     if(predictor)particle(posPred)=dist+particle(posPred)
  else                       ! default for now - KW
     lostParticles=lostParticles+1
     particle(gr_ptBlk)=NONEXISTENT
  end if

  return
end subroutine gr_ptOneFaceBC


subroutine reflect_velocity(pos, vel, axis, face, vel_reflected)

  implicit none
#include "constants.h"
#include "Flash.h"

  ! Input/output
  real, dimension(MDIM), intent(in) :: pos, vel
  integer, intent(in) :: axis, face
  real, dimension(MDIM), intent(out) :: vel_reflected

  ! aux variables
  real :: theta, phi
  real :: sin_theta, cos_theta, sin_phi, cos_phi
  real, dimension(MDIM) :: sph_vel, sph_vel_reflected
  real, dimension(MDIM) :: r_hat, theta_hat, phi_hat
  real, dimension(MDIM) :: x_hat, y_hat, z_hat
  real :: vr, vt, vp


  if (NDIM == 3) then
    theta = pos(JAXIS)
    phi   = pos(KAXIS)

    sin_theta = sin(theta)
    cos_theta = cos(theta)
    sin_phi   = sin(phi)
    cos_phi   = cos(phi)

    r_hat = (/ sin_theta*cos_phi,&
               sin_theta*sin_phi,&
               cos_theta /)

    theta_hat = (/ cos_theta*cos_phi,&
                   cos_theta*sin_phi,&
                   -sin_theta /)
    phi_hat   = (/ -sin_phi,&
                    cos_phi,&
                    0.0d0 /)
    vr = dot_product(vel, r_hat)
    vt = dot_product(vel, theta_hat)
    vp = dot_product(vel, phi_hat)

    sph_vel = (/ vr, vt, vp /)

    ! Benny debugging
    if ((face == 1) .and. (vt > 0.0)) then
      print *, "wrong face 1"
      print *, "pos1", pos
      print *, "vel1", vel
      print *, "wrongflip1", sph_vel 
    end if
    if ((face == 2) .and. (vt < 0.0)) then
      print *, "wrong face 2"
      print *, "pos2", pos
      print *, "vel2", vel
      print *, "wrongflip2", sph_vel 
    end if

    sph_vel_reflected = sph_vel
    sph_vel_reflected(axis) = -sph_vel_reflected(axis)

    ! Benny debugging
!    if (axis == 2) then
!      print *, "theta crossing", face
!      print *, "flipping", sph_vel
!      print *, "flipped", sph_vel_reflected
!    end if

    ! Convert it back to Cartesian cooridinates
    x_hat = (/ 1.0d0, 0.0d0, 0.0d0 /)
    y_hat = (/ 0.0d0, 1.0d0, 0.0d0 /)
    z_hat = (/ 0.0d0, 0.0d0, 1.0d0 /)

    !vel_reflected(IAXIS) = dot_product(sph_vel_reflected, x_hat)
    !vel_reflected(JAXIS) = dot_product(sph_vel_reflected, y_hat)
    !vel_reflected(KAXIS) = dot_product(sph_vel_reflected, z_hat)

    ! Convert velocity (vr, vt, vp) back to (vx, vy, vz)
    vel_reflected = sph_vel_reflected(IAXIS)*r_hat +&
                    sph_vel_reflected(JAXIS)*theta_hat +&
                    sph_vel_reflected(KAXIS)*phi_hat

  else
    call Driver_abortFlash("reflect_velocity: &
            1D/2D applications not implemented yet.")
  end if

end subroutine reflect_velocity
