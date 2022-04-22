module spherical

  implicit none
  contains

subroutine get_cartesian_position(sph_pos, cart_pos)

  implicit none
#include "constants.h"

  ! Input/output
  real, dimension(MDIM), intent(in) :: sph_pos
  real, dimension(MDIM), intent(out) :: cart_pos

  ! aux variables
  real :: pos_x, pos_y, pos_z
  real :: pos_r, pos_t, pos_p

  pos_r = sph_pos(IAXIS)
  pos_t = sph_pos(JAXIS) ! p_t = 1.0 in 1D
  pos_p = sph_pos(KAXIS) ! p_p = 1.0 in 1D/2D

  pos_x = pos_r * sin(pos_t) * cos(pos_p)
  pos_y = pos_r * sin(pos_t) * sin(pos_p)
  pos_z = pos_r * cos(pos_t)

  cart_pos = (/ pos_x, pos_y, pos_z /)

end subroutine get_cartesian_position


subroutine get_spherical_position(cart_pos, sph_pos, debug)
  use Particles_data, only : pt_quadrant_I, pt_quadrant_II,&
                             pt_quadrant_III, pt_quadrant_IV
  implicit none
#include "constants.h"

  ! Input/output
  real, dimension(MDIM), intent(in) :: cart_pos
  real, dimension(MDIM), intent(out) :: sph_pos
  logical, optional :: debug

  ! aux variables
  real :: pos_x, pos_y, pos_z
  real :: pos_r, pos_t, pos_p
  integer :: quadrant

  pos_x = cart_pos(IAXIS)
  pos_y = cart_pos(JAXIS)
  pos_z = cart_pos(KAXIS)

  pos_r = sqrt(pos_x*pos_x + pos_y*pos_y + pos_z*pos_z)
  pos_t = acos(pos_z/pos_r)

  if (pos_x /= 0.0d0) then
    pos_p = atan(pos_y/pos_x)
  else
    if (pos_y > 0.0d0) pos_p = 0.5*PI
    if (pos_y < 0.0d0) pos_p = 1.5*PI
  end if
  ! Use the following safe version that takes into
  ! account different quadrants
  ! + PI to offset it to [0, 2PI]
  !pos_p = atan2(pos_y, pos_x) + PI

  ! Adjusting for different quadrants
  ! Quadrant II or III
  call get_quadrant(pos_x, pos_y, quadrant)

  if ((quadrant == pt_quadrant_II) .or.&
      (quadrant == pt_quadrant_III)) then
    pos_p = pos_p + PI
  else if (quadrant == pt_quadrant_IV) then
    pos_p = pos_p + 2.0d0*PI
  end if

  if (present(debug)) then
    print *, "pos_x/y", pos_x, pos_y
    print *, "pos_r", pos_r
    print *, "pos_t", pos_t
    print *, "pos_p", pos_p
  end if

  if (pos_p < 0.0d0) then
    print *, "NegGetSphPos"
    print *, "Pos_xyz", pos_x, pos_y, pos_z
    print *, "Pos_rtp", pos_r, pos_t, pos_p
    print *, "quadrant", quadrant
  end if

  sph_pos = (/ pos_r, pos_t, pos_p /)

end subroutine get_spherical_position


subroutine get_quadrant(posx, posy, quadrant)
  use Particles_data, only : pt_quadrant_I, pt_quadrant_II,&
                             pt_quadrant_III, pt_quadrant_IV
  implicit none

  real, intent(in) :: posx, posy
  integer, intent(out) :: quadrant

  quadrant = pt_quadrant_I

  if ((posx < 0.0d0) .and. (posy >= 0.0d0)) then
    ! Putting pos_y = 0 case and add PI
    quadrant = pt_quadrant_II
  else if ((posx < 0.0d0) .and. (posy < 0.0d0)) then
    quadrant = pt_quadrant_III
  else if ((posx > 0.0d0) .and. (posy < 0.0d0)) then
    quadrant = pt_quadrant_IV
  end if

end subroutine


subroutine get_spherical_velocity(pos, vel, sph_vel)

  implicit none
#include "constants.h"
#include "Flash.h"

  ! Input/output
  real, dimension(MDIM), intent(in) :: pos, vel
  real, dimension(MDIM), intent(out) :: sph_vel

  ! aux variables
  real :: theta, phi
  real :: sin_theta, cos_theta, sin_phi, cos_phi
  real, dimension(MDIM) :: r_hat, theta_hat, phi_hat
  real :: vr, vt, vp


  if (NDIM <= 3) then
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
  else
    call Driver_abortFlash("get_spherical_velocity: &
            1D/2D applications not implemented yet.")
  end if

end subroutine get_spherical_velocity

! Helper subroutine to convert (vr, vt, vp) into (vx, vy, vz)
subroutine get_cartesian_velocity(pos, vel, cart_vel)

  implicit none
#include "constants.h"
#include "Flash.h"

  ! Input/output
  real, dimension(MDIM), intent(in) :: pos, vel
  real, dimension(MDIM), intent(out) :: cart_vel

  ! aux variables
  real :: theta, phi
  real :: sin_theta, cos_theta, sin_phi, cos_phi
  real, dimension(MDIM) :: r_hat, theta_hat, phi_hat
  real :: vr, vt, vp


  if (NDIM <= 3) then
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

    ! add the three components
    cart_vel = vel(IAXIS)*r_hat
    cart_vel = cart_vel + vel(JAXIS)*theta_hat
    cart_vel = cart_vel + vel(KAXIS)*phi_hat
  else
    call Driver_abortFlash("get_spherical_velocity: &
            1D/2D applications not implemented yet.")
  end if

end subroutine get_cartesian_velocity

! subroutine to compute surface area of spherical cells
subroutine get_spherical_cell_face_area(bndBox, deltaCell,&
                                        axis, faceID, area)
  implicit none
#include "constants.h"
#include "Flash.h"

  ! Input/Output
  real, dimension(LOW:HIGH, MDIM), intent(in) :: bndBox
  real, dimension(MDIM), intent(in) :: deltaCell
  integer, intent(in) :: axis
  integer, dimension(MDIM) :: faceID
  real, intent(out) :: area

  ! aux variables
  real :: r_in, r_out, theta_in, theta_out, phi_in, phi_out
  real :: r_0, theta_0

  ! radial, faceID(IAXIS) is the face index
  if (axis == IAXIS) then

    r_0 = bndBox(LOW, IAXIS) ! left boundary
    r_0 = r_0 + (faceID(IAXIS) - NGUARD - 1)*deltaCell(IAXIS)

    area = 4.0*PI*r_0*r_0

    if (NDIM >= 2) then

      theta_in = bndBox(LOW, JAXIS) ! left boundary
      theta_in = theta_in + (faceID(JAXIS) - NGUARD - 1)*deltaCell(JAXIS)
      theta_out = theta_in + deltaCell(JAXIS)

      area = 2.0*PI*r_0*r_0*(cos(theta_in) - cos(theta_out))

      if (NDIM == 3) then
        phi_in = bndBox(LOW, KAXIS) ! left boundary
        phi_in = phi_in + (faceID(KAXIS) - NGUARD - 1)*deltaCell(KAXIS)
        phi_out = phi_in + deltaCell(KAXIS)

        area = r_0*r_0*(phi_out - phi_in)*(cos(theta_in) - cos(theta_out))
      end if
    end if

  ! theta, faceID(JAXIS) is the face index
  else if (axis == JAXIS) then

    theta_0 = bndBox(LOW, JAXIS)
    theta_0 = theta_0 + (faceID(JAXIS) - NGUARD - 1)*deltaCell(JAXIS)

    r_in = bndBox(LOW, IAXIS)
    r_in = r_in + (faceID(IAXIS) - NGUARD - 1)*deltaCell(IAXIS)
    r_out = r_in + deltaCell(IAXIS)

    if (NDIM >= 2) then
      ! 2D result
      area = PI * sin(theta_0) * (r_out*r_out - r_in*r_in)

      if (NDIM == 3) then
        phi_in = bndBox(LOW, KAXIS) ! left boundary
        phi_in = phi_in + (faceID(KAXIS) - NGUARD - 1)*deltaCell(KAXIS)
        phi_out = phi_in + deltaCell(KAXIS)

        area = 0.5 * sin(theta_0) * (phi_out - phi_in) *&
                     (r_out*r_out - r_in*r_in)
      end if

    end if
    if (NDIM == 1) then
      call Driver_abortFlash("get_spherical_cell_face_area: 1D but axis == JAXIS.")
    end if

  ! phi
  else if (axis == KAXIS) then
    r_in = bndBox(LOW, IAXIS)
    r_in = r_in + (faceID(IAXIS) - NGUARD - 1)*deltaCell(IAXIS)
    r_out = r_in + deltaCell(IAXIS)

    theta_in = bndBox(LOW, JAXIS) ! left boundary
    theta_in = theta_in + (faceID(JAXIS) - NGUARD - 1)*deltaCell(JAXIS)
    theta_out = theta_in + deltaCell(JAXIS)    


    if (NDIM == 3) then
      area = 0.5 * (theta_out - theta_in) * (r_out*r_out - r_in*r_in)
    else
      call Driver_abortFlash("get_spherical_cell_face_area: 1/2D but axis == KAXIS.")
    end if

  end if

end subroutine get_spherical_cell_face_area

end module spherical
