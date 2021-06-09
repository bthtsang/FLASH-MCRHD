module new_mcp

  implicit none
  contains

subroutine sample_cell_position(bndBox, deltaCell, cellID, newxyz)
  use Grid_data, only: gr_geometry
  use random, only : randnozero
  implicit none

#include "constants.h"
#include "Flash.h"

  ! Input/Output
  real, dimension(LOW:HIGH, MDIM), intent(in) :: bndBox
  real, dimension(MDIM), intent(in) :: deltaCell
  integer, dimension(MDIM), intent(in) :: cellID
  real, dimension(MDIM), intent(out) :: newxyz

  ! aux variables
  real :: r_in, r_out, r_rand
  real :: theta_in, theta_out, theta_rand
  real :: cos_theta_in, cos_theta_out, cos_theta_rand
  real :: phi_in, phi_out, phi_rand
  real :: xi
  real :: small_frac = 1.0e-15

  newxyz = 1.0d0 ! initialization 

  if (gr_geometry == SPHERICAL) then
    xi = randnozero()
    r_in  = bndBox(LOW, IAXIS) + (cellID(IAXIS) - NGUARD - 1) * deltaCell(IAXIS)
    r_out = r_in + deltaCell(IAXIS)
    ! resample for safety
    if (xi <= small_frac*(r_in/(r_out - r_in))**3) xi = randnozero()
    r_rand = (r_in**3 + xi*(r_out - r_in)**3 )**(1.0/3.0)
    newxyz(IAXIS) = r_rand

    if (r_rand == r_in) then
      print *, "sampled on edge", xi, r_rand 
      print *, "bndBox", bndBox(:, IAXIS)
      print *, "deltaCel", deltaCell
      print *, "rin/out", r_in, r_out
    end if

    if (NDIM >= 2) then
      theta_in  = bndBox(LOW, JAXIS) + (cellID(JAXIS) - NGUARD - 1)*deltaCell(JAXIS)
      theta_out = theta_in + deltaCell(JAXIS)
      cos_theta_in  = cos(theta_in)
      cos_theta_out = cos(theta_out)

      cos_theta_rand = cos_theta_out + randnozero() * (cos_theta_in - cos_theta_out)
      theta_rand = acos(cos_theta_rand)
      newxyz(JAXIS) = theta_rand

      if (NDIM == 3) then
        phi_in  = bndBox(LOW, KAXIS) + (cellID(KAXIS) - NGUARD - 1)*deltaCell(KAXIS)
        phi_out = phi_in + deltaCell(KAXIS)
        phi_rand = phi_in + (phi_out - phi_in) * randnozero()
        newxyz(KAXIS) = phi_rand
      end if
    end if

  else if (gr_geometry == CARTESIAN) then

    newxyz(IAXIS) = bndBox(LOW, IAXIS) + (cellID(IAXIS) - NGUARD - 1 + randnozero()) * deltaCell(IAXIS)
    if (NDIM >= 2) then
      newxyz(JAXIS) = bndBox(LOW, JAXIS) + (cellID(JAXIS) - NGUARD - 1 + randnozero()) * deltaCell(JAXIS)

      if (NDIM == 3) then
        newxyz(KAXIS) = bndBox(LOW, KAXIS) + (cellID(KAXIS) - NGUARD - 1 + randnozero()) * deltaCell(KAXIS)
      end if
    end if

  end if

end subroutine sample_cell_position


subroutine sample_blk_position(bndBox, newxyz)
  use Grid_data, only: gr_geometry
  use random, only : rand
  implicit none

#include "constants.h"
#include "Flash.h"

  ! Input/Output
  real, dimension(LOW:HIGH, MDIM), intent(in) :: bndBox
  real, dimension(MDIM), intent(out) :: newxyz

  ! aux variables
  real :: r_in, r_out, r_rand
  real :: theta_in, theta_out, theta_rand
  real :: cos_theta_in, cos_theta_out, cos_theta_rand
  real :: phi_in, phi_out, phi_rand
  real :: dx, dy, dz

  newxyz = 1.0d0 ! initialization 

  if (gr_geometry == SPHERICAL) then
    r_in  = bndBox(LOW, IAXIS) 
    r_out = bndBox(HIGH, IAXIS)
    r_rand = (r_in**3 + rand()*(r_out - r_in)**3 )**(1.0/3.0)
    newxyz(IAXIS) = r_rand

    if (NDIM >= 2) then
      theta_in  = bndBox(LOW, JAXIS)
      theta_out = bndBox(HIGH, JAXIS)
      cos_theta_in  = cos(theta_in)
      cos_theta_out = cos(theta_out)

      cos_theta_rand = cos_theta_in - rand() * (cos_theta_in - cos_theta_out)
      theta_rand = acos(cos_theta_rand)
      newxyz(JAXIS) = theta_rand

      if (NDIM <= 3) then
        phi_in  = bndBox(LOW, KAXIS)
        phi_out = bndBox(HIGH, KAXIS)
        phi_rand = phi_in + (phi_out - phi_in) * rand()
        newxyz(KAXIS) = phi_rand
      end if
    end if

  else if (gr_geometry == CARTESIAN) then

    dx = bndBox(HIGH, IAXIS) - bndBox(LOW, IAXIS)
    newxyz(IAXIS) = bndBox(LOW, IAXIS) + rand() * dx
    if (NDIM >= 2) then
      dy = bndBox(HIGH, JAXIS) - bndBox(LOW, JAXIS)
      newxyz(JAXIS) = bndBox(LOW, JAXIS) + rand() * dy 

      if (NDIM == 3) then
        dz = bndBox(HIGH, KAXIS) - bndBox(LOW, KAXIS)
        newxyz(KAXIS) = bndBox(LOW, KAXIS) + rand() * dz 
      end if
    end if

  end if

end subroutine sample_blk_position


subroutine sample_iso_velocity(velvec)
  use Simulation_data, only : clight
  use random, only : rand
  implicit none

#include "constants.h"
#include "Flash.h"

  ! Input/Output
  real, dimension(MDIM), intent(out) :: velvec

  ! aux variables
  real :: mu
  real :: theta, sintheta, costheta
  real :: phi, sinphi, cosphi

  ! Initialization 
  velvec = 1.0d0

  mu = (2.0d0 * rand()) - 1.d0

  theta = ACOS(mu)
  sintheta = SIN(theta)
  costheta = COS(theta)

  phi = 2.0d0 * PI * rand()
  sinphi = SIN(phi)
  cosphi = COS(phi)

  velvec(IAXIS) = clight * sintheta * sinphi
  velvec(JAXIS) = clight * sintheta * cosphi
  velvec(KAXIS) = clight * costheta

end subroutine sample_iso_velocity


subroutine sample_therm_face_velocity(r_hat, theta_prime, velvec)
  use Simulation_data, only : clight
  use random, only : rand
  implicit none

#include "constants.h"
#include "Flash.h"

  ! Input/Output
  real, dimension(MDIM), intent(in)  :: r_hat
  real, intent(in) :: theta_prime
  real, dimension(MDIM), intent(out) :: velvec

  ! aux variables
  real :: xi, mu
  real :: theta, sintheta, costheta
  real :: phi, sinphi, cosphi

  real, dimension(MDIM) :: z_hat, k_hat
  real, dimension(MDIM) :: vel, vel_rot
  real :: costheta_prime

  ! Initialization 
  velvec = 1.0d0

  ! Sample thermal face velocity
  xi = rand() ! Planck surface has theta in [0, 90], so mu in [0, 1]
  mu = sqrt(xi)

  theta = ACOS(mu)
  sintheta = SIN(theta)
  costheta = COS(theta)

  phi = 2.0d0 * PI * rand()
  sinphi = SIN(phi)
  cosphi = COS(phi)

  ! This velvec is in 
  vel(IAXIS) = clight * sintheta * sinphi
  vel(JAXIS) = clight * sintheta * cosphi
  vel(KAXIS) = clight * costheta

  ! Perform Rodrigues' rotation to transform the temp vel vector
  ! to the r_hat direction
  z_hat = (/ 0.0d0, 0.0d0, 1.0d0 /)
  k_hat = cross_product(z_hat, r_hat) / sin(theta_prime)

  ! Rodrigues' rotation formula
  vel_rot = vel*cos(theta_prime) + cross_product(k_hat, vel)*sin(theta_prime)&
            + k_hat*dot_product(k_hat, vel)*(1.0d0 - cos(theta_prime))

  velvec = vel_rot

end subroutine sample_therm_face_velocity


! Cartesian version of thermal face velocity sampling
subroutine sample_cart_therm_face_velocity(side, axis, velvec)
  use Simulation_data, only : clight
  use random, only : rand
  implicit none

#include "constants.h"
#include "Flash.h"

  ! Input/Output
  integer, intent(in) :: side
  integer, intent(in) :: axis
  real, dimension(MDIM), intent(out) :: velvec

  ! aux variables
  real :: v_forward, v_2, v_3
  real :: xi, mu
  real :: theta, sintheta, costheta
  real :: phi, sinphi, cosphi

  ! Initialization
  velvec = 1.0d0

  xi = rand() ! Planck surface has theta in [0, 90], so mu in [0, 1]
  mu = sqrt(xi) ! outward pointing direction should take sqrt(xi)

  theta = ACOS(mu)
  sintheta = SIN(theta)
  costheta = COS(theta)

  phi = 2.0d0 * PI * rand()
  sinphi = SIN(phi)
  cosphi = COS(phi)

  ! Setting forward and other components
  v_forward = clight * costheta
  v_2       = clight * sintheta * sinphi
  v_3       = clight * sintheta * cosphi

  ! Assigning velocity components
  velvec(axis) = v_forward
  ! The order of the non-forward direction does not matter
  if (axis .EQ. IAXIS) then
    velvec(JAXIS) = v_2
    velvec(KAXIS) = v_3
  else if (axis .EQ. JAXIS) then
    velvec(IAXIS) = v_2
    velvec(KAXIS) = v_3
  else if (axis .EQ. KAXIS) then
    velvec(IAXIS) = v_2
    velvec(JAXIS) = v_3
  else
    print *, "Error: Unknown axis for vel. sampling"
  end if

end subroutine sample_cart_therm_face_velocity


function cross_product(va, vb)
  implicit none
#include "constants.h"

  real, dimension(MDIM) :: cross_product 
  real, dimension(MDIM), intent(in) :: va, vb

  cross_product(1) = va(2)*vb(3) - va(3)*vb(2)
  cross_product(2) = va(3)*vb(1) - va(1)*vb(3)
  cross_product(3) = va(1)*vb(2) - va(2)*vb(1)

end function cross_product


function f_planck(x)
  implicit none
#include "constants.h"
  real, intent(in) :: x
  real :: f_planck

  if (x <= 20.0) then
    f_planck = (15.0 * (x**3)) / ((PI ** 4.0) * (EXP(x) - 1.0))
  else
    f_planck = 0.0d0
  end if

  return

end function f_planck


subroutine sample_time(dtfull, dtrand)
  use random, only : rand
  implicit none
  real, intent(in)  :: dtfull
  real, intent(out) :: dtrand

  dtrand = dtfull * rand()

end subroutine sample_time


subroutine sample_energy(solnData, cellID, mode, fixed_temp, eps)
  use Particles_data, only : ev2erg, pt_is_mcp_grey, pt_grey_eps,&
                             pt_energy_min_eV, pt_energy_max_eV,&
                             pt_samp_uniform, pt_samp_Tgas, pt_samp_fixedT
  use Simulation_data, only : kB
  use random, only : rand

  implicit none
#include "constants.h"

  ! Input/output
  real, pointer :: solnData(:,:,:,:)
  integer, dimension(MDIM) :: cellID
  integer, intent(in) :: mode
  real, intent(in)  :: fixed_temp
  real, intent(out) :: eps

  ! aux variables
  real :: temp, x, f_x
  real, parameter :: f_max = 0.218886


  if (pt_is_mcp_grey) then
    eps = pt_grey_eps
  else ! Thermal emission
    ! split into different sampling modes
    if (mode == pt_samp_uniform) then
      x = pt_energy_min_eV + rand() * (pt_energy_max_eV - pt_energy_min_eV)
      eps = x * ev2erg
    else ! sampling planck function
      select case(mode)
        case (pt_samp_Tgas) 
          temp = solnData(TEMP_VAR, cellID(IAXIS), cellID(JAXIS), cellID(KAXIS))
        case (pt_samp_fixedT)
          temp = fixed_temp
      end select
      do 
        x = pt_energy_min_eV + rand() * (pt_energy_max_eV - pt_energy_min_eV)
        x = (x * ev2erg) / (kB * temp)
        f_x = f_planck(x)
        if (f_max * rand() .LT. f_x) exit
      end do
      eps = (x * kB * temp) ! in erg
    end if
  end if

end subroutine sample_energy

end module new_mcp
