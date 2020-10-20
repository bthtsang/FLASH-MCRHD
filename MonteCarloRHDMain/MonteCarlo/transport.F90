module transport

  implicit none
  contains


subroutine transport_mcps(dtOld, dtNew, particles, p_count, maxcount, ind)
  use Grid_interface, only : Grid_getDeltas,&
                             Grid_outsideBoundBox, Grid_getBlkType,&
                             Grid_sortParticles, Grid_moveParticles,&
                             Grid_getBlkIDFromPos,&
                             Grid_getTileIterator, Grid_releaseTileIterator
  use Grid_data, only : gr_geometry
  use Grid_tile, only : Grid_tile_t
  use Grid_iterator, ONLY : Grid_iterator_t
  use Particles_interface, only : Particles_getGlobalNum
  use Particles_data, only : pt_typeInfo, pt_maxPerProc, pt_indexList,&
                             pt_indexCount, pt_meshMe, pt_numLocal
  use Particles_data, only : pt_is_eff_scattering, pt_max_rt_iterations,&
                             pt_is_scat_elastic, pt_is_escat_elastic,&
                             pt_is_scat_iso, pt_is_escat_iso,&
                             pt_STAY_ID, pt_SCAT_ID, pt_ESCAT_ID,&
                             pt_ABS_ID, pt_CROSS_ID, pt_smlpush,&
                             pt_PION_ESCAT_ID,&
                             pt_quadrant_I, pt_quadrant_II,&
                             pt_quadrant_III, pt_quadrant_IV,&
                             pt_is_veldp,&
                             pt_is_photoionization, pt_is_es_photoionization,&
                             pt_is_caseB, pt_is_veldp
  use Paramesh_comm_data, only : amr_mpi_meshComm
  use Driver_interface, only : Driver_abortFlash
  use spherical, only : get_cartesian_position, get_spherical_velocity,&
                        get_quadrant
  use opacity, only : calc_abs_opac, calc_sca_opac, get_fleck
  use scattering, only : scatter_mcp
  use ionization, only : calc_pi_opac, deposit_ionizing_radiation,&
                         eff_scatter_ionizing_mcp,&
                         sanitize_fully_ionized_cell
  use relativity, only : compute_dshift
  use Simulation_data, only : clight
  use gr_interface, only : gr_xyzToBlock
  use Timers_interface, only : Timers_start, Timers_stop

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Particles.h"
#include "mpif.h"

  ! Input/output 
  real, intent(in)  :: dtOld, dtNew
  integer, intent(in) :: p_count, maxcount, ind
  real,dimension(NPART_PROPS,maxcount), intent(inout) :: particles

  ! aux variables
  integer :: globalNumParticles
  integer :: nthreads, ierr
  integer :: num_sinks, num_done_local, num_done_global
  integer, dimension(:), allocatable :: num_done_list
  integer :: numpass
  real :: frac_done

  integer :: i, mcp_begin, mcp_end, tid
  integer, dimension(MAXBLOCKS,NPART_TYPES) :: particlesPerBlk

  integer :: currentBlk, currentProc
  real, dimension(MDIM) :: currentPos
  integer, dimension(MDIM) :: cellID, xcellID
  real :: dvol
  real, dimension(LOW:HIGH, MDIM) :: bndBox
  real, dimension(MDIM) :: deltaCell

  real :: dt
  real, pointer, dimension(:,:,:,:) :: solnVec => NULL()
  real :: k_a, k_s
  real :: fleck, min_time, min_dist
  real :: k_a_cmf, k_s_cmf
  real :: k_ion, k_ion_ea, k_ion_es, fleckp, nH1, N_H1
  integer :: mcp_fate
  logical :: is_empty_cell_event

  logical :: isoutside
  integer, dimension(MDIM) :: neghdir
  integer,dimension(BLKNO:PROCNO) :: neghID
  real, dimension(MDIM) :: newPos
  integer, dimension(MDIM) :: chk_cellID
  integer, dimension(MDIM) :: new_cellID
  integer :: blkType
  logical :: is_crossed
  logical :: is_crossproc

  real :: dshift

  ! Debug
  real, dimension(MDIM) :: currentVel, currentCartPos
  real, dimension(MDIM) :: or_hat, on_hat, nr_hat, nn_hat
  real, dimension(MDIM) :: sph_vel
  real, dimension(MDIM) :: cart_pos
  real :: theta_dist, phi_dist, closest_dist
  integer :: theta_dir, phi_dir
  real :: pos_x, pos_y, pos_z, pos_r, pos_t, pos_p
  integer :: quadrant
  integer :: ii
  real, dimension(LOW:HIGH, MDIM) :: new_bndBox
  real, dimension(MDIM) :: new_deltaCell
  integer :: actualBlkID, actualProcID
  real, dimension(MDIM) :: pos_for_negh, pos_after_period, pos_after_push
  real, dimension(MDIM) :: pos_before_adv, pos_after_adv
  integer, dimension(MDIM) :: xcellID_original
  integer :: crossproc_tag
  integer :: n_it_current

  ! variables for flash5 grid interface
  type(Grid_iterator_t) :: itor, newItor
  type(Grid_tile_t)    :: tileDesc, newTileDesc
  real, allocatable :: cellVolumes(:,:,:)
  real :: cellVols(1:1, 1:1, 1:1)
  integer :: icid(MDIM)
  integer, dimension(MDIM) :: local_cellID
  integer :: ansproc, ansblk

  call Timers_start("MCP Transport")

  ! Obtaining the global number of particles (MCPs and sinks)
  call Particles_getGlobalNum(globalNumParticles)

  ! Local variable for holding number of finished MCP globally
  num_done_global = 0

  ! For hybrid MPI/OpenMP applications
  !$omp parallel
  !$ nthreads = omp_get_num_threads()
  !$omp end parallel

  ! For MPI-only applications
  nthreads = 1

  allocate(num_done_list(nthreads), stat=ierr)
  ! Initialization for when a processor holds no MCP
  mcp_begin = 1
  mcp_end   = 0

  if (globalNumParticles > 0) then
    frac_done = real(num_done_global) / real(globalNumParticles)
    if (pt_meshMe .EQ. 0) call progress(frac_done)
  end if

  numpass = 0

  do while (num_done_global /= globalNumParticles)
    ! Status reporting
    !if (pt_meshMe .EQ. 0) print *, "RT:", num_done_global, globalNumParticles
    num_sinks = sum(pt_typeInfo(PART_LOCAL,:)) - pt_typeInfo(PART_LOCAL,ind)

    ! Offsetting with the number of active particles
    ! RHS is total number - photons, i.e. active particles
    num_done_list = 0
    num_done_list(1) = num_done_list(1) + num_sinks

    mcp_begin = pt_typeInfo(PART_TYPE_BEGIN,ind)
    mcp_end   = pt_typeInfo(PART_TYPE_BEGIN,ind) + pt_typeInfo(PART_LOCAL,ind) - 1

    call Timers_start("MCP - whileLoop")
    do i = mcp_begin, mcp_end

      ! Hybrid MPI/OpenMP applications
      !tid = omp_get_thread_num()
      ! MPI-only applications
      tid = 0

      n_it_current = 0
      ! Perform actual radiation transport
      do while ((particles(TREM_PART_PROP, i) .GT. 0.0) .and. (int(particles(ISCP_PART_PROP, i)) .NE. 1))
        currentBlk = int(particles(BLK_PART_PROP,i))
        currentProc = int(particles(PROC_PART_PROP,i))

        ! Get block properties
        ! HACK - this is paramesh specific
        itor%curBlk = currentBlk
        call itor%currentTile(tileDesc)
        call tileDesc%getDataPtr(solnVec, CENTER)
        icid  = tileDesc%cid
        
        currentPos = particles(POSX_PART_PROP:POSZ_PART_PROP, i)

        if (gr_geometry == SPHERICAL) then
          call get_cartesian_position(currentPos, currentCartPos)
        else if (gr_geometry == CARTESIAN) then
          currentCartPos = currentPos
        end if
        currentVel = particles(VELX_PART_PROP:VELZ_PART_PROP, i)
        or_hat = currentCartPos / sqrt(dot_product(currentCartPos, currentCartPos))
        on_hat = currentVel / clight

        call tileDesc%boundBox(bndBox)
        call tileDesc%deltas(deltaCell)

        call get_cellID(bndBox, deltaCell, currentPos, local_cellID)
        ! adding to offset for FLASH5 indexing
        cellID = local_cellID + icid - 1

        n_it_current = n_it_current + 1
        
        call Grid_getCellVolumes(tileDesc%level, &
                              cellID, cellID, &
                              cellVols)
        ! lbound(cellVols)==1 and ubound(cellVols)==1
        dvol = cellVols(1,1,1)

        ! Compute the dshift term for opacity calculation, 
        ! the particles array is not modified
        ! VEL*_VAR/VEL*_PART_PROP are vx/vy/vz regardless of 
        ! coordinate systems.
        dshift = 1.0d0
        if (pt_is_veldp) then
          call compute_dshift(cellID, solnVec, .false.,&
                              particles(:,i), dshift)
        end if

        ! Opacities are evaluated in the comoving frame
        call calc_abs_opac(cellID, solnVec, &
                particles(ENER_PART_PROP,i), k_a)
        call calc_sca_opac(cellID, solnVec, &
                particles(ENER_PART_PROP,i), k_s)
        call get_fleck(cellID, solnVec, fleck)
        k_a_cmf = k_a
        k_s_cmf = k_s

        ! Convert opacities into the lab frame for determine_fate
        k_a = dshift*k_a
        k_s = dshift*k_s

        k_ion = 0.0d0
        fleckp = 1.0d0
        nH1 = 0.0d0
#ifdef FLASH_MCRHD_IONIZATION
        fleckp = solnVec(FLEP_VAR, cellID(1), cellID(2), cellID(3))
        call calc_pi_opac(cellID, solnVec,&
                particles(ENER_PART_PROP,i), k_ion, nH1)
#endif
        N_H1 = nH1*dvol

        !if (fleckp < 0.5) then
        !  print *, "before fate", fleckp, solnVec(FLEP_VAR, cellID(1), cellID(2), cellID(3))
        !  print *, "H1_spec", solnVec(H1_SPEC, cellID(1), cellID(2), cellID(3))
        !  print *, "now H1", N_H1, nH1, dvol
        !end if
        if (fleckp < 1.0e-10) then
          print *, "fpDENS", solnVec(DENS_VAR, cellID(1), cellID(2), cellID(3))
          print *, "fpTEMP", solnVec(TEMP_VAR, cellID(1), cellID(2), cellID(3))
          print *, "fpbndBoxX", bndBox(:,1)
          print *, "fpbndBoxY", bndBox(:,2)
          print *, "fpbndBoxZ", bndBox(:,3)
        end if

        call determine_fate(bndBox, deltaCell, particles(:, i), local_cellID,&
                            k_a, k_s, fleck, k_ion, fleckp, N_H1,&
                            mcp_fate, xcellID, min_time, min_dist,&
                            is_empty_cell_event)

        !if (min_time < 1.0d-20) then
        !  print *, "smldt_pos", particles(POSX_PART_PROP:POSZ_PART_PROP,i)
        !  print *, "fate,ka,ks", mcp_fate, k_a, k_s
        !  print *, "fleck", fleck, fleckp
        !  print *, "vel*", particles(VELX_PART_PROP:VELZ_PART_PROP,i)
        !  print *, "cellID", cellID
        !  print *, "xcellID", xcellID
        !end if

        if (numpass > pt_max_rt_iterations) then
          print *, "many iterations...", mcp_fate, min_time, min_dist
          print *, "current it", numpass, i
          print *, "currentBlk", currentBlk
          print *, "currentProc", currentProc
          print *, "trem", particles(TREM_PART_PROP,i)
          print *, "tag", particles(TAG_PART_PROP,i)
          print *, "pos", particles(POSX_PART_PROP:POSZ_PART_PROP,i)
          print *, "vel", particles(VELX_PART_PROP:VELZ_PART_PROP,i)
          print *, "now cellID", cellID
          print *, "target xcellID", xcellID
          !print *, "ka/ks", k_a, k_s, fleck
          !print *, "k_ion", k_ion, fleckp
          !print *, "nH1, N_H1", nH1, N_H1
          !print *, "dvol", dvol
          !print *, "rho", solnVec(DENS_VAR, cellID(1), cellID(2), cellID(3))
          !print *, "temp", solnVec(TEMP_VAR, cellID(1), cellID(2), cellID(3))


          call get_spherical_velocity(particles(POSX_PART_PROP:POSZ_PART_PROP,i),&
                                      particles(VELX_PART_PROP:VELZ_PART_PROP,i),&
                                      sph_vel)
          print *, "sph vel", sph_vel
          !call Driver_abortFlash("Too many RT subcycles! Aborted.")
        end if

        ! Perform actions according to mcp_fate
        select case (mcp_fate)

          case (pt_STAY_ID)
            ! advance MCP to end of time step
            dt = particles(TREM_PART_PROP,i)
            !print *, "before", particles(POSX_PART_PROP:POSZ_PART_PROP, i)
            call advance_mcp(particles(:,i), dt)
            !print *, "after", particles(POSX_PART_PROP:POSZ_PART_PROP, i)

            !currentPos = particles(POSX_PART_PROP:POSZ_PART_PROP, i)
            !call get_cartesian_position(currentPos, currentCartPos)
            !currentVel = particles(VELX_PART_PROP:VELZ_PART_PROP, i)
            !nr_hat = currentCartPos / sqrt(dot_product(currentCartPos, currentCartPos))
            !nn_hat = currentVel / clight

            !print *, "fate", mcp_fate
            !print *, "or_hat", i, or_hat
            !print *, "nr_hat", i, nr_hat
            !print *, "on_hat", i, on_hat
            !print *, "nn_hat", i, nn_hat
            !print *, "******************"

            call deposit_energy_momentum(solnVec, cellID, particles(:,i),&
                                         mcp_fate, dt, dtNew,&
                                         fleck, k_a_cmf, k_s_cmf, dvol)
            call deposit_ionizing_radiation(solnVec, cellID, particles(:,i),&
                                         mcp_fate, dt, dtNew,&
                                         fleckp, k_ion, dvol,&
                                         is_empty_cell_event)

            ! Place in census
            particles(TREM_PART_PROP,i) = -1.0d0

            !currentPos = particles(POSX_PART_PROP:POSZ_PART_PROP, i)
            !call get_cartesian_position(currentPos, currentCartPos)
            !currentVel = particles(VELX_PART_PROP:VELZ_PART_PROP, i)
            !nr_hat = currentCartPos / sqrt(dot_product(currentCartPos, currentCartPos))
            !nn_hat = currentVel / clight

            !print *, "fate", mcp_fate
            !print *, "or_hat", i, or_hat
            !print *, "nr_hat", i, nr_hat
            !print *, "on_hat", i, on_hat
            !print *, "nn_hat", i, nn_hat
            !print *, "=================="

          case (pt_SCAT_ID:pt_PION_ESCAT_ID) ! Physical and effective scattering
            ! advance to scattering location
            call advance_mcp(particles(:, i), min_time)
            call deposit_energy_momentum(solnVec, cellID, particles(:,i),&
                                         mcp_fate, min_time, dtNew,&
                                         fleck, k_a_cmf, k_s_cmf, dvol)
            call deposit_ionizing_radiation(solnVec, cellID, particles(:,i),&
                                         mcp_fate, min_time, dtNew,&
                                         fleckp, k_ion, dvol,&
                                         is_empty_cell_event)
            if (mcp_fate == pt_SCAT_ID) then
              call scatter_mcp(solnVec, cellID, particles(:,i),&
                               pt_is_scat_elastic, pt_is_scat_iso)
            else if (mcp_fate == pt_ESCAT_ID) then
              call scatter_mcp(solnVec, cellID, particles(:,i),&
                               pt_is_escat_elastic, pt_is_escat_iso)
            else if (mcp_fate == pt_PION_ESCAT_ID) then
              call eff_scatter_ionizing_mcp(solnVec, cellID, particles(:,i),&
                               pt_is_caseB)
            end if

            ! Update remaining time step
            particles(TREM_PART_PROP, i) = particles(TREM_PART_PROP, i)&
                                         - min_time

          case (pt_ABS_ID) 
            ! advance MCP to absorption location
            call advance_mcp(particles(:,i), min_time)
            call deposit_energy_momentum(solnVec, cellID, particles(:,i),&
                                         mcp_fate, min_time, dtNew,&
                                         fleck, k_a_cmf, k_s_cmf, dvol)
            call deposit_ionizing_radiation(solnVec, cellID, particles(:,i),&
                                         mcp_fate, min_time, dtNew,&
                                         fleckp, k_ion, dvol,&
                                         is_empty_cell_event)

            if (.not. is_empty_cell_event) then
              ! Case where the particle is exhausted
              if (particles(NUMP_PART_PROP,i) .eq. 0.0d0) then
                particles(ISAB_PART_PROP, i) = 1.0
                particles(TREM_PART_PROP, i) = -1.0
              else 
                call Driver_abortFlash("transport_mcps: absorbed MCP does not&
                                        have zero final weight.")
              end if
            else
              ! neutral hydrogen is fully consumed, make sure
              call sanitize_fully_ionized_cell(cellID, solnVec)

              ! Update remaining time step, do not set isab to 1.
              particles(TREM_PART_PROP, i) = particles(TREM_PART_PROP, i)&
                                           - min_time
            end if


!            if (particles(NUMP_PART_PROP,i) .eq. 0.0d0) then
!              particles(ISAB_PART_PROP, i) = 1.0
!              particles(TREM_PART_PROP, i) = -1.0
!            else 
              ! maybe we should let it pass if the ionizing photons
              ! clear the neutral fraction of the current cell.
              ! dont set isab to 1, update trem by -= min_time
              ! pick up from here!


!              print *, "NUMP remaining", particles(NUMP_PART_PROP,i)
!              print *, "Event:", fleck, k_a, k_s, fleckp, k_ion, min_time
!              call Driver_abortFlash("transport_mcps: absorbed MCP does not &
!                                      have zero final weight.")
!            end if

          case (pt_CROSS_ID)
            ! advance MCP over to the next cell
            pos_before_adv = particles(POSX_PART_PROP:POSZ_PART_PROP, i)
            call advance_mcp(particles(:, i), min_time)
            pos_after_adv = particles(POSX_PART_PROP:POSZ_PART_PROP, i)
            xcellID_original = xcellID

            ! At this point, the MCP must be near a boundary
            ! Make sure to push it across the boundary
            call sanitize_boundary_mcp(local_cellID, xcellID, &
                                       bndBox, deltaCell, particles(:,i))
            pos_after_push  = particles(POSX_PART_PROP:POSZ_PART_PROP, i)

            ! Obtain the neighbor direction before BC check to preserve the
            ! correct neghdir
            newPos = particles(POSX_PART_PROP:POSZ_PART_PROP, i)
            pos_for_negh = newPos
            call Grid_outsideBoundBox(newPos, bndBox, isoutside, neghdir)

            ! Take care of the periodic BC in the phi direction when needed
            ! xcellID updated if periodic BC is triggered.
            !call sanitize_phi_periodic_BCs(cellID, xcellID, &
            !                           bndBox, deltaCell, particles(:,i))

            ! Make sure the xcellID is updated.
            ! For boundary crossing, 4/13 are corrected to 12/5.
            !call sanitize_xblock_cellID(cellID, xcellID)
            pos_after_period = particles(POSX_PART_PROP:POSZ_PART_PROP, i)

            call deposit_energy_momentum(solnVec, cellID, particles(:,i),&
                                         mcp_fate, min_time, dtNew,&
                                         fleck, k_a_cmf, k_s_cmf, dvol)
            call deposit_ionizing_radiation(solnVec, cellID, particles(:,i),&
                                         mcp_fate, min_time, dtNew,&
                                         fleckp, k_ion, dvol,&
                                         is_empty_cell_event)

            ! Get the updated cell if MCP is not moving out of block, or
            ! if isoutside is false
            newPos = particles(POSX_PART_PROP:POSZ_PART_PROP, i)
            call get_cellID(bndBox, deltaCell, newPos, new_cellID)

            !print *, "newcell", new_cellID
            crossproc_tag = int(particles(TAG_PART_PROP,i))

            is_crossproc = .false.
            if (isoutside) then
              !call Grid_getBlkIDFromPos(currentBlk, newPos, neghdir, neghID)
              ! no comm needed if using BITTREE
              call Grid_getBlkIDFromPos(newPos, neghID(BLKNO), neghID(PROCNO))

              if (neghID(PROCNO) /= currentProc) then
                ! Make sure the phi values are legal for the new_cellID 
                !call sanitize_phi_periodic_BCs(cellID, xcellID, &
                !                           bndBox, deltaCell, particles(:,i))
                particles(ISCP_PART_PROP,i) = 1.0d0
                is_crossproc = .true.
              else
                particles(BLK_PART_PROP, i) = neghID(BLKNO)
                call Grid_getBlkType(neghID(BLKNO), blkType)

                ! For checking the boundaries of the to-be new block
                newItor%curBlk = neghID(BLKNO)
                call newItor%currentTile(newTileDesc)
                call newTileDesc%boundBox(new_bndBox)
                call newTileDesc%deltas(new_deltaCell)

                ! Make sure the phi values are legal for the new_cellID 
                call sanitize_phi_periodic_BCs(cellID, xcellID, &
                                           bndBox, deltaCell, particles(:,i))
                ! Do the same check for Cartesian coordinate
                call sanitize_cart_periodic_BCs(cellID, xcellID, &
                                           bndBox, deltaCell, particles(:,i))

                ! xcellID is updated to legal values also
                call sanitize_xblock_cellID(local_cellID, xcellID)

                ! The newPos here should be free of negative phi
                newPos = particles(POSX_PART_PROP:POSZ_PART_PROP, i)
                ! Update the new_cellID for new bndBox and deltas
                call get_cellID(new_bndBox, new_deltaCell, newPos, chk_cellID)
                call get_cellID(new_bndBox, new_deltaCell, newPos, new_cellID)

                ! If the new position is still not within the new block,
                ! adjust it to be a cross-processor MCP.
                ! This could happen if a periodic BC sends a MCP farther
                ! than one block away
!                do ii = 1, NDIM
!                  if ((newPos(ii) < new_bndBox(LOW,ii)) .or.&
!                      (newPos(ii) > new_bndBox(HIGH,ii))) then
!                    print *, "StillWrongBlk", ii
!                    print *, ">cellID", cellID
!                    print *, ">xcellID", xcellID
!                    print *, ">xcellID_original", xcellID_original
!                    print *, ">chk_cellID", chk_cellID
!                    print *, ">new_cellID", new_cellID
!                    print *, "> vel", particles(VELX_PART_PROP:VELZ_PART_PROP, i)

!                    print *, ">pos_before_adv", pos_before_adv
!                    print *, ">pos_after_adv", pos_after_adv
!                    print *, ">pos_after_push", pos_after_push
!                    print *, ">pos_after_period", pos_after_period
!                    print *, ">pos_for_negh", pos_for_negh
!                    print *, ">neghdir", neghdir
!                    print *, ">WrongBCurBlk", currentBlk
!                    print *, ">WrongBlkProcID", neghID(BLKNO), neghID(PROCNO)
!                    print *, ">WrongBlkTag", particles(TAG_PART_PROP,i)
!                    print *, ">WrongBlkPos", newPos
!                    print *, ">new_bndBox(LOW:HIGH,ii)", new_bndBox(LOW,ii), new_bndBox(HIGH,ii)
!                    print *, ">bndBox(LOW:HIGH,ii)", bndBox(LOW,ii), bndBox(HIGH,ii)
!                    particles(ISCP_PART_PROP,i) = 1.0d0
!                    is_crossproc = .true.
!                  end if

!                  if ((chk_cellID(ii) <= NGUARD) .or.&
!                      (chk_cellID(ii) > NGUARD + NXB)) then
!                    print *, "StillWrongCell", ii
!                    print *, "WrongCCurBlk", currentBlk
!                    print *, "WrongCell", chk_cellID
!                    print *, "xcellID_original", xcellID_original
!                    print *, "WrongCBlkProcID", neghID(BLKNO), neghID(PROCNO)
!                    print *, "WrongCellTag", particles(TAG_PART_PROP,i)
!                    print *, "WrongCPos", newPos
!                    print *, "pos_after_adv", pos_after_adv
!                    print *, "pos_after_push", pos_after_push
!                    print *, "pos_after_period", pos_after_period
!                    print *, "pos_for_negh", pos_for_negh
!                    print *, "WrongCbnd", new_bndBox(LOW,ii), new_bndBox(HIGH,ii)
!                    print *, "PrevCbnd", bndBox(LOW,ii), bndBox(HIGH,ii)
!                    particles(ISCP_PART_PROP,i) = 1.0d0
!                    is_crossproc = .true.
!                  end if
!                end do                
              end if
            end if ! if not (isoutside)


            ! In case of negative crossing, I think it may not be
            ! set with crossproc
            !if (newPos(3) < 0.0d0) then
            !if (particles(POSZ_PART_PROP,i) < 0.0d0) then
            if ((particles(POSZ_PART_PROP,i) < 0.0d0) .and.&
                (.not. is_crossproc)) then
              print *, "Neg phi", newPos
              print *, "cellID", cellID
              print *, "xcellID", xcellID
              print *, "is_crossProc", is_crossproc
              print *, "isoutside", isoutside
              print *, "neghdir", neghdir
              print *, "neghID(BLK,PROC)", neghID(BLKNO), neghID(PROCNO)
              print *, "currentBlk,Proc", currentBlk, currentProc
            end if
            !print *, "cross_tag", particles(TAG_PART_PROP, i)
            !print *, "cross_time", min_time
            !print *, "cross_cellID", cellID
            !print *, "cross_xcellID", xcellID
            !print *, "cross_old", currentPos(3), particles(POSZ_PART_PROP, i), is_crossproc
            !print *, "cross_id", particles(PROC_PART_PROP,i), particles(BLK_PART_PROP, i)
            !print *, "cross_future", neghID(PROCNO),  neghID(BLKNO)
            !print *, "cross_bnd", bndBox(:,3)
            !print *, "cross_cell", cellID

            is_crossed = all((new_cellID == xcellID), 1)

            if (numpass > pt_max_rt_iterations) then
              print *, "CROSSB", numpass, i
              print *, "new_cellID", new_cellID
              print *, "new proc", particles(PROC_PART_PROP, i)
              print *, "new blk", particles(BLK_PART_PROP, i)
              print *, "pos_before_adv", pos_before_adv
              print *, "pos_after_adv", pos_after_adv
              print *, "pos_after_push", pos_after_push
              print *, "pos_after_period", pos_after_period
              print *, "pos_for_negh", pos_for_negh
              print *, "actualPos", particles(POSX_PART_PROP:POSZ_PART_PROP, i)
              print *, "actualVel", particles(VELX_PART_PROP:VELZ_PART_PROP, i)
              print *, "is_crossProc", is_crossproc
              print *, "isoutside", isoutside
              print *, "neghdir", neghdir
            end if
        
            if ((.NOT. is_crossed) .and. (.not. is_crossproc)) then
            !if ((.NOT. is_crossed) .and. (.not. is_crossproc)) then
              print *, "Old cell", local_cellID 
              print *, "Target cell", xcellID
              print *, "Current cell", new_cellID
              print *, "Old position", currentPos
              print *, "Current position", newPos
              print *, ">pos_before_adv", pos_before_adv
              print *, ">pos_after_adv", pos_after_adv
              print *, ">pos_after_push", pos_after_push
              print *, ">pos_after_period", pos_after_period
              print *, ">pos_for_negh", pos_for_negh
              print *, "outside?", isoutside
              print *, "neghDir", neghDir
              print *, "procID", neghID(PROCNO), currentProc, particles(PROC_PART_PROP, i)
              print *, "blkID", neghID(BLKNO), currentBlk, particles(BLK_PART_PROP, i)
              print *, "bndbox-x", bndBox(LOW,1) + deltaCell(1)*(cellID(1) - NGUARD - 1),&
                                   bndBox(LOW,1) + deltaCell(1)*(cellID(1) - NGUARD) 
              print *, "bndbox-y", bndBox(LOW,2) + deltaCell(2)*(cellID(2) - NGUARD - 1),&
                                   bndBox(LOW,2) + deltaCell(2)*(cellID(2) - NGUARD) 
              print *, "bndbox-z", bndBox(LOW,3) + deltaCell(3)*(cellID(3) - NGUARD - 1),&
                                   bndBox(LOW,3) + deltaCell(3)*(cellID(3) - NGUARD) 
              print *, "new-bndbox-z", new_bndBox(LOW,3) , new_bndBox(HIGH,3)
              print *, "deltas", deltaCell
              print *, "smlpush", pt_smlpush
              print *, "min_time", min_time
              print *, "min_dist", min_dist
              print *, "t_remain", particles(TREM_PART_PROP, i), "of", dtNew

              print *, "Re-advance MCP"
              call get_cartesian_position(currentPos, cart_pos)
              print *, "OriginalSph", currentPos
              print *, "OriginalCart", cart_pos
              print *, "Originalratios", cart_pos(2)/cart_pos(1)
              call advance_mcp_cartesian(cart_pos, currentVel, min_time)
              print *, "advancedCart", cart_pos
              print *, "advancedratios", cart_pos(2)/cart_pos(1)

              pos_x = cart_pos(IAXIS)
              pos_y = cart_pos(JAXIS)
              pos_z = cart_pos(KAXIS)

              pos_r = sqrt(pos_x*pos_x + pos_y*pos_y + pos_z*pos_z)
              pos_t = acos(pos_z/pos_r)
              pos_p = atan(pos_y/pos_x)

              call get_quadrant(pos_x, pos_y, quadrant)

              if ((quadrant == pt_quadrant_II) .or.&
                  (quadrant == pt_quadrant_III)) then
                pos_p = pos_p + PI
              else if (quadrant == pt_quadrant_IV) then
                pos_p = pos_p + 2.0d0*PI
              end if
              print *, "quadrant", quadrant
              print *, "cart2sph", pos_r, pos_t, pos_p

              call get_spherical_velocity(currentPos, currentVel, sph_vel)
              print *, "sph_vel_old", sph_vel
              call get_spherical_velocity(newPos, currentVel, sph_vel)
              print *, "sph_vel_new", sph_vel

              call calc_distance_to_theta_wall(bndBox, deltaCell,&
                          particles(:, i), cellID, theta_dist, theta_dir, .true.)
              call calc_distance_to_phi_wall(bndBox, deltaCell,&
                          particles(:, i), cellID, phi_dist, phi_dir, .true.)
              call distance_to_closest_wall(bndBox, deltaCell, particles(:, i), cellID,&
                                xcellID, closest_dist, .true.)

              call Driver_abortFlash("transport: cell boundary not crossed.")
            end if

            particles(TREM_PART_PROP, i) = particles(TREM_PART_PROP, i)&
                                         - min_time
        end select
      
        !currentPos = particles(POSX_PART_PROP:POSZ_PART_PROP, i)
        !call get_cartesian_position(currentPos, currentCartPos)
        !currentVel = particles(VELX_PART_PROP:VELZ_PART_PROP, i)
        !nr_hat = currentCartPos / sqrt(dot_product(currentCartPos, currentCartPos))
        !nn_hat = currentVel / clight

        !print *, "fate", mcp_fate
        !print *, "or_hat", i, or_hat
        !print *, "nr_hat", i, nr_hat
        !print *, "on_hat", i, on_hat
        !print *, "nn_hat", i, nr_hat
        !print *, "=================="

        call tileDesc%releaseDataPtr(solnVec, CENTER)

        !deallocate(cellVolumes)
      end do ! end of trem/iscp while loop for one MCP

      ! finished MCPs will not enter the above while loop,
      ! they get accounted for here.
      if (particles(TREM_PART_PROP, i) <= 0.0d0) then
         num_done_list(tid+1) = num_done_list(tid+1) + 1 ! avoid false sharing
      end if

    end do ! Done one pass of MCPs
    call Timers_stop("MCP - whileLoop")

    ! Update the processor and block IDs
    call Timers_start("MCP - moveParticles")
    call Grid_moveParticles(particles, NPART_PROPS, pt_maxPerProc,&
                  pt_numLocal, pt_indexList, pt_indexCount,.FALSE.)
    call Timers_stop("MCP - moveParticles")

    call Timers_start("MCP - sortParticles")
#ifdef TYPE_PART_PROP
    call Grid_sortParticles(particles, NPART_PROPS, pt_numLocal, NPART_TYPES,&
                            pt_maxPerProc, particlesPerBlk, BLK_PART_PROP,&
                            TYPE_PART_PROP)
#else
    call Grid_sortParticles(particles, NPART_PROPS, pt_numLocal, NPART_TYPES,&
                            pt_maxPerProc, particlesPerBlk, BLK_PART_PROP)
#endif
    call Timers_stop("MCP - sortParticles")
    call pt_updateTypeDS(particlesPerBlk)

    num_done_local = sum(num_done_list)

    ! resetting iscp attribute for all MCPs
    particles(ISCP_PART_PROP, mcp_begin:mcp_end) = 0.0d0

    call MPI_AllReduce(num_done_local, num_done_global, 1, MPI_INTEGER,&
                       MPI_SUM, amr_mpi_meshComm, ierr)
    call Particles_getGlobalNum(globalNumParticles)

    if (globalNumParticles > 0) then
      frac_done = real(num_done_global) / real(globalNumParticles)
    else
      frac_done = 1.0 ! Ran out of global particles

      if ((num_done_global /= globalNumParticles) .and.&
         (pt_meshMe == 0)) then
        call Driver_abortFlash("transport_mcps: globalnumparticles = 0&
                                but numbers not match.")
      end if
    end if

    if (pt_meshMe .EQ. 0) call progress(frac_done)

    numpass = numpass + 1

  end do ! while loop for all particles are done
  deallocate(num_done_list)

  if (pt_meshMe == 0) then
    print *, "RT Finished in", numpass, "subcycles."
    print *, "RT time step", dtNew
  end if

  ! Marking MCPs to be removed if marked 'absorbed'.
  do i = mcp_begin, mcp_end
    if (int(particles(ISAB_PART_PROP, i)) == 1) then
      particles(BLK_PART_PROP,i) = NONEXISTENT
    end if
  end do

  ! Actually removing the 'absorbed' MCPs
#ifdef TYPE_PART_PROP
  CALL Grid_sortParticles(particles, NPART_PROPS, pt_numLocal,&
                          NPART_TYPES, pt_maxPerProc, particlesPerBlk,&
                          BLK_PART_PROP, TYPE_PART_PROP)
#else
  CALL Grid_sortParticles(particles, NPART_PROPS, pt_numLocal,&
                          NPART_TYPES, pt_maxPerProc, particlesPerBlk,&
                          BLK_PART_PROP)
#endif
  call pt_updateTypeDS(particlesPerBlk)

  call Timers_stop("MCP Transport")

  return

end subroutine transport_mcps


subroutine determine_fate(bndBox, deltaCell, particle, cellID, k_a, k_s,&
                          fleck, k_ion, fleckp, N_H1, mcp_fate, xcellID,&
                          min_time, min_dist, is_empty_cell_event)
  use Particles_data, only : pt_STAY_ID, pt_SCAT_ID, pt_ESCAT_ID,&
                             pt_ABS_ID, pt_CROSS_ID, pt_PION_ESCAT_ID,&
                             pt_is_eff_scattering, pt_is_es_photoionization
  use Simulation_data, only : clight
  use random, only : rand
  implicit none
#include "constants.h"
#include "Flash.h"

  ! Input/output
  real, dimension(LOW:HIGH, MDIM), intent(in) :: bndBox
  real, dimension(MDIM), intent(in) :: deltaCell
  real, dimension(NPART_PROPS), intent(in) :: particle
  integer, dimension(MDIM) :: cellID
  real, intent(in) :: k_a, k_s, fleck, k_ion, fleckp, N_H1
  integer, intent(out) :: mcp_fate
  integer, dimension(MDIM), intent(out) :: xcellID
  real, intent(out) :: min_time, min_dist
  logical, intent(out) :: is_empty_cell_event

  ! aux variables
  real :: k_ion_ea, k_ion_es
  real :: k_a_sample, k_s_sample, ion_s_frac
  real :: trem, ini_weight, now_weight
  real, dimension(4) :: d_i
  real :: xi, sca_prob
  integer :: min_id
  logical :: is_empty_cell

  ! Getting MCP information 
  trem = particle(TREM_PART_PROP)
  ini_weight = particle(NPIN_PART_PROP)
  now_weight = particle(NUMP_PART_PROP)

  ! Initialization of d_i
  d_i = huge(1.0d0)
  is_empty_cell = .false.
  is_empty_cell_event = .false.

  ! distance to end of time step
  d_i(1) = trem * clight

  k_a_sample = 0.0
  k_s_sample = 0.0

  if (pt_is_eff_scattering) then
    k_a_sample = fleck * k_a
    k_s_sample = k_s + (1.0 - fleck) * k_a
  else 
    k_a_sample = k_a
    k_s_sample = k_s
  end if

  if ((pt_is_es_photoionization) .and. (fleckp /= 1.0d0)) then
    k_ion_ea   = fleckp * k_ion
    k_ion_es   = (1.0 - fleckp) * k_ion
    k_a_sample = k_a_sample + k_ion_ea
    k_s_sample = k_s_sample + k_ion_es
    ion_s_frac = 0.0
    if (k_ion_es > 0.0) ion_s_frac = k_ion_es / k_s_sample
  else
    k_a_sample = k_a_sample + k_ion ! directly adding photoionization
    ion_s_frac = 0.0
  end if

  ! Sample distances to scattering and absorption
  ! d_i = (d_stay, d_coll, d_boundary, d_absorption)
  call calc_distance_to_collision(k_s_sample, d_i(2))
  !call calc_distance_to_absorption(k_a_sample, ini_weight, now_weight, d_i(4))
  call calc_distance_to_absorption_with_ionization(k_a_sample, N_H1,&
                                      ini_weight, now_weight, d_i(4),&
                                      is_empty_cell)

  ! distance to the closest cell boundary
  call distance_to_closest_wall(bndBox, deltaCell, particle, cellID,& 
                                xcellID, d_i(3))

  ! Determine the MCP's fate
  min_dist = minval(d_i, 1)
  min_time = min_dist / clight
  min_id = minloc(d_i, 1)

  xi = 1.0
  if (min_id == 1) then
    mcp_fate = pt_STAY_ID
  else if (min_id == 2) then
    xi = rand()

    if (xi <= ion_s_frac) then
      mcp_fate = pt_PION_ESCAT_ID ! Effective ionizing scattering
    else ! physical or effective scattering 
      ! Scattering or effective scattering
      xi = rand() ! another xi for type of scattering
      if (pt_is_eff_scattering) then
        sca_prob = k_s / (k_s + (1.0d0 - fleck) * k_a)
        if (xi <= sca_prob) then
          mcp_fate = pt_SCAT_ID ! Physical scattering
        else 
          mcp_fate = pt_ESCAT_ID ! Effective scattering
        end if
      else
        mcp_fate = pt_SCAT_ID ! Physical scattering
      end if

    end if
  else if (min_id == 3) then ! Boundary crossing 
    mcp_fate = pt_CROSS_ID
  else if (min_id == 4) then
    mcp_fate = pt_ABS_ID ! Physical or effective absorption
    if (is_empty_cell) is_empty_cell_event = .true.
  end if

end subroutine determine_fate


subroutine get_cellID(bndBox, deltaCell, pos, cellID)

  implicit none
#include "Flash.h"
#include "constants.h"

  ! Input/Output
  real, dimension(LOW:HIGH, MDIM), intent(in) :: bndBox
  real, dimension(MDIM), intent(in) :: pos, deltaCell
  integer, dimension(MDIM), intent(out) :: cellID

  ! aux variables
  real :: dx_block_i, dy_block_i, dz_block_i, xp, yp, zp
  integer :: ip, jp, kp

  ! This works for both Cartesian and Spherical coord. 
  dx_block_i = 1.0/deltaCell(1)
  xp = (pos(IAXIS) - bndBox(LOW,IAXIS)) * dx_block_i
  ip = floor(xp) + 1

  if (NDIM >= 2) then
    dy_block_i = 1.0/deltaCell(2)
    yp = (pos(JAXIS) - bndBox(LOW,JAXIS)) * dy_block_i
    jp = floor(yp) + 1
  else
    jp = 1
  endif

  if (NDIM == 3) then
    dz_block_i = 1.0/deltaCell(3)
    zp = (pos(KAXIS) - bndBox(LOW,KAXIS)) * dz_block_i
    kp = floor(zp) + 1
  else
    kp = 1
  endif

  cellID = (/ ip, jp, kp /)

end subroutine get_cellID

subroutine calc_distance_to_collision(k_opac, d_coll)

  USE random, ONLY : randnozero
  implicit none
  real, intent(in)  :: k_opac
  real, intent(out) :: d_coll

  real :: xi

  xi = randnozero()
  if (xi .eq. 1.0d0) xi = 0.9999999 ! since log(1) is zero
  if (xi .eq. 0.0d0) xi = 1.0d-10

  if (k_opac .EQ. 0.0d0) then
    d_coll = huge(1.0d0)
  else
    d_coll = -log(xi) / k_opac
  end if

end subroutine calc_distance_to_collision


subroutine calc_distance_to_absorption(k_abs, ini_num, now_num, d_abs)
  use Particles_data, only : pt_abs_threshold

  implicit none
  real, intent(in)  :: k_abs, ini_num, now_num
  real, intent(out) :: d_abs

  if (k_abs .EQ. 0.0d0) then
    d_abs = huge(1.0d0)
  else
    d_abs = log(now_num / (pt_abs_threshold * ini_num)) / k_abs
  end if

end subroutine calc_distance_to_absorption


subroutine calc_distance_to_absorption_with_ionization(k_abs, N_H1,&
                                  ini_num, now_num, d_ab, is_cell_empty)
  use Particles_data, only : pt_is_photoionization, pt_abs_threshold

  implicit none
  ! Input/Output
  real, intent(in)  :: k_abs, N_H1, ini_num, now_num
  real, intent(out) :: d_ab
  logical, intent(out) :: is_cell_empty

  ! aux variable
  real :: d_exp, d_ex

  is_cell_empty = .false.
  if (k_abs == 0.0d0) then
    d_ab = huge(1.0d0)
  else
    if (.not. pt_is_photoionization) then
      d_ab = log(now_num / (pt_abs_threshold * ini_num)) / k_abs
    else 
      if (N_H1 >= now_num) then
        ! d_ab set by exp(-tau)
        d_ab = log(now_num / (pt_abs_threshold * ini_num)) / k_abs
      else
        d_exp = log(now_num / (pt_abs_threshold * ini_num)) / k_abs
        ! d_ex is set by available number of neutral hydrogen
        d_ex  = -log(1.0 - (N_H1 / now_num)) / k_abs

        if (d_exp < d_ex) then
          d_ab = d_exp
        else
          d_ab = d_ex
          is_cell_empty = .true.
        end if
      end if
    end if
  end if

end subroutine calc_distance_to_absorption_with_ionization


subroutine advance_mcp_cartesian(pos, vel, dt)
  implicit none 
#include "constants.h"
#include "Flash.h"

  real, dimension(MDIM), intent(inout) :: pos
  real, dimension(MDIM), intent(in) :: vel 
  real, intent(in) :: dt

  pos(IAXIS) = pos(IAXIS) + dt*vel(IAXIS)

  if (NDIM >= 2) then
    pos(JAXIS) = pos(JAXIS) + dt*vel(JAXIS)

    if (NDIM == 3) then
      pos(KAXIS) = pos(KAXIS) + dt*vel(KAXIS)
    end if
  end if

end subroutine advance_mcp_cartesian


subroutine advance_mcp(particle, dt)
  use Grid_data, only : gr_geometry
  use Simulation_data, only : clight
  use Driver_interface, only : Driver_abortFlash
  use spherical, only : get_cartesian_position,&
                        get_spherical_position

  implicit none
#include "constants.h"
#include "Flash.h"

  ! Input/output
  real, dimension(NPART_PROPS), intent(inout) :: particle
  real, intent(in) :: dt

  ! aux variables 
  real, dimension(MDIM) :: pos, cart_pos, r_hat, sph_pos_new
  real, dimension(MDIM) :: vel, n_hat
  real :: mu, d_traveled, r_old, r_new

  real, dimension(MDIM) :: old_cart_pos

  ! Grabbing the position and velocity attributes
  pos = particle(POSX_PART_PROP:POSZ_PART_PROP)
  vel = particle(VELX_PART_PROP:VELZ_PART_PROP)

  if (gr_geometry == CARTESIAN) then
    call advance_mcp_cartesian(pos, vel, dt)
  else if (gr_geometry == SPHERICAL) then

    r_old = pos(IAXIS) ! Spherical coord
    call get_cartesian_position(pos, cart_pos)
    old_cart_pos = cart_pos
    r_hat = cart_pos / r_old
    n_hat = vel / clight

    mu = dot_product(r_hat, n_hat) ! directional cosine

    if (NDIM <= 2) then
      d_traveled = clight * dt

      ! cosine law
      r_new = sqrt(r_old*r_old + d_traveled*d_traveled +&
                   2.0d0*r_old*d_traveled*mu)
      pos(1) = r_new
      call Driver_abortFlash("advance_mcp:&
            1D/2D spherical RT not implemented yet.")

    else if (NDIM == 3) then
      call advance_mcp_cartesian(cart_pos, vel, dt)
      call get_spherical_position(cart_pos, sph_pos_new)
     
      d_traveled = clight * dt
      r_new = sqrt(r_old*r_old + d_traveled*d_traveled +&
                   2.0d0*r_old*d_traveled*mu)
      sph_pos_new(1) = r_new
      !print *, "advmcp-old", r_hat
      !print *, "n_hat", n_hat
      if (sph_pos_new(1) == pos(1)) then
        print *, "not advancing" 
        print *, "radiuscorr", pos(1), sph_pos_new(1)
        print *, "dt", dt
        print *, "vel", vel
        print *, "oldcart", old_cart_pos
        print *, "newcart", cart_pos
        print *, "rold, rnew", r_old, r_new
        call Driver_abortFlash("advance_mcp:&
                                MCP not advanced in sph. coord.!")
      end if
      !print *, "advmcp-new", cart_pos / sqrt(dot_product(cart_pos, cart_pos))

      pos = sph_pos_new
    end if 
  end if

  ! Actually update the position
  particle(POSX_PART_PROP:POSZ_PART_PROP) = pos

end subroutine advance_mcp


subroutine sanitize_boundary_mcp(cellID, xcellID,&
                                 bndBox, deltaCell, particle)
  use Particles_data, only : pt_smlpush
  use Grid_data, only: gr_geometry
  use Grid_iterator, ONLY : Grid_iterator_t
  use Grid_tile,        ONLY : Grid_tile_t

  implicit none
#include "constants.h"
#include "Flash.h"

  ! Input/output
  integer, dimension(MDIM), intent(in) :: cellID
  integer, dimension(MDIM), intent(in) :: xcellID
  real, dimension(LOW:HIGH, MDIM), intent(in) :: bndBox
  real, dimension(MDIM), intent(in) :: deltaCell
  real, dimension(NPART_PROPS), intent(inout) :: particle

  ! aux variables
  real, dimension(MDIM) :: now_pos
  integer, dimension(MDIM) :: delta_cellID
  integer :: currentBlk, ii
  integer, dimension(2,MDIM) :: faces, onBoundary

  ! variables for flash5 grid interface
  type(Grid_iterator_t) :: itor
  type(Grid_tile_t)    :: tileDesc

  ! Gathering particle information
  ! This now_pos should be near a cell boundary
  now_pos = particle(POSX_PART_PROP:POSZ_PART_PROP)
  currentBlk = particle(BLK_PART_PROP)

  ! Check for periodic BC, for phi direction
  ! HACK - this is paramesh specific
  itor%curBlk = currentBlk
  call itor%currentTile(tileDesc)
  call tileDesc%faceBCs(faces, onBoundary)

  delta_cellID = xcellID - cellID

  ! Sanitizing
  ! If the MCP is moving up/down a cell, give it a little push
  do ii = IAXIS, KAXIS
    if (delta_cellID(ii) == 1) then
      now_pos(ii) = now_pos(ii) + deltaCell(ii)*pt_smlpush

      ! Sanitize the phi direction in case of periodic BCs.
      if (ii == KAXIS) then
        if ((gr_geometry == SPHERICAL) .and.&
            (onBoundary(HIGH, KAXIS) == PERIODIC) .and.&
            (xcellID(KAXIS) == NGUARD+NZB+1)) then
          if (now_pos(ii) < 0.5*PI) then
            now_pos(ii) = now_pos(ii) + 2.0*PI
          end if
        end if
      end if

    else if (delta_cellID(ii) == -1) then
      now_pos(ii) = now_pos(ii) - deltaCell(ii)*pt_smlpush

      ! Sanitize the phi direction in case of periodic BCs.
      if (ii == KAXIS) then
        if ((gr_geometry == SPHERICAL) .and.&
            (onBoundary(LOW, KAXIS) == PERIODIC) .and.&
            (xcellID(KAXIS) == NGUARD)) then
          if (now_pos(ii) > 1.5*PI) then
            now_pos(ii) = now_pos(ii) - 2.0*PI
          end if
        end if
      end if
    end if

  end do

  ! Update to the shifted position
  particle(POSX_PART_PROP:POSZ_PART_PROP) = now_pos

end subroutine sanitize_boundary_mcp


! Subroutine to push MCP across cell boundaries
subroutine sanitize_boundary_mcp_wrong(cellID, xcellID,&
                                 bndBox, deltaCell, particle)
  use Particles_data, only : pt_smlpush
  implicit none
#include "constants.h"
#include "Flash.h"

  ! Input/output
  integer, dimension(MDIM), intent(in) :: cellID
  integer, dimension(MDIM), intent(in) :: xcellID
  real, dimension(LOW:HIGH, MDIM), intent(in) :: bndBox
  real, dimension(MDIM), intent(in) :: deltaCell
  real, dimension(NPART_PROPS), intent(inout) :: particle

  ! aux variables
  real, dimension(MDIM) :: now_pos
  integer, dimension(MDIM) :: delta_cellID
  integer :: currentBlk, ii

  ! Gathering particle information
  ! This now_pos should be near a cell boundary
  now_pos = particle(POSX_PART_PROP:POSZ_PART_PROP)
  currentBlk = particle(BLK_PART_PROP)

  delta_cellID = xcellID - cellID

  ! Sanitizing
  ! If the MCP is moving up/down a cell, give it a little push
  do ii = IAXIS, KAXIS
    if (delta_cellID(ii) == 1) then
      now_pos(ii) = now_pos(ii) + deltaCell(ii)*pt_smlpush
    else if (delta_cellID(ii) == -1) then
      now_pos(ii) = now_pos(ii) - deltaCell(ii)*pt_smlpush

      ! Make sure the phi position is legal
      !if ((ii == KAXIS) .and. (now_pos(ii) < 0.0d0)) then
      !  now_pos(ii) = now_pos(ii) + 2.0*PI
      !end if
    end if
  end do

  ! Update to the shifted position
  particle(POSX_PART_PROP:POSZ_PART_PROP) = now_pos

end subroutine sanitize_boundary_mcp_wrong


! Subroutine to apply the periodic boundary condition 
! in the phi direction
subroutine sanitize_phi_periodic_BCs(cellID, xcellID,&
                                 bndBox, deltaCell, particle)
  use Particles_data, only : pt_smlpush
  use Grid_data, only: gr_geometry
  use Simulation_data, only : sim_zMin, sim_zMax
  use Grid_iterator, ONLY : Grid_iterator_t
  use Grid_tile,        ONLY : Grid_tile_t
  implicit none
#include "constants.h"
#include "Flash.h"

  ! Input/output
  integer, dimension(MDIM), intent(in) :: cellID
  integer, dimension(MDIM), intent(in) :: xcellID
  real, dimension(LOW:HIGH, MDIM), intent(in) :: bndBox
  real, dimension(MDIM), intent(in) :: deltaCell
  real, dimension(NPART_PROPS), intent(inout) :: particle

  ! aux variables
  real, dimension(MDIM) :: now_pos
  integer, dimension(MDIM) :: delta_cellID
  integer :: currentBlk, ii
  real:: zp
  integer, dimension(2,MDIM) :: faces, onBoundary

  ! variables for flash5 grid interface
  type(Grid_iterator_t) :: itor
  type(Grid_tile_t)    :: tileDesc

  ! Gathering particle information
  ! This now_pos should be near a cell boundary
  now_pos = particle(POSX_PART_PROP:POSZ_PART_PROP)
  currentBlk = particle(BLK_PART_PROP)

  ! Check for periodic BC, for phi direction
  ! HACK - this is paramesh specific
  itor%curBlk = currentBlk
  call itor%currentTile(tileDesc)
  call tileDesc%faceBCs(faces, onBoundary)

  delta_cellID = xcellID - cellID

  do ii = IAXIS, KAXIS
    if (delta_cellID(ii) == 1) then
      ! Sanitize the phi direction in case of periodic BCs.
      if (ii == KAXIS) then
        ! When MCP goes beyond 12th to 13th phi cell
        if ((gr_geometry == SPHERICAL) .and.&
            (onBoundary(HIGH, KAXIS) == PERIODIC) .and.&
            (xcellID(KAXIS) == NGUARD+NZB+1)) then
          !if (now_pos(ii) < 0.5*PI) then
          !  now_pos(ii) = now_pos(ii) + 2.0*PI
          !end if
          if (now_pos(ii) > sim_zMax) then
            now_pos(ii) = now_pos(ii) - sim_zMax
          end if
          !zp = (now_pos(ii) - bndBox(LOW, KAXIS)) / deltaCell(KAXIS)
          !xcellID(ii) = floor(zp) + 1 + NGUARD
        end if
      end if

    else if (delta_cellID(ii) == -1) then
      ! Sanitize the phi direction in case of periodic BCs.
      if (ii == KAXIS) then
        if ((gr_geometry == SPHERICAL) .and.&
            (onBoundary(LOW, KAXIS) == PERIODIC) .and.&
            (xcellID(KAXIS) == NGUARD)) then
          !if (now_pos(ii) > 1.5*PI) then
          !  now_pos(ii) = now_pos(ii) - 2.0*PI
          !end if
          if (now_pos(ii) < sim_zMin) then
            now_pos(ii) = now_pos(ii) + sim_zMax
          end if
          !zp = (now_pos(ii) - bndBox(LOW, KAXIS)) / deltaCell(KAXIS)
          !xcellID(ii) = floor(zp) + 1 + NGUARD
        end if
      end if
    end if
  end do

  ! Update to the shifted position
  particle(POSX_PART_PROP:POSZ_PART_PROP) = now_pos

end subroutine sanitize_phi_periodic_BCs


! Subroutine to apply the periodic boundary condition 
! in all the xyz directions in Cartesian coordinate
subroutine sanitize_cart_periodic_BCs(cellID, xcellID,&
                                 bndBox, deltaCell, particle)
  use Particles_data, only : pt_smlpush
  use Grid_data, only: gr_geometry
  use Simulation_data, only : sim_xMin, sim_xMax,&
                              sim_yMin, sim_yMax,&
                              sim_zMin, sim_zMax
  use Grid_iterator, ONLY : Grid_iterator_t
  use Grid_tile,        ONLY : Grid_tile_t
  implicit none
#include "constants.h"
#include "Flash.h"

  ! Input/output
  integer, dimension(MDIM), intent(in) :: cellID
  integer, dimension(MDIM), intent(in) :: xcellID
  real, dimension(LOW:HIGH, MDIM), intent(in) :: bndBox
  real, dimension(MDIM), intent(in) :: deltaCell
  real, dimension(NPART_PROPS), intent(inout) :: particle

  ! aux variables
  real, dimension(MDIM) :: now_pos
  integer, dimension(MDIM) :: delta_cellID
  integer :: currentBlk, ii
  real:: zp
  integer, dimension(2,MDIM) :: faces, onBoundary
  real, dimension(MDIM) :: domain_min, domain_max

  ! variables for flash5 grid interface
  type(Grid_iterator_t) :: itor
  type(Grid_tile_t)    :: tileDesc

  ! Gathering particle information
  ! This now_pos should be near a cell boundary
  now_pos = particle(POSX_PART_PROP:POSZ_PART_PROP)
  currentBlk = particle(BLK_PART_PROP)

  ! Check for periodic BC, for phi direction
  ! HACK - this is paramesh specific
  itor%curBlk = currentBlk
  call itor%currentTile(tileDesc)
  call tileDesc%faceBCs(faces, onBoundary)

  delta_cellID = xcellID - cellID

  domain_min(1) = sim_xMin
  domain_min(2) = sim_yMin
  domain_min(3) = sim_zMin
  domain_max(1) = sim_xMax
  domain_max(2) = sim_yMax
  domain_max(3) = sim_zMax

  do ii = IAXIS, KAXIS
    if (delta_cellID(ii) == 1) then
      ! For Cartesian we need to take care of x/y/z directions
      if (gr_geometry == CARTESIAN) then
        if (now_pos(ii) > domain_max(ii)) then
          now_pos(ii) = now_pos(ii) - domain_max(ii)
        end if
      end if

    else if (delta_cellID(ii) == -1) then

      ! For Cartesian we need to take care of x/y/z directions
      if (gr_geometry == CARTESIAN) then
        if (now_pos(ii) < domain_min(ii)) then
          now_pos(ii) = now_pos(ii) + domain_max(ii)
        end if
      end if

    end if
  end do

  ! Update to the shifted position
  particle(POSX_PART_PROP:POSZ_PART_PROP) = now_pos

end subroutine sanitize_cart_periodic_BCs


! Subroutine to modify xcellID to take into account
! cross-block movement of MCPs.
subroutine sanitize_xblock_cellID(cellID, xcellID)
  implicit none
#include "constants.h"
#include "Flash.h"

  ! Input/output
  integer, dimension(MDIM), intent(in) :: cellID
  integer, dimension(MDIM), intent(inout) :: xcellID

  ! aux variables
  integer, dimension(MDIM) :: delta_cellID
  integer :: ii
  integer, dimension(MDIM), parameter :: niib = (/ NXB, NYB, NZB /)

  ! change of cellID
  delta_cellID = xcellID - cellID

  do ii = IAXIS, KAXIS
    if (delta_cellID(ii) == 1) then
      ! Hitting the 9th cell
      if (xcellID(ii) == niib(ii)+1) then
        xcellID(ii) = xcellID(ii) - niib(ii)
      end if

    else if (delta_cellID(ii) == -1) then
      ! Hitting the 0th cell
      if (xcellID(ii) == 0) then
        xcellID(ii) = xcellID(ii) + niib(ii)
      end if
    end if
  end do

end subroutine sanitize_xblock_cellID



! Check whether a cellID array is valid in the NDIM directions
subroutine check_cellID(cellID, is_legal)

  implicit none
#include "Flash.h"

  ! Input/output
  integer, dimension(MDIM), intent(in) :: cellID
  logical, intent(out) :: is_legal
  ! aux variables


  is_legal = .false.




end subroutine


subroutine calc_distance_to_cartesian_wall(bndBox, deltaCell, axis,&
                                           particle, d_cart, cart_dir)
  use Simulation_data, only : clight
  use Particles_data, only : pt_smlpush
  use Driver_interface, only : Driver_abortFlash
  implicit none
#include "constants.h"
#include "Flash.h"

  ! Input/output
  real, dimension(NPART_PROPS), intent(in) :: particle
  real, dimension(LOW:HIGH, MDIM), intent(in) :: bndBox
  real, dimension(MDIM), intent(in) :: deltaCell
  integer, intent(in) :: axis
  real, intent(out) :: d_cart
  integer, intent(out):: cart_dir

  ! aux variables
  real, dimension(MDIM) :: pos, vel
  real :: dx, xp, x_inner, x_outer, min_time
  integer :: id, ii
  real, dimension(2) :: soln
  real :: pos_expected

  ! Getting MCP position
  pos = particle(POSX_PART_PROP:POSZ_PART_PROP)
  vel = particle(VELX_PART_PROP:VELZ_PART_PROP)

  dx = 1.0d0 / deltaCell(axis)
  xp = (pos(axis) - bndBox(LOW,axis)) * dx ! Offset from block edge
  id = floor(xp) 

  ! Current cell boundaries
  x_inner = bndBox(LOW,axis) + id*deltaCell(axis)
  x_outer = x_inner + deltaCell(axis)

  soln(1) = (x_inner - pos(axis)) / vel(axis) ! lower bnd
  soln(2) = (x_outer - pos(axis)) / vel(axis) ! outer bnd

  ! Select the minimum positive time 
  min_time = huge(1.0d0)
  cart_dir = 0
  do ii = 1, 2
    if ((soln(ii) > 0.0d0) .and. (soln(ii) < min_time)) then
      min_time = soln(ii)
      if (ii == 1) then
        cart_dir = -1
      else if (ii == 2) then
        cart_dir = 1
      end if
    end if
  end do

  if (min_time < 0.0d0) then
    call Driver_abortFlash("calc_dist_to_cartesian_wall: negative time resulted.")
  end if

  if ((cart_dir /= 1) .and. (cart_dir /= -1)) then
    !print *, "MCP on-boundary along axis", axis
    !print *, "pos/vel", pos(axis), vel(axis)
    !print *, "soln", soln
    !call Driver_abortFlash("calc_dist_to_cartesian_wall: unknown cart_dir.")

    ! Go through soln for zeros
    ! Offer a small fraction of cell crossing time
    print *, "Try resolving on-boundary MCP"
    do ii = 1, 2
      if (soln(ii) == 0.0) then
        min_time = pt_smlpush*deltaCell(axis)/abs(vel(axis))
        pos_expected = pos(axis) + vel(axis)*min_time
        if (vel(axis) > 0.0) then
          cart_dir = 1
        else if (vel(axis) < 0.0) then
          cart_dir = -1
        end if

        print *, "MCP on-boundary along axis", axis
        print *, "pos/vel", pos(axis), vel(axis)
        print *, "soln", soln
        print *, "Proposed cart_dir", cart_dir
        print *, "Proposed min_time", min_time
        print *, "Next expected position =", pos_expected
      end if
    end do
  end if

  ! Second check
  if ((cart_dir /= 1) .and. (cart_dir /= -1)) then
    print *, "MCP on-boundary along axis", axis
    print *, "pos/vel", pos(axis), vel(axis)
    print *, "soln", soln
    call Driver_abortFlash("calc_dist_to_cartesian_wall: unknown cart_dir.")
  end if

  d_cart = min_time * clight

end subroutine calc_distance_to_cartesian_wall



! Solving for the intersection with the next spherical wall
! The quadratic system is a_sys*x^2 + b_sys*x + c_sys = 0
subroutine calc_distance_to_spherical_wall(bndBox, deltaCell,&
                                           particle, cellID, d_r, r_dir)
  use Simulation_data, only : clight
  use Driver_interface, only : Driver_abortFlash
  use spherical, only : get_cartesian_position
  implicit none
#include "constants.h"
#include "Flash.h"

  ! Input/output
  real, dimension(NPART_PROPS), intent(in) :: particle
  integer, dimension(MDIM), intent(in) :: cellID
  real, dimension(LOW:HIGH, MDIM), intent(in) :: bndBox
  real, dimension(MDIM), intent(in) :: deltaCell
  real, intent(out) :: d_r
  integer, intent(out):: r_dir

  ! aux variables
  real, dimension(MDIM) :: sph_pos, vel, n_hat
  real, dimension(MDIM) :: cart_pos
  real :: r_in, r_out
  real :: a_sys, b_sys, c_sys
  real, dimension(4) :: soln
  integer :: ii
  real :: min_time

  ! Gathering MCP information
  sph_pos = particle(POSX_PART_PROP:POSZ_PART_PROP)
  call get_cartesian_position(sph_pos, cart_pos) 
  vel = particle(VELX_PART_PROP:VELZ_PART_PROP)
  n_hat = vel / clight
 
  ! Spherical grid information
  r_in  = bndBox(LOW, IAXIS) +&
            (cellID(IAXIS) - NGUARD - 1) * deltaCell(IAXIS)
  r_out = r_in + deltaCell(IAXIS)

  a_sys = clight*clight
  b_sys = 2.0d0 * clight * dot_product(n_hat, cart_pos)

  ! Initialization
  soln = -1.0d0

  ! First system
  c_sys = dot_product(cart_pos, cart_pos) - r_in**2
  call quadratic(a_sys, b_sys, c_sys, soln(1:2))

  ! Second system
  c_sys = dot_product(cart_pos, cart_pos) - r_out**2
  call quadratic(a_sys, b_sys, c_sys, soln(3:4))
 
  ! Select the minimum positive time 
  min_time = huge(1.0d0)
  r_dir = 0
  do ii = 1, 4
    if ((soln(ii) > 0.0d0) .and. (soln(ii) < min_time)) then
      min_time = soln(ii)
      if (ii <= 2) then
        r_dir = -1
      else 
        r_dir = 1
      end if
    end if 
  end do

  if (min_time < 0.0d0) then
    call Driver_abortFlash("calc_dist_to_spherical_wall: negative time resulted.")
  end if

  if ((r_dir /= 1) .and. (r_dir /= -1)) then
    print *, "soln", soln
    call Driver_abortFlash("calc_dist_to_spherical_wall: unknown r_dir.")
  end if

  d_r = min_time * clight

end subroutine calc_distance_to_spherical_wall


! Re-express vector in the new, flipped xyz coordinate:
! theta -> pi - theta
! phi   -> -phi
subroutine flip_vector(r_original, r_flipped)

  implicit none
#include "constants.h"
#include "Flash.h"

  ! Input/output
  real, dimension(MDIM), intent(in) :: r_original
  real, dimension(MDIM), intent(out) :: r_flipped

  ! aux variables
  real :: theta_o, phi_o, theta_f, phi_f

  ! Input vectors are in spherical coordinates
  theta_o = r_original(JAXIS)
  phi_o   = r_original(KAXIS)

  theta_f = PI - theta_o
  phi_f   = 2.0*PI - phi_o

  r_flipped(1) = r_original(1) ! copy radius
  r_flipped(2) = theta_f ! update angles
  r_flipped(3) =  phi_f

end subroutine flip_vector


subroutine setup_theta_wall_linear_system(theta, sph_pos, n_hat,&
                                          a_sys, b_sys, c_sys)
  use Simulation_data, only : clight
  use spherical, only : get_cartesian_position
  implicit none
#include "constants.h"
#include "Flash.h"

  ! Input/output
  real, intent(in) :: theta
  real, dimension(MDIM), intent(in) :: sph_pos, n_hat
  real, intent(out) :: a_sys, b_sys, c_sys

  ! aux variables
  real :: theta_up
  real, dimension(MDIM) :: sph_pos_up, n_hat_up
  real, dimension(MDIM) :: cart_pos_up


  !real :: theta_new
  !real, dimension(MDIM) :: cart_pos_new, n_hat_new
  !real, dimension(MDIM), parameter :: x_hat = (/1.0d0, 0.0d0, 0.0d0/) 
  !real, dimension(MDIM) :: old_sph_vel, new_sph_vel, new_sph_pos

  ! Flip position is it is below 
  theta_up   = theta
  sph_pos_up = sph_pos
  n_hat_up   = n_hat

  if (theta > 0.5*PI) then
    theta_up = PI - theta
    call flip_vector(sph_pos, sph_pos_up)
 
    ! flipping y and z axes
    n_hat_up(JAXIS:KAXIS) = -n_hat(JAXIS:KAXIS)
  end if

  ! Norminal case
  call get_cartesian_position(sph_pos_up, cart_pos_up)

  a_sys = clight*clight * (-n_hat_up(KAXIS)**2 + cos(theta_up)**2)
  b_sys = 2.0d0 * clight * &
             ( cos(theta_up)**2 * dot_product(n_hat_up, cart_pos_up)&
             - n_hat_up(KAXIS)*cart_pos_up(KAXIS))
  c_sys = cos(theta_up)**2 * dot_product(cart_pos_up, cart_pos_up)&
             - cart_pos_up(KAXIS)**2

  ! Mid-plane crossing
  if (theta == 0.5*PI) then
    a_sys = 0.0d0
    b_sys = n_hat(KAXIS)*clight
    c_sys = cart_pos_up(KAXIS)
  end if

  ! Dummy system for non-boundary crossing
  if ((theta == 0.0d0) .or. (theta == PI)) then
    a_sys = 1.0d0
    b_sys = 0.0d0
    c_sys = 1.0d0
  end if

  ! If theta > pi/2, modify system by 
  ! (1) flipping theta to (pi - theta)
  ! (2) rotating n_hat around x-axis by pi
!  if (theta > 0.5*PI) then
!    theta_new = PI - theta_new

    ! Rodrigues' rotation formula
!    n_hat_new    = -n_hat_new +&
!                 2.0d0*x_hat*dot_product(x_hat, n_hat_new)
!    cart_pos_new = -cart_pos_new +&
!                 2.0d0*x_hat*dot_product(x_hat, cart_pos_new)
    !n_hat_new(KAXIS) = -n_hat_new(KAXIS)
    !cart_pos_new(KAXIS) = -cart_pos_new(KAXIS)
!    call get_spherical_position(cart_pos_new, new_sph_pos)
!    call get_spherical_velocity(new_sph_pos, n_hat_new*clight, new_sph_vel)

    !print *, "old_sph_vel", old_sph_vel
    !print *, "new_sph_vel", new_sph_vel
!  end if

  !a_sys = clight*clight * (-n_hat_new(KAXIS)**2 + cos(theta_new)**2)
  !b_sys = 2.0d0 * clight * &
  !           ( cos(theta_new)**2 * dot_product(n_hat_new, cart_pos_new)&
  !           - n_hat_new(KAXIS)*cart_pos_new(KAXIS))
  !c_sys = cos(theta_new)**2 * dot_product(cart_pos_new, cart_pos_new)&
  !           - cart_pos_new(KAXIS)**2

end subroutine setup_theta_wall_linear_system


subroutine calc_distance_to_theta_wall(bndBox, deltaCell,&
                                       particle, cellID, d_t, t_dir, chk_stat)
  use Simulation_data, only : clight
  use Driver_interface, only : Driver_abortFlash
  use spherical, only : get_cartesian_position, get_spherical_velocity
  implicit none
#include "constants.h"
#include "Flash.h"

  ! Input/output
  real, dimension(NPART_PROPS), intent(in) :: particle
  integer, dimension(MDIM), intent(in) :: cellID
  real, dimension(LOW:HIGH, MDIM), intent(in) :: bndBox
  real, dimension(MDIM), intent(in) :: deltaCell
  real, intent(out) :: d_t
  integer, intent(out) :: t_dir
  logical, intent(in), optional :: chk_stat

  ! aux variables
  real, dimension(MDIM) :: sph_pos, vel, n_hat
  real, dimension(MDIM) :: cart_pos
  real :: theta_in, theta_out
  real :: a_sys, b_sys, c_sys
  real, dimension(4) :: soln
  integer :: ii
  real :: min_time
  real :: time_to_midplane

  real, dimension(MDIM) :: sph_vel

  ! Gathering MCP information
  sph_pos = particle(POSX_PART_PROP:POSZ_PART_PROP)
  call get_cartesian_position(sph_pos, cart_pos)
  vel = particle(VELX_PART_PROP:VELZ_PART_PROP)
  n_hat = vel / clight
  !print *, "b4,r_hat", cart_pos / sqrt(dot_product(cart_pos, cart_pos))
  !print *, "b4,n_hat", n_hat

  ! Spherical grid information
  theta_in  = bndBox(LOW, JAXIS) +&
                (cellID(JAXIS) - NGUARD - 1) * deltaCell(JAXIS)
  theta_out = theta_in + deltaCell(JAXIS)

  ! Initialization
  soln = -1.0d0

  call setup_theta_wall_linear_system(theta_in, sph_pos, n_hat,&
                                          a_sys, b_sys, c_sys)
  call quadratic(a_sys, b_sys, c_sys, soln(1:2))

  call setup_theta_wall_linear_system(theta_out, sph_pos, n_hat,&
                                          a_sys, b_sys, c_sys)
  call quadratic(a_sys, b_sys, c_sys, soln(3:4))

  ! Select the minimum positive time 
  min_time = huge(1.0d0)
  t_dir = 0
  do ii = 1, 4
    if ((soln(ii) > 0.0d0) .and. (soln(ii) < min_time)) then
      min_time = soln(ii)

      if (ii <= 2) then
        t_dir = -1
      else
        t_dir = 1
      end if
    end if
  end do

  if (min_time < 0.0d0) then
    call Driver_abortFlash("calc_dist_to_theta_wall: negative time resulted.")
  end if
  if ((t_dir /= 1) .and. (t_dir /= -1)) then
    !call Driver_abortFlash("calc_dist_to_theta_wall: unknown t_dir.")
    d_t = huge(1.0d0)
  else 
    d_t = min_time * clight

    !print *, "cross-theta", soln
    !print *, "ct,r_hat", cart_pos / sqrt(dot_product(cart_pos, cart_pos))
    !print *, "ct,n_hat", n_hat

    ! Check direction
    call get_spherical_velocity(sph_pos, vel, sph_vel)

    !if (((sph_vel(2) > 0) .and. (t_dir == -1)) .or. &
    !    ((sph_vel(2) < 0) .and. (t_dir == 1))) then
    !  print *, "Wrong theta-crossing"
    !  print *, "theta_in/out", theta_in, theta_out
    !  print *, "cart pos", cart_pos
    !  print *, "sph_pos", sph_pos
    !  print *, "vel", vel
    !  print *, "sph_vel", sph_vel
    !  print *, "t_dir", t_dir
    !  print *, "soln", soln
    !end if

  end if

  if (present(chk_stat)) then
    print *, "ThetaWall", soln
    print *, "ThetaVal", theta_in, theta_out, 0.5*PI
    print *, "ThetaCheck", (theta_in > 0.5*PI), (theta_out > 0.5*PI)
  end if

end subroutine calc_distance_to_theta_wall


subroutine calc_distance_to_phi_wall(bndBox, deltaCell,&
                                     particle, cellID, d_p, p_dir, chk_stat)
  use Simulation_data, only : clight
  use Driver_interface, only : Driver_abortFlash
  use spherical, only : get_cartesian_position
  implicit none
#include "constants.h"
#include "Flash.h"

  ! Input/output
  real, dimension(NPART_PROPS), intent(in) :: particle
  integer, dimension(MDIM), intent(in) :: cellID
  real, dimension(LOW:HIGH, MDIM), intent(in) :: bndBox
  real, dimension(MDIM), intent(in) :: deltaCell
  real, intent(out) :: d_p
  integer, intent(out) :: p_dir
  logical, intent(in), optional:: chk_stat

  ! aux variables
  real, dimension(MDIM) :: sph_pos, vel, n_hat
  real, dimension(MDIM) :: cart_pos
  real :: phi_in, phi_out
  real :: sin_phi_in, sin_phi_out
  real :: cos_phi_in, cos_phi_out
  real :: x0, y0
  real, dimension(2) :: soln_t
  integer :: ii
  real :: min_time

  ! Gathering MCP information
  sph_pos = particle(POSX_PART_PROP:POSZ_PART_PROP)
  call get_cartesian_position(sph_pos, cart_pos)
  vel = particle(VELX_PART_PROP:VELZ_PART_PROP)
  n_hat = vel / clight

  ! Spherical grid information
  phi_in  = bndBox(LOW, KAXIS) +&
              (cellID(KAXIS) - NGUARD - 1) * deltaCell(KAXIS)
  sin_phi_in = sin(phi_in)
  cos_phi_in = cos(phi_in)

  phi_out = phi_in + deltaCell(KAXIS)
  sin_phi_out = sin(phi_out)
  cos_phi_out = cos(phi_out)

  ! Initialization
  soln_t = -1.0d0

  x0 = cart_pos(IAXIS)
  y0 = cart_pos(JAXIS)

  ! Inner phi wall
  soln_t(1) = (x0*sin_phi_in - y0*cos_phi_in)&
              / clight&
              / (n_hat(JAXIS)*cos_phi_in - n_hat(IAXIS)*sin_phi_in)
  soln_t(2) = (x0*sin_phi_out - y0*cos_phi_out)&
              / clight&
              / (n_hat(JAXIS)*cos_phi_out - n_hat(IAXIS)*sin_phi_out)

  !soln_t(1) = (y0 - x0*sin_phi_in) &
  !            / (n_hat(IAXIS)*sin_phi_in - n_hat(JAXIS))
  !soln_t(2) = (y0 - x0*sin_phi_out) &
  !            / (n_hat(IAXIS)*sin_phi_out - n_hat(JAXIS))

  ! Select the minimum positive time 
  min_time = huge(1.0d0)
  p_dir = 0
  do ii = 1, 2
    if ((soln_t(ii) > 0.0d0) .and. (soln_t(ii) < min_time)) then
      min_time = soln_t(ii)

      if (ii == 1) then
        p_dir = -1
      else if (ii == 2) then
        p_dir = 1
      end if
    end if
  end do

  if (min_time < 0.0d0) then
    call Driver_abortFlash("calc_dist_to_phi_wall: negative time resulted.")
  end if
  if ((p_dir /= 1) .and. (p_dir /= -1)) then
    d_p = huge(1.0d0)
  else
    d_p = min_time * clight
  end if

  if (present(chk_stat)) then
    print *, "PhiWall", soln_t
  end if

end subroutine calc_distance_to_phi_wall


subroutine quadratic(a_sys, b_sys, c_sys, soln)

  implicit none
#include "constants.h"
#include "Flash.h"

  ! Input/output
  real, intent(in) :: a_sys, b_sys, c_sys
  real, dimension(2), intent(out) :: soln

  ! aux variables
  real :: Delta_sys

  ! Compute determinant
  Delta_sys = b_sys*b_sys - 4.0*a_sys*c_sys

  ! Initialization of soln array
  soln = -1.0
  
  if (a_sys == 0.0d0) then
    soln = -c_sys / b_sys
  else if (Delta_sys < 0.0d0) then
    soln = huge(1.0d0)
  else if (Delta_sys == 0.0d0) then
    soln = -b_sys / (2.0d0 * a_sys)
  else 
    soln(1) = (-b_sys + sqrt(Delta_sys)) / (2.0d0 * a_sys)
    soln(2) = (-b_sys - sqrt(Delta_sys)) / (2.0d0 * a_sys)
  end if

end subroutine quadratic 


subroutine distance_to_closest_wall(bndBox, deltaCell,&
                      particle, cellID, xcellID, min_dist, chk_stat)
  use Particles_data, only : pt_smlpush
  use Grid_data, only: gr_geometry

  implicit none
#include "constants.h"
#include "Flash.h"
#include "Particles.h"

  ! Input/output
  real, dimension(LOW:HIGH, MDIM), intent(in) :: bndBox
  real, dimension(MDIM), intent(in) :: deltaCell
  real, dimension(NPART_PROPS), intent(in) :: particle
  integer, dimension(MDIM), intent(in) :: cellID
  integer, dimension(MDIM), intent(out) :: xcellID
  real, intent(out) :: min_dist
  logical, intent(in), optional :: chk_stat

  ! aux variables
  real, dimension(MDIM) :: distances
  integer, dimension(MDIM) :: dir
  integer :: xdim, xdir

  distances = huge(1.0d0) ! Initialize to huge distances

  if (gr_geometry == SPHERICAL) then
    ! Gather distances to walls
    call calc_distance_to_spherical_wall(bndBox, deltaCell,&
                            particle, cellID, distances(IAXIS), dir(IAXIS))
    if (NDIM > 1) then
      call calc_distance_to_theta_wall(bndBox, deltaCell,&
                              particle, cellID, distances(JAXIS), dir(JAXIS))
      if (NDIM == 3) then
        call calc_distance_to_phi_wall(bndBox, deltaCell,&
                              particle, cellID, distances(KAXIS), dir(KAXIS))
      end if
    end if
  else 
    ! Gather distances to Cartesian walls
    call calc_distance_to_cartesian_wall(bndBox, deltaCell, IAXIS,&
                                         particle, distances(IAXIS), dir(IAXIS))
    if (NDIM > 1) then
      call calc_distance_to_cartesian_wall(bndBox, deltaCell, JAXIS,&
                                           particle, distances(JAXIS), dir(JAXIS))
      if (NDIM == 3) then
        call calc_distance_to_cartesian_wall(bndBox, deltaCell, KAXIS,&
                                             particle, distances(KAXIS), dir(KAXIS))
      end if
    end if
  end if

  min_dist = minval(distances, 1)
  ! Make sure MCP gets across the cell boundary
  !min_dist = min_dist * (1.0d0 + pt_smlpush)

  ! Bookkeeping 
  xdim     = minloc(distances, 1)
  xdir     = dir(xdim)
  xcellID = cellID
  xcellID(xdim) = xcellID(xdim) + xdir

  if (present(chk_stat)) then
    print *, "WallDist", distances
    print *, "WallDir", dir
    print *, "WallID", xdim, xdir
  end if

end subroutine distance_to_closest_wall


subroutine deposit_energy_momentum(solnVec, cellID, particle,&
                                   mcp_fate, dt, dtNew,&
                                   fleck, k_a, k_s, dvol)
  use Grid_data, only : gr_geometry
  use Simulation_data, only : clight
  use Particles_data, only : pt_ABS_ID, pt_is_corrdl,&
                             pt_is_deposit_urad, pt_is_deposit_energy,&
                             pt_is_deposit_momentum, &
                             pt_is_photoionization, pt_is_veldp
  use rhd, only : cellAddVar
  use relativity, only : transform_lab_to_comoving,&
                         transform_comoving_to_lab
  use Driver_interface, only : Driver_abortFlash
  implicit none 
#include "constants.h"
#include "Flash.h"

  ! Input/output
  real, pointer :: solnVec(:,:,:,:)
  integer, dimension(MDIM), intent(in) :: cellID
  real, dimension(NPART_PROPS), intent(inout) :: particle
  integer, intent(in) :: mcp_fate
  real, intent(in) :: dt, dtNew, fleck, k_a, k_s, dvol

  ! aux variables
  real :: dshift
  real :: init_weight, old_weight, new_weight
  real, dimension(MDIM) :: mcp_pos, mcp_vel
  real :: k_ea, k_es, kappa_a, kappa_s
  real :: dl, dl_corr, dtau
  real :: mcp_eps, mcp_energy
  real :: delta_e, rho, u_rad
  real :: radflux
  real :: theta, phi
  real, dimension(MDIM) :: x_hat, y_hat, z_hat
  real, dimension(MDIM) :: n_hat
  real :: a_x, a_y, a_z

  rho = solnVec(DENS_VAR, cellID(IAXIS), cellID(JAXIS), cellID(KAXIS))

  dshift = 1.0d0
  if (pt_is_veldp) then 
    call transform_lab_to_comoving(cellID, solnVec, particle, dshift)
  end if

  old_weight  = particle(NUMP_PART_PROP)  ! frame invariant
  init_weight = particle(NPIN_PART_PROP)  ! frame invariant
  mcp_eps     = particle(ENER_PART_PROP)
  mcp_pos     = particle(POSX_PART_PROP:POSZ_PART_PROP)

  ! Inputs of k_a and k_s are assumed to be in the comoving frame
  k_ea = fleck * k_a
  k_es = (1.0d0 - fleck) * k_a

  dl = clight * dt * dshift ! convert dt to comoving frame
  ! Correction for absorption during the integrated path 
  dl_corr = 1.0d0
  if (pt_is_corrdl .and. (k_ea > 0.0)) then
    dl_corr = (1.0d0 - exp(-dtau)) / dtau
  end if
  dl = dl * dl_corr

  dtau = k_ea * dl

  if (pt_is_deposit_energy) then
    new_weight = old_weight * exp(-dtau)

    if ((mcp_fate == pt_ABS_ID) .and. (.not. pt_is_photoionization)) then
      new_weight = 0.0d0
    end if

    ! no dshift needed, mcp_eps already in comoving frame
    delta_e = (mcp_eps * (old_weight - new_weight)) / (rho * dvol)
    call cellAddVar(solnVec, cellID, ABSE_VAR, delta_e)

    ! Update MCP weights
    particle(NUMP_PART_PROP) = new_weight
    particle(NUM0_PART_PROP) = new_weight
  end if

  if (pt_is_deposit_urad) then
    mcp_energy = mcp_eps * old_weight
    u_rad = mcp_energy * dl / (dvol * clight * dtNew)
    call cellAddVar(solnVec, cellID, URAD_VAR, u_rad)
  end if

  if (pt_is_deposit_momentum) then
    mcp_vel = particle(VELX_PART_PROP:VELZ_PART_PROP)
    n_hat = mcp_vel / clight

    if (gr_geometry == CARTESIAN) then
      x_hat = (/ 1.0d0, 0.0d0, 0.0d0 /)
      y_hat = (/ 0.0d0, 1.0d0, 0.0d0 /)
      z_hat = (/ 0.0d0, 0.0d0, 1.0d0 /)
    else if (gr_geometry == SPHERICAL) then
      theta = mcp_pos(JAXIS)
      phi   = mcp_pos(KAXIS)

      x_hat = (/ sin(theta)*cos(phi),&
                 sin(theta)*sin(phi),&
                 cos(theta) /)
      y_hat = (/ cos(theta)*cos(phi),&
                 cos(theta)*sin(phi),&
                -sin(theta) /)
      z_hat = (/-sin(phi),&
                 cos(phi),&
                 0.0d0 /)
    end if

    kappa_a = k_ea / rho
    kappa_s = (k_es + k_s) / rho

    radflux = mcp_energy * dl / (clight * dvol * dtNew)

    a_x = radflux * dot_product(n_hat, x_hat)
    call cellAddVar(solnVec, cellID, ABMX_VAR, a_x*kappa_a)
    call cellAddVar(solnVec, cellID, SCMX_VAR, a_x*kappa_s)

    if (NDIM >= 2) then
      a_y = radflux * dot_product(n_hat, y_hat)
      call cellAddVar(solnVec, cellID, ABMY_VAR, a_y*kappa_a)
      call cellAddVar(solnVec, cellID, SCMY_VAR, a_y*kappa_s)
 
      if (NDIM == 3) then
        a_z = radflux * dot_product(n_hat, z_hat)
        call cellAddVar(solnVec, cellID, ABMZ_VAR, a_z*kappa_a)
        call cellAddVar(solnVec, cellID, SCMZ_VAR, a_z*kappa_s)
      end if
    end if

  end if

  if (pt_is_veldp) then 
    call transform_comoving_to_lab(cellID, solnVec, particle, dshift)
  end if

end subroutine deposit_energy_momentum


! Support subroutine for showing a progress bar
subroutine progress(frac_done)
  implicit none
  real, intent(in) :: frac_done

  logical, save :: first_call = .true.
  integer, save :: last_j

  ! aux variables
  integer :: j, k
  character(len=17) :: bar

  ! Resetting bar
  bar="???% |          |"

  if (first_call) then
    last_j = -1
    first_call = .false.
  end if

  ! Convert fraction to tens of percent
  j = int(frac_done*100)

  write(unit=bar(1:3),fmt="(i3)") j

  do k = 1, j/10
    bar(6+k:6+k) = "*"
  end do

  ! print the progress bar.
  if (j /= last_j) write(unit=6,fmt="(a3,a17)") 'RT:', bar
  last_j = j

  return
end subroutine progress

end module transport
