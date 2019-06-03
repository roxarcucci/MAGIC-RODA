  !    Copyright (C) 2006 Imperial College London and others.
  !
  !    Please see the AUTHORS file in the main source directory for a full list
  !    of copyright holders.
  !
  !    Prof. C Pain
  !    Applied Modelling and Computation Group
  !    Department of Earth Science and Engineering
  !    Imperial College London
  !
  !    amcgsoftware@imperial.ac.uk
  !
  !    This library is free software; you can redistribute it and/or
  !    modify it under the terms of the GNU Lesser General Public
  !    License as published by the Free Software Foundation,
  !    version 2.1 of the License.
  !
  !    This library is distributed in the hope that it will be useful,
  !    but WITHOUT ANY WARRANTY; without even the implied warranty of
  !    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  !    Lesser General Public License for more details.
  !
  !    You should have received a copy of the GNU Lesser General Public
  !    License along with this library; if not, write to the Free Software
  !    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
  !    USA
  
#include "fdebug.h"
program rmse_coor
  use spud
  use fields
  use state_module
  use write_state_module
  use timeloop_utilities
  use global_parameters, only: option_path_len, current_time, dt, FIELD_NAME_LEN
  use FLDebug
  use snapsvd_module
  use vtk_interfaces
  use memory_diagnostics
  use populate_state_module
  use python_state
  use vtk_interfaces
  use quadrature
  use diagnostic_fields_wrapper
  use field_options
  use momentum_equation
  use solvers 
  use momentum_cg 
 ! use momentum_equation_reduced
!  use checkpoint

  implicit none
#ifdef HAVE_PETSC
#include "finclude/petsc.h"
#endif
  integer, dimension(:), allocatable :: indices,indices_tmp
  type(state_type), dimension(:), allocatable :: state_full, state_pod 
  integer :: timestep
  integer :: ierr
  integer :: total_dumps
  character(len = OPTION_PATH_LEN) :: simulation_name
  integer :: numphase
#ifdef HAVE_MPI
  call mpi_init(ierr)
#endif

#ifdef HAVE_PETSC
  call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
#endif

  call python_init()
  call read_command_line()

  
  !call calculate_rsme_coor_mp()
    call calculate_rsme_coor_standard_velocity()
    call calculate_rsme_coor_standard_pressure()
    call calculate_rsme_coor_standard_density()
   call calculate_rsme_coor_standard_solidc()
  call deallocate(state_full)
  call deallocate(state_pod)
  call deallocate_transform_cache()

  call print_references(0)
#ifdef HAVE_MEMORY_STATS
  call print_current_memory_stats(0)
#endif

#ifdef HAVE_MPI
  call mpi_finalize(ierr)
#endif

contains
  subroutine calculate_rsme_coor_mp() 
    integer :: snapshots, u_nodes, p_nodes, nsvd
    integer :: i,dump_no, d, dim, kk
    integer :: stat 
    type(vector_field), pointer :: u_full,u_pod
     type(scalar_field), pointer :: vf_full,vf_pod
     real ::  nodesum !RMSE
     integer  :: ntime
    integer dump_no_tmp, int_dump_period,total_timestep
    real current_time,finish_time,dt
    real :: nodeval_full(2),nodeval_pod(2),rmsetemp,coor(15000),fullnodesum,podnodesum,coor_xy,coor_x,coor_y , rmse(15000)
     call get_option("/timestepping/current_time", current_time)
  call get_option("/timestepping/finish_time", finish_time)       
  call get_option("/timestepping/timestep", dt)
  call get_option('/simulation_name', simulation_name)
  
  ntime=int((finish_time-current_time)/dt)
  total_timestep=int((finish_time-current_time)/dt)
      if(have_option("/io/dump_period_in_timesteps/constant")) then
       call get_option("/io/dump_period_in_timesteps/constant", int_dump_period)
    else
       int_dump_period = 1
    endif
   ! call get_option('/reduced_model/pod_basis_formation/pod_basis_count', nsvd)
    call get_option('/simulation_name',simulation_name)
   ! call read_input_states_fromvtu(state_full, state_pod) 
   
    call read_input_states (state_full, state_pod) 
      do kk =1, total_dumps!,ntime,int_dump_period
      
       ! call print_state(state_full(1))
      !  vf_full=>extract_scalar_field (state_full(kk), "phase1::PhaseVolumeFraction") ! 1 donates phase 1.
      !  vf_pod=>extract_scalar_field (state_pod(kk), "phase1::PhaseVolumeFraction") 
              vf_full=>extract_scalar_field (state_full(kk), "Pressure") ! 1 donates phase 1.
        vf_pod=>extract_scalar_field (state_pod(kk), "Pressure") 
              u_nodes=node_count(vf_full)  
                            nodesum=0
                            fullnodesum=0
                            podnodesum=0
                        do i=1,u_nodes
                            nodeval_full=node_val(vf_full,i)
                            nodeval_pod=node_val(vf_pod,i)
                            !print *,'nodeval_full(2)nodeval_pod', nodeval_full(2),nodeval_pod                                      
                            rmsetemp=(nodeval_full(2)-nodeval_pod(2))**2
                            nodesum=nodesum+rmsetemp
                            fullnodesum= fullnodesum+nodeval_full(2) !coor
                            podnodesum=podnodesum+nodeval_pod(2)!coor
                        enddo
                            rmse(kk)=SQRT(nodesum/u_nodes) 
                           ! rmseone=SQRT(nodesum/u_nodes) 
                      coor_xy=0 
                      coor_x=0
                      coor_y=0
                do i=1,u_nodes
                    nodeval_full=node_val(vf_full,i)
                    nodeval_pod=node_val(vf_pod,i)
                   coor_xy=coor_xy+(nodeval_full(2)-fullnodesum/u_nodes)*(nodeval_pod(2)-podnodesum/u_nodes)
                   coor_x=coor_x+(nodeval_full(2)-fullnodesum/u_nodes)*(nodeval_full(2)-fullnodesum/u_nodes)
                   coor_y=coor_y+(nodeval_pod(2)-podnodesum/u_nodes)*(nodeval_pod(2)-podnodesum/u_nodes)                        
                enddo
                  coor(kk)=coor_xy/(SQRT(coor_x)*SQRT(coor_y))
                   coor(kk)=abs(coor(kk))
                ! coorone=coor_xy/(SQRT(coor_x)*SQRT(coor_y))
               !  coorone=abs(coor(kk))
             !print *, 'rmse'
    enddo !do kk =1, ntime,int_dump_period
 
      open(unit=61,file='coorelation')
  do i=1, total_dumps
    ! write(61,*) ((i-1)*4+100)*2, coor(i)
     write(61,*) i, coor(i)! depends on what you want to draw the numer
   enddo
   close(61) 
  open(unit=60,file='rsme')
   do i=1, total_dumps
    ! write(60,*) ((i-1)*4+100)*2, rmse(i)
     write(60,*) i, rmse(i)
   enddo 
   close(60)
  end subroutine calculate_rsme_coor_mp
  
    subroutine calculate_rsme_coor_standard_velocity() 
    integer :: snapshots, u_nodes, p_nodes, nsvd
    integer :: i,dump_no, d, dim, kk
    integer :: stat 
    type(vector_field), pointer :: u_full,u_pod
    type(scalar_field), pointer :: vf_full,vf_pod
    real ::  nodesum, rmse_initial !RMSE
    integer  :: ntime
    integer dump_no_tmp, int_dump_period,total_timestep
    real current_time,finish_time,dt
    real :: nodeval_full(2),nodeval_pod(2),rmsetemp,coor(15000),fullnodesum,podnodesum,coor_xy,coor_x,coor_y , rmse(15000), nrmse(15000),y_max,y_min, y_sum
    call get_option("/timestepping/current_time", current_time)
    call get_option("/timestepping/finish_time", finish_time)       
    call get_option("/timestepping/timestep", dt)
    call get_option('/simulation_name', simulation_name)
  
    ntime=int((finish_time-current_time)/dt)
    total_timestep=int((finish_time-current_time)/dt)
    if(have_option("/io/dump_period_in_timesteps/constant")) then
       call get_option("/io/dump_period_in_timesteps/constant", int_dump_period)
    else
       int_dump_period = 1
    endif
    ! call get_option('/reduced_model/pod_basis_formation/pod_basis_count', nsvd)
    call get_option('/simulation_name',simulation_name)
    call read_input_states (state_full, state_pod)
    ! call put_checkpoint_vtu_to_states(state_full, state_pod)
     call print_state(state_full(1))
     print *, '--------------------------------------------------------------------------------'
     call print_state(state_pod(1))
      do kk =1, total_dumps-1!,ntime,int_dump_period
          !  stop 121
      !  vf_full=>extract_scalar_field (state_full(kk), "phase1::PhaseVolumeFraction") ! 1 donates phase 1.
      !  vf_pod=>extract_scalar_field (state_pod(kk), "phase1::PhaseVolumeFraction") 
              u_full=>extract_vector_field (state_full(kk+1), "Velocity") ! 1 donates phase 1.
              !u_pod=>extract_vector_field (state_pod(kk), "::Velocity") 
              u_pod=>extract_vector_field (state_pod(kk), "Velocity") 
              u_nodes=node_count(u_full)  
                            nodesum=0
                            fullnodesum=0
                            podnodesum=0
                        do i=1,u_nodes
                            nodeval_full=node_val(u_full,i)
                            nodeval_pod=node_val(u_pod,i)
                            !print *,'nodeval_full(2)nodeval_pod', nodeval_full(2),nodeval_pod                                      
                            rmsetemp=(nodeval_full(1)-nodeval_pod(1))**2
                            nodesum=nodesum+rmsetemp
                            fullnodesum= fullnodesum+nodeval_full(1) !coor
                            podnodesum=podnodesum+nodeval_pod(1)!coor
                        enddo
                            rmse(kk)=SQRT(nodesum/u_nodes) 
                           ! rmseone=SQRT(nodesum/u_nodes) 
                      coor_xy=0 
                      coor_x=0
                      coor_y=0
                do i=1,u_nodes
                    nodeval_full=node_val(u_full,i)
                    nodeval_pod=node_val(u_pod,i)
                   coor_xy=coor_xy+(nodeval_full(1)-fullnodesum/u_nodes)*(nodeval_pod(1)-podnodesum/u_nodes)
                   coor_x=coor_x+(nodeval_full(1)-fullnodesum/u_nodes)*(nodeval_full(1)-fullnodesum/u_nodes)
                   coor_y=coor_y+(nodeval_pod(1)-podnodesum/u_nodes)*(nodeval_pod(1)-podnodesum/u_nodes)                        
                enddo
                  coor(kk)=coor_xy/(SQRT(coor_x)*SQRT(coor_y))
                   coor(kk)=abs(coor(kk))
                ! coorone=coor_xy/(SQRT(coor_x)*SQRT(coor_y))
               !  coorone=abs(coor(kk))
             !print *, 'rmse'
    enddo !do kk =1, ntime,int_dump_period
       
  open(unit=61,file='coorelation_velocity')
  do i=1, total_dumps
    ! write(61,*) ((i-1)*4+100)*2, coor(i)
     write(61,*) i, coor(i)! depends on what you want to draw the numer
   enddo
   close(61) 
  open(unit=60,file='rmse_velocity')
   do i=1, total_dumps
    ! write(60,*) ((i-1)*4+100)*2, rmse(i)
     write(60,*) i, rmse(i)
   enddo 
   close(60)
   open(unit=62,file='nrmse_velocity')
   do i=1, total_dumps
    ! write(60,*) ((i-1)*4+100)*2, rmse(i)
     write(62,*) i, nrmse(i)
   enddo 
   close(62)

  ! open(unit=63,file='rmse_initial')
  ! write(63,*) rmse_initial  ! depends on what you want to draw the number 
  ! close(63)

 end subroutine calculate_rsme_coor_standard_velocity

  subroutine calculate_rsme_coor_standard_pressure() 
    integer :: snapshots, u_nodes, p_nodes, nsvd
    integer :: i,dump_no, d, dim, kk
    integer :: stat 
    type(vector_field), pointer :: u_full,u_pod
    type(scalar_field), pointer :: vf_full,vf_pod
    real ::  nodesum, rmse_initial !RMSE
    integer  :: ntime
    integer dump_no_tmp, int_dump_period,total_timestep
    real current_time,finish_time,dt
    real :: nodeval_full(2),nodeval_pod(2),rmsetemp,coor(15000),fullnodesum,podnodesum,coor_xy,coor_x,coor_y , rmse(15000), nrmse(15000),y_max,y_min, y_sum
    call get_option("/timestepping/current_time", current_time)
    call get_option("/timestepping/finish_time", finish_time)       
    call get_option("/timestepping/timestep", dt)
    call get_option('/simulation_name', simulation_name)
  
    ntime=int((finish_time-current_time)/dt)
    total_timestep=int((finish_time-current_time)/dt)
    if(have_option("/io/dump_period_in_timesteps/constant")) then
       call get_option("/io/dump_period_in_timesteps/constant", int_dump_period)
    else
       int_dump_period = 1
    endif
    ! call get_option('/reduced_model/pod_basis_formation/pod_basis_count', nsvd)
    call get_option('/simulation_name',simulation_name)  
    
      call read_input_states (state_full, state_pod) 
    !  call put_checkpoint_vtu_to_states(state_full, state_pod)
      do kk =1, total_dumps-1!,ntime,int_dump_period
        call print_state(state_full(1))
        call print_state(state_pod(1))
      !  vf_full=>extract_scalar_field (state_full(kk), "phase1::PhaseVolumeFraction") ! 1 donates phase 1.
      !  vf_pod=>extract_scalar_field (state_pod(kk), "phase1::PhaseVolumeFraction") 
              vf_full=>extract_scalar_field (state_full(kk+1), "Pressure") ! 1 donates phase 1.
           !   vf_pod=>extract_scalar_field (state_pod(kk), "::Pressure") 
              vf_pod=>extract_scalar_field (state_pod(kk), "Pressure") 
              u_nodes=node_count(vf_full)  
                            nodesum=0
                            fullnodesum=0
                            podnodesum=0
                            y_max=0
                            y_min=0
                        do i=1,u_nodes
                            nodeval_full=node_val(vf_full,i)
                            nodeval_pod=node_val(vf_pod,i)
                            if (nodeval_full(1)>y_max) then
                              y_max=nodeval_full(1)
                            endif
                             if (nodeval_full(1)<y_min) then
                              y_min=nodeval_full(1)
                            endif
                                               !print *,'nodeval_full(2)nodeval_pod', nodeval_full(2),nodeval_pod                                      
                            rmsetemp=(nodeval_full(1)-nodeval_pod(1))**2
                            nodesum=nodesum+rmsetemp
                            fullnodesum= fullnodesum+nodeval_full(1) !coor
                            podnodesum=podnodesum+nodeval_pod(1)!coor
                        enddo
                            rmse(kk)=SQRT(nodesum/u_nodes) 
                            nrmse(kk)=rmse(kk)/(y_max-y_min)
                           ! rmseone=SQRT(nodesum/u_nodes) 
                            rmse_initial=SQRT(nodesum)
                      coor_xy=0 
                      coor_x=0
                      coor_y=0
                do i=1,u_nodes
                    nodeval_full=node_val(vf_full,i)
                    nodeval_pod=node_val(vf_pod,i)
                   coor_xy=coor_xy+(nodeval_full(1)-fullnodesum/u_nodes)*(nodeval_pod(1)-podnodesum/u_nodes)
                   coor_x=coor_x+(nodeval_full(1)-fullnodesum/u_nodes)*(nodeval_full(1)-fullnodesum/u_nodes)
                   coor_y=coor_y+(nodeval_pod(1)-podnodesum/u_nodes)*(nodeval_pod(1)-podnodesum/u_nodes)                        
                enddo
                  coor(kk)=coor_xy/(SQRT(coor_x)*SQRT(coor_y))
                   coor(kk)=abs(coor(kk))
                ! coorone=coor_xy/(SQRT(coor_x)*SQRT(coor_y))
               !  coorone=abs(coor(kk))
             !print *, 'rmse'
    enddo !do kk =1, ntime,int_dump_period
       
  open(unit=61,file='coorelation_pressure')
  do i=1, total_dumps
    ! write(61,*) ((i-1)*4+100)*2, coor(i)
     write(61,*) i, coor(i)! depends on what you want to draw the numer
   enddo
   close(61) 
  open(unit=60,file='rmse_pressure')
   do i=1, total_dumps
    ! write(60,*) ((i-1)*4+100)*2, rmse(i)
     write(60,*) i, rmse(i)
   enddo 
   close(60)
   open(unit=62,file='nrmse_pressure')
   do i=1, total_dumps
    ! write(60,*) ((i-1)*4+100)*2, rmse(i)
     write(62,*) i, nrmse(i)
   enddo 
   close(62)

  ! open(unit=63,file='rmse_initial')
  ! write(63,*) rmse_initial  ! depends on what you want to draw the number 
  ! close(63)

 end subroutine calculate_rsme_coor_standard_pressure

 subroutine calculate_rsme_coor_standard_density() 
    integer :: snapshots, u_nodes, p_nodes, nsvd
    integer :: i,dump_no, d, dim, kk
    integer :: stat 
    type(vector_field), pointer :: u_full,u_pod
    type(scalar_field), pointer :: vf_full,vf_pod
    real ::  nodesum, rmse_initial !RMSE
    integer  :: ntime
    integer dump_no_tmp, int_dump_period,total_timestep
    real current_time,finish_time,dt
    real :: nodeval_full(2),nodeval_pod(2),rmsetemp,coor(15000),fullnodesum,podnodesum,coor_xy,coor_x,coor_y , rmse(15000), nrmse(15000),y_max,y_min, y_sum
    call get_option("/timestepping/current_time", current_time)
    call get_option("/timestepping/finish_time", finish_time)       
    call get_option("/timestepping/timestep", dt)
    call get_option('/simulation_name', simulation_name)
  
    ntime=int((finish_time-current_time)/dt)
    total_timestep=int((finish_time-current_time)/dt)
    if(have_option("/io/dump_period_in_timesteps/constant")) then
       call get_option("/io/dump_period_in_timesteps/constant", int_dump_period)
    else
       int_dump_period = 1
    endif
    ! call get_option('/reduced_model/pod_basis_formation/pod_basis_count', nsvd)
    call get_option('/simulation_name',simulation_name) 
    
     call read_input_states (state_full, state_pod)
    ! call put_checkpoint_vtu_to_states(state_full, state_pod) 
      do kk =1, total_dumps-1!,ntime,int_dump_period
      !   call print_state(state_full(1))
     ! call print_state(state_pod(1))
    ! stop 232
      !  vf_full=>extract_scalar_field (state_full(kk), "phase1::PhaseVolumeFraction") ! 1 donates phase 1.
      !  vf_pod=>extract_scalar_field (state_pod(kk), "phase1::PhaseVolumeFraction") 
            vf_full=>extract_scalar_field (state_full(kk+1), "Density") ! 1 donates phase 1.
           !   vf_pod=>extract_scalar_field (state_pod(kk), "::Density") 
            vf_pod=>extract_scalar_field (state_pod(kk), "Density") 
              u_nodes=node_count(vf_full)  
                            nodesum=0
                            fullnodesum=0
                            podnodesum=0
                            y_max=0
                            y_min=0
                        do i=1,u_nodes
                            nodeval_full=node_val(vf_full,i)
                            nodeval_pod=node_val(vf_pod,i)
                            if (nodeval_full(1)>y_max) then
                              y_max=nodeval_full(1)
                            endif
                             if (nodeval_full(1)<y_min) then
                              y_min=nodeval_full(1)
                            endif
                                               !print *,'nodeval_full(2)nodeval_pod', nodeval_full(2),nodeval_pod                                      
                            rmsetemp=(nodeval_full(1)-nodeval_pod(1))**2
                            nodesum=nodesum+rmsetemp
                            fullnodesum= fullnodesum+nodeval_full(1) !coor
                            podnodesum=podnodesum+nodeval_pod(1)!coor
                        enddo
                            rmse(kk)=SQRT(nodesum/u_nodes) 
                            nrmse(kk)=rmse(kk)/(y_max-y_min)
                           ! rmseone=SQRT(nodesum/u_nodes) 
                            rmse_initial=SQRT(nodesum)
                      coor_xy=0 
                      coor_x=0
                      coor_y=0
                do i=1,u_nodes
                    nodeval_full=node_val(vf_full,i)
                    nodeval_pod=node_val(vf_pod,i)
                   coor_xy=coor_xy+(nodeval_full(1)-fullnodesum/u_nodes)*(nodeval_pod(1)-podnodesum/u_nodes)
                   coor_x=coor_x+(nodeval_full(1)-fullnodesum/u_nodes)*(nodeval_full(1)-fullnodesum/u_nodes)
                   coor_y=coor_y+(nodeval_pod(1)-podnodesum/u_nodes)*(nodeval_pod(1)-podnodesum/u_nodes)                        
                enddo
                  coor(kk)=coor_xy/(SQRT(coor_x)*SQRT(coor_y))
                   coor(kk)=abs(coor(kk))
                ! coorone=coor_xy/(SQRT(coor_x)*SQRT(coor_y))
               !  coorone=abs(coor(kk))
             !print *, 'rmse'
    enddo !do kk =1, ntime,int_dump_period
       
  open(unit=61,file='coorelation_density')
  do i=1, total_dumps
    ! write(61,*) ((i-1)*4+100)*2, coor(i)
     write(61,*) i, coor(i)! depends on what you want to draw the numer
   enddo
   close(61) 
  open(unit=60,file='rmse_density')
   do i=1, total_dumps
    ! write(60,*) ((i-1)*4+100)*2, rmse(i)
     write(60,*) i, rmse(i)
   enddo 
   close(60)
   open(unit=62,file='nrmse_density')
   do i=1, total_dumps
    ! write(60,*) ((i-1)*4+100)*2, rmse(i)
     write(62,*) i, nrmse(i)
   enddo 
   close(62)

  ! open(unit=63,file='rmse_initial')
  ! write(63,*) rmse_initial  ! depends on what you want to draw the number 
  ! close(63)

 end subroutine calculate_rsme_coor_standard_density

  subroutine calculate_rsme_coor_standard_solidc() 
    integer :: snapshots, u_nodes, p_nodes, nsvd
    integer :: i,dump_no, d, dim, kk
    integer :: stat 
    type(vector_field), pointer :: u_full,u_pod
    type(scalar_field), pointer :: vf_full,vf_pod
    real ::  nodesum, rmse_initial !RMSE
    integer  :: ntime
    integer dump_no_tmp, int_dump_period,total_timestep
    real current_time,finish_time,dt
    real :: nodeval_full(2),nodeval_pod(2),rmsetemp,coor(15000),fullnodesum,podnodesum,coor_xy,coor_x,coor_y , rmse(15000), nrmse(15000),y_max,y_min, y_sum
    call get_option("/timestepping/current_time", current_time)
    call get_option("/timestepping/finish_time", finish_time)       
    call get_option("/timestepping/timestep", dt)
    call get_option('/simulation_name', simulation_name)
  
    ntime=int((finish_time-current_time)/dt)
    total_timestep=int((finish_time-current_time)/dt)
    if(have_option("/io/dump_period_in_timesteps/constant")) then
       call get_option("/io/dump_period_in_timesteps/constant", int_dump_period)
    else
       int_dump_period = 1
    endif
    ! call get_option('/reduced_model/pod_basis_formation/pod_basis_count', nsvd)
    call get_option('/simulation_name',simulation_name) 
    
     call read_input_states (state_full, state_pod) 
   !   call put_checkpoint_vtu_to_states(state_full, state_pod)
      do kk =1, total_dumps-1!,ntime,int_dump_period
      !  call print_state(state_full(1))
      !  vf_full=>extract_scalar_field (state_full(kk), "phase1::PhaseVolumeFraction") ! 1 donates phase 1.
      !  vf_pod=>extract_scalar_field (state_pod(kk), "phase1::PhaseVolumeFraction") 
              u_full=>extract_vector_field (state_full(kk+1), "Velocity") ! 1 donates phase 1.
              ! u_pod=>extract_vector_field (state_pod(kk), "::Velocity") !nomean-bad
              u_pod=>extract_vector_field (state_pod(kk), "Velocity") 
              u_nodes=node_count(u_full)  
                            nodesum=0
                            fullnodesum=0
                            podnodesum=0
                            y_max=0
                            y_min=0
                        do i=1,u_nodes
                            nodeval_full=node_val(u_full,i)
                            nodeval_pod=node_val(u_pod,i)
                            if (nodeval_full(2)>y_max) then
                              y_max=nodeval_full(2)
                            endif
                             if (nodeval_full(2)<y_min) then
                              y_min=nodeval_full(2)
                            endif
                                               !print *,'nodeval_full(2)nodeval_pod', nodeval_full(2),nodeval_pod                                      
                            rmsetemp=(nodeval_full(2)-nodeval_pod(2))**2
                            nodesum=nodesum+rmsetemp
                            fullnodesum= fullnodesum+nodeval_full(2) !coor
                            podnodesum=podnodesum+nodeval_pod(2)!coor
                        enddo
                            rmse(kk)=SQRT(nodesum/u_nodes) 
                            nrmse(kk)=rmse(kk)/(y_max-y_min)
                           ! rmseone=SQRT(nodesum/u_nodes) 
                            rmse_initial=SQRT(nodesum)
                      coor_xy=0 
                      coor_x=0
                      coor_y=0
                do i=1,u_nodes
                    nodeval_full=node_val(u_full,i)
                    nodeval_pod=node_val(u_pod,i)
                   coor_xy=coor_xy+(nodeval_full(2)-fullnodesum/u_nodes)*(nodeval_pod(2)-podnodesum/u_nodes)
                   coor_x=coor_x+(nodeval_full(2)-fullnodesum/u_nodes)*(nodeval_full(2)-fullnodesum/u_nodes)
                   coor_y=coor_y+(nodeval_pod(2)-podnodesum/u_nodes)*(nodeval_pod(2)-podnodesum/u_nodes)                        
                enddo
                  coor(kk)=coor_xy/(SQRT(coor_x)*SQRT(coor_y))
                   coor(kk)=abs(coor(kk))
                ! coorone=coor_xy/(SQRT(coor_x)*SQRT(coor_y))
               !  coorone=abs(coor(kk))
             !print *, 'rmse'
    enddo !do kk =1, ntime,int_dump_period
       
  open(unit=61,file='coorelation_solidc')
  do i=1, total_dumps
    ! write(61,*) ((i-1)*4+100)*2, coor(i)
     write(61,*) i, coor(i)! depends on what you want to draw the numer
   enddo
   close(61) 
  open(unit=60,file='rmse_solidc')
   do i=1, total_dumps
    ! write(60,*) ((i-1)*4+100)*2, rmse(i)
     write(60,*) i, rmse(i)
   enddo 
   close(60)
   open(unit=62,file='nrmse_solidc')
   do i=1, total_dumps
    ! write(60,*) ((i-1)*4+100)*2, rmse(i)
     write(62,*) i, nrmse(i)
   enddo 
   close(62)

  ! open(unit=63,file='rmse_initial')
  ! write(63,*) rmse_initial  ! depends on what you want to draw the number 
  ! close(63)

 end subroutine calculate_rsme_coor_standard_solidc

subroutine read_input_states(state_full, state_pod)
 
    type(state_type), intent(out), dimension(:), allocatable :: state_full,state_pod
    character(len=1024) :: filename

    integer :: dump_sampling_period, quadrature_degree
    integer :: i,j,k,stable_dumps
    real :: pr

    call get_option('/reduced_model/pod_basis_formation/dump_sampling_period',dump_sampling_period)
    call get_option('/geometry/quadrature/degree', quadrature_degree)

    ewrite(3,*) 'dump_sampling_period',dump_sampling_period
 
    total_dumps=count_dumps(dump_sampling_period)-1
    allocate(state_full(total_dumps)) 
    allocate(state_pod(total_dumps)) 

    do i=1, total_dumps/2  
      
       write(filename, '(a, i0, a)') trim(simulation_name)//'_', (i)*dump_sampling_period,".vtu"   
       call vtk_read_state(filename, state_full(i), quadrature_degree) 
     !  write(filename, '(a, i0, a)') trim(simulation_name)//'_POD_', (i)*dump_sampling_period,".vtu"  !coupling
       ! write(filename, '(a, i0, a)') trim(simulation_name)//'_POD_final_', (i)*dump_sampling_period,".vtu"   !multi
      write(filename, '(a, i0, a)') trim(simulation_name)//'_POD_', (i)*dump_sampling_period*2,".vtu" ! 
       call vtk_read_state(filename, state_pod(i), quadrature_degree) 
     end do
  ! do i=1, 332
   !  pr=100-(real(100*i-50)/332)
   !  print *, pr
  ! enddo

  end subroutine read_input_states

 

 function count_dumps(dump_sampling_period) result (count)
    !! Work out how many dumps we're going to read in.
    integer :: count,dump_sampling_period

    logical :: exists
    !      character(len=FILE_NAME_LEN) :: filename
    character(len=1024) :: filename
    dump_sampling_period=1
    count=1

    do 
       !! Note that this won't work in parallel. Have to look for the pvtu in that case.
   
       write(filename, '(a, i0, a)') trim(simulation_name)//'_', (count)*dump_sampling_period,".vtu" 
       inquire(file=trim(filename), exist=exists)
       if (.not. exists) then
          if(dump_sampling_period.eq.1) then
             count=count-2
          else
             count=count-1
          endif
          exit
       end if

       count=count+1
    end do

    if (count==0) then
       FLExit("No .vtu files found!")
    end if
    count=count-1
  end function count_dumps

  subroutine read_command_line()
    implicit none
    ! Read the input filename.
    character(len=1024) :: argument, filename
    integer :: status, argn, level

    call set_global_debug_level(0)

    argn=1
    do 

       call get_command_argument(argn, value=argument, status=status)
       argn=argn+1

       if (status/=0) then
          call usage
          stop
       end if

       if (argument=="-v") then
          call get_command_argument(argn, value=argument, status=status)
          argn=argn+1

          if (status/=0) then
             call usage
             stop
          end if

          read(argument, "(i1)", err=666) level
          call set_global_debug_level(level)

       else

          ! Must be the filename
          filename=argument

       end if

       if (argn>=command_argument_count()) exit
    end do

    call load_options(filename)

    return

666 call usage
    stop

  end subroutine read_command_line

  subroutine usage
    implicit none

    write (0,*) "usage: form_pod_basis [-v n] <options_file>"
    write (0,*) ""
    write (0,*) "-v n sets the verbosity of debugging"
  end subroutine usage 

end program rmse_coor
