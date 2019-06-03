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
program non_intrusive_high_materials
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
  !use rbf_interp
  !use reduced_model_runtime
  !use momentum_equation_reduced
  !USE IFPORT
  !Use DFLib
  !use checkpoint





  implicit none
#ifdef HAVE_PETSC
!#include "finclude/petsc.h"
#endif
  integer, dimension(:), allocatable :: indices,indices_tmp
  type(state_type), dimension(:,:,:), allocatable :: state ,state_p  
  integer :: timestep, dimen
  integer :: ierr
  integer :: numphase,total_dumps
  character(len = OPTION_PATH_LEN) :: simulation_name
  integer :: no_smolyak_nodes, m, sparse_max,snapshots, u_nodes, statei,p_nodes,v_nodes
  type(state_type), dimension(:), pointer :: state_show=> null()
  type(state_type), dimension(:,:,:), allocatable :: pod_state 
  type(state_type), dimension(:,:,:), allocatable :: pod_state_p, pod_state_v
  real, dimension(:,:,:,:), allocatable ::  interpolate_velocity
    real, dimension(:,:,:,:), allocatable :: snapmatrix_pressure, snapmatrix_volumefraction
    real, dimension(:,:,:), allocatable ::  interpolate_pressure,interpolate_volumefraction
    real, dimension(:,:,:,:), allocatable :: leftsvd_velocity 
    real, dimension(:,:,:), allocatable :: leftsvd_pressure,  leftsvd_volumefraction
    real, dimension(:,:,:), allocatable :: svdval_velocity 
    real, dimension(:,:), allocatable :: svdval_pressure, svdval_pressure_deim,svdval_volumefraction
    real, dimension(:,:,:), allocatable :: snapmean_velocity  
    real, dimension(:,:), allocatable :: snapmean_pressure, snapmean_volumefraction 
    real, dimension(:,:,:), allocatable :: pod_coef_all_obv
    integer :: nsvd
    real, dimension(:,:), allocatable :: input_parameters 
    integer :: j,k
#ifdef HAVE_MPI
  call mpi_init(ierr)
#endif

#ifdef HAVE_PETSC
 ! call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
#endif

  call python_init()
  call read_command_line()
  !m = 4 
  !
 ! call get_option( '/reduced_model/variables_to_interpolate/dimension_number', m)
 ! call get_option(&            '/reduced_model/smolyak_order', sparse_max) 

  if (m .eq. 4 .and. sparse_max .eq. 1) then
     no_smolyak_nodes=9!41
  else if (m .eq. 4 .and. sparse_max .eq. 2) then   
     no_smolyak_nodes=41
  else if (m .eq. 8 .and. sparse_max .eq. 1) then 
     no_smolyak_nodes=17
   endif
 ! sparse_max =1
  

  
  
 ! call interpolate_snapshots_main() 
  call interpolate_channel_22d_main()
  do statei=1, size(state,1)
  call deallocate(state(statei,:,:))
  call deallocate(state_p(statei,:,:))
  enddo
    deallocate(pod_coef_all_obv)
    deallocate(snapmean_velocity)
    deallocate(snapmean_pressure)  
    deallocate(snapmean_volumefraction )  
 
  call deallocate_transform_cache()

  !call print_references(0)
#ifdef HAVE_MEMORY_STATS
  call print_current_memory_stats(0)
#endif

#ifdef HAVE_MPI
  call mpi_finalize(ierr)
#endif

contains 
   subroutine interpolate_channel_22d_main()
    !!< Matrices containing the snapshots for arpack
  
    real, dimension(:,:,:,:,:), allocatable :: snapmatrix_velocity
    !real, dimension(:,:,:,:), allocatable ::  interpolate_velocity
    !real, dimension(:,:,:,:), allocatable :: snapmatrix_pressure, snapmatrix_volumefraction
    !real, dimension(:,:,:), allocatable ::  interpolate_pressure,interpolate_volumefraction
    
    type(vector_field), pointer :: velocity_field
    type(scalar_field), pointer :: pressure_field, volumefraction_field
    integer :: i, dump_no, d, j,k , i2,i3,ph,p
    integer :: stat    
    
  
    integer :: udim,deim_number_tmp  
    type(vector_field), pointer :: velocityudim 
    type(vector_field), pointer :: vfield
    character(len=1024) :: phase_name
    character(len = OPTION_PATH_LEN) :: simulation_name1
     call get_option('/reduced_model/pod_basis_formation/pod_basis_count', nsvd)
    call get_option('/simulation_name',simulation_name)
    !nsvd =12
    numphase=option_count("/material_phase")

    m=22
    sparse_max=1
    no_smolyak_nodes=23
   allocate(input_parameters(no_smolyak_nodes, m))
 
   open(unit=61, file="input_permeability")   
     read(61,*)((input_parameters(j,k),k=1,m), j=1,no_smolyak_nodes)  
    close(61) 

    call initialise_write_state
    ! Read state from .flml file
   
    ! call populate_multi_state(state_show)
    ! call print_state(state_show(1))
    
     
    call read_input_states(state,state_p)  
     print *, ' retrieve_nodesvalue_from_states_22d begin'
    call retrieve_nodesvalue_from_states_22d(state,  snapmatrix_velocity, snapmatrix_pressure, snapmatrix_volumefraction)
     print *, ' interpolation begin'
     
      call interpolate_snapmatrix_velocity_rbf_22d( snapmatrix_velocity, interpolate_velocity) !allocate(snapmatrix_velocity(no_smolyak_nodes,u_nodes,snapshots,dim,numphase))
     call interpolate_snapmatrix_pressure_rbf_22d( snapmatrix_pressure, interpolate_pressure) !allocate(snapmatrix_pressure(no_smolyak_nodes,p_nodes,snapshots,numphase))
      call interpolate_snapmatrix_volumefraction_rbf_22d( snapmatrix_volumefraction, interpolate_volumefraction)
    print *, ' interpolation over'
     
     
    !   visiualize the vtus 
    !   print *, 'size: 2, 1, 1152, 1', size(snapmean_pressure%val), size(pressure_field%val )  
    !   pressure_field%val=snapmean_pressure%val  
    !   volumefraction_field%val(:)=volumefraction_field%val(:) 
    !   volumefraction_field%val(1,2,:)=1-volumefraction_field%val(1,1,:) 
    if(.false.) then
    do i=1, snapshots
     do ph=1, 1!numphase
     pressure_field => extract_scalar_field( state_show(ph), "Pressure" )
     volumefraction_field=> extract_scalar_field(state_show(ph), "PhaseVolumeFraction" )
     velocity_field=>extract_vector_field(state_show(ph),"Velocity")
    !  print *, 'velocity_field%val', size(velocity_field%val,1), size(velocity_field%val,2)
      velocity_field%val(1,:)=interpolate_velocity(:,i,1,ph)
      velocity_field%val(2,:)=interpolate_velocity(:,i,2,ph)
      call insert(state_show(ph), velocity_field, "Velocity")
     print *, 'size(pressure_field%val', size(pressure_field%val), size(interpolate_pressure,1)
     
     pressure_field%val=interpolate_pressure(:,i,ph)
     
     call insert(state_show(ph), pressure_field, "pressure")
     volumefraction_field%val=interpolate_volumefraction(:,i,ph)
     call insert(state_show(ph), volumefraction_field, "PhaseVolumeFraction")
     dump_no=i
     call set_option('/simulation_name', trim(simulation_name)//'_interp') 
     call write_state(dump_no, state_show)
    enddo
    enddo
   endif

   if(.true.) then
     do i=1, snapshots
     do ph=1, 1!numphase
     	pressure_field => extract_scalar_field(state(no_smolyak_nodes,i,ph), "Pressure" )
    	volumefraction_field => extract_scalar_field(state(no_smolyak_nodes,i,ph), "PhaseVolumeFraction" )
     	velocity_field=>extract_vector_field(state(no_smolyak_nodes,i,ph), "phase"//trim(int2str(ph))//"::Velocity")
    	!  print *, 'velocity_field%val', size(velocity_field%val,1), size(velocity_field%val,2)
      	velocity_field%val(1,:)=interpolate_velocity(:,i,1,ph)
      	velocity_field%val(2,:)=interpolate_velocity(:,i,2,ph)
      	call insert(state(no_smolyak_nodes,i,ph), velocity_field, "phase"//trim(int2str(ph))//"::Velocity")
     	print *, 'size(pressure_field%val', size(pressure_field%val), size(interpolate_pressure,1)
     
    	 pressure_field%val=interpolate_pressure(:,i,ph)
     
    	 call insert(state(no_smolyak_nodes,i,ph), pressure_field, "pressure")
    	 volumefraction_field%val=interpolate_volumefraction(:,i,ph)
     	call insert(state(no_smolyak_nodes,i,ph), volumefraction_field, "PhaseVolumeFraction")
     	dump_no=i
     	call set_option('/simulation_name', trim(simulation_name)//'_interp') 
     	!call write_state(dump_no, state(no_smolyak_nodes,i,ph))
     enddo
     enddo
       print *, 'write_state(dump_no, state(no_smolyak_nodes,:,1))' 
      !call write_state(dump_no, state(no_smolyak_nodes,:,1))

     do ph=1, numphase
           call get_option("/material_phase["//int2str(ph-1)//"]/name", phase_name)
     do i=1,snapshots 
       dump_no=i 
      call vtk_write_state(filename=trim(simulation_name)//"_interp"//trim(phase_name)//"Velocity", index=dump_no, model="Mesh", state=state(no_smolyak_nodes,:,ph))
      call vtk_write_state(filename=trim(simulation_name)//"_interp"//trim(phase_name)//"Pressure", index=dump_no, model="PressureMesh", state=state(no_smolyak_nodes,:,ph))
      call vtk_write_state(filename=trim(simulation_name)//"_interp"//trim(phase_name)//"Volumefraction", index=dump_no, model="PressureMesh", state=state(no_smolyak_nodes,:,ph))
             
    enddo
    enddo !numphase       
   endif 
    !!!!!!!form podbasis for interpolated snapshots
    allocate(snapmean_velocity(u_nodes,dimen,numphase))
    allocate(snapmean_pressure(p_nodes,numphase))  
    allocate(snapmean_volumefraction(v_nodes,numphase))  
     snapmean_velocity=0
     snapmean_pressure=0
    snapmean_volumefraction=0 
     do p = 1, numphase
     do i=1, snapshots
       do d=1, dimen
          snapmean_velocity(:,d,p)= snapmean_velocity(:,d,p)+interpolate_velocity(:,i,d,p)
       enddo
       snapmean_pressure(:,p)=snapmean_pressure(:,p)+interpolate_pressure(:,i,p)
       snapmean_volumefraction(:,p)=snapmean_volumefraction(:,p)+interpolate_volumefraction(:,i,p)
    end do

    do d=1,dimen
       snapmean_velocity(:,d,p)=snapmean_velocity(:,d,p)/snapshots
    enddo
       snapmean_pressure(:,p)=snapmean_pressure(:,p)/snapshots
       snapmean_volumefraction(:,p)=snapmean_volumefraction(:,p)/snapshots
       snapmean_velocity=0
     snapmean_pressure=0
    snapmean_volumefraction=0
    do i=1,snapshots
       do d=1,dimen
          interpolate_velocity(:,i,d,p)=interpolate_velocity(:,i,d,p)-snapmean_velocity(:,d,p)
       enddo
       interpolate_pressure(:,i,p)=interpolate_pressure(:,i,p)-snapmean_pressure(:,p)
       interpolate_volumefraction(:,i,p)=interpolate_volumefraction(:,i,p)-snapmean_volumefraction(:,p)
    enddo
    enddo !end numphase  
    print *, 'begin forming basis functions'
    call form_svd_mp(interpolate_velocity, interpolate_pressure, interpolate_volumefraction,&
         & leftsvd_velocity, leftsvd_pressure, leftsvd_volumefraction, svdval_velocity, svdval_pressure, svdval_volumefraction, snapshots)
    print *, 'begin forming podstate'    
    call form_podstate_velocity_22d(state(1,:,:),pod_state,leftsvd_velocity,leftsvd_pressure, snapmean_velocity, snapmean_pressure)
    call form_podstate_pressure(state_p(1,:,:),pod_state_p,leftsvd_velocity,leftsvd_pressure, snapmean_velocity, snapmean_pressure)
    call form_podstate_volumefraction(state_p(1,:,:),pod_state_v,leftsvd_volumefraction, snapmean_volumefraction)
    ! enddo
    call add_option("/reduced_model/execute_reduced_model",stat)
    call get_option('/simulation_name',simulation_name)
    do ph=1, numphase
           call get_option("/material_phase["//int2str(ph-1)//"]/name", phase_name)
    do i=1,nsvd 
        dump_no=i 
      call vtk_write_state(filename=trim(simulation_name)//"_"//trim(phase_name)//"_PODBasisVelocity", index=dump_no, model="VelocityMesh", state=pod_state(i,:,ph))
      call vtk_write_state(filename=trim(simulation_name)//"_"//trim(phase_name)//"_PODBasisPressure", index=dump_no, model="PressureMesh", state=pod_state_p(i,:,ph))
      call vtk_write_state(filename=trim(simulation_name)//"_"//trim(phase_name)//"_PODBasisVolumefraction", index=dump_no, model="VolumefractionMesh", state=pod_state_v(i,:,ph))
             
    enddo
    enddo !numphase 
  
     call set_option('/simulation_name', trim(simulation_name)//'_POD') 
     call get_option('/simulation_name',simulation_name) 
     call write_options(trim(simulation_name)//".mpml")

      call calculate_pod_coef_all_timesteps()
     call non_intrusive_22d()

    deallocate(snapmatrix_velocity)
    deallocate(snapmatrix_pressure) 
    deallocate(snapmatrix_volumefraction) 
    deallocate(interpolate_pressure)
    deallocate(interpolate_velocity)
    deallocate(interpolate_volumefraction)
       deallocate(input_parameters)
    
  end subroutine interpolate_channel_22d_main

  subroutine interpolate_snapshots_main()
    !!< Matrices containing the snapshots for arpack
  
    real, dimension(:,:,:,:,:), allocatable :: snapmatrix_velocity
    !real, dimension(:,:,:,:), allocatable ::  interpolate_velocity
    !real, dimension(:,:,:,:), allocatable :: snapmatrix_pressure, snapmatrix_volumefraction
    !real, dimension(:,:,:), allocatable ::  interpolate_pressure,interpolate_volumefraction
    
    type(vector_field), pointer :: velocity_field
    type(scalar_field), pointer :: pressure_field, volumefraction_field
    integer :: i, dump_no, d, j,k ,m,i2,i3,ph,p
    integer :: stat    
    
  
    integer :: udim,deim_number_tmp  
    type(vector_field), pointer :: velocityudim 
    type(vector_field), pointer :: vfield
    character(len=1024) :: phase_name
    character(len = OPTION_PATH_LEN) :: simulation_name1
   ! call get_option('/reduced_model/pod_basis_formation/pod_basis_count', nsvd)
    call get_option('/simulation_name',simulation_name)
    nsvd =12
    numphase=option_count("/material_phase")
    call initialise_write_state
    ! Read state from .flml file
      call populate_state(state_show)
   ! call populate_multi_state(state_show)
    call print_state(state_show(1))
    
    call put_vtu_to_states(state,state_p)  
    call retrieve_nodesvalue_from_states(state,  snapmatrix_velocity, snapmatrix_pressure, snapmatrix_volumefraction)
    if (.false.) then
    if(m .eq.8 .and. sparse_max .eq.1) then
      call interpolate_snapmatrix_velocity_8dim_level1( snapmatrix_velocity, interpolate_velocity) !allocate(snapmatrix_velocity(no_smolyak_nodes,u_nodes,snapshots,dim,numphase))
      call interpolate_snapmatrix_pressure_8dim_level1( snapmatrix_pressure, interpolate_pressure) !allocate(snapmatrix_pressure(no_smolyak_nodes,p_nodes,snapshots,numphase))
      call interpolate_snapmatrix_volumefraction_8dim_level1( snapmatrix_volumefraction, interpolate_volumefraction)
    else if (m .eq.4 .and. sparse_max .eq.1) then
      call interpolate_snapmatrix_velocity_4dim_level1( snapmatrix_velocity, interpolate_velocity) !allocate(snapmatrix_velocity(no_smolyak_nodes,u_nodes,snapshots,dim,numphase))
      call interpolate_snapmatrix_pressure_4dim_level1( snapmatrix_pressure, interpolate_pressure) !allocate(snapmatrix_pressure(no_smolyak_nodes,p_nodes,snapshots,numphase))
      call interpolate_snapmatrix_volumefraction_4dim_level1( snapmatrix_volumefraction, interpolate_volumefraction)
    else if  (m .eq.4 .and. sparse_max .eq.2) then
       call interpolate_snapmatrix_velocity_4dim ( snapmatrix_velocity, interpolate_velocity) !allocate(snapmatrix_velocity(no_smolyak_nodes,u_nodes,snapshots,dim,numphase))
      call interpolate_snapmatrix_pressure_4dim ( snapmatrix_pressure, interpolate_pressure) !allocate(snapmatrix_pressure(no_smolyak_nodes,p_nodes,snapshots,numphase))
      call interpolate_snapmatrix_volumefraction_4dim ( snapmatrix_volumefraction, interpolate_volumefraction)
    endif 
    endif
     if (.false.) then
       call interpolate_snapmatrix_velocity( snapmatrix_velocity, interpolate_velocity) !allocate(snapmatrix_velocity(no_smolyak_nodes,u_nodes,snapshots,dim,numphase))
       call interpolate_snapmatrix_pressure( snapmatrix_pressure, interpolate_pressure) !allocate(snapmatrix_pressure(no_smolyak_nodes,p_nodes,snapshots,numphase))
       call interpolate_snapmatrix_volumefraction( snapmatrix_volumefraction, interpolate_volumefraction)
    endif

    if (.true.) then !CMAME-material-property
      call interpolate_snapmatrix_velocity_rbf_8d( snapmatrix_velocity, interpolate_velocity) !allocate(snapmatrix_velocity(no_smolyak_nodes,u_nodes,snapshots,dim,numphase))
      call interpolate_snapmatrix_pressure_rbf_8d( snapmatrix_pressure, interpolate_pressure) !allocate(snapmatrix_pressure(no_smolyak_nodes,p_nodes,snapshots,numphase))
      call interpolate_snapmatrix_volumefraction_rbf_8d( snapmatrix_volumefraction, interpolate_volumefraction)
    endif

    !   visiualize the vtus 
    !   print *, 'size: 2, 1, 1152, 1', size(snapmean_pressure%val), size(pressure_field%val )  
    !   pressure_field%val=snapmean_pressure%val  
    !   volumefraction_field%val(:)=volumefraction_field%val(:) 
    !   volumefraction_field%val(1,2,:)=1-volumefraction_field%val(1,1,:) 
    if(.true.) then
    do i=1, snapshots
     do ph=1, numphase
     pressure_field => extract_scalar_field( state_show(ph), "Pressure" )
     volumefraction_field=> extract_scalar_field(state_show(ph), "PhaseVolumeFraction" )
     velocity_field=>extract_vector_field(state_show(ph),"Velocity")
    !  print *, 'velocity_field%val', size(velocity_field%val,1), size(velocity_field%val,2)
      velocity_field%val(1,:)=interpolate_velocity(:,i,1,ph)
      velocity_field%val(2,:)=interpolate_velocity(:,i,2,ph)
      call insert(state_show(ph), velocity_field, "Velocity")
   !  print *, 'size(pressure_field%val', size(pressure_field%val), size(interpolate_pressure,1)
     pressure_field%val=interpolate_pressure(:,i,ph)
     
     call insert(state_show(ph), pressure_field, "pressure")
     volumefraction_field%val=interpolate_volumefraction(:,i,ph)
     call insert(state_show(ph), volumefraction_field, "PhaseVolumeFraction")
     dump_no=i
     call set_option('/simulation_name', trim(simulation_name)//'_interp') 
     call write_state(dump_no, state_show)
    enddo
    enddo
   endif
    
    !!!!!!!form podbasis for interpolated snapshots
    allocate(snapmean_velocity(u_nodes,dimen,numphase))
    allocate(snapmean_pressure(p_nodes,numphase))  
    allocate(snapmean_volumefraction(v_nodes,numphase))  
     snapmean_velocity=0
     snapmean_pressure=0
    snapmean_volumefraction=0 
     do p = 1, numphase
     do i=1, snapshots
       do d=1, dimen
          snapmean_velocity(:,d,p)= snapmean_velocity(:,d,p)+interpolate_velocity(:,i,d,p)
       enddo
       snapmean_pressure(:,p)=snapmean_pressure(:,p)+interpolate_pressure(:,i,p)
       snapmean_volumefraction(:,p)=snapmean_volumefraction(:,p)+interpolate_volumefraction(:,i,p)
    end do

    do d=1,dimen
       snapmean_velocity(:,d,p)=snapmean_velocity(:,d,p)/snapshots
    enddo
       snapmean_pressure(:,p)=snapmean_pressure(:,p)/snapshots
       snapmean_volumefraction(:,p)=snapmean_volumefraction(:,p)/snapshots

    do i=1,snapshots
       do d=1,dimen
          interpolate_velocity(:,i,d,p)=interpolate_velocity(:,i,d,p)-snapmean_velocity(:,d,p)
       enddo
       interpolate_pressure(:,i,p)=interpolate_pressure(:,i,p)-snapmean_pressure(:,p)
       interpolate_volumefraction(:,i,p)=interpolate_volumefraction(:,i,p)-snapmean_volumefraction(:,p)
    enddo
    enddo !end numphase  

    call form_svd_mp(interpolate_velocity, interpolate_pressure, interpolate_volumefraction,&
         & leftsvd_velocity, leftsvd_pressure, leftsvd_volumefraction, svdval_velocity, svdval_pressure, svdval_volumefraction, snapshots)
        
    call form_podstate_velocity(state(1,:,:),pod_state,leftsvd_velocity,leftsvd_pressure, snapmean_velocity, snapmean_pressure)
    call form_podstate_pressure(state_p(1,:,:),pod_state_p,leftsvd_velocity,leftsvd_pressure, snapmean_velocity, snapmean_pressure)
    call form_podstate_volumefraction(state_p(1,:,:),pod_state_v,leftsvd_volumefraction, snapmean_volumefraction)
    ! enddo
    call add_option("/reduced_model/execute_reduced_model",stat)
    call get_option('/simulation_name',simulation_name)
    do ph=1, numphase
           call get_option("/material_phase["//int2str(ph-1)//"]/name", phase_name)
    do i=1,nsvd 
        dump_no=i 
      call vtk_write_state(filename=trim(simulation_name)//"_"//trim(phase_name)//"_PODBasisVelocity", index=dump_no, model="VelocityMesh", state=pod_state(i,:,ph))
      call vtk_write_state(filename=trim(simulation_name)//"_"//trim(phase_name)//"_PODBasisPressure", index=dump_no, model="PressureMesh", state=pod_state_p(i,:,ph))
      call vtk_write_state(filename=trim(simulation_name)//"_"//trim(phase_name)//"_PODBasisVolumefraction", index=dump_no, model="VolumefractionMesh", state=pod_state_v(i,:,ph))
             
    enddo
    enddo !numphase 
  
     call set_option('/simulation_name', trim(simulation_name)//'_POD') 
     call get_option('/simulation_name',simulation_name) 
     call write_options(trim(simulation_name)//".mpml")

      call calculate_pod_coef_all_timesteps()
     call non_intrusive()

    deallocate(snapmatrix_velocity)
    deallocate(snapmatrix_pressure) 
    deallocate(snapmatrix_volumefraction) 
    deallocate(interpolate_pressure)
    deallocate(interpolate_velocity)
    deallocate(interpolate_volumefraction)
    
  end subroutine interpolate_snapshots_main 
 
   subroutine populate_multi_state(state)
        implicit none
        type(state_type), dimension(:), pointer, intent(inout) :: state
        !Local variables
        integer :: i, nphase, npres
        character( len = option_path_len ) :: option_path
        nphase = option_count("/material_phase")
        npres = option_count("/material_phase/scalar_field::Pressure/prognostic")
        !Adjust nphase to not account for the extra phases added by the wells
        nphase = nphase/npres

        !Prepare some specific modifications prior to populating state
        !If the extra mesh have not been created, create them here
        if (.not. have_option("/geometry/mesh::VelocityMesh_Continuous")) then
            call copy_option("/geometry/mesh::VelocityMesh", "/geometry/mesh::VelocityMesh_Continuous")
            call set_option("/geometry/mesh::VelocityMesh_Continuous/from_mesh/mesh_continuity", "continuous")
        end if

        if (.not. have_option("/geometry/mesh::PressureMesh_Continuous")) then
            call copy_option("/geometry/mesh::PressureMesh", "/geometry/mesh::PressureMesh_Continuous")
            call set_option("/geometry/mesh::PressureMesh_Continuous/from_mesh/mesh_continuity", "continuous")
        end if

        if (.not. have_option("/geometry/mesh::PressureMesh_Discontinuous")) then
            call copy_option("/geometry/mesh::PressureMesh", "/geometry/mesh::PressureMesh_Discontinuous")
            call set_option("/geometry/mesh::PressureMesh_Discontinuous/from_mesh/mesh_continuity", "discontinuous")
        end if

        if (.not. have_option("/geometry/mesh::P0DG")) then
            call copy_option("/geometry/mesh::PressureMesh", "/geometry/mesh::P0DG")
            call set_option("/geometry/mesh::P0DG/from_mesh/mesh_shape/polynomial_degree", 0)
            call set_option("/geometry/mesh::P0DG/from_mesh/mesh_continuity", "discontinuous")
        end if


        if (have_option('/mesh_adaptivity/hr_adaptivity/adapt_mesh_within_FPI')) then
            !Create necessary backups for storing the saturation (in a way that it is also adapted)
            do i = 1, nphase
                option_path = "/material_phase["// int2str( i - 1 )//"]/scalar_field::Saturation_bak"
                call copy_option("/material_phase["// int2str( i - 1 )//"]/scalar_field::PhaseVolumeFraction",&
                     trim(option_path))
                !Make sure the field is not shown
                if (.not.have_option(trim(option_path)//"/prognostic/output/exclude_from_vtu")) then
                    !Don't know how to set exclude_from_vtu to true from the spud options, hence,
                    !since Porous_media HAS to be true I copy it to obtain the same effect
                    call copy_option("/geometry/mesh::VelocityMesh/from_mesh/mesh_shape/Porous_media",&
                     trim(option_path)//"/prognostic/output/exclude_from_vtu")
                end if
                !Make sure that this field is not the objective of adaptivity
                if (have_option(trim(option_path)//"/prognostic/adaptivity_options")) then
                    call delete_option(trim(option_path)//"/prognostic/adaptivity_options")
                end if
            end do
        end if


        if (have_option('/physical_parameters/black-oil_PVT_table')) then

        !Maybe no need for a field in state? and just calculate using required conditions????

            !Create necessary memory to store the mass fraction of one of the pseudo-components (in a way that it is also adapted)
            if (nphase /= 3)then
                FLAbort('Black-Oil modelling requires three phases. Phase 1 Aqua, phase 2 liquid, phase 3 vapour')
            end if

            option_path = "/material_phase["// int2str( nphase -1 )//"]/scalar_field::VapourMassFraction"
            call copy_option("/material_phase["// int2str( nphase - 1 )//"]/scalar_field::PhaseVolumeFraction",&
                 trim(option_path))
            !Make sure the field is not shown
            if (.not.have_option(trim(option_path)//"/prognostic/output/exclude_from_vtu")) then
                !Don't know how to set exclude_from_vtu to true from the spud options, hence,
                !since Porous_media HAS to be true I copy it to obtain the same effect
                call copy_option("/geometry/mesh::VelocityMesh/from_mesh/mesh_shape/Porous_media",&
                 trim(option_path)//"/prognostic/output/exclude_from_vtu")
            end if
            !Make sure that this field is not the objective of adaptivity
            if (have_option(trim(option_path)//"/prognostic/adaptivity_options")) then
                call delete_option(trim(option_path)//"/prognostic/adaptivity_options")
            end if
        end if


        

        !Call fluidity to populate state
        call populate_state(state)

    end subroutine populate_multi_state


  subroutine calculate_pod_coef_all_timesteps()    
     type(vector_field), pointer :: snapmean_velocity1
      type(scalar_field), pointer :: snapmean_pressure1, snapmean_volumefraction1
      type(vector_field), pointer :: PODBasis_velocity 
      type(scalar_field), pointer :: PODBasis_pressure, podbasis_volumefraction
    integer :: i, j, k,s 
     if (have_option("/reduced_model/Velocity_Interpolation")) then
           allocate(pod_coef_all_obv(snapshots,(dimen+2)*nsvd,numphase)) 
           ! call read_pod_basis_interp(POD_state, state_show)
              snapmean_velocity=0
     snapmean_pressure=0
    snapmean_volumefraction=0
            do s=1, snapshots
                 do k=1, numphase
                             !pressure_field=>extract_scalar_field(packed_state,"FEPressure")
                    	     !volumefraction_field=>extract_tensor_field(packed_state,"PackedPhaseVolumeFraction")
                             !velocity_field=>extract_tensor_field(packed_state,"PackedVelocity") 
                   do i=1,size(POD_state,1)
                            ! PODBasis_velocity=>extract_vector_field(POD_state(i,1,k), "PODVelocity") !!podnum,:,numphase
                            ! snapmean_velocity1=>extract_vector_field(POD_state(i,1,k),"SnapmeanVelocity")
                            ! snapmean_pressure1=>extract_Scalar_field(POD_state(i,2,k),"SnapmeanPressure")
                             !PODBasis_pressure=>extract_scalar_field(POD_state(i,2,k), "PODPressure")
                            ! PODBasis_volumefraction=>extract_scalar_field(POD_state(i,3,k), "PODVolumefraction")
                           !  snapmean_volumefraction1=>extract_Scalar_field(POD_state(i,3,k),"SnapmeanVolumefraction")
 
                             !romdim=velocity%dim
                             do j=1,dimen!velocity%dim
                             !  print *, 'PODBasis_velocity%val(j,:), u%val(j,:)',PODBasis_velocity%val(j,:), u%val(j,:),snapmean_velocity%val(j,:)
                             !allocate(snapmean_velocity(u_nodes,dimen,numphase)) allocate(snapmean_pressure(p_nodes,numphase))
   				! allocate(leftsvd_velocity(u_nodes,nsvd,dim,numphase))    allocate(leftsvd_pressure(p_nodes,nsvd,numphase))
                             pod_coef_all_obv(s,i+nsvd*(j-1),k)=&
                                        dot_product(leftsvd_velocity(:,i,j,k),(interpolate_velocity(:,s,j,k)-snapmean_velocity(:,j,k)))   
                                                                     !(:,i,1,ph)(u_nodes,snapshots,dim,numphase))
                             enddo 
                             pod_coef_all_obv(s,i+nsvd*dimen,k)= &
                                        dot_product(leftsvd_pressure(:,i,k),(interpolate_pressure(:,s,k)-snapmean_pressure(:,k)))
                             
                             pod_coef_all_obv(s,i+nsvd*(dimen+1),k)= &
                                        dot_product(leftsvd_volumefraction(:,i,k),(interpolate_volumefraction(:,s,k)-snapmean_volumefraction(:,k)))
                     enddo   
                       if((s .eq. 1).and.(k.eq.1)) then                       
                    	open(11,file="coef_pod_all_obv"//"_"//trim(simulation_name))  
                           write(11,*)(pod_coef_all_obv(s,1:nsvd*(dimen+1),k)) !no volumefraction
                           print *, 'total_timestep----timestep===',simulation_name, snapshots, s 
                   	close(11) 
                    else 
                   	open(11,file="coef_pod_all_obv"//"_"//trim(simulation_name), position='append',ACTION='WRITE')  
                         write(11,*)(pod_coef_all_obv(s,1:nsvd*(dimen+1),k)) !no volumefraction
                         print *, 'total_timestep----timestep===', simulation_name, snapshots, s                         
                   	close(11)
                     endif
                   ! stop 12121
                  enddo !numphase
              enddo
         else 
             allocate(pod_coef_all_obv(snapshots,2*nsvd,numphase)) 
           ! call read_pod_basis_interp(POD_state, state_show)
             !snapmean_pressure=0
            ! snapmean_volumefraction=0
            do s=1, snapshots
                 do k=1, numphase 
                   do i=1,size(POD_state,1)                            
                             pod_coef_all_obv(s,i,k)=dot_product(leftsvd_pressure(:,i,k),(interpolate_pressure(:,s,k)))!-snapmean_pressure(:,k)
                             
                             pod_coef_all_obv(s,i+nsvd,k)= dot_product(leftsvd_volumefraction(:,i,k),(interpolate_volumefraction(:,s,k)))!-snapmean_volumefraction(:,k)
                     enddo   
                       if((s .eq. 1).and.(k.eq.1)) then                       
                    	open(11,file="coef_pod_all_obv"//"_"//trim(simulation_name))  
                           write(11,*)(pod_coef_all_obv(s,:,k)) 
                           print *, 'total_timestep----timestep===', simulation_name,snapshots, s 
                   	close(11) 
                    else 
                   	open(11,file="coef_pod_all_obv"//"_"//trim(simulation_name), position='append',ACTION='WRITE')  
                         write(11,*)(pod_coef_all_obv(s,:,k))
                         print *, 'total_timestep----timestep===', simulation_name, snapshots, s                         
                   	close(11)
                     endif
                   ! stop 12121
                  enddo !numphase
              enddo

         endif 
      
     
  end subroutine calculate_pod_coef_all_timesteps

  subroutine non_intrusive()
    type(scalar_field), pointer :: pressure_field, volumefraction_field
    real, dimension(:), allocatable ::  pod_coef_orig,pod_coef_new,pod_coef_old
    integer :: i, ph,j,k,s,dump_no 
    integer ::  total_timestep
    real, allocatable :: w(:)
    real, allocatable :: wm(:,:)
    real :: r0
    real , allocatable :: xd(:,:),xd_lq(:,:), pod_coef_all(:,:), smoylak_all(:,:) 
    real , allocatable :: xi(:,:)
    real , allocatable :: zd(:),fd(:)
    real , allocatable :: ze(:),temp_coef(:)
    real , allocatable :: zi(:),minin(:),maxin(:)
    real , allocatable :: zpi(:)
    real , allocatable :: fe(:,:) !exact values of function produced from fluidity
    real ( kind = 8 ), allocatable :: fi(:) 
    real :: maxall,minall
    integer :: m,ll,kk,u_nonods,l,nd, podnum,ni
    numphase=option_count("/material_phase")
    total_timestep=total_dumps
    if (have_option("/reduced_model/Velocity_Interpolation")) then 
        m = (dimen+1)*nsvd*numphase
        
      else 
        m = 2*nsvd*numphase
    endif
    podnum =m   
     allocate(pod_coef_all(total_timestep,podnum))    
     allocate(pod_coef_new(m)) 
     allocate(pod_coef_old(m)) 
    allocate(w(m))
   if (have_option("/reduced_model/Velocity_Interpolation")) then 
      open(unit=61,file="coef_pod_all_obv"//"_"//trim(simulation_name))   
      read(61,*)((pod_coef_all(j,k),k=1,podnum),j=1,total_timestep)
      close(61)      
   else
     open(unit=61,file="coef_pod_all_obv"//"_"//trim(simulation_name))   
     read(61,*)((pod_coef_all(j,k),k=1,podnum),j=1,total_timestep)      
     close(61)      
    endif
    do j=1, total_timestep-1
     ! pod_coef_all(j,:)=pod_coef_all(j+1,:)
    enddo
   Do i=1,m   ! cannot use minval because some of the vector is null.          
          minall=pod_coef_all(1,1)
          maxall=pod_coef_all(1,1)      
     do j=1,total_timestep
           if(minall>pod_coef_all(j,i)) then
              minall=pod_coef_all(j,i)          
           endif
           if(maxall<pod_coef_all(j,i)) then
              maxall=pod_coef_all(j,i)
           endif
      enddo
   ENDDO 
  
   
   ni = 1
   nd = total_timestep-1
   allocate ( xi(1:m,1:ni) )
   allocate ( fi(1:ni) ) 
   allocate (fd(nd)) 
   allocate (wm(nd,m))
   allocate (xd(podnum,nd)) 
     do j=1,m
         do k=1,nd
           xd(j,k)=pod_coef_all(k,j)
         enddo
       enddo
     r0 = ( maxall - minall ) / real ( nd, kind = 8 )
     do l=1, m
        do k=1,nd
          fd(k)=pod_coef_all(k+1,l)  ! target function value is the l^th next timestep's pod_coefficient  
         ! print *,'fd', fd(k)           
        enddo 
      call rbf_weight ( m, nd, xd, r0, phi4, fd, wm(:,l) )
     enddo
     
    pod_coef_old=pod_coef_all(1,:)
    DO s=1, total_timestep!=total_dumps     
         xi(:,1)= pod_coef_old !pod_coef_all_obv(kk,:)pod_coef_all(s,:) ! 
       ! print *,'xiiiiiiiiiii=',  xi
     do l=1, m   
           
      call rbf_interp_nd ( m, nd, xd, r0, phi4, wm(:,l), ni, xi, fi )
      !PHI4 evaluates the gaussian radial basis function,, best
      ! PHI3 evaluates the thin-plate spline radial basis function. terrible
      !PHI2 evaluates the inverse multiquadric radial basis function.,! 
      ! PHI1 evaluates the multiquadric radial basis function. terrible.! 
       pod_coef_new(l)=fi(1)  
     enddo  !do l=1, m      
    ! visulisation
      !print *, 'pod_coef_new', pod_coef_new
     do k=1, numphase 
      pressure_field => extract_scalar_field(state_show(k), "Pressure" )
      volumefraction_field=> extract_scalar_field(state_show(k), "PhaseVolumeFraction" )  
       pressure_field%val(:)=0
       volumefraction_field%val(:)=0  
           do i=1,nsvd  
              !print *, 'difference', pod_coef_new(i+(k-1)*(m/numphase))-pod_coef_all_obv(s,i,k)
               pressure_field%val(:)=pressure_field%val(:)+pod_coef_all_obv(s,i,k)*leftsvd_pressure(:,i,k)
               volumefraction_field%val(:)=volumefraction_field%val(:)+pod_coef_all_obv(s,i+nsvd,k)*leftsvd_volumefraction(:,i,k)  
              ! pressure_field%val(:)=pressure_field%val(:)+pod_coef_new(i+(k-1)*(m/numphase))*leftsvd_pressure(:,i,k)
              ! volumefraction_field%val(:)=volumefraction_field%val(:)+pod_coef_new(i+nsvd+(k-1)*(m/numphase))*leftsvd_volumefraction(:,i,k) 
           enddo 
      snapmean_velocity=0
     snapmean_pressure=0
    snapmean_volumefraction=0
     pressure_field%val=pressure_field%val(:)+ snapmean_pressure(:,k)   
     call insert(state_show(k), pressure_field, "pressure")
     volumefraction_field%val=volumefraction_field%val+snapmean_volumefraction(:,k)
     call insert(state_show(k), volumefraction_field, "PhaseVolumeFraction")
      enddo !do k=1, numphase 
     dump_no=s
      call set_option('/simulation_name', trim(simulation_name)//'_final') 
     call write_state(dump_no, state_show)

      pod_coef_old=pod_coef_new ! iteration  
    ENDDO ! Do s=1, total_timestep
     deallocate(pod_coef_new) 
    deallocate(pod_coef_old) 
     deallocate ( fi )
     deallocate ( xi ) 
     deallocate(fd) 
     deallocate(wm)
     deallocate(xd)
     deallocate(w)
    end subroutine non_intrusive

 
  subroutine non_intrusive_22d()
    type(scalar_field), pointer :: pressure_field, volumefraction_field
    real, dimension(:), allocatable ::  pod_coef_orig,pod_coef_new,pod_coef_old
    integer :: i, ph,j,k,s,dump_no 
    integer ::  total_timestep
    real, allocatable :: w(:)
    real, allocatable :: wm(:,:)
    real :: r0
    real , allocatable :: xd(:,:),xd_lq(:,:), pod_coef_all(:,:), smoylak_all(:,:) 
    real , allocatable :: xi(:,:)
    real , allocatable :: zd(:),fd(:)
    real , allocatable :: ze(:),temp_coef(:)
    real , allocatable :: zi(:),minin(:),maxin(:)
    real , allocatable :: zpi(:)
    real , allocatable :: fe(:,:) !exact values of function produced from fluidity
    real ( kind = 8 ), allocatable :: fi(:) 
    real :: maxall,minall
    integer :: m,ll,kk,u_nonods,l,nd, podnum,ni
    numphase=option_count("/material_phase")
    total_timestep=total_dumps
    if (have_option("/reduced_model/Velocity_Interpolation")) then 
        m = (dimen+1)*nsvd*numphase 
    else 
        m = 2*nsvd*numphase
    endif
    podnum =m   
     allocate(pod_coef_all(total_timestep,podnum))    
     allocate(pod_coef_new(m)) 
     allocate(pod_coef_old(m)) 
    allocate(w(m))
   if (have_option("/reduced_model/Velocity_Interpolation")) then 
      open(unit=61,file="coef_pod_all_obv"//"_"//trim(simulation_name))   
      read(61,*)((pod_coef_all(j,k),k=1,podnum),j=1,total_timestep)
      close(61)      
   else
     open(unit=61,file="coef_pod_all_obv"//"_"//trim(simulation_name))   
     read(61,*)((pod_coef_all(j,k),k=1,podnum),j=1,total_timestep)      
     close(61)      
    endif
    do j=1, total_timestep-1
     ! pod_coef_all(j,:)=pod_coef_all(j+1,:)
    enddo
   Do i=1,m   ! cannot use minval because some of the vector is null.          
          minall=pod_coef_all(1,1)
          maxall=pod_coef_all(1,1)      
     do j=1,total_timestep
           if(minall>pod_coef_all(j,i)) then
              minall=pod_coef_all(j,i)          
           endif
           if(maxall<pod_coef_all(j,i)) then
              maxall=pod_coef_all(j,i)
           endif
      enddo
   ENDDO 
   ni = 1
   nd = total_timestep-1
   allocate ( xi(1:m,1:ni) )
   allocate ( fi(1:ni) ) 
   allocate (fd(nd)) 
   allocate (wm(nd,m))
   allocate (xd(podnum,nd)) 
     do j=1,m
         do k=1,nd
           xd(j,k)=pod_coef_all(k,j)
         enddo
       enddo
     r0 = ( maxall - minall ) / real ( nd, kind = 8 )
     do l=1, m
        do k=1,nd
          fd(k)=pod_coef_all(k+1,l)  ! target function value is the l^th next timestep's pod_coefficient  
         ! print *,'fd', fd(k)           
        enddo 
      call rbf_weight ( m, nd, xd, r0, phi4, fd, wm(:,l) )
     enddo
     
    pod_coef_old=pod_coef_all(1,:)
    DO s=1, total_timestep!=total_dumps     
         xi(:,1)= pod_coef_old !pod_coef_all_obv(kk,:)pod_coef_all(s,:) ! 
       ! print *,'xiiiiiiiiiii=',  xi
     do l=1, m   
           
      call rbf_interp_nd ( m, nd, xd, r0, phi4, wm(:,l), ni, xi, fi )
      ! PHI4 evaluates the gaussian radial basis function, best
      ! PHI3 evaluates the thin-plate spline radial basis function. terrible
      ! PHI2 evaluates the inverse multiquadric radial basis function ! 
      ! PHI1 evaluates the multiquadric radial basis function. terrible ! 
       pod_coef_new(l)=fi(1)  
     enddo  !do l=1, m      
    ! visulisation
      !print *, 'pod_coef_new', pod_coef_new
     do k=1, numphase 
      pressure_field => extract_scalar_field(state_show(k), "Pressure" )
      volumefraction_field=> extract_scalar_field(state_show(k), "PhaseVolumeFraction" )  
       pressure_field%val(:)=0
       volumefraction_field%val(:)=0  
           do i=1,nsvd  
              !print *, 'difference', pod_coef_new(i+(k-1)*(m/numphase))-pod_coef_all_obv(s,i,k)
               pressure_field%val(:)=pressure_field%val(:)+pod_coef_all_obv(s,i,k)*leftsvd_pressure(:,i,k)
               volumefraction_field%val(:)=volumefraction_field%val(:)+pod_coef_all_obv(s,i+nsvd,k)*leftsvd_volumefraction(:,i,k)  
              ! pressure_field%val(:)=pressure_field%val(:)+pod_coef_new(i+(k-1)*(m/numphase))*leftsvd_pressure(:,i,k)
              ! volumefraction_field%val(:)=volumefraction_field%val(:)+pod_coef_new(i+nsvd+(k-1)*(m/numphase))*leftsvd_volumefraction(:,i,k) 
           enddo 
      snapmean_velocity=0
      snapmean_pressure=0
      snapmean_volumefraction=0

    
     pressure_field%val=pressure_field%val(:)+ snapmean_pressure(:,k)   
     call insert(state_show(k), pressure_field, "pressure")
     volumefraction_field%val=volumefraction_field%val+snapmean_volumefraction(:,k)
     call insert(state_show(k), volumefraction_field, "PhaseVolumeFraction")
      enddo !do k=1, numphase 
     dump_no=s
      call set_option('/simulation_name', trim(simulation_name)//'_final') 
     call write_state(dump_no, state_show)
      
      pod_coef_old=pod_coef_new ! iteration  
    ENDDO ! Do s=1, total_timestep

    
   ! do ph=1, numphase
      !     call get_option("/material_phase["//int2str(ph-1)//"]/name", phase_name)
  !  do i=1,nsvd 
    !    dump_no=i 
    !  call vtk_write_state(filename=trim(simulation_name)//"_"//trim(phase_name)//"_PODBasisVelocity", index=dump_no, model="VelocityMesh", state=pod_state(i,:,ph))
    !  call vtk_write_state(filename=trim(simulation_name)//"_"//trim(phase_name)//"_PODBasisPressure", index=dump_no, model="PressureMesh", state=pod_state_p(i,:,ph))
    !  call vtk_write_state(filename=trim(simulation_name)//"_"//trim(phase_name)//"_PODBasisVolumefraction", index=dump_no, model="VolumefractionMesh", state=pod_state_v(i,:,ph))
             
   ! enddo
  !  enddo !numphase 
  
   
     deallocate(pod_coef_new) 
     deallocate(pod_coef_old) 
     deallocate ( fi )
     deallocate ( xi ) 
     deallocate(fd) 
     deallocate(wm)
     deallocate(xd)
     deallocate(w)
    end subroutine non_intrusive_22d

  subroutine put_vtu_to_states(state,state_p)
    !!< Read the input states from the vtu dumps of the forward run.
    type(state_type), intent(out), dimension(:,:,:), allocatable :: state
    type(state_type), intent(out), dimension(:,:,:), allocatable :: state_p
    character(len=1024) :: filename
     type(scalar_field) ::pressure
    integer :: dump_sampling_period, quadrature_degree
    integer :: i,j,k,stable_dumps
    type(mesh_type) ::pressuremesh
    integer :: numphase
    character(len=1024) :: phase_name
    
   ! call get_option('/reduced_model/pod_basis_formation/dump_sampling_period',dump_sampling_period)
    call get_option('/geometry/quadrature/degree', quadrature_degree)
    dump_sampling_period =1
    ewrite(3,*) 'dump_sampling_period',dump_sampling_period
    
    !substract gyre_0.vtu
    total_dumps=count_dumps(dump_sampling_period)
    !print *,  total_dumps, no_smolyak_nodes
   ! stop 333
    numphase=option_count("/material_phase")
    allocate(state(no_smolyak_nodes,total_dumps,numphase))
    allocate(state_p(no_smolyak_nodes,total_dumps,numphase))
    !no_smolyak_node=22
    !    vtu--->state
    do k=1, no_smolyak_nodes
     do j=1, numphase
         call get_option("/material_phase["//int2str(j-1)//"]/name", phase_name) ! not'but" 
         ! call get_option("/material_phase["//int2str(m)//"]/name", mat_name)
         print *,'phase_name=',  int2str(j), phase_name
     do i=1, total_dumps
         ! if(numphase .eq. 1) then
       !      write(filename, '(a, i0, a, a)') '/data/dxiao/Channel_model_23_cases/Case_'//trim(int2str(k))//'/'//simulation_name//'_'//trim(phase_name)//'_VelocityMesh_', (i)*dump_sampling_period,'_checkpoint',".vtu" 
             write(filename, '(a, i0, a, a)') trim(int2str(k))//'_'//trim(phase_name)//'_PressureMesh_',  (i)*dump_sampling_period,'_checkpoint',".vtu" 
             print *, filename
             stop 333
             call vtk_read_state(filename, state(k,i,j), quadrature_degree)
             print *,'filename=', filename  
             write(filename, '(a, i0, a,a)')  trim(int2str(k))//'_'//trim(phase_name)//'_PressureMesh_', (i)*dump_sampling_period,'_checkpoint',".vtu" 
             call vtk_read_state(filename, state_p(k,i,j), quadrature_degree) 
              print *,'filename=', filename
        !  else
          !   write(filename, '(a, i0, a,a)')  trim(k)//'_'//trim(phase_name)//'_VelocityMesh_', (i)*dump_sampling_period,'_checkpoint',".vtu"
           !  call vtk_read_state(filename, state(k,i,j), quadrature_degree)
           !  write(filename, '(a, i0, a,a)')  trim(k)//'_'//trim(phase_name)//'_PressureMesh_', (i)*dump_sampling_period,'_checkpoint',".vtu" 
           !  call vtk_read_state(filename, state_p(k,i,j), quadrature_degree)
        !  endif           
          pressuremesh =extract_mesh(state_p(k,i,j), "Mesh") 
        call insert(state(k,i,j), pressuremesh, name="PressureMesh") 
    enddo
    enddo  
    enddo
   end subroutine put_vtu_to_states

subroutine read_input_states(state,state_p)
 
     type(state_type), intent(inout), dimension(:,:,:), allocatable :: state 
     type(state_type), intent(out), dimension(:,:,:), allocatable :: state_p
    character(len=1024) :: filename
      
    integer :: dump_sampling_period, quadrature_degree
    integer :: i,j,k,stable_dumps
    type(mesh_type) ::pressuremesh
    integer :: numphase
    character(len=1024) :: phase_name
    type(scalar_field), pointer :: pressure, volumefraction
   ! call get_option('/reduced_model/pod_basis_formation/dump_sampling_period',dump_sampling_period)
    call get_option('/geometry/quadrature/degree', quadrature_degree)
    dump_sampling_period =1
    ewrite(3,*) 'dump_sampling_period',dump_sampling_period
    
    !substract gyre_0.vtu
    total_dumps=count_vtu(dump_sampling_period)
    !print *,  total_dumps, no_smolyak_nodes
   ! stop 333
    numphase=option_count("/material_phase")
    allocate(state(no_smolyak_nodes,total_dumps,numphase))
    allocate(state_p(no_smolyak_nodes,total_dumps,numphase))
    
    do k=1, no_smolyak_nodes
     do j=1, numphase
         call get_option("/material_phase["//int2str(j-1)//"]/name", phase_name) ! not'but"  
       do i=1, total_dumps    
         write(filename, '(a, i0, a)') trim(int2str(k))//'_', (i)*dump_sampling_period,".vtu"  
            call vtk_read_state(filename, state(k,i,j), quadrature_degree) 
           write(filename, '(a, i0, a,a)')  trim(int2str(k))//'_phase2'//'_PressureMesh_', (i)*dump_sampling_period,'_checkpoint',".vtu" 
             call vtk_read_state(filename, state_p(k,i,j), quadrature_degree) 
              print *,'filename=', filename 
              pressuremesh =extract_mesh(state_p(k,i,j), "Mesh") 
              call insert(state(k,i,j), pressuremesh, name="PressureMesh") 
              print *, filename
              pressure => extract_scalar_field(state_p(k,i,j), "Density")
              volumefraction=> extract_scalar_field(state_p(k,i,j), "PhaseVolumeFraction")
              call insert(state(k,i,j), pressure, name="Pressure")
              call insert(state(k,i,j), volumefraction, name="PhaseVolumeFraction")
             ! call vtk_read_state(filename, state_p(k,i,j), quadrature_degree)  
     end do
    enddo 
   enddo
     
       
  end subroutine read_input_states

    subroutine read_pod_basis_interp(POD_state, state)
    
    !! Read the podbasis from the set of vtu files.
   
    character(len=1024) :: simulation_name, filename, phase_name 
    integer :: dump_period, quadrature_degree, numphase
    integer :: i,j,POD_num,k,nfield
    ! type(state_type), dimension(:,:,:), allocatable :: POD_state
    type(state_type), dimension(:,:,:), allocatable :: POD_state,POD_state_p,POD_state_v!POD_state(podnum,:,numphase)
    type(state_type), dimension(:) :: state
    type(vector_field) :: podVelocity, newpodVelocity 
    type(scalar_field) :: podPressure, newpodPressure, podVolumefraction, newpodVolumefraction, podTemperature, newpodTemperature ,snapmean_pressure, snapmean_volumefraction
    type(mesh_type) :: VelocityMesh, PressureMesh,TemperatureMesh
    
    type(scalar_field), pointer :: pres
    integer :: pod_pnodes,pod_unodes,p_nodes,u_nodes
    !type(mesh_type) ,pointer ::pmesh
    call get_option('/simulation_name', simulation_name)
    print *, 'simulation_name', simulation_name
    ! if (have_option("/reduced_model/execute_reduced_model")) then
   ! If(have_option("/reduced_model/Smolyak")) then
    if (have_option("/reduced_model/Smolyak") .or. have_option("/reduced_model/RBF_interpolation")) then
        simulation_name(len_trim(simulation_name)-3:len_trim(simulation_name))="" !delete _POD in the name
    endif
    print *, 'simulation_name', simulation_name
    call get_option('/geometry/quadrature/degree', quadrature_degree)
    call get_option(&
         "/reduced_model/pod_basis_formation/pod_basis_count", POD_num) 
    numphase=option_count("/material_phase")
    !total_dumps=POD_num!count_dumps(simulation_name)
    nfield = vector_field_count( state(1) )+scalar_field_count( state(1))+tensor_field_count(state(1))
    
   ! allocate(pod_state(POD_num,nfield,numphase))
    allocate(pod_state_p(POD_num,nfield,numphase))
    allocate(pod_state_v(POD_num,nfield,numphase))
    VelocityMesh=extract_velocity_mesh(state)
    PressureMesh=extract_pressure_mesh(state)  
    !flag=1
  
    do k =1, numphase
       call get_option("/material_phase["//int2str(k-1)//"]/name", phase_name)
    do i=1, POD_num
       !! Note that this won't work in parallel. Have to look for the pvtu in that case.  
       write(filename, '(a, i0, a)') trim(simulation_name)//"_"//trim(phase_name)//'_PODBasisVelocity_', i,".vtu"        
       call vtk_read_state(filename, POD_state(i,1,k), quadrature_degree)

       write(filename, '(a, i0, a)') trim(simulation_name)//"_"//trim(phase_name)//'_PODBasisPressure_', i,".vtu"         
       call vtk_read_state(filename, POD_state_p(i,2,k), quadrature_degree) 
       !call print_state(POD_state_p(1,2,1))
       !!!no velocity mean, because it pushs to pod-state directly rather than transfer from pod-state_p to pod-state
       snapmean_pressure=extract_scalar_field(POD_state_p(i,2,k),"SnapmeanPressure")
       call insert(pod_state(i,2,k), snapmean_pressure, name="SnapmeanPressure")

       write(filename, '(a, i0, a)') trim(simulation_name)//"_"//trim(phase_name)//'_PODBasisVolumefraction_', i,".vtu"         
       call vtk_read_state(filename, POD_state_v(i,1,k), quadrature_degree)      
       !call print_state(POD_state_v(1,1,1))
       snapmean_volumefraction=extract_scalar_field(POD_state_v(i,1,k),"SnapmeanVolumefraction")
       call insert(pod_state(i,3,k), snapmean_volumefraction, name="SnapmeanVolumefraction") 
 
       PODVelocity=extract_vector_field(POD_state(i,1,k),"PODVelocity")
       call allocate(newpodVelocity, podVelocity%dim, VelocityMesh, "PODVelocity")
       call remap_field(from_field=podVelocity, to_field=newpodVelocity)
       call insert(POD_state(i,1,k), newpodVelocity, "PODVelocity")
       call deallocate(newpodVelocity)
       !print *, 'ggggggggggggggggggggggggggggggggggg'
    
       PODPressure=extract_scalar_field(POD_state_p(i,2,k),"PODPressure")   !!!!!!!!!!!!!!!!!!! 
       call allocate(newpodPressure, PressureMesh, "PODPressure")      
       call remap_field(from_field=podPressure, to_field=newpodPressure)     
       call insert(POD_state(i,2,k), newpodPressure, "PODPressure")      
       call deallocate(newpodPressure)

       podVolumefraction=extract_scalar_field(POD_state_v(i,1,k),"PODVolumefraction")   !!!!!!!!!!!!!!!!!!! 
       call allocate(newpodVolumefraction, PressureMesh, "PODVolumefraction")      
       call remap_field(from_field=podVolumefraction, to_field=newpodVolumefraction)     
       call insert(POD_state(i,3,k), newpodVolumefraction, "PODVolumefraction")      
       call deallocate(newpodVolumefraction)

      ! print *,'remap_field(from_field=podPressure, to_field=newpodPressure)  pass '
      if(.false.)then
       if(have_option('/material_phase::ocean/scalar_field::Temperature'))then
          TemperatureMesh=extract_mesh(state,"CoordinateMesh")
          PODTemperature=extract_scalar_field(POD_state(i,1,k),"PODTemperature")
          call allocate(newpodTemperature, TemperatureMesh, "PODTemperaturecsr_mult")
          call remap_field(from_field=podTemperature, to_field=newpodTemperature)
          call insert(POD_state(i,1,k), newpodTemperature, "PODTemperature")
          call deallocate(newpodTemperature)
       endif
      endif
     end do ! pod_num
    enddo  ! numphase
   !  deallocate(pod_state)
     deallocate(pod_state_p)
     deallocate(pod_state_v) 

  end subroutine read_pod_basis_interp

   subroutine retrieve_nodesvalue_from_states(state, snapmatrix_velocity, snapmatrix_pressure, snapmatrix_volumefraction)
             
   !!< 
    type(state_type), intent(in), dimension(:,:,:) :: state
    real, dimension(:,:,:,:,:), allocatable :: snapmatrix_velocity
    real, dimension(:,:,:,:), allocatable :: snapmatrix_pressure,snapmatrix_volumefraction
     
    ! integer :: numphase
    integer ::  i, d, p,k
    type(vector_field), pointer :: velocity
    type(scalar_field), pointer :: pressure, volumefraction  
    
    velocity => extract_vector_field(state(1,1,1), "Velocity")
    pressure => extract_scalar_field(state_p(1,1,1), "Pressure") 
    volumefraction=> extract_scalar_field(state_p(1,1,1), "PhaseVolumeFraction")
    dimen=velocity%dim
    p_nodes=node_count(pressure)
    u_nodes=node_count(velocity)
    v_nodes=node_count(volumefraction)
    snapshots=size(state(1,:,1))
    allocate(snapmatrix_velocity(no_smolyak_nodes,u_nodes,snapshots,dimen,numphase))
    allocate(snapmatrix_pressure(no_smolyak_nodes,p_nodes,snapshots,numphase)) 
    allocate(snapmatrix_volumefraction(no_smolyak_nodes,v_nodes,snapshots,numphase))
    
    snapmatrix_velocity=0.0
    snapmatrix_pressure=0.0  
    snapmatrix_volumefraction=0.0 
   do k=1, no_smolyak_nodes
   do p = 1, numphase
    do i = 1, snapshots
       velocity => extract_vector_field(state(k,i,p), "Velocity")
       pressure => extract_scalar_field(state_p(k,i,1), "Pressure")
       volumefraction=> extract_scalar_field(state_p(k,i,p), "PhaseVolumeFraction")
       do d = 1, dimen 
           snapmatrix_velocity(k,:,i,d,p)=velocity%val(d,:)
       end do 
           snapmatrix_pressure(k,:,i,p)=pressure%val
           snapmatrix_volumefraction(k,:,i,p)=volumefraction%val!(p,:)
    end do
   enddo 
   enddo !end numphase  
    
  end subroutine retrieve_nodesvalue_from_states
 
   subroutine retrieve_nodesvalue_from_states_22d(state, snapmatrix_velocity, snapmatrix_pressure, snapmatrix_volumefraction)
        
    type(state_type), intent(in), dimension(:,:,:) :: state
    real, dimension(:,:,:,:,:), allocatable :: snapmatrix_velocity
    real, dimension(:,:,:,:), allocatable :: snapmatrix_pressure,snapmatrix_volumefraction
    
    integer :: i, d, p, k
    type(vector_field), pointer :: velocity
    type(scalar_field), pointer :: pressure, volumefraction  
    
    velocity => extract_vector_field(state(1,1,1), "phase1::Velocity")
    pressure => extract_scalar_field(state(1,1,1), "Pressure") 
    volumefraction=> extract_scalar_field(state(1,1,1), "PhaseVolumeFraction")
    dimen=velocity%dim
    p_nodes=node_count(pressure)
    u_nodes=node_count(velocity)
    v_nodes=node_count(volumefraction)
    snapshots=size(state(1,:,1))
    allocate(snapmatrix_velocity(no_smolyak_nodes,u_nodes,snapshots,dimen,numphase))
    allocate(snapmatrix_pressure(no_smolyak_nodes,p_nodes,snapshots,numphase)) 
    allocate(snapmatrix_volumefraction(no_smolyak_nodes,v_nodes,snapshots,numphase))
    
    snapmatrix_velocity=0.0
    snapmatrix_pressure=0.0  
    snapmatrix_volumefraction=0.0 
   do k=1, no_smolyak_nodes
   do p = 1, numphase
    do i = 1, snapshots
       velocity => extract_vector_field(state(k,i,p), "phase1::Velocity")
       pressure => extract_scalar_field(state(k,i,1), "Pressure")
       volumefraction=> extract_scalar_field(state(k,i,p), "PhaseVolumeFraction")
       do d = 1, dimen 
           snapmatrix_velocity(k,:,i,d,p)=velocity%val(d,:)
       end do 
           snapmatrix_pressure(k,:,i,p)=pressure%val
           snapmatrix_volumefraction(k,:,i,p)=1-volumefraction%val!(p,:)
    end do
   enddo 
   enddo !end numphase  
    
  end subroutine retrieve_nodesvalue_from_states_22d

  subroutine form_svd_mp(snapmatrix_velocity, snapmatrix_pressure,snapmatrix_volumefraction,&
       & leftsvd_velocity, leftsvd_pressure, leftsvd_volumefraction, svdval_velocity, svdval_pressure, svdval_volumefraction, snapshots)
    
    real, dimension(:,:,:,:), intent(in) :: snapmatrix_velocity
    real, dimension(:,:,:), intent(in) :: snapmatrix_pressure, snapmatrix_volumefraction
    real, dimension(:,:,:,:), allocatable, intent(out) :: leftsvd_velocity
    real, dimension(:,:,:), allocatable, intent(out) :: leftsvd_pressure, leftsvd_volumefraction
    real, dimension(:,:,:), allocatable, intent(out) :: svdval_velocity
    real, dimension(:,:), allocatable, intent(out) :: svdval_pressure,svdval_volumefraction
    integer i, d, dim ,nsvd, snapshots, p_nodes, u_nodes, v_nodes, p

    call get_option(&
         '/reduced_model/pod_basis_formation/pod_basis_count', nsvd)

    dim=size(snapmatrix_velocity,3)
    p_nodes=size(snapmatrix_pressure,1)
    u_nodes=size(snapmatrix_velocity,1)
    v_nodes=size(snapmatrix_volumefraction,1)
    allocate(leftsvd_velocity(u_nodes,nsvd,dim,numphase))
    allocate(leftsvd_pressure(p_nodes,nsvd,numphase))
    allocate(svdval_velocity(nsvd,dim,numphase))
    allocate(svdval_pressure(nsvd,numphase))
    allocate(leftsvd_volumefraction(v_nodes,nsvd,numphase)) 
    allocate(svdval_volumefraction(nsvd,numphase))
    do p=1,numphase
    do d=1,dim
       call snapsvd(u_nodes,snapshots,snapmatrix_velocity(:,:,d,p),&
            nsvd,nsvd,leftsvd_velocity(:,:,d,p),svdval_velocity(:,d,p))
    end do

    call snapsvd(p_nodes,snapshots,snapmatrix_pressure(:,:,p),nsvd,nsvd,leftsvd_pressure(:,:,p),svdval_pressure(:,p))
    call snapsvd(v_nodes,snapshots,snapmatrix_volumefraction(:,:,p),nsvd,nsvd,leftsvd_volumefraction(:,:,p),svdval_volumefraction(:,p))
    enddo !numphase

  end subroutine form_svd_mp 
   
 subroutine form_podstate_velocity(state, pod_state, leftsvd_u, leftsvd_p, snapmean_u, snapmean_p)
   
    type(state_type), intent(in), dimension(:,:) :: state    
    type(state_type), intent(out), dimension(:,:,:), allocatable :: pod_state
    real, intent(in), dimension(:,:,:,:) :: leftsvd_u
    real, intent(in), dimension(:,:,:) :: leftsvd_p
    real, intent(in), dimension(:,:,:) :: snapmean_u
    real, intent(in), dimension(:,:) :: snapmean_p 

    type(mesh_type), pointer :: pod_xmesh, pod_umesh, pod_pmesh, pmesh, pod_mesh
    type(element_type) :: pod_xshape, pod_ushape, pod_pshape
    type(vector_field), pointer :: pod_positions, velocity
    type(scalar_field), pointer :: pressure

    type(vector_field) :: pod_velocity
    type(scalar_field) :: pod_pressure
    type(vector_field) :: snapmean_velocity
    type(scalar_field) :: snapmean_pressure

    real, dimension(:), pointer :: x_ptr,y_ptr,z_ptr
    real, dimension(:), allocatable :: x,y,z   
    
    character(len=1024) :: filename
    character(len = FIELD_NAME_LEN) :: field_name

    integer :: dump_sampling_period, quadrature_degree,nonods
    integer :: i,j,k,nod,total_dumps,POD_num,stat,f,d,p
    logical :: all_meshes_same
    integer :: pod_pnodes,pod_unodes

    call get_option(&
         "/reduced_model/pod_basis_formation/pod_basis_count", POD_num) 

    allocate(pod_state(POD_num,1,numphase))
    call nullify(pod_state)
    
    do p=1,numphase
    do i = 1,POD_num

      ! pod_mesh => extract_mesh(state(1), "Mesh")

       all_meshes_same = .true.

        print *, 'aaaaaaaaaaaaaaaaaaaaaaaa'
        pod_xmesh => extract_mesh(state(1,p), "Mesh")
        print *, 'bbbbbbbbbbbbbbbbbbbbbbbb'
        pod_umesh =>extract_mesh(state(1,p), "Mesh")
        pod_positions => extract_vector_field(state(1,p), "Coordinate")

 
       call insert(pod_state(i,1,p), pod_umesh, "CoordinateMesh")
       call insert(pod_state(i,1,p), pod_umesh, "VelocityMesh")
    
       call insert(pod_state(i,1,p), pod_positions, "Coordinate")

       velocity => extract_vector_field(state(1,p), "Velocity")

       call allocate(pod_velocity, velocity%dim, pod_umesh, "PODVelocity")
       call zero(pod_velocity)
       do d=1,velocity%dim
          call set_all(pod_velocity, d, leftsvd_u(:,i,d,p))
       end do
       call insert(pod_state(i,1,p), pod_velocity, name="PODVelocity")
       call deallocate(pod_velocity) 
       !!insert snapmean data into state
     
       call allocate(snapmean_velocity, velocity%dim, pod_umesh, "SnapmeanVelocity")
       call zero(snapmean_velocity)
       do d=1,velocity%dim
          call set_all(snapmean_velocity, d, snapmean_u(:,d,p))
       end do
       call insert(pod_state(i,1,p), snapmean_velocity, name="SnapmeanVelocity")
       call deallocate(snapmean_velocity)

   
  !test the podnodes        

      velocity => extract_vector_field(pod_state(i,1,p), "PODVelocity")
   
      pod_unodes=node_count(velocity)
     
     ! print*, "podddddddddddddddddddddddvelocity", pod_unodes
   
    enddo  
     enddo !numphase
  end subroutine form_podstate_velocity

  subroutine form_podstate_pressure(state_p, pod_state_p, leftsvd_u, leftsvd_p, snapmean_u, snapmean_p)
   
    type(state_type), intent(in), dimension(:,:) :: state_p
   ! type(state_type), pointer, dimension(:):: state
    type(state_type), intent(out), dimension(:,:,:), allocatable :: pod_state_p 
    real, intent(in), dimension(:,:,:,:) :: leftsvd_u
    real, intent(in), dimension(:,:,:) :: leftsvd_p
    real, intent(in), dimension(:,:,:) :: snapmean_u
    real, intent(in), dimension(:,:) :: snapmean_p 

    type(mesh_type), pointer :: pod_xmesh, pod_umesh, pod_pmesh, pmesh, pod_mesh
    type(element_type) :: pod_xshape, pod_ushape, pod_pshape
    type(vector_field), pointer :: pod_positions, velocity
    type(scalar_field), pointer :: pressure

    type(vector_field) :: pod_velocity
    type(scalar_field) :: pod_pressure
    type(vector_field) :: snapmean_velocity
    type(scalar_field) :: snapmean_pressure

    real, dimension(:), pointer :: x_ptr,y_ptr,z_ptr
    real, dimension(:), allocatable :: x,y,z   

    character(len=1024) :: filename
    character(len = FIELD_NAME_LEN) :: field_name

    integer :: dump_sampling_period, quadrature_degree,nonods
    integer :: i,j,k,nod,total_dumps,POD_num,stat,f,d,p
    logical :: all_meshes_same
    integer :: pod_pnodes,pod_unodes  

    call get_option(&
         "/reduced_model/pod_basis_formation/pod_basis_count", POD_num) 

    allocate(pod_state_p(POD_num,1,numphase))  
    call nullify(pod_state_p)
    do p=1,numphase
    do i = 1,POD_num

      ! pod_mesh => extract_mesh(state(1), "Mesh")

       all_meshes_same = .true.

     
        pod_xmesh => extract_mesh(state_p(1,p), "Mesh")
       
        pod_pmesh =>extract_mesh(state_p(1,p), "Mesh") 
 
       pod_positions => extract_vector_field(state_p(1,p), "Coordinate")

 
       call insert(pod_state_p(i,1,p), pod_xmesh, "CoordinateMesh") 
       call insert(pod_state_p(i,1,p), pod_pmesh, "PressureMesh")
       call insert(pod_state_p(i,1,p), pod_positions, "Coordinate")
     

       call allocate(pod_pressure, pod_pmesh, "PODPressure")    
       call zero(pod_pressure)
       call set_all(pod_pressure, leftsvd_p(:,i,p))
       call insert(pod_state_p(i,1,p), pod_pressure, name="PODPressure")
     
       call deallocate(pod_pressure)

       !!insert snapmean data into state 
       call allocate(snapmean_pressure, pod_pmesh, "SnapmeanPressure")
       call zero(snapmean_pressure)
       call set_all(snapmean_pressure, snapmean_p(:,p))
       call insert(pod_state_p(i,1,p), snapmean_pressure, name="SnapmeanPressure") 
       call deallocate(snapmean_pressure)
  !test the podnodes        

       
      pressure => extract_scalar_field(pod_state_p(i,1,p), "PODPressure")       
      pod_pnodes=node_count(pressure)  
   
   !   print*, "podddddddddddddddddddddddpressure", pod_pnodes
   
    enddo
    enddo     ! numphase  
  end subroutine form_podstate_pressure

  subroutine form_podstate_volumefraction(state_p,pod_state_v,leftsvd_v, snapmean_v)
    type(state_type), intent(in), dimension(:,:) :: state_p
    ! type(state_type), pointer, dimension(:):: state
    type(state_type), intent(out), dimension(:,:,:), allocatable :: pod_state_v     
    real, intent(in), dimension(:,:,:) ::  leftsvd_v     
    real, intent(in), dimension(:,:) :: snapmean_v
    type(mesh_type), pointer :: pod_xmesh, pod_umesh, pod_pmesh, pod_vmesh, pmesh, pod_mesh
    type(element_type) :: pod_xshape, pod_ushape, pod_pshape
    type(vector_field), pointer :: pod_positions, velocity
    type(scalar_field), pointer :: pressure,volumefraction

    type(vector_field) :: pod_velocity
    type(scalar_field) :: pod_pressure, pod_volumefraction
    type(vector_field) :: snapmean_velocity
    type(scalar_field) :: snapmean_pressure, snapmean_volumefraction

    real, dimension(:), pointer :: x_ptr,y_ptr,z_ptr
    real, dimension(:), allocatable :: x,y,z   

    character(len=1024) :: filename
    character(len = FIELD_NAME_LEN) :: field_name

    integer :: dump_sampling_period, quadrature_degree,nonods
    integer :: i,j,k,nod,total_dumps,POD_num,stat,f,d,p
    logical :: all_meshes_same
    integer :: pod_pnodes, pod_unodes , pod_vnodes, nstates 

    call get_option(&
         "/reduced_model/pod_basis_formation/pod_basis_count", POD_num) 
    nstates=option_count("/material_phase")
    allocate(pod_state_v(POD_num, 1,numphase))  
    call nullify(pod_state_v)
    do p=1,numphase
    do i = 1,POD_num

      ! pod_mesh => extract_mesh(state(1), "Mesh")

       all_meshes_same = .true. 
       pod_xmesh => extract_mesh(state_p(1,p), "Mesh") 
       pod_pmesh =>extract_mesh(state_p(1,p), "Mesh")  
       pod_positions => extract_vector_field(state_p(1,p), "Coordinate") 
       call insert(pod_state_v(i,1,p), pod_xmesh, "CoordinateMesh") 
       call insert(pod_state_v(i,1,p), pod_pmesh, "VolumefractionMesh")
       call insert(pod_state_v(i,1,p), pod_positions, "Coordinate") 
       call allocate(pod_volumefraction, pod_pmesh, "PODVolumefraction")    
       call zero(pod_volumefraction)
       call set_all(pod_volumefraction, leftsvd_v(:,i,p))
       call insert(pod_state_v(i,1,p), pod_volumefraction, name="PODVolumefraction")
     
       call deallocate(pod_volumefraction)

       !!insert snapmean data into state 
       call allocate(snapmean_volumefraction, pod_pmesh, "SnapmeanVolumefraction")
       call zero(snapmean_volumefraction)
       call set_all(snapmean_volumefraction, snapmean_v(:,p))
       call insert(pod_state_v(i,1,p), snapmean_volumefraction, name="SnapmeanVolumefraction") 
       call deallocate(snapmean_volumefraction)
  !test the podnodes  
      volumefraction => extract_scalar_field(pod_state_v(i,1,p), "PODVolumefraction")  
      volumefraction => extract_scalar_field(pod_state_v(i,1,p), "SnapmeanVolumefraction")      
     ! pod_vnodes=node_count(volumefraction)  
      !print*, "podddddddddddddddddddddddvolumefraction", volumefraction%val!pod_vnodes
   
    enddo
    enddo     ! numphase  
  end subroutine form_podstate_volumefraction
 subroutine form_podstate_velocity_22d(state, pod_state, leftsvd_u, leftsvd_p, snapmean_u, snapmean_p)
   
    type(state_type), intent(in), dimension(:,:) :: state    
    type(state_type), intent(out), dimension(:,:,:), allocatable :: pod_state
    real, intent(in), dimension(:,:,:,:) :: leftsvd_u
    real, intent(in), dimension(:,:,:) :: leftsvd_p
    real, intent(in), dimension(:,:,:) :: snapmean_u
    real, intent(in), dimension(:,:) :: snapmean_p 

    type(mesh_type), pointer :: pod_xmesh, pod_umesh, pod_pmesh, pmesh, pod_mesh
    type(element_type) :: pod_xshape, pod_ushape, pod_pshape
    type(vector_field), pointer :: pod_positions, velocity
    type(scalar_field), pointer :: pressure

    type(vector_field) :: pod_velocity
    type(scalar_field) :: pod_pressure
    type(vector_field) :: snapmean_velocity
    type(scalar_field) :: snapmean_pressure

    real, dimension(:), pointer :: x_ptr,y_ptr,z_ptr
    real, dimension(:), allocatable :: x,y,z   
    
    character(len=1024) :: filename
    character(len = FIELD_NAME_LEN) :: field_name

    integer :: dump_sampling_period, quadrature_degree,nonods
    integer :: i,j,k,nod,total_dumps,POD_num,stat,f,d,p
    logical :: all_meshes_same
    integer :: pod_pnodes,pod_unodes

    call get_option(&
         "/reduced_model/pod_basis_formation/pod_basis_count", POD_num) 

    allocate(pod_state(POD_num,1,numphase))
    call nullify(pod_state)
    
    do p=1,numphase
    do i = 1,POD_num

      ! pod_mesh => extract_mesh(state(1), "Mesh")

       all_meshes_same = .true.

        print *, 'aaaaaaaaaaaaaaaaaaaaaaaa'
        pod_xmesh => extract_mesh(state(1,p), "Mesh")
        print *, 'bbbbbbbbbbbbbbbbbbbbbbbb'
        pod_umesh =>extract_mesh(state(1,p), "Mesh")
        pod_positions => extract_vector_field(state(1,p), "Coordinate")

 
       call insert(pod_state(i,1,p), pod_umesh, "CoordinateMesh")
       call insert(pod_state(i,1,p), pod_umesh, "VelocityMesh")
    
       call insert(pod_state(i,1,p), pod_positions, "Coordinate")

       velocity => extract_vector_field(state(1,p), "phase"//int2str(p)//"::Velocity")

       call allocate(pod_velocity, velocity%dim, pod_umesh, "PODVelocity")
       call zero(pod_velocity)
       do d=1,velocity%dim
          call set_all(pod_velocity, d, leftsvd_u(:,i,d,p))
       end do
       call insert(pod_state(i,1,p), pod_velocity, name="PODVelocity")
       call deallocate(pod_velocity) 
       !!insert snapmean data into state
     
       call allocate(snapmean_velocity, velocity%dim, pod_umesh, "SnapmeanVelocity")
       call zero(snapmean_velocity)
       do d=1,velocity%dim
          call set_all(snapmean_velocity, d, snapmean_u(:,d,p))
       end do
       call insert(pod_state(i,1,p), snapmean_velocity, name="SnapmeanVelocity")
       call deallocate(snapmean_velocity)

   
  !test the podnodes        

      velocity => extract_vector_field(pod_state(i,1,p), "PODVelocity")
   
      pod_unodes=node_count(velocity)
     
     ! print*, "podddddddddddddddddddddddvelocity", pod_unodes
   
    enddo  
     enddo !numphase
  end subroutine form_podstate_velocity_22d

  subroutine interpolate_snapmatrix_velocity( snapmatrix_velocity, interpolate_velocity)  
    
! top-bottom between 1-2, pod_coef_all_obv is used for storing some podbasis to avoid to much modifications.
 
!  Parameters:
!
!    Local, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) SPARSE_MAX, the maximum sparse grid level to try.
!
!  Local Parameters:
!
!    Local, real ( kind = 8 ) A(M), B(M), the upper and lower variable limits 
!    in each dimension.
!
!    Local, real ( kind = 8 ) APP_ERROR, the averaged Euclidean norm of the 
!    difference between the sparse interpolant and the exact function at 
!    the interpolation points.
!
!    Local, integer ( kind = 4 ) C(0:L_MAX), the sparse grid coefficient vector.
!    Results at level L are weighted by C(L).
!
!    Local, integer ( kind = 4 ) IND(M), the 1D indices defining a Lagrange grid.
!    Each entry is a 1d "level" that specifies the order of a 
!    Clenshaw Curtis 1D grid.
!
!    Local, integer ( kind = 4 ) L, the current Lagrange grid level.
!
!    Local, integer ( kind = 4 ) L_MAX, the current sparse grid level.
!
!    Local, integer ( kind = 4 ) MORE, used to control the enumeration of all the
!    Lagrange grids at a current grid level.
!
!    Local, integer ( kind = 4 ) ND, the number of points used in a Lagrange grid.
!
!    Local, integer ( kind = 4 ) ND_TOTAL, the total number of points used in all the
!    Lagrange interpolants for a given sparse interpolant points that occur
!    in more than one Lagrange grid are counted multiple times.
!
!    Local, integer ( kind = 4 ) NI, the number of interpolant evaluation points.
!
!    Local, integer ( kind = 4 ) SPARSE_MIN, the minimum sparse grid level to try.
!
!    Local, real ( kind = 8 ) XD(M,ND), the data points for a Lagrange grid.
!
!    Local, real ( kind = 8 ) XI(M,NI), the interpolant evaluation points.
!
!    Local, real ( kind = 8 ) ZD(ND), the data values for a Lagrange grid.
!
!    Local, real ( kind = 8 ) ZE(NI), the exact function values at XI.
!
!    Local, real ( kind = 8 ) ZI(NI), the sparse interpolant values at XI.
!
!    Local, real ( kind = 8 ) ZPI(NI), one set of Lagrange interpolant values at XI.
!
  implicit none
  !type(state_type), intent(in), dimension(:,:,:) :: state
  real, dimension(:,:,:,:,:), allocatable :: snapmatrix_velocity ! allocate(snapmatrix_velocity(no_smolyak_nodes,u_nodes,snapshots,dim,numphase))
  real, dimension(:,:,:,:), intent(out),allocatable ::  interpolate_velocity 
  real ( kind = 8 ), allocatable :: a(:)
  real ( kind = 8 ) app_error
  real ( kind = 8 ), allocatable :: b(:)
  integer ( kind = 4 ), allocatable :: c(:)
  integer ( kind = 4 ) h
  integer ( kind = 4 ), allocatable :: ind(:)
  integer ( kind = 4 ) l
  integer ( kind = 4 ) l_max
  integer ( kind = 4 ) l_min
   
  logical more
  integer ( kind = 4 ) nd
  integer ( kind = 4 ) nd_total
  integer ( kind = 4 ) ni
  real ( kind = 8 ) r8vec_norm_affine
  integer ( kind = 4 ) seed
 ! integer ( kind = 4 ) sparse_max
  integer ( kind = 4 ) sparse_min
  integer ( kind = 4 ) t
  integer ( kind = 4 ), allocatable :: w(:)
  real ( kind = 8 ), allocatable :: xd(:,:)
  real ( kind = 8 ), allocatable :: xi(:,:)
  real ( kind = 8 ), allocatable :: zd(:)
  real ( kind = 8 ), allocatable :: ze(:)
  real ( kind = 8 ), allocatable :: zi(:)
  real ( kind = 8 ), allocatable :: zpi(:)
  integer :: i, ii, total, total_timestep,num_pod,j,k,jj,kk, smolyak_total_points,u_nodesi,snapshotsi,dimi,numphasei
   
  allocate (ind(1:m))

  ! Define the region.
  allocate ( a(1:m) )
  allocate ( b(1:m) )

  a(1:m) = 1
  b(1:m) = 2 
   
  !total_timestep=864
  !num_pod=12*36
  !smolyak_total_points=41!2**sparse_max+1
  !print *, '  smolyak_total_points=2**sparse_max+1=', smolyak_total_points 
  !  stop 2323
  !  Define the interpolation evaluation information.
  ni = 1
  seed = 123456789
  allocate ( xi(m,ni) )
  !call r8mat_uniform_abvec ( m, ni, a, b, seed, xi )
 !  print * , 'xi', xi
  print *, 'm,ni', m, ni 
 ! stop 345
  xi(1,1)= 1 
  xi(2,1)= 2
  xi(3,1)= 2
  xi(4,1)= 2!1.7 !2.0000000000000000   1.5000000000000000   2.0000000000000000   1.5000000000000000,21
    call get_option(&
         '/reduced_model/variables_to_interpolate/dimension_one', xi(1,1))
   call get_option(&
         '/reduced_model/variables_to_interpolate/dimension_two', xi(2,1))
 call get_option(&
         '/reduced_model/variables_to_interpolate/dimension_three', xi(3,1))
 call get_option(&
         '/reduced_model/variables_to_interpolate/dimension_four', xi(4,1))
  
  write ( *, '(a,i6)' ) '  Number of interpolation points is NI = ', ni
  allocate ( ze(1:ni) )
 ! call f_sinr ( m, ni, xi, ze )
  print *, 'actual function value is, will interpolate this value later on,then do comparison', xi,ze
!
!  Compute a sequence of sparse grid interpolants of increasing level.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '   L     Order    ApproxError'
  write ( *, '(a)' ) ''
  allocate (interpolate_velocity(u_nodes,snapshots,dimen,numphase))
  allocate ( zpi(1:ni) )
  allocate ( zi(1:ni) )
  !smolynode, u_nodes,snapshots,dim,numphase
  Do u_nodesi=1, u_nodes  
   Do snapshotsi=1, snapshots  
   do dimi=1, size(snapmatrix_velocity,4)
   do numphasei=1, size(snapmatrix_velocity,5)
  total = 0
  sparse_min = 0

  ! do l_max = sparse_min, sparse_max
    do l_max = sparse_max, sparse_max
    allocate ( c(0:l_max) )
    allocate ( w(0:l_max) )
    call smolyak_coefficients ( l_max, m, c, w )

    zi(1:ni) = 0.0D+00
    nd_total = 0
    
    l_min = max ( l_max + 1 - m, 0 )

    do l = l_min, l_max

      more = .false.
      
      do
     ! print *,'l_max,l,' ,l_max,l   
!  Define the next product grid.
!
        call comp_next ( l, m, ind, more, h, t )
!
!  Count the grid, find its points, evaluate the data there.
!
        call lagrange_interp_nd_size2 ( m, ind, nd )
        allocate ( xd(1:m,1:nd) )
        call lagrange_interp_nd_grid2 ( m, ind, a, b, nd, xd )        
         !  print *,'xd', xd, 点1 坐标为： xd(:,1)
        allocate ( zd(1:nd) )

       ! print *,'---------------------nd',nd  
        do ii=1, nd           
              if(xd(1,ii) .eq. 1.5 .and. xd(2,ii) .eq. 1.5   .and. xd(3,ii) .eq. 1.5 .and. xd(4,ii) .eq. 1.5) then 
!( 1.5000000000000000    1.5000000000000000   1.5000000000000000  1.5000000000000000)
           zd(ii)=snapmatrix_velocity(1,u_nodesi,snapshotsi,dimi,numphasei)
         else if(xd(1,ii) .eq. 1 .and. xd(2,ii) .eq. 1.5   .and. xd(3,ii) .eq. 1.5 .and. xd(4,ii) .eq. 1.5) then 
!(1.0000000000000000   1.5000000000000000   1.5000000000000000   1.5000000000000000)
           zd(ii)=snapmatrix_velocity(2,u_nodesi,snapshotsi,dimi,numphasei)
         else if(xd(1,ii) .eq. 2 .and. xd(2,ii) .eq. 1.5   .and. xd(3,ii) .eq. 1.5 .and. xd(4,ii) .eq. 1.5) then
 !( 2.0000000000000000   1.5000000000000000   1.5000000000000000   1.5000000000000000)
           zd(ii)=snapmatrix_velocity(3,u_nodesi,snapshotsi,dimi,numphasei)
         else if(xd(1,ii) .eq. 1.5 .and. xd(2,ii) .eq. 1   .and. xd(3,ii) .eq. 1.5 .and. xd(4,ii) .eq. 1.5) then 
!(1.5000000000000000   1.0000000000000000   1.5000000000000000   1.5000000000000000 )
           zd(ii)=snapmatrix_velocity(4,u_nodesi,snapshotsi,dimi,numphasei)
         else if(xd(1,ii) .eq. 1.5 .and. xd(2,ii) .eq. 2   .and. xd(3,ii) .eq. 1.5 .and. xd(4,ii) .eq. 1.5) then 
!(1.5000000000000000   2.0000000000000000   1.5000000000000000   1.5000000000000000 )
           zd(ii)=snapmatrix_velocity(5,u_nodesi,snapshotsi,dimi,numphasei)
        else if(xd(1,ii) .eq. 1.5 .and. xd(2,ii) .eq. 1.5   .and. xd(3,ii) .eq. 1 .and. xd(4,ii) .eq. 1.5) then 
!(1.5000000000000000   1.5000000000000000   1.0000000000000000   1.5000000000000000  )
           zd(ii)=snapmatrix_velocity(6,u_nodesi,snapshotsi,dimi,numphasei)
         else if(xd(1,ii) .eq. 1.5 .and. xd(2,ii) .eq. 1.5   .and. xd(3,ii) .eq. 2 .and. xd(4,ii) .eq. 1.5) then 
!(1.5000000000000000   1.5000000000000000   2.0000000000000000   1.5000000000000000  )
           zd(ii)=snapmatrix_velocity(7,u_nodesi,snapshotsi,dimi,numphasei)
         else if(xd(1,ii) .eq. 1.5 .and. xd(2,ii) .eq. 1.5   .and. xd(3,ii) .eq. 1.5 .and. xd(4,ii) .eq. 1) then
 !(1.5000000000000000   1.5000000000000000   1.5000000000000000   1.0000000000000000   )
           zd(ii)=snapmatrix_velocity(8,u_nodesi,snapshotsi,dimi,numphasei)
         else if(xd(1,ii) .eq. 1.5 .and. xd(2,ii) .eq. 1.5   .and. xd(3,ii) .eq. 1.5 .and. xd(4,ii) .eq. 2) then 
!(  1.5000000000000000   1.5000000000000000   1.5000000000000000   2.0000000000000000  )
           zd(ii)=snapmatrix_velocity(9,u_nodesi,snapshotsi,dimi,numphasei) 
         else if(int(xd(1,ii)*1000)/1000 .eq. 1.146 .and. xd(2,ii) .eq. 1.5 .and. xd(3,ii) .eq. 1.5 .and. xd(4,ii) .eq. 1.5) then 
               !(1.1464466094067263   1.5000000000000000   1.5000000000000000   1.5000000000000000 .and. xd(1,ii) <= 1.15 ..and. xd(2:4,ii) .eq. 1.5  )
           zd(ii)=snapmatrix_velocity(10,u_nodesi,snapshotsi,dimi,numphasei)
         elseif(int(xd(1,ii)*1000)/1000 .eq. 1.853 .and. xd(2,ii) .eq. 1.5 .and. xd(3,ii) .eq. 1.5 .and. xd(4,ii) .eq. 1.5) then 
               !(1.8535533905932737   1.5000000000000000   1.5000000000000000   1.5000000000000000.and. xd(3,ii) .eq. 1.5 .and. xd(4,ii) .eq. 1.5 )
           zd(ii)=snapmatrix_velocity(11,u_nodesi,snapshotsi,dimi,numphasei)
         else if(xd(1,ii) .eq. 1 .and. xd(2,ii) .eq. 1   .and. xd(3,ii) .eq. 1.5 .and. xd(4,ii) .eq. 1.5) then !(  1.0000000000000000   1.0000000000000000   1.5000000000000000   1.5000000000000000 )
           zd(ii)=snapmatrix_velocity(12,u_nodesi,snapshotsi,dimi,numphasei)
         else if(xd(1,ii) .eq. 2 .and. xd(2,ii) .eq. 1   .and. xd(3,ii) .eq. 1.5 .and. xd(4,ii) .eq. 1.5) then !( 2.0000000000000000   1.0000000000000000   1.5000000000000000   1.5000000000000000)
           zd(ii)=snapmatrix_velocity(13,u_nodesi,snapshotsi,dimi,numphasei)
          else if(xd(1,ii) .eq. 1 .and. xd(2,ii) .eq. 2   .and. xd(3,ii) .eq. 1.5 .and. xd(4,ii) .eq. 1.5) then !(1.0000000000000000   2.0000000000000000   1.5000000000000000   1.5000000000000000)
           zd(ii)=snapmatrix_velocity(14,u_nodesi,snapshotsi,dimi,numphasei)
         else if(xd(1,ii) .eq. 2 .and. xd(2,ii) .eq. 2   .and. xd(3,ii) .eq. 1.5 .and. xd(4,ii) .eq. 1.5) then !(2.0000000000000000   2.0000000000000000   1.5000000000000000   1.5000000000000000)
           zd(ii)=snapmatrix_velocity(15,u_nodesi,snapshotsi,dimi,numphasei)
         else if(xd(1,ii) .eq. 1.5 .and. int(xd(2,ii)*1000)/1000 .eq. 1.146 .and. xd(3,ii) .eq. 1.5 .and. xd(4,ii) .eq. 1.5) then 
            !(1.5000000000000000   1.1464466094067263   1.5000000000000000   1.5000000000000000 xd(2,ii) >=1.14  .and. xd(2,ii) <=1.15)
           zd(ii)=snapmatrix_velocity(16,u_nodesi,snapshotsi,dimi,numphasei)
         else if(xd(1,ii) .eq. 1.5 .and. int(xd(2,ii)*1000)/1000 .eq. 1.853 .and. xd(3,ii) .eq. 1.5 .and. xd(4,ii) .eq. 1.5) then
          !(1.5000000000000000   1.8535533905932737   1.5000000000000000   1.5000000000000000)
           zd(ii)=snapmatrix_velocity(17,u_nodesi,snapshotsi,dimi,numphasei)
           else if(xd(1,ii) .eq. 1 .and. xd(2,ii) .eq. 1.5 .and. xd(3,ii) .eq. 1 .and. xd(4,ii) .eq. 1.5) then !( 1.0000000000000000   1.5000000000000000   1.0000000000000000  1.5000000000000000)
           zd(ii)=snapmatrix_velocity(18,u_nodesi,snapshotsi,dimi,numphasei)
         else if(xd(1,ii) .eq. 2 .and. xd(2,ii) .eq. 1.5   .and. xd(3,ii) .eq. 1 .and. xd(4,ii) .eq. 1.5) then !(2.0000000000000000   1.5000000000000000   1.0000000000000000   1.5000000000000000)
           zd(ii)=snapmatrix_velocity(19,u_nodesi,snapshotsi,dimi,numphasei)
         else if(xd(1,ii) .eq. 1 .and. xd(2,ii) .eq. 1.5   .and. xd(3,ii) .eq. 2 .and. xd(4,ii) .eq. 1.5) then !(1.0000000000000000   1.5000000000000000   2.0000000000000000   1.5000000000000000)
           zd(ii)=snapmatrix_velocity(20,u_nodesi,snapshotsi,dimi,numphasei)
         else if(xd(1,ii) .eq. 2 .and. xd(2,ii) .eq. 1.5   .and. xd(3,ii) .eq. 2 .and. xd(4,ii) .eq. 1.5) then !(2.0000000000000000   1.5000000000000000   2.0000000000000000   1.5000000000000000)
           zd(ii)=snapmatrix_velocity(21,u_nodesi,snapshotsi,dimi,numphasei) 
         else if(xd(1,ii) .eq. 1.5 .and. xd(2,ii) .eq. 1   .and. xd(3,ii) .eq. 1 .and. xd(4,ii) .eq. 1.5) then !(1.5000000000000000   1.0000000000000000   1.0000000000000000   1.5000000000000000)
           zd(ii)=snapmatrix_velocity(22,u_nodesi,snapshotsi,dimi,numphasei)
         else if(xd(1,ii) .eq. 1.5 .and. xd(2,ii) .eq. 2   .and. xd(3,ii) .eq. 1 .and. xd(4,ii) .eq. 1.5) then !(1.5000000000000000   2.0000000000000000   1.0000000000000000   1.5000000000000000)
           zd(ii)=snapmatrix_velocity(23,u_nodesi,snapshotsi,dimi,numphasei)
         else if(xd(1,ii) .eq. 1.5 .and. xd(2,ii) .eq. 1   .and. xd(3,ii) .eq. 2 .and. xd(4,ii) .eq. 1.5) then !(1.5000000000000000   1.0000000000000000   2.0000000000000000   1.5000000000000000)
           zd(ii)=snapmatrix_velocity(24,u_nodesi,snapshotsi,dimi,numphasei)
         else if(xd(1,ii) .eq. 1.5 .and. xd(2,ii) .eq. 2   .and. xd(3,ii) .eq. 2 .and. xd(4,ii) .eq. 1.5) then !(1.5000000000000000   2.0000000000000000   2.0000000000000000   1.5000000000000000)
           zd(ii)=snapmatrix_velocity(25,u_nodesi,snapshotsi,dimi,numphasei)
          else if(xd(1,ii) .eq. 1.5 .and. xd(2,ii) .eq. 1.5 .and. int(xd(3,ii)*1000)/1000 .eq. 1.146 .and. xd(4,ii) .eq. 1.5) then 
       !(1.5000000000000000   1.5000000000000000  1.1464466094067263   1.5000000000000000)
           zd(ii)=snapmatrix_velocity(26,u_nodesi,snapshotsi,dimi,numphasei)
         else if(xd(1,ii) .eq. 1.5  .and. xd(2,ii) .eq. 1.5  .and. int(xd(3,ii)*1000)/1000 .eq. 1.853 .and. xd(4,ii) .eq. 1.5) then 
       !( 1.5000000000000000   1.5000000000000000   1.8535533905932737   1.5000000000000000 )
           zd(ii)=snapmatrix_velocity(27,u_nodesi,snapshotsi,dimi,numphasei)
         else if(xd(1,ii) .eq. 1 .and. xd(2,ii) .eq. 1.5   .and. xd(3,ii) .eq. 1.5 .and. xd(4,ii) .eq. 1) then !(1.0000000000000000   1.5000000000000000   1.5000000000000000   1.0000000000000000)
           zd(ii)=snapmatrix_velocity(28,u_nodesi,snapshotsi,dimi,numphasei)
         else if(xd(1,ii) .eq. 2 .and. xd(2,ii) .eq. 1.5   .and. xd(3,ii) .eq. 1.5 .and. xd(4,ii) .eq. 1) then !(2.0000000000000000   1.5000000000000000   1.5000000000000000   1.0000000000000000)
           zd(ii)=snapmatrix_velocity(29,u_nodesi,snapshotsi,dimi,numphasei)
           else if(xd(1,ii) .eq. 1 .and. xd(2,ii) .eq. 1.5   .and. xd(3,ii) .eq. 1.5 .and. xd(4,ii) .eq. 2) then !(1.0000000000000000   1.5000000000000000   1.5000000000000000   2.0000000000000000)
           zd(ii)=snapmatrix_velocity(30,u_nodesi,snapshotsi,dimi,numphasei)
         else if(xd(1,ii) .eq. 2 .and. xd(2,ii) .eq. 1.5   .and. xd(3,ii) .eq. 1.5 .and. xd(4,ii) .eq. 2) then !(2.0000000000000000   1.5000000000000000   1.5000000000000000   2.0000000000000000)
           zd(ii)=snapmatrix_velocity(31,u_nodesi,snapshotsi,dimi,numphasei)
         else if(xd(1,ii) .eq. 1.5 .and. xd(2,ii) .eq. 1   .and. xd(3,ii) .eq. 1.5 .and. xd(4,ii) .eq. 1) then !(1.5000000000000000   1.0000000000000000   1.5000000000000000   1.0000000000000000)
           zd(ii)=snapmatrix_velocity(32,u_nodesi,snapshotsi,dimi,numphasei)
         else if(xd(1,ii) .eq. 1.5 .and. xd(2,ii) .eq. 2  .and. xd(3,ii) .eq. 1.5 .and. xd(4,ii) .eq. 1) then !(1.5000000000000000   2.0000000000000000   1.5000000000000000   1.0000000000000000)
           zd(ii)=snapmatrix_velocity(33,u_nodesi,snapshotsi,dimi,numphasei) 
         else if(xd(1,ii) .eq. 1.5 .and. xd(2,ii) .eq. 1   .and. xd(3,ii) .eq. 1.5 .and. xd(4,ii) .eq. 2) then !(1.5000000000000000   1.0000000000000000   1.5000000000000000   2.0000000000000000)
           zd(ii)=snapmatrix_velocity(34,u_nodesi,snapshotsi,dimi,numphasei)
         else if(xd(1,ii) .eq. 1.5 .and. xd(2,ii) .eq. 2   .and. xd(3,ii) .eq. 1.5 .and. xd(4,ii) .eq. 2) then !(1.5000000000000000   2.0000000000000000   1.5000000000000000   2.0000000000000000)
           zd(ii)=snapmatrix_velocity(35,u_nodesi,snapshotsi,dimi,numphasei)
         else if(xd(1,ii) .eq. 1.5 .and. xd(2,ii) .eq. 1.5   .and. xd(3,ii) .eq. 1 .and. xd(4,ii) .eq. 1) then !(1.5000000000000000   1.5000000000000000   1.0000000000000000   1.0000000000000000)
           zd(ii)=snapmatrix_velocity(36,u_nodesi,snapshotsi,dimi,numphasei)
         else if(xd(1,ii) .eq. 1.5 .and. xd(2,ii) .eq. 1.5   .and. xd(3,ii) .eq. 2 .and. xd(4,ii) .eq. 1) then !( 1.5000000000000000   1.5000000000000000   2.0000000000000000   1.0000000000000000)
           zd(ii)=snapmatrix_velocity(37,u_nodesi,snapshotsi,dimi,numphasei)
          else if(xd(1,ii) .eq. 1.5 .and. xd(2,ii) .eq. 1.5   .and. xd(3,ii) .eq. 1 .and. xd(4,ii) .eq. 2) then !(1.5000000000000000   1.5000000000000000   1.0000000000000000   2.0000000000000000)
           zd(ii)=snapmatrix_velocity(38,u_nodesi,snapshotsi,dimi,numphasei)
         else if(xd(1,ii) .eq. 1.5 .and. xd(2,ii) .eq. 1.5   .and. xd(3,ii) .eq. 2 .and. xd(4,ii) .eq. 2) then !(1.5000000000000000   1.5000000000000000   2.0000000000000000   2.0000000000000000)
           zd(ii)=snapmatrix_velocity(39,u_nodesi,snapshotsi,dimi,numphasei)
         else if(xd(1,ii) .eq. 1.5 .and. xd(2,ii) .eq. 1.5 .and.xd(3,ii) .eq. 1.5 .and.xd(4,ii) >= 1.14 .and. xd(4,ii) <= 1.15) then 
      !(1.5000000000000000   1.5000000000000000   1.5000000000000000   1.1464466094067263 )
           zd(ii)=snapmatrix_velocity(40,u_nodesi,snapshotsi,dimi,numphasei)
         else if(xd(1,ii) .eq. 1.5.and. xd(2,ii) .eq. 1.5 .and.xd(3,ii) .eq. 1.5.and. xd(4,ii) >= 1.85 .and. xd(4,ii) <= 1.855) then 
     !(1.5000000000000000   1.5000000000000000   1.5000000000000000   1.8535533905932737     )
           zd(ii)=snapmatrix_velocity(41,u_nodesi,snapshotsi,dimi,numphasei)
         endif
        enddo
          !call f_sinr ( m, nd, xd, zd )       
           
          ! print *, 'zd',  zd(:)
            
!  Use the grid to evaluate the interpolant.
        call lagrange_interp_nd_value2 ( m, ind, a, b, nd, zd, ni, xi, zpi )
 
!  Weighted the interpolant values and add to the sparse grid interpolant.
 
        nd_total = nd_total + nd
        zi(1:ni) = zi(1:ni) + c(l) * zpi(1:ni)
        !  print *, ' nd_total,level',  nd_total, l, c(l), nd
        !  print *, 'zizizizizizizic(l)zpi', zi(1)!,c(l),zpi(1)
        deallocate ( xd )
        deallocate ( zd )     
       
       if ( .not. more ) then
          exit
        end if    
 
      end do !DO

    end do !do l = l_min, l_max   
    deallocate ( c )
    deallocate ( w )
    
   end do ! do l_max = sparse_max, sparse_max
   interpolate_velocity(u_nodesi,snapshotsi,dimi,numphasei)=zi(ni)
   enddo !Do u_nodesi  
  enddo !Do snapshotsi
  enddo !dimi
  enddo !numphasei

 !  print *, ' total', total , nd_total
 ! open(unit=66,file='podbasis_interp') 
 !  do j=1, total_timestep 
 !     write(66,*)  pod_coef_all_obv_41i(j,1:num_pod)
 ! enddo
 ! close(66)
 
  deallocate ( a )
  deallocate ( b )
  deallocate ( ind )
  deallocate ( xi )
  deallocate ( ze )
  deallocate ( zi )
  deallocate ( zpi )
  end subroutine interpolate_snapmatrix_velocity

 subroutine interpolate_snapmatrix_pressure( snapmatrix_pressure, interpolate_pressure)  
    
! top-bottom between 1-2, pod_coef_all_obv is used for storing some podbasis to avoid to much modifications.
 
!  Parameters:
!
!    Local, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) SPARSE_MAX, the maximum sparse grid level to try.
!
!  Local Parameters:
!
!    Local, real ( kind = 8 ) A(M), B(M), the upper and lower variable limits 
!    in each dimension.
!
!    Local, real ( kind = 8 ) APP_ERROR, the averaged Euclidean norm of the 
!    difference between the sparse interpolant and the exact function at 
!    the interpolation points.
!
!    Local, integer ( kind = 4 ) C(0:L_MAX), the sparse grid coefficient vector.
!    Results at level L are weighted by C(L).
!
!    Local, integer ( kind = 4 ) IND(M), the 1D indices defining a Lagrange grid.
!    Each entry is a 1d "level" that specifies the order of a 
!    Clenshaw Curtis 1D grid.
!
!    Local, integer ( kind = 4 ) L, the current Lagrange grid level.
!
!    Local, integer ( kind = 4 ) L_MAX, the current sparse grid level.
!
!    Local, integer ( kind = 4 ) MORE, used to control the enumeration of all the
!    Lagrange grids at a current grid level.
!
!    Local, integer ( kind = 4 ) ND, the number of points used in a Lagrange grid.
!
!    Local, integer ( kind = 4 ) ND_TOTAL, the total number of points used in all the
!    Lagrange interpolants for a given sparse interpolant points that occur
!    in more than one Lagrange grid are counted multiple times.
!
!    Local, integer ( kind = 4 ) NI, the number of interpolant evaluation points.
!
!    Local, integer ( kind = 4 ) SPARSE_MIN, the minimum sparse grid level to try.
!
!    Local, real ( kind = 8 ) XD(M,ND), the data points for a Lagrange grid.
!
!    Local, real ( kind = 8 ) XI(M,NI), the interpolant evaluation points.
!
!    Local, real ( kind = 8 ) ZD(ND), the data values for a Lagrange grid.
!
!    Local, real ( kind = 8 ) ZE(NI), the exact function values at XI.
!
!    Local, real ( kind = 8 ) ZI(NI), the sparse interpolant values at XI.
!
!    Local, real ( kind = 8 ) ZPI(NI), one set of Lagrange interpolant values at XI.
!
  implicit none
  !type(state_type), intent(in), dimension(:,:,:) :: state
  real, dimension(:,:,:,:), allocatable :: snapmatrix_pressure ! allocate(snapmatrix_pressure(no_smolyak_nodes,p_nodes,snapshots,numphase))
  real, dimension(:,:,:), intent(out), allocatable ::  interpolate_pressure 
  real ( kind = 8 ), allocatable :: a(:)
  real ( kind = 8 ) app_error
  real ( kind = 8 ), allocatable :: b(:)
  integer ( kind = 4 ), allocatable :: c(:)
  integer ( kind = 4 ) h
  integer ( kind = 4 ), allocatable :: ind(:)
  integer ( kind = 4 ) l
  integer ( kind = 4 ) l_max
  integer ( kind = 4 ) l_min
  !integer ( kind = 4 ) m
  logical more
  integer ( kind = 4 ) nd
  integer ( kind = 4 ) nd_total
  integer ( kind = 4 ) ni
  real ( kind = 8 ) r8vec_norm_affine
  integer ( kind = 4 ) seed
 ! integer ( kind = 4 ) sparse_max
  integer ( kind = 4 ) sparse_min
  integer ( kind = 4 ) t
  integer ( kind = 4 ), allocatable :: w(:)
  real ( kind = 8 ), allocatable :: xd(:,:)
  real ( kind = 8 ), allocatable :: xi(:,:)
  real ( kind = 8 ), allocatable :: zd(:)
  real ( kind = 8 ), allocatable :: ze(:)
  real ( kind = 8 ), allocatable :: zi(:)
  real ( kind = 8 ), allocatable :: zpi(:)
  integer :: i, ii, total, total_timestep,num_pod,j,k,jj,kk, smolyak_total_points,u_nodesi,snapshotsi,dimi,numphasei
   
  allocate (ind(1:m))

  ! Define the region.
  allocate ( a(1:m) )
  allocate ( b(1:m) )

  a(1:m) = 1
  b(1:m) = 2 
   
  !total_timestep=864
  !num_pod=12*36
  !smolyak_total_points=41!2**sparse_max+1
  !print *, '  smolyak_total_points=2**sparse_max+1=', smolyak_total_points 
  !  stop 2323
  !  Define the interpolation evaluation information.
  ni = 1
  seed = 123456789
  allocate ( xi(m,ni) )
  !call r8mat_uniform_abvec ( m, ni, a, b, seed, xi )
!  print * , 'xi', xi
 
  xi(1,1)= 1 
  xi(2,1)= 2
  xi(3,1)= 2
  xi(4,1)=2!1.7 !2.0000000000000000   1.5000000000000000   2.0000000000000000   1.5000000000000000,21
    call get_option(&
         '/reduced_model/variables_to_interpolate/dimension_one', xi(1,1))
   call get_option(&
         '/reduced_model/variables_to_interpolate/dimension_two', xi(2,1))
 call get_option(&
         '/reduced_model/variables_to_interpolate/dimension_three', xi(3,1))
 call get_option(&
         '/reduced_model/variables_to_interpolate/dimension_four', xi(4,1))
  
  write ( *, '(a,i6)' ) '  Number of interpolation points is NI = ', ni
  allocate ( ze(1:ni) )
 ! call f_sinr ( m, ni, xi, ze )
  print *, 'actual function value is, will interpolate this value later on,then do comparison', xi,ze
!
!  Compute a sequence of sparse grid interpolants of increasing level.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '   L     Order    ApproxError'
  write ( *, '(a)' ) ''
  allocate(interpolate_pressure(p_nodes,snapshots,numphase)) 
  allocate ( zpi(1:ni) )
  allocate ( zi(1:ni) )
  !smolynode, u_nodes,snapshots,dim,numphase
  Do u_nodesi=1, p_nodes  
   Do snapshotsi=1, snapshots  
   do numphasei=1, numphase
  total = 0
  sparse_min = 0

  ! do l_max = sparse_min, sparse_max
    do l_max = sparse_max, sparse_max
    allocate ( c(0:l_max) )
    allocate ( w(0:l_max) )
    call smolyak_coefficients ( l_max, m, c, w )

    zi(1:ni) = 0.0D+00
    nd_total = 0
    
    l_min = max ( l_max + 1 - m, 0 )

    do l = l_min, l_max

      more = .false.
      
      do
     ! print *,'l_max,l,' ,l_max,l   
!  Define the next product grid.
!
        call comp_next ( l, m, ind, more, h, t )
!
!  Count the grid, find its points, evaluate the data there.
!
        call lagrange_interp_nd_size2 ( m, ind, nd )
        allocate ( xd(1:m,1:nd) )
        call lagrange_interp_nd_grid2 ( m, ind, a, b, nd, xd )        
         !  print *,'xd', xd, 点1 坐标为： xd(:,1)
        allocate ( zd(1:nd) )

       ! print *,'---------------------nd',nd  
        do ii=1, nd           
              if(xd(1,ii) .eq. 1.5 .and. xd(2,ii) .eq. 1.5   .and. xd(3,ii) .eq. 1.5 .and. xd(4,ii) .eq. 1.5) then 
!( 1.5000000000000000    1.5000000000000000   1.5000000000000000  1.5000000000000000)
           zd(ii)=snapmatrix_pressure(1,u_nodesi,snapshotsi,numphasei)
         else if(xd(1,ii) .eq. 1 .and. xd(2,ii) .eq. 1.5   .and. xd(3,ii) .eq. 1.5 .and. xd(4,ii) .eq. 1.5) then 
!(1.0000000000000000   1.5000000000000000   1.5000000000000000   1.5000000000000000)
           zd(ii)=snapmatrix_pressure(2,u_nodesi,snapshotsi,numphasei)
         else if(xd(1,ii) .eq. 2 .and. xd(2,ii) .eq. 1.5   .and. xd(3,ii) .eq. 1.5 .and. xd(4,ii) .eq. 1.5) then
 !( 2.0000000000000000   1.5000000000000000   1.5000000000000000   1.5000000000000000)
           zd(ii)=snapmatrix_pressure(3,u_nodesi,snapshotsi,numphasei)
         else if(xd(1,ii) .eq. 1.5 .and. xd(2,ii) .eq. 1   .and. xd(3,ii) .eq. 1.5 .and. xd(4,ii) .eq. 1.5) then 
!(1.5000000000000000   1.0000000000000000   1.5000000000000000   1.5000000000000000 )
           zd(ii)=snapmatrix_pressure(4,u_nodesi,snapshotsi,numphasei)
         else if(xd(1,ii) .eq. 1.5 .and. xd(2,ii) .eq. 2   .and. xd(3,ii) .eq. 1.5 .and. xd(4,ii) .eq. 1.5) then 
!(1.5000000000000000   2.0000000000000000   1.5000000000000000   1.5000000000000000 )
           zd(ii)=snapmatrix_pressure(5,u_nodesi,snapshotsi,numphasei)
        else if(xd(1,ii) .eq. 1.5 .and. xd(2,ii) .eq. 1.5   .and. xd(3,ii) .eq. 1 .and. xd(4,ii) .eq. 1.5) then 
!(1.5000000000000000   1.5000000000000000   1.0000000000000000   1.5000000000000000  )
           zd(ii)=snapmatrix_pressure(6,u_nodesi,snapshotsi,numphasei)
         else if(xd(1,ii) .eq. 1.5 .and. xd(2,ii) .eq. 1.5   .and. xd(3,ii) .eq. 2 .and. xd(4,ii) .eq. 1.5) then 
!(1.5000000000000000   1.5000000000000000   2.0000000000000000   1.5000000000000000  )
           zd(ii)=snapmatrix_pressure(7,u_nodesi,snapshotsi,numphasei)
         else if(xd(1,ii) .eq. 1.5 .and. xd(2,ii) .eq. 1.5   .and. xd(3,ii) .eq. 1.5 .and. xd(4,ii) .eq. 1) then
 !(1.5000000000000000   1.5000000000000000   1.5000000000000000   1.0000000000000000   )
           zd(ii)=snapmatrix_pressure(8,u_nodesi,snapshotsi,numphasei)
         else if(xd(1,ii) .eq. 1.5 .and. xd(2,ii) .eq. 1.5   .and. xd(3,ii) .eq. 1.5 .and. xd(4,ii) .eq. 2) then 
!(  1.5000000000000000   1.5000000000000000   1.5000000000000000   2.0000000000000000  )
           zd(ii)=snapmatrix_pressure(9,u_nodesi,snapshotsi,numphasei) 
         else if(int(xd(1,ii)*1000)/1000 .eq. 1.146 .and. xd(2,ii) .eq. 1.5 .and. xd(3,ii) .eq. 1.5 .and. xd(4,ii) .eq. 1.5) then 
               !(1.1464466094067263   1.5000000000000000   1.5000000000000000   1.5000000000000000 .and. xd(1,ii) <= 1.15 ..and. xd(2:4,ii) .eq. 1.5  )
           zd(ii)=snapmatrix_pressure(10,u_nodesi,snapshotsi,numphasei)
         elseif(int(xd(1,ii)*1000)/1000 .eq. 1.853 .and. xd(2,ii) .eq. 1.5 .and. xd(3,ii) .eq. 1.5 .and. xd(4,ii) .eq. 1.5) then 
               !(1.8535533905932737   1.5000000000000000   1.5000000000000000   1.5000000000000000.and. xd(3,ii) .eq. 1.5 .and. xd(4,ii) .eq. 1.5 )
           zd(ii)=snapmatrix_pressure(11,u_nodesi,snapshotsi,numphasei)
         else if(xd(1,ii) .eq. 1 .and. xd(2,ii) .eq. 1   .and. xd(3,ii) .eq. 1.5 .and. xd(4,ii) .eq. 1.5) then !(  1.0000000000000000   1.0000000000000000   1.5000000000000000   1.5000000000000000 )
           zd(ii)=snapmatrix_pressure(12,u_nodesi,snapshotsi,numphasei)
         else if(xd(1,ii) .eq. 2 .and. xd(2,ii) .eq. 1   .and. xd(3,ii) .eq. 1.5 .and. xd(4,ii) .eq. 1.5) then !( 2.0000000000000000   1.0000000000000000   1.5000000000000000   1.5000000000000000)
           zd(ii)=snapmatrix_pressure(13,u_nodesi,snapshotsi,numphasei)
          else if(xd(1,ii) .eq. 1 .and. xd(2,ii) .eq. 2   .and. xd(3,ii) .eq. 1.5 .and. xd(4,ii) .eq. 1.5) then !(1.0000000000000000   2.0000000000000000   1.5000000000000000   1.5000000000000000)
           zd(ii)=snapmatrix_pressure(14,u_nodesi,snapshotsi,numphasei)
         else if(xd(1,ii) .eq. 2 .and. xd(2,ii) .eq. 2   .and. xd(3,ii) .eq. 1.5 .and. xd(4,ii) .eq. 1.5) then !(2.0000000000000000   2.0000000000000000   1.5000000000000000   1.5000000000000000)
           zd(ii)=snapmatrix_pressure(15,u_nodesi,snapshotsi,numphasei)
         else if(xd(1,ii) .eq. 1.5 .and. int(xd(2,ii)*1000)/1000 .eq. 1.146 .and. xd(3,ii) .eq. 1.5 .and. xd(4,ii) .eq. 1.5) then 
            !(1.5000000000000000   1.1464466094067263   1.5000000000000000   1.5000000000000000 xd(2,ii) >=1.14  .and. xd(2,ii) <=1.15)
           zd(ii)=snapmatrix_pressure(16,u_nodesi,snapshotsi,numphasei)
         else if(xd(1,ii) .eq. 1.5 .and. int(xd(2,ii)*1000)/1000 .eq. 1.853 .and. xd(3,ii) .eq. 1.5 .and. xd(4,ii) .eq. 1.5) then
          !(1.5000000000000000   1.8535533905932737   1.5000000000000000   1.5000000000000000)
           zd(ii)=snapmatrix_pressure(17,u_nodesi,snapshotsi,numphasei)
           else if(xd(1,ii) .eq. 1 .and. xd(2,ii) .eq. 1.5 .and. xd(3,ii) .eq. 1 .and. xd(4,ii) .eq. 1.5) then !( 1.0000000000000000   1.5000000000000000   1.0000000000000000  1.5000000000000000)
           zd(ii)=snapmatrix_pressure(18,u_nodesi,snapshotsi,numphasei)
         else if(xd(1,ii) .eq. 2 .and. xd(2,ii) .eq. 1.5   .and. xd(3,ii) .eq. 1 .and. xd(4,ii) .eq. 1.5) then !(2.0000000000000000   1.5000000000000000   1.0000000000000000   1.5000000000000000)
           zd(ii)=snapmatrix_pressure(19,u_nodesi,snapshotsi,numphasei)
         else if(xd(1,ii) .eq. 1 .and. xd(2,ii) .eq. 1.5   .and. xd(3,ii) .eq. 2 .and. xd(4,ii) .eq. 1.5) then !(1.0000000000000000   1.5000000000000000   2.0000000000000000   1.5000000000000000)
           zd(ii)=snapmatrix_pressure(20,u_nodesi,snapshotsi,numphasei)
         else if(xd(1,ii) .eq. 2 .and. xd(2,ii) .eq. 1.5   .and. xd(3,ii) .eq. 2 .and. xd(4,ii) .eq. 1.5) then !(2.0000000000000000   1.5000000000000000   2.0000000000000000   1.5000000000000000)
           zd(ii)=snapmatrix_pressure(21,u_nodesi,snapshotsi,numphasei) 
         else if(xd(1,ii) .eq. 1.5 .and. xd(2,ii) .eq. 1   .and. xd(3,ii) .eq. 1 .and. xd(4,ii) .eq. 1.5) then !(1.5000000000000000   1.0000000000000000   1.0000000000000000   1.5000000000000000)
           zd(ii)=snapmatrix_pressure(22,u_nodesi,snapshotsi,numphasei)
         else if(xd(1,ii) .eq. 1.5 .and. xd(2,ii) .eq. 2   .and. xd(3,ii) .eq. 1 .and. xd(4,ii) .eq. 1.5) then !(1.5000000000000000   2.0000000000000000   1.0000000000000000   1.5000000000000000)
           zd(ii)=snapmatrix_pressure(23,u_nodesi,snapshotsi,numphasei)
         else if(xd(1,ii) .eq. 1.5 .and. xd(2,ii) .eq. 1   .and. xd(3,ii) .eq. 2 .and. xd(4,ii) .eq. 1.5) then !(1.5000000000000000   1.0000000000000000   2.0000000000000000   1.5000000000000000)
           zd(ii)=snapmatrix_pressure(24,u_nodesi,snapshotsi,numphasei)
         else if(xd(1,ii) .eq. 1.5 .and. xd(2,ii) .eq. 2   .and. xd(3,ii) .eq. 2 .and. xd(4,ii) .eq. 1.5) then !(1.5000000000000000   2.0000000000000000   2.0000000000000000   1.5000000000000000)
           zd(ii)=snapmatrix_pressure(25,u_nodesi,snapshotsi,numphasei)
          else if(xd(1,ii) .eq. 1.5 .and. xd(2,ii) .eq. 1.5 .and. int(xd(3,ii)*1000)/1000 .eq. 1.146 .and. xd(4,ii) .eq. 1.5) then 
       !(1.5000000000000000   1.5000000000000000  1.1464466094067263   1.5000000000000000)
           zd(ii)=snapmatrix_pressure(26,u_nodesi,snapshotsi,numphasei)
         else if(xd(1,ii) .eq. 1.5  .and. xd(2,ii) .eq. 1.5  .and. int(xd(3,ii)*1000)/1000 .eq. 1.853 .and. xd(4,ii) .eq. 1.5) then 
       !( 1.5000000000000000   1.5000000000000000   1.8535533905932737   1.5000000000000000 )
           zd(ii)=snapmatrix_pressure(27,u_nodesi,snapshotsi,numphasei)
         else if(xd(1,ii) .eq. 1 .and. xd(2,ii) .eq. 1.5   .and. xd(3,ii) .eq. 1.5 .and. xd(4,ii) .eq. 1) then !(1.0000000000000000   1.5000000000000000   1.5000000000000000   1.0000000000000000)
           zd(ii)=snapmatrix_pressure(28,u_nodesi,snapshotsi,numphasei)
         else if(xd(1,ii) .eq. 2 .and. xd(2,ii) .eq. 1.5   .and. xd(3,ii) .eq. 1.5 .and. xd(4,ii) .eq. 1) then !(2.0000000000000000   1.5000000000000000   1.5000000000000000   1.0000000000000000)
           zd(ii)=snapmatrix_pressure(29,u_nodesi,snapshotsi,numphasei)
           else if(xd(1,ii) .eq. 1 .and. xd(2,ii) .eq. 1.5   .and. xd(3,ii) .eq. 1.5 .and. xd(4,ii) .eq. 2) then !(1.0000000000000000   1.5000000000000000   1.5000000000000000   2.0000000000000000)
           zd(ii)=snapmatrix_pressure(30,u_nodesi,snapshotsi,numphasei)
         else if(xd(1,ii) .eq. 2 .and. xd(2,ii) .eq. 1.5   .and. xd(3,ii) .eq. 1.5 .and. xd(4,ii) .eq. 2) then !(2.0000000000000000   1.5000000000000000   1.5000000000000000   2.0000000000000000)
           zd(ii)=snapmatrix_pressure(31,u_nodesi,snapshotsi,numphasei)
         else if(xd(1,ii) .eq. 1.5 .and. xd(2,ii) .eq. 1   .and. xd(3,ii) .eq. 1.5 .and. xd(4,ii) .eq. 1) then !(1.5000000000000000   1.0000000000000000   1.5000000000000000   1.0000000000000000)
           zd(ii)=snapmatrix_pressure(32,u_nodesi,snapshotsi,numphasei)
         else if(xd(1,ii) .eq. 1.5 .and. xd(2,ii) .eq. 2  .and. xd(3,ii) .eq. 1.5 .and. xd(4,ii) .eq. 1) then !(1.5000000000000000   2.0000000000000000   1.5000000000000000   1.0000000000000000)
           zd(ii)=snapmatrix_pressure(33,u_nodesi,snapshotsi,numphasei) 
         else if(xd(1,ii) .eq. 1.5 .and. xd(2,ii) .eq. 1   .and. xd(3,ii) .eq. 1.5 .and. xd(4,ii) .eq. 2) then !(1.5000000000000000   1.0000000000000000   1.5000000000000000   2.0000000000000000)
           zd(ii)=snapmatrix_pressure(34,u_nodesi,snapshotsi,numphasei)
         else if(xd(1,ii) .eq. 1.5 .and. xd(2,ii) .eq. 2   .and. xd(3,ii) .eq. 1.5 .and. xd(4,ii) .eq. 2) then !(1.5000000000000000   2.0000000000000000   1.5000000000000000   2.0000000000000000)
           zd(ii)=snapmatrix_pressure(35,u_nodesi,snapshotsi,numphasei)
         else if(xd(1,ii) .eq. 1.5 .and. xd(2,ii) .eq. 1.5   .and. xd(3,ii) .eq. 1 .and. xd(4,ii) .eq. 1) then !(1.5000000000000000   1.5000000000000000   1.0000000000000000   1.0000000000000000)
           zd(ii)=snapmatrix_pressure(36,u_nodesi,snapshotsi,numphasei)
         else if(xd(1,ii) .eq. 1.5 .and. xd(2,ii) .eq. 1.5   .and. xd(3,ii) .eq. 2 .and. xd(4,ii) .eq. 1) then !( 1.5000000000000000   1.5000000000000000   2.0000000000000000   1.0000000000000000)
           zd(ii)=snapmatrix_pressure(37,u_nodesi,snapshotsi,numphasei)
          else if(xd(1,ii) .eq. 1.5 .and. xd(2,ii) .eq. 1.5   .and. xd(3,ii) .eq. 1 .and. xd(4,ii) .eq. 2) then !(1.5000000000000000   1.5000000000000000   1.0000000000000000   2.0000000000000000)
           zd(ii)=snapmatrix_pressure(38,u_nodesi,snapshotsi,numphasei)
         else if(xd(1,ii) .eq. 1.5 .and. xd(2,ii) .eq. 1.5   .and. xd(3,ii) .eq. 2 .and. xd(4,ii) .eq. 2) then !(1.5000000000000000   1.5000000000000000   2.0000000000000000   2.0000000000000000)
           zd(ii)=snapmatrix_pressure(39,u_nodesi,snapshotsi,numphasei)
         else if(xd(1,ii) .eq. 1.5 .and. xd(2,ii) .eq. 1.5 .and.xd(3,ii) .eq. 1.5 .and.xd(4,ii) >= 1.14 .and. xd(4,ii) <= 1.15) then 
      !(1.5000000000000000   1.5000000000000000   1.5000000000000000   1.1464466094067263 )
           zd(ii)=snapmatrix_pressure(40,u_nodesi,snapshotsi,numphasei)
         else if(xd(1,ii) .eq. 1.5.and. xd(2,ii) .eq. 1.5 .and.xd(3,ii) .eq. 1.5.and. xd(4,ii) >= 1.85 .and. xd(4,ii) <= 1.855) then 
     !(1.5000000000000000   1.5000000000000000   1.5000000000000000   1.8535533905932737     )
           zd(ii)=snapmatrix_pressure(41,u_nodesi,snapshotsi,numphasei)
         endif
        enddo
          !call f_sinr ( m, nd, xd, zd )       
           
          ! print *, 'zd',  zd(:)
            
!  Use the grid to evaluate the interpolant.
        call lagrange_interp_nd_value2 ( m, ind, a, b, nd, zd, ni, xi, zpi )
 
!  Weighted the interpolant values and add to the sparse grid interpolant.
 
        nd_total = nd_total + nd
        zi(1:ni) = zi(1:ni) + c(l) * zpi(1:ni)
        !  print *, ' nd_total,level',  nd_total, l, c(l), nd
        !  print *, 'zizizizizizizic(l)zpi', zi(1)!,c(l),zpi(1)
        deallocate ( xd )
        deallocate ( zd )     
       
       if ( .not. more ) then
          exit
        end if    
 
      end do !DO

    end do !do l = l_min, l_max   
    deallocate ( c )
    deallocate ( w )
    
   end do ! do l_max = sparse_max, sparse_max
   interpolate_pressure(u_nodesi,snapshotsi,numphasei)=zi(ni)
   enddo !Do u_nodesi  
  enddo !Do snapshotsi 
  enddo !numphasei

 !  print *, ' total', total , nd_total
 ! open(unit=66,file='podbasis_interp') 
 !  do j=1, total_timestep 
 !     write(66,*)  pod_coef_all_obv_41i(j,1:num_pod)
 ! enddo
 ! close(66)
 
  deallocate ( a )
  deallocate ( b )
  deallocate ( ind )
  deallocate ( xi )
  deallocate ( ze )
  deallocate ( zi )
  deallocate ( zpi )
  end subroutine interpolate_snapmatrix_pressure

 subroutine interpolate_snapmatrix_volumefraction( snapmatrix_volumefraction, interpolate_volumefraction)  
    
! top-bottom between 1-2, pod_coef_all_obv is used for storing some podbasis to avoid to much modifications.
 
!  Parameters:
!
!    Local, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) SPARSE_MAX, the maximum sparse grid level to try.
!
!  Local Parameters:
!
!    Local, real ( kind = 8 ) A(M), B(M), the upper and lower variable limits 
!    in each dimension.
!
!    Local, real ( kind = 8 ) APP_ERROR, the averaged Euclidean norm of the 
!    difference between the sparse interpolant and the exact function at 
!    the interpolation points.
!
!    Local, integer ( kind = 4 ) C(0:L_MAX), the sparse grid coefficient vector.
!    Results at level L are weighted by C(L).
!
!    Local, integer ( kind = 4 ) IND(M), the 1D indices defining a Lagrange grid.
!    Each entry is a 1d "level" that specifies the order of a 
!    Clenshaw Curtis 1D grid.
!
!    Local, integer ( kind = 4 ) L, the current Lagrange grid level.
!
!    Local, integer ( kind = 4 ) L_MAX, the current sparse grid level.
!
!    Local, integer ( kind = 4 ) MORE, used to control the enumeration of all the
!    Lagrange grids at a current grid level.
!
!    Local, integer ( kind = 4 ) ND, the number of points used in a Lagrange grid.
!
!    Local, integer ( kind = 4 ) ND_TOTAL, the total number of points used in all the
!    Lagrange interpolants for a given sparse interpolant points that occur
!    in more than one Lagrange grid are counted multiple times.
!
!    Local, integer ( kind = 4 ) NI, the number of interpolant evaluation points.
!
!    Local, integer ( kind = 4 ) SPARSE_MIN, the minimum sparse grid level to try.
!
!    Local, real ( kind = 8 ) XD(M,ND), the data points for a Lagrange grid.
!
!    Local, real ( kind = 8 ) XI(M,NI), the interpolant evaluation points.
!
!    Local, real ( kind = 8 ) ZD(ND), the data values for a Lagrange grid.
!
!    Local, real ( kind = 8 ) ZE(NI), the exact function values at XI.
!
!    Local, real ( kind = 8 ) ZI(NI), the sparse interpolant values at XI.
!
!    Local, real ( kind = 8 ) ZPI(NI), one set of Lagrange interpolant values at XI.
!
  implicit none
  !type(state_type), intent(in), dimension(:,:,:) :: state
  real, dimension(:,:,:,:), allocatable :: snapmatrix_volumefraction ! allocate(snapmatrix_pressure(no_smolyak_nodes,p_nodes,snapshots,numphase))
  real, dimension(:,:,:), intent(out), allocatable ::  interpolate_volumefraction 
  real ( kind = 8 ), allocatable :: a(:)
  real ( kind = 8 ) app_error
  real ( kind = 8 ), allocatable :: b(:)
  integer ( kind = 4 ), allocatable :: c(:)
  integer ( kind = 4 ) h
  integer ( kind = 4 ), allocatable :: ind(:)
  integer ( kind = 4 ) l
  integer ( kind = 4 ) l_max
  integer ( kind = 4 ) l_min
  !integer ( kind = 4 ) m
  logical more
  integer ( kind = 4 ) nd
  integer ( kind = 4 ) nd_total
  integer ( kind = 4 ) ni
  real ( kind = 8 ) r8vec_norm_affine
  integer ( kind = 4 ) seed
  !integer ( kind = 4 ) sparse_max
  integer ( kind = 4 ) sparse_min
  integer ( kind = 4 ) t
  integer ( kind = 4 ), allocatable :: w(:)
  real ( kind = 8 ), allocatable :: xd(:,:)
  real ( kind = 8 ), allocatable :: xi(:,:)
  real ( kind = 8 ), allocatable :: zd(:)
  real ( kind = 8 ), allocatable :: ze(:)
  real ( kind = 8 ), allocatable :: zi(:)
  real ( kind = 8 ), allocatable :: zpi(:)
  integer :: i, ii, total, total_timestep,num_pod,j,k,jj,kk, smolyak_total_points,u_nodesi,snapshotsi,dimi,numphasei
   
  allocate (ind(1:m))

  ! Define the region.
  allocate ( a(1:m) )
  allocate ( b(1:m) )

  a(1:m) = 1
  b(1:m) = 2 
   
  !total_timestep=864
  !num_pod=12*36
  !smolyak_total_points=41!2**sparse_max+1
  !print *, '  smolyak_total_points=2**sparse_max+1=', smolyak_total_points 
  
  !  Define the interpolation evaluation information.
  ni = 1
  seed = 123456789
  allocate ( xi(m,ni) )
  !call r8mat_uniform_abvec ( m, ni, a, b, seed, xi )
!  print * , 'xi', xi
 
  xi(1,1)= 1 
  xi(2,1)= 2
  xi(3,1)= 2
  xi(4,1)= 2!1.7 !2.0000000000000000   1.5000000000000000   2.0000000000000000   1.5500000000000000 
    call get_option(&
         '/reduced_model/variables_to_interpolate/dimension_one', xi(1,1))
   call get_option(&
         '/reduced_model/variables_to_interpolate/dimension_two', xi(2,1))
 call get_option(&
         '/reduced_model/variables_to_interpolate/dimension_three', xi(3,1))
 call get_option(&
         '/reduced_model/variables_to_interpolate/dimension_four', xi(4,1))
   
  write ( *, '(a,i6)' ) '  Number of interpolation points is NI = ', ni
  allocate ( ze(1:ni) )
 ! call f_sinr ( m, ni, xi, ze )
  print *, 'actual function value is, will interpolate this value later on,then do comparison', xi,ze
!
!  Compute a sequence of sparse grid interpolants of increasing level.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '   L     Order    ApproxError'
  write ( *, '(a)' ) ''

  allocate ( zpi(1:ni) )
  allocate ( zi(1:ni) )
  allocate(interpolate_volumefraction(p_nodes,snapshots,numphase))
  !smolynode, u_nodes,snapshots,dim,numphase
  Do u_nodesi=1, p_nodes  
   Do snapshotsi=1, snapshots  
   do numphasei=1, numphase
  total = 0
  sparse_min = 0

  ! do l_max = sparse_min, sparse_max
    do l_max = sparse_max, sparse_max
    allocate ( c(0:l_max) )
    allocate ( w(0:l_max) )
    call smolyak_coefficients ( l_max, m, c, w )

    zi(1:ni) = 0.0D+00
    nd_total = 0
    
    l_min = max ( l_max + 1 - m, 0 )

    do l = l_min, l_max

      more = .false.
      
      do
     ! print *,'l_max,l,' ,l_max,l   
!  Define the next product grid.
!
        call comp_next ( l, m, ind, more, h, t )
!
!  Count the grid, find its points, evaluate the data there.
!
        call lagrange_interp_nd_size2 ( m, ind, nd )
        allocate ( xd(1:m,1:nd) )
        call lagrange_interp_nd_grid2 ( m, ind, a, b, nd, xd )        
         !  print *,'xd', xd, 点1 坐标为： xd(:,1)
        allocate ( zd(1:nd) )

       ! print *,'---------------------nd',nd  
        do ii=1, nd           
              if(xd(1,ii) .eq. 1.5 .and. xd(2,ii) .eq. 1.5   .and. xd(3,ii) .eq. 1.5 .and. xd(4,ii) .eq. 1.5) then 
!( 1.5000000000000000    1.5000000000000000   1.5000000000000000  1.5000000000000000)
           zd(ii)=snapmatrix_volumefraction(1,u_nodesi,snapshotsi,numphasei)
         else if(xd(1,ii) .eq. 1 .and. xd(2,ii) .eq. 1.5   .and. xd(3,ii) .eq. 1.5 .and. xd(4,ii) .eq. 1.5) then 
!(1.0000000000000000   1.5000000000000000   1.5000000000000000   1.5000000000000000)
           zd(ii)=snapmatrix_volumefraction(2,u_nodesi,snapshotsi,numphasei)
         else if(xd(1,ii) .eq. 2 .and. xd(2,ii) .eq. 1.5   .and. xd(3,ii) .eq. 1.5 .and. xd(4,ii) .eq. 1.5) then
 !( 2.0000000000000000   1.5000000000000000   1.5000000000000000   1.5000000000000000)
           zd(ii)=snapmatrix_volumefraction(3,u_nodesi,snapshotsi,numphasei)
         else if(xd(1,ii) .eq. 1.5 .and. xd(2,ii) .eq. 1   .and. xd(3,ii) .eq. 1.5 .and. xd(4,ii) .eq. 1.5) then 
!(1.5000000000000000   1.0000000000000000   1.5000000000000000   1.5000000000000000 )
           zd(ii)=snapmatrix_volumefraction(4,u_nodesi,snapshotsi,numphasei)
         else if(xd(1,ii) .eq. 1.5 .and. xd(2,ii) .eq. 2   .and. xd(3,ii) .eq. 1.5 .and. xd(4,ii) .eq. 1.5) then 
!(1.5000000000000000   2.0000000000000000   1.5000000000000000   1.5000000000000000 )
           zd(ii)=snapmatrix_volumefraction(5,u_nodesi,snapshotsi,numphasei)
        else if(xd(1,ii) .eq. 1.5 .and. xd(2,ii) .eq. 1.5   .and. xd(3,ii) .eq. 1 .and. xd(4,ii) .eq. 1.5) then 
!(1.5000000000000000   1.5000000000000000   1.0000000000000000   1.5000000000000000  )
           zd(ii)=snapmatrix_volumefraction(6,u_nodesi,snapshotsi,numphasei)
         else if(xd(1,ii) .eq. 1.5 .and. xd(2,ii) .eq. 1.5   .and. xd(3,ii) .eq. 2 .and. xd(4,ii) .eq. 1.5) then 
!(1.5000000000000000   1.5000000000000000   2.0000000000000000   1.5000000000000000  )
           zd(ii)=snapmatrix_volumefraction(7,u_nodesi,snapshotsi,numphasei)
         else if(xd(1,ii) .eq. 1.5 .and. xd(2,ii) .eq. 1.5   .and. xd(3,ii) .eq. 1.5 .and. xd(4,ii) .eq. 1) then
 !(1.5000000000000000   1.5000000000000000   1.5000000000000000   1.0000000000000000   )
           zd(ii)=snapmatrix_volumefraction(8,u_nodesi,snapshotsi,numphasei)
         else if(xd(1,ii) .eq. 1.5 .and. xd(2,ii) .eq. 1.5   .and. xd(3,ii) .eq. 1.5 .and. xd(4,ii) .eq. 2) then 
!(  1.5000000000000000   1.5000000000000000   1.5000000000000000   2.0000000000000000  )
           zd(ii)=snapmatrix_volumefraction(9,u_nodesi,snapshotsi,numphasei) 
         else if(int(xd(1,ii)*1000)/1000 .eq. 1.146 .and. xd(2,ii) .eq. 1.5 .and. xd(3,ii) .eq. 1.5 .and. xd(4,ii) .eq. 1.5) then 
               !(1.1464466094067263   1.5000000000000000   1.5000000000000000   1.5000000000000000 .and. xd(1,ii) <= 1.15 ..and. xd(2:4,ii) .eq. 1.5  )
           zd(ii)=snapmatrix_volumefraction(10,u_nodesi,snapshotsi,numphasei)
         elseif(int(xd(1,ii)*1000)/1000 .eq. 1.853 .and. xd(2,ii) .eq. 1.5 .and. xd(3,ii) .eq. 1.5 .and. xd(4,ii) .eq. 1.5) then 
               !(1.8535533905932737   1.5000000000000000   1.5000000000000000   1.5000000000000000.and. xd(3,ii) .eq. 1.5 .and. xd(4,ii) .eq. 1.5 )
           zd(ii)=snapmatrix_volumefraction(11,u_nodesi,snapshotsi,numphasei)
         else if(xd(1,ii) .eq. 1 .and. xd(2,ii) .eq. 1   .and. xd(3,ii) .eq. 1.5 .and. xd(4,ii) .eq. 1.5) then !(  1.0000000000000000   1.0000000000000000   1.5000000000000000   1.5000000000000000 )
           zd(ii)=snapmatrix_volumefraction(12,u_nodesi,snapshotsi,numphasei)
         else if(xd(1,ii) .eq. 2 .and. xd(2,ii) .eq. 1   .and. xd(3,ii) .eq. 1.5 .and. xd(4,ii) .eq. 1.5) then !( 2.0000000000000000   1.0000000000000000   1.5000000000000000   1.5000000000000000)
           zd(ii)=snapmatrix_volumefraction(13,u_nodesi,snapshotsi,numphasei)
          else if(xd(1,ii) .eq. 1 .and. xd(2,ii) .eq. 2   .and. xd(3,ii) .eq. 1.5 .and. xd(4,ii) .eq. 1.5) then !(1.0000000000000000   2.0000000000000000   1.5000000000000000   1.5000000000000000)
           zd(ii)=snapmatrix_volumefraction(14,u_nodesi,snapshotsi,numphasei)
         else if(xd(1,ii) .eq. 2 .and. xd(2,ii) .eq. 2   .and. xd(3,ii) .eq. 1.5 .and. xd(4,ii) .eq. 1.5) then !(2.0000000000000000   2.0000000000000000   1.5000000000000000   1.5000000000000000)
           zd(ii)=snapmatrix_volumefraction(15,u_nodesi,snapshotsi,numphasei)
         else if(xd(1,ii) .eq. 1.5 .and. int(xd(2,ii)*1000)/1000 .eq. 1.146 .and. xd(3,ii) .eq. 1.5 .and. xd(4,ii) .eq. 1.5) then 
            !(1.5000000000000000   1.1464466094067263   1.5000000000000000   1.5000000000000000 xd(2,ii) >=1.14  .and. xd(2,ii) <=1.15)
           zd(ii)=snapmatrix_volumefraction(16,u_nodesi,snapshotsi,numphasei)
         else if(xd(1,ii) .eq. 1.5 .and. int(xd(2,ii)*1000)/1000 .eq. 1.853 .and. xd(3,ii) .eq. 1.5 .and. xd(4,ii) .eq. 1.5) then
          !(1.5000000000000000   1.8535533905932737   1.5000000000000000   1.5000000000000000)
           zd(ii)=snapmatrix_volumefraction(17,u_nodesi,snapshotsi,numphasei)
           else if(xd(1,ii) .eq. 1 .and. xd(2,ii) .eq. 1.5 .and. xd(3,ii) .eq. 1 .and. xd(4,ii) .eq. 1.5) then !( 1.0000000000000000   1.5000000000000000   1.0000000000000000  1.5000000000000000)
           zd(ii)=snapmatrix_volumefraction(18,u_nodesi,snapshotsi,numphasei)
         else if(xd(1,ii) .eq. 2 .and. xd(2,ii) .eq. 1.5   .and. xd(3,ii) .eq. 1 .and. xd(4,ii) .eq. 1.5) then !(2.0000000000000000   1.5000000000000000   1.0000000000000000   1.5000000000000000)
           zd(ii)=snapmatrix_volumefraction(19,u_nodesi,snapshotsi,numphasei)
         else if(xd(1,ii) .eq. 1 .and. xd(2,ii) .eq. 1.5   .and. xd(3,ii) .eq. 2 .and. xd(4,ii) .eq. 1.5) then !(1.0000000000000000   1.5000000000000000   2.0000000000000000   1.5000000000000000)
           zd(ii)=snapmatrix_volumefraction(20,u_nodesi,snapshotsi,numphasei)
         else if(xd(1,ii) .eq. 2 .and. xd(2,ii) .eq. 1.5   .and. xd(3,ii) .eq. 2 .and. xd(4,ii) .eq. 1.5) then !(2.0000000000000000   1.5000000000000000   2.0000000000000000   1.5000000000000000)
           zd(ii)=snapmatrix_volumefraction(21,u_nodesi,snapshotsi,numphasei) 
         else if(xd(1,ii) .eq. 1.5 .and. xd(2,ii) .eq. 1   .and. xd(3,ii) .eq. 1 .and. xd(4,ii) .eq. 1.5) then !(1.5000000000000000   1.0000000000000000   1.0000000000000000   1.5000000000000000)
           zd(ii)=snapmatrix_volumefraction(22,u_nodesi,snapshotsi,numphasei)
         else if(xd(1,ii) .eq. 1.5 .and. xd(2,ii) .eq. 2   .and. xd(3,ii) .eq. 1 .and. xd(4,ii) .eq. 1.5) then !(1.5000000000000000   2.0000000000000000   1.0000000000000000   1.5000000000000000)
           zd(ii)=snapmatrix_volumefraction(23,u_nodesi,snapshotsi,numphasei)
         else if(xd(1,ii) .eq. 1.5 .and. xd(2,ii) .eq. 1   .and. xd(3,ii) .eq. 2 .and. xd(4,ii) .eq. 1.5) then !(1.5000000000000000   1.0000000000000000   2.0000000000000000   1.5000000000000000)
           zd(ii)=snapmatrix_volumefraction(24,u_nodesi,snapshotsi,numphasei)
         else if(xd(1,ii) .eq. 1.5 .and. xd(2,ii) .eq. 2   .and. xd(3,ii) .eq. 2 .and. xd(4,ii) .eq. 1.5) then !(1.5000000000000000   2.0000000000000000   2.0000000000000000   1.5000000000000000)
           zd(ii)=snapmatrix_volumefraction(25,u_nodesi,snapshotsi,numphasei)
          else if(xd(1,ii) .eq. 1.5 .and. xd(2,ii) .eq. 1.5 .and. int(xd(3,ii)*1000)/1000 .eq. 1.146 .and. xd(4,ii) .eq. 1.5) then 
       !(1.5000000000000000   1.5000000000000000  1.1464466094067263   1.5000000000000000)
           zd(ii)=snapmatrix_volumefraction(26,u_nodesi,snapshotsi,numphasei)
         else if(xd(1,ii) .eq. 1.5  .and. xd(2,ii) .eq. 1.5  .and. int(xd(3,ii)*1000)/1000 .eq. 1.853 .and. xd(4,ii) .eq. 1.5) then 
       !( 1.5000000000000000   1.5000000000000000   1.8535533905932737   1.5000000000000000 )
           zd(ii)=snapmatrix_volumefraction(27,u_nodesi,snapshotsi,numphasei)
         else if(xd(1,ii) .eq. 1 .and. xd(2,ii) .eq. 1.5   .and. xd(3,ii) .eq. 1.5 .and. xd(4,ii) .eq. 1) then !(1.0000000000000000   1.5000000000000000   1.5000000000000000   1.0000000000000000)
           zd(ii)=snapmatrix_volumefraction(28,u_nodesi,snapshotsi,numphasei)
         else if(xd(1,ii) .eq. 2 .and. xd(2,ii) .eq. 1.5   .and. xd(3,ii) .eq. 1.5 .and. xd(4,ii) .eq. 1) then !(2.0000000000000000   1.5000000000000000   1.5000000000000000   1.0000000000000000)
           zd(ii)=snapmatrix_volumefraction(29,u_nodesi,snapshotsi,numphasei)
           else if(xd(1,ii) .eq. 1 .and. xd(2,ii) .eq. 1.5   .and. xd(3,ii) .eq. 1.5 .and. xd(4,ii) .eq. 2) then !(1.0000000000000000   1.5000000000000000   1.5000000000000000   2.0000000000000000)
           zd(ii)=snapmatrix_volumefraction(30,u_nodesi,snapshotsi,numphasei)
         else if(xd(1,ii) .eq. 2 .and. xd(2,ii) .eq. 1.5   .and. xd(3,ii) .eq. 1.5 .and. xd(4,ii) .eq. 2) then !(2.0000000000000000   1.5000000000000000   1.5000000000000000   2.0000000000000000)
           zd(ii)=snapmatrix_volumefraction(31,u_nodesi,snapshotsi,numphasei)
         else if(xd(1,ii) .eq. 1.5 .and. xd(2,ii) .eq. 1   .and. xd(3,ii) .eq. 1.5 .and. xd(4,ii) .eq. 1) then !(1.5000000000000000   1.0000000000000000   1.5000000000000000   1.0000000000000000)
           zd(ii)=snapmatrix_volumefraction(32,u_nodesi,snapshotsi,numphasei)
         else if(xd(1,ii) .eq. 1.5 .and. xd(2,ii) .eq. 2  .and. xd(3,ii) .eq. 1.5 .and. xd(4,ii) .eq. 1) then !(1.5000000000000000   2.0000000000000000   1.5000000000000000   1.0000000000000000)
           zd(ii)=snapmatrix_volumefraction(33,u_nodesi,snapshotsi,numphasei) 
         else if(xd(1,ii) .eq. 1.5 .and. xd(2,ii) .eq. 1   .and. xd(3,ii) .eq. 1.5 .and. xd(4,ii) .eq. 2) then !(1.5000000000000000   1.0000000000000000   1.5000000000000000   2.0000000000000000)
           zd(ii)=snapmatrix_volumefraction(34,u_nodesi,snapshotsi,numphasei)
         else if(xd(1,ii) .eq. 1.5 .and. xd(2,ii) .eq. 2   .and. xd(3,ii) .eq. 1.5 .and. xd(4,ii) .eq. 2) then !(1.5000000000000000   2.0000000000000000   1.5000000000000000   2.0000000000000000)
           zd(ii)=snapmatrix_volumefraction(35,u_nodesi,snapshotsi,numphasei)
         else if(xd(1,ii) .eq. 1.5 .and. xd(2,ii) .eq. 1.5   .and. xd(3,ii) .eq. 1 .and. xd(4,ii) .eq. 1) then !(1.5000000000000000   1.5000000000000000   1.0000000000000000   1.0000000000000000)
           zd(ii)=snapmatrix_volumefraction(36,u_nodesi,snapshotsi,numphasei)
         else if(xd(1,ii) .eq. 1.5 .and. xd(2,ii) .eq. 1.5   .and. xd(3,ii) .eq. 2 .and. xd(4,ii) .eq. 1) then !( 1.5000000000000000   1.5000000000000000   2.0000000000000000   1.0000000000000000)
           zd(ii)=snapmatrix_volumefraction(37,u_nodesi,snapshotsi,numphasei)
          else if(xd(1,ii) .eq. 1.5 .and. xd(2,ii) .eq. 1.5   .and. xd(3,ii) .eq. 1 .and. xd(4,ii) .eq. 2) then !(1.5000000000000000   1.5000000000000000   1.0000000000000000   2.0000000000000000)
           zd(ii)=snapmatrix_volumefraction(38,u_nodesi,snapshotsi,numphasei)
         else if(xd(1,ii) .eq. 1.5 .and. xd(2,ii) .eq. 1.5   .and. xd(3,ii) .eq. 2 .and. xd(4,ii) .eq. 2) then !(1.5000000000000000   1.5000000000000000   2.0000000000000000   2.0000000000000000)
           zd(ii)=snapmatrix_volumefraction(39,u_nodesi,snapshotsi,numphasei)
         else if(xd(1,ii) .eq. 1.5 .and. xd(2,ii) .eq. 1.5 .and.xd(3,ii) .eq. 1.5 .and.xd(4,ii) >= 1.14 .and. xd(4,ii) <= 1.15) then 
      !(1.5000000000000000   1.5000000000000000   1.5000000000000000   1.1464466094067263 )
           zd(ii)=snapmatrix_volumefraction(40,u_nodesi,snapshotsi,numphasei)
         else if(xd(1,ii) .eq. 1.5.and. xd(2,ii) .eq. 1.5 .and.xd(3,ii) .eq. 1.5.and. xd(4,ii) >= 1.85 .and. xd(4,ii) <= 1.855) then 
     !(1.5000000000000000   1.5000000000000000   1.5000000000000000   1.8535533905932737     )
           zd(ii)=snapmatrix_volumefraction(41,u_nodesi,snapshotsi,numphasei)
         endif
        enddo
          !call f_sinr ( m, nd, xd, zd )       
           
          ! print *, 'zd',  zd(:)
            
!  Use the grid to evaluate the interpolant.
        call lagrange_interp_nd_value2 ( m, ind, a, b, nd, zd, ni, xi, zpi )
 
!  Weighted the interpolant values and add to the sparse grid interpolant.
 
        nd_total = nd_total + nd
        zi(1:ni) = zi(1:ni) + c(l) * zpi(1:ni)
        !  print *, ' nd_total,level',  nd_total, l, c(l), nd
        !  print *, 'zizizizizizizic(l)zpi', zi(1)!,c(l),zpi(1)
        deallocate ( xd )
        deallocate ( zd )     
       
       if ( .not. more ) then
          exit
        end if    
 
      end do !DO

    end do !do l = l_min, l_max   
    deallocate ( c )
    deallocate ( w )
    
   end do ! do l_max = sparse_max, sparse_max
   interpolate_volumefraction(u_nodesi,snapshotsi,numphasei)=zi(ni)
   enddo !Do u_nodesi  
  enddo !Do snapshotsi 
  enddo !numphasei

 !  print *, ' total', total , nd_total
 ! open(unit=66,file='podbasis_interp') 
 !  do j=1, total_timestep 
 !     write(66,*)  pod_coef_all_obv_41i(j,1:num_pod)
 ! enddo
 ! close(66)
 
  deallocate ( a )
  deallocate ( b )
  deallocate ( ind )
  deallocate ( xi )
  deallocate ( ze )
  deallocate ( zi )
  deallocate ( zpi )
  end subroutine interpolate_snapmatrix_volumefraction



 subroutine interpolate_snapmatrix_velocity_rbf_8d( snapmatrix_velocity, interpolate_velocity)     ! CMAME-material-property
  implicit none
  
  real, dimension(:,:,:,:,:), allocatable :: snapmatrix_velocity ! allocate(snapmatrix_velocity(no_smolyak_nodes,u_nodes,snapshots,dim,numphase))
  real, dimension(:,:,:,:), intent(out),allocatable ::  interpolate_velocity 
 
  integer ( kind = 4 ) nd
  
  integer ( kind = 4 ) ni 
  integer ( kind = 4 ), allocatable :: w(:) 
   
  integer :: i, ii, total, total_timestep,num_pod,j,k,jj,kk, smolyak_total_points,u_nodesi,snapshotsi,dimi,numphasei
  real, allocatable :: wm(:,:)
  real :: r0
  real , allocatable :: xd(:,:) 
  real , allocatable :: xi(:,:)
  real , allocatable ::  fd(:) 
  real ( kind = 8 ), allocatable :: fi(:) 
  ni = 1
   
  !call r8mat_uniform_abvec ( m, ni, a, b, seed, xi )
 !  print * , 'xi', xi
  print *, 'm,ni', m, ni 
 ! stop 345
 
   ! call get_option(&
    !     '/reduced_model/variables_to_interpolate/dimension_one', xi(1,1))
  ! call get_option(&
   !      '/reduced_model/variables_to_interpolate/dimension_two', xi(2,1))
  ! call get_option(&
  !        '/reduced_model/variables_to_interpolate/dimension_three', xi(3,1))
  ! call get_option(&
  !       '/reduced_model/variables_to_interpolate/dimension_four', xi(4,1))
   
   
  allocate (interpolate_velocity(u_nodes,snapshots,dimen,numphase))  
   
   ni = 1
   nd = no_smolyak_nodes
   allocate ( xi(1:m,1:ni) )
   allocate ( fi(1:ni) ) 
   allocate (fd(nd)) 
   allocate (wm(nd,1))
   allocate (xd(m,nd))  
  !smolynode, u_nodes,snapshots,dim,numphase 
  Do u_nodesi=1, u_nodes  
   Do snapshotsi=1, snapshots  
   do dimi=1, size(snapmatrix_velocity,4)
   do numphasei=1, size(snapmatrix_velocity,5)  
    xd(1,1)=0.3! 
    xd(2,1)=0.3! 
    xd(3,1)=0.3
    xd(4,1)=0.3
    xd(5,1)=0.3
    xd(6,1)=0.3
    xd(7,1)=0.3
    xd(8,1)=0.3
     fd(1)=snapmatrix_velocity(1,u_nodesi,snapshotsi,dimi,numphasei)
 
    xd(1,2)=0.1! 
    xd(2,2)=0.3! 
    xd(3,2)=0.3
    xd(4,2)=0.3
    xd(5,2)=0.3
    xd(6,2)=0.3
    xd(7,2)=0.3
    xd(8,2)=0.3
    
      fd(2)=snapmatrix_velocity(2,u_nodesi,snapshotsi,dimi,numphasei)

    xd(1,3)=0.5! 
    xd(2,3)=0.3! 
    xd(3,3)=0.3
    xd(4,3)=0.3
    xd(5,3)=0.3
    xd(6,3)=0.3
    xd(7,3)=0.3
    xd(8,3)=0.3
  
      fd(3)=snapmatrix_velocity(3,u_nodesi,snapshotsi,dimi,numphasei)

    xd(1,4)=0.3   
    xd(2,4)=0.1   
    xd(3,4)=0.3
    xd(4,4)=0.3
    xd(5,4)=0.3
    xd(6,4)=0.3
    xd(7,4)=0.3
    xd(8,4)=0.3
      
       fd(4)=snapmatrix_velocity(4,u_nodesi,snapshotsi,dimi,numphasei)
     
    xd(1,5)=0.3! 
    xd(2,5)=0.5! 
    xd(3,5)=0.3
    xd(4,5)=0.3
    xd(5,5)=0.3
    xd(6,5)=0.3
    xd(7,5)=0.3
    xd(8,5)=0.3
       fd(5)=snapmatrix_velocity(5,u_nodesi,snapshotsi,dimi,numphasei)

    xd(1,6)=0.3! 
    xd(2,6)=0.3! 
    xd(3,6)=0.1
    xd(4,6)=0.3
    xd(5,6)=0.3
    xd(6,6)=0.3
    xd(7,6)=0.3
    xd(8,6)=0.3
      fd(6)=snapmatrix_velocity(6,u_nodesi,snapshotsi,dimi,numphasei)

    xd(1,7)=0.3 
    xd(2,7)=0.3  
    xd(3,7)=0.5
    xd(4,7)=0.3
    xd(5,7)=0.3
    xd(6,7)=0.3
    xd(7,7)=0.3
    xd(8,7)=0.3
        fd(7)=snapmatrix_velocity(7,u_nodesi,snapshotsi,dimi,numphasei)
    xd(1,8)=0.3! 
    xd(2,8)=0.3! 
    xd(3,8)=0.3
    xd(4,8)=0.1
    xd(5,8)=0.3
    xd(6,8)=0.3
    xd(7,8)=0.3
    xd(8,8)=0.3
        fd(8)=snapmatrix_velocity(8,u_nodesi,snapshotsi,dimi,numphasei)

    xd(1,9)=0.3! 
    xd(2,9)=0.3! 
    xd(3,9)=0.3
    xd(4,9)=0.5
    xd(5,9)=0.3
    xd(6,9)=0.3
    xd(7,9)=0.3
    xd(8,9)=0.3
        fd(9)=snapmatrix_velocity(9,u_nodesi,snapshotsi,dimi,numphasei)
    xd(1,10)=0.3 
    xd(2,10)=0.3! 
    xd(3,10)=0.3
    xd(4,10)=0.3
    xd(5,10)=0.1
    xd(6,10)=0.3
    xd(7,10)=0.3
    xd(8,10)=0.3
        fd(10)=snapmatrix_velocity(10,u_nodesi,snapshotsi,dimi,numphasei)
    xd(1,11)=0.3 
    xd(2,11)=0.3! 
    xd(3,11)=0.3
    xd(4,11)=0.3
    xd(5,11)=0.5
    xd(6,11)=0.3
    xd(7,11)=0.3
    xd(8,11)=0.3
        fd(11)=snapmatrix_velocity(11,u_nodesi,snapshotsi,dimi,numphasei)
    xd(1,12)=0.3 
    xd(2,12)=0.3! 
    xd(3,12)=0.3
    xd(4,12)=0.3
    xd(5,12)=0.3
    xd(6,12)=0.1
    xd(7,12)=0.3
    xd(8,12)=0.3
        fd(12)=snapmatrix_velocity(12,u_nodesi,snapshotsi,dimi,numphasei)
    xd(1,13)=0.1! 
    xd(2,13)=0.3! 
    xd(3,13)=0.3
    xd(4,13)=0.3
    xd(5,13)=0.3
    xd(6,13)=0.5
    xd(7,13)=0.3
    xd(8,13)=0.3
      fd(13)=snapmatrix_velocity(13,u_nodesi,snapshotsi,dimi,numphasei)
    xd(1,14)=0.3 
    xd(2,14)=0.3! 
    xd(3,14)=0.3
    xd(4,14)=0.3
    xd(5,14)=0.3
    xd(6,14)=0.3
    xd(7,14)=0.1
    xd(8,14)=0.3
      fd(14)=snapmatrix_velocity(14,u_nodesi,snapshotsi,dimi,numphasei)

    xd(1,15)=0.3! 
    xd(2,15)=0.3! 
    xd(3,15)=0.3
    xd(4,15)=0.3
    xd(5,15)=0.3
    xd(6,15)=0.3
    xd(7,15)=0.5
    xd(8,15)=0.3
      fd(15)=snapmatrix_velocity(15,u_nodesi,snapshotsi,dimi,numphasei)
    xd(1,16)=0.3 
    xd(2,16)=0.3! 
    xd(3,16)=0.3
    xd(4,16)=0.3
    xd(5,16)=0.3
    xd(6,16)=0.3
    xd(7,16)=0.3
    xd(8,16)=0.1
      fd(16)=snapmatrix_velocity(16,u_nodesi,snapshotsi,dimi,numphasei)
    xd(1,17)=0.3 
    xd(2,17)=0.3! 
    xd(3,17)=0.3
    xd(4,17)=0.3
    xd(5,17)=0.3
    xd(6,17)=0.3
    xd(7,17)=0.3
    xd(8,17)=0.5
      fd(17)=snapmatrix_velocity(17,u_nodesi,snapshotsi,dimi,numphasei)
      
    
      r0 = ( 1000- 50 ) / real ( nd, kind = 8 )
      call rbf_weight ( m, nd, xd, r0, phi4, fd, wm(:,1) )
     ! enddo 
        xi(1,1)=0.1 
        xi(2,1)=0.3 
        xi(3,1)=0.3 
        xi(4,1)=0.3
        xi(5,1)=0.1 
        xi(6,1)=0.3
        xi(7,1)=0.3 
        xi(8,1)=0.3
          
      !   xi(:,1)= pod_coef_old !pod_coef_all_obv(kk,:)pod_coef_all(s,:) !  
      call rbf_interp_nd ( m, nd, xd, r0, phi4, wm(:,1), ni, xi, fi )
     ! fi(1)=fi(1)*0.0000000000009869233/1000
     
      !PHI4 evaluates the gaussian radial basis function, best
      ! PHI3 evaluates the thin-plate spline radial basis function. terrible
      !PHI2 evaluates the inverse multiquadric radial basis function.,! 
      ! PHI1 evaluates the multiquadric radial basis function. terrible.! 
      ! pod_coef_new(l)=fi(1) 
     interpolate_velocity(u_nodesi,snapshotsi,dimi,numphasei)=fi(1) !zi(ni)
     !enddo  !do l=1, m     
   
 
   enddo !Do u_nodesi  
  enddo !Do snapshotsi
  enddo !dimi
  enddo !numphasei 
  
  deallocate ( fi )
  deallocate ( xi ) 
  deallocate(fd) 
  deallocate(wm)
  deallocate(xd)
  !deallocate(w)
  end subroutine interpolate_snapmatrix_velocity_rbf_8d

 subroutine interpolate_snapmatrix_pressure_rbf_8d( snapmatrix_pressure, interpolate_pressure)  
    
  implicit none
  !type(state_type), intent(in), dimension(:,:,:) :: state
  real, dimension(:,:,:,:), allocatable :: snapmatrix_pressure ! allocate(snapmatrix_pressure(no_smolyak_nodes,p_nodes,snapshots,numphase))
  real, dimension(:,:,:), intent(out), allocatable ::  interpolate_pressure 
   
  integer ( kind = 4 ) h 
  integer ( kind = 4 ) nd
  
  integer ( kind = 4 ) ni
  real :: r0
  integer ( kind = 4 ) t
  integer ( kind = 4 ), allocatable :: w(:)
  real ( kind = 8 ), allocatable :: xd(:,:)
  real ( kind = 8 ), allocatable :: xi(:,:) 
   real , allocatable ::  fd(:) 
   real ( kind = 8 ), allocatable :: fi(:) 
 real, allocatable :: wm(:,:)
  integer :: i, ii, total, total_timestep,num_pod,j,k,jj,kk, smolyak_total_points,u_nodesi,snapshotsi,dimi,numphasei
    ni = 1
   nd = no_smolyak_nodes
   allocate ( xi(1:m,1) )
   allocate ( fi(1:ni) ) 
   allocate (fd(nd)) 
   allocate (wm(nd,1))
   allocate (xd(m,nd))   
   
  
  allocate(interpolate_pressure(p_nodes,snapshots,numphase)) 
   
  !smolynode, u_nodes,snapshots,dim,numphase
  Do u_nodesi=1, p_nodes  
   Do snapshotsi=1, snapshots  
   do numphasei=1, numphase
     
    xd(1,1)=0.3! 
    xd(2,1)=0.3! 
    xd(3,1)=0.3
    xd(4,1)=0.3
    xd(5,1)=0.3
    xd(6,1)=0.3
    xd(7,1)=0.3
    xd(8,1)=0.3
     fd(1)=snapmatrix_pressure(1,u_nodesi,snapshotsi,numphasei) 

    xd(1,2)=0.1! 
    xd(2,2)=0.3! 
    xd(3,2)=0.3
    xd(4,2)=0.3
    xd(5,2)=0.3
    xd(6,2)=0.3
    xd(7,2)=0.3
    xd(8,2)=0.3
    
      fd(2)=snapmatrix_pressure(2,u_nodesi,snapshotsi,numphasei)

    xd(1,3)=0.5! 
    xd(2,3)=0.3! 
    xd(3,3)=0.3
    xd(4,3)=0.3
    xd(5,3)=0.3
    xd(6,3)=0.3
    xd(7,3)=0.3
    xd(8,3)=0.3
  
      fd(3)=snapmatrix_pressure(3,u_nodesi,snapshotsi,numphasei)

    xd(1,4)=0.3   
    xd(2,4)=0.1   
    xd(3,4)=0.3
    xd(4,4)=0.3
    xd(5,4)=0.3
    xd(6,4)=0.3
    xd(7,4)=0.3
    xd(8,4)=0.3
      
       fd(4)=snapmatrix_pressure(4,u_nodesi,snapshotsi,numphasei)
     
    xd(1,5)=0.3! 
    xd(2,5)=0.5! 
    xd(3,5)=0.3
    xd(4,5)=0.3
    xd(5,5)=0.3
    xd(6,5)=0.3
    xd(7,5)=0.3
    xd(8,5)=0.3
       fd(5)=snapmatrix_pressure(5,u_nodesi,snapshotsi,numphasei)

    xd(1,6)=0.3! 
    xd(2,6)=0.3! 
    xd(3,6)=0.1
    xd(4,6)=0.3
    xd(5,6)=0.3
    xd(6,6)=0.3
    xd(7,6)=0.3
    xd(8,6)=0.3
      fd(6)=snapmatrix_pressure(6,u_nodesi,snapshotsi,numphasei)

    xd(1,7)=0.3 
    xd(2,7)=0.3  
    xd(3,7)=0.5
    xd(4,7)=0.3
    xd(5,7)=0.3
    xd(6,7)=0.3
    xd(7,7)=0.3
    xd(8,7)=0.3
        fd(7)=snapmatrix_pressure(7,u_nodesi,snapshotsi,numphasei)
    xd(1,8)=0.3! 
    xd(2,8)=0.3! 
    xd(3,8)=0.3
    xd(4,8)=0.1
    xd(5,8)=0.3
    xd(6,8)=0.3
    xd(7,8)=0.3
    xd(8,8)=0.3
        fd(8)=snapmatrix_pressure(8,u_nodesi,snapshotsi,numphasei)

    xd(1,9)=0.3! 
    xd(2,9)=0.3! 
    xd(3,9)=0.3
    xd(4,9)=0.5
    xd(5,9)=0.3
    xd(6,9)=0.3
    xd(7,9)=0.3
    xd(8,9)=0.3
        fd(9)=snapmatrix_pressure(9,u_nodesi,snapshotsi,numphasei)
    xd(1,10)=0.3 
    xd(2,10)=0.3! 
    xd(3,10)=0.3
    xd(4,10)=0.3
    xd(5,10)=0.1
    xd(6,10)=0.3
    xd(7,10)=0.3
    xd(8,10)=0.3
        fd(10)=snapmatrix_pressure(10,u_nodesi,snapshotsi,numphasei)
    xd(1,11)=0.3 
    xd(2,11)=0.3! 
    xd(3,11)=0.3
    xd(4,11)=0.3
    xd(5,11)=0.5
    xd(6,11)=0.3
    xd(7,11)=0.3
    xd(8,11)=0.3
        fd(11)=snapmatrix_pressure(11,u_nodesi,snapshotsi,numphasei)
    xd(1,12)=0.3 
    xd(2,12)=0.3! 
    xd(3,12)=0.3
    xd(4,12)=0.3
    xd(5,12)=0.3
    xd(6,12)=0.1
    xd(7,12)=0.3
    xd(8,12)=0.3
        fd(12)=snapmatrix_pressure(12,u_nodesi,snapshotsi,numphasei)
    xd(1,13)=0.1! 
    xd(2,13)=0.3! 
    xd(3,13)=0.3
    xd(4,13)=0.3
    xd(5,13)=0.3
    xd(6,13)=0.5
    xd(7,13)=0.3
    xd(8,13)=0.3
        fd(13)=snapmatrix_pressure(13,u_nodesi,snapshotsi,numphasei)
    xd(1,14)=0.3 
    xd(2,14)=0.3! 
    xd(3,14)=0.3
    xd(4,14)=0.3
    xd(5,14)=0.3
    xd(6,14)=0.3
    xd(7,14)=0.1
    xd(8,14)=0.3
        fd(14)=snapmatrix_pressure(14,u_nodesi,snapshotsi,numphasei)

    xd(1,15)=0.3! 
    xd(2,15)=0.3! 
    xd(3,15)=0.3
    xd(4,15)=0.3
    xd(5,15)=0.3
    xd(6,15)=0.3
    xd(7,15)=0.5
    xd(8,15)=0.3
      fd(15)=snapmatrix_pressure(15,u_nodesi,snapshotsi,numphasei)
    xd(1,16)=0.3 
    xd(2,16)=0.3! 
    xd(3,16)=0.3
    xd(4,16)=0.3
    xd(5,16)=0.3
    xd(6,16)=0.3
    xd(7,16)=0.3
    xd(8,16)=0.1
        fd(16)=snapmatrix_pressure(16,u_nodesi,snapshotsi,numphasei)
    xd(1,17)=0.3 
    xd(2,17)=0.3! 
    xd(3,17)=0.3
    xd(4,17)=0.3
    xd(5,17)=0.3
    xd(6,17)=0.3
    xd(7,17)=0.3
    xd(8,17)=0.5
        fd(17)=snapmatrix_pressure(17,u_nodesi,snapshotsi,numphasei)
       
   
     ! print * ,'xd',xd
     ! r0 = ( 9.869233E-013 - 3.2897443e-13 ) / real ( nd, kind = 8 )
     ! r0 = ( 1000- 50 ) / real ( nd, kind = 8 )
      r0=1
      call rbf_weight ( m, nd, xd, r0, phi4, fd, wm(:,1) )
     ! enddo 
        xi(1,1)=0.1 
        xi(2,1)=0.3 
        xi(3,1)=0.3 
        xi(4,1)=0.3
        xi(5,1)=0.1 
        xi(6,1)=0.3
        xi(7,1)=0.3 
        xi(8,1)=0.3
        
         !xi(4,1)= 2!1.7 !2.0000000000000000   1.5000000000000000   2.0000000000000000   1.5000000000000000,21
      !   xi(:,1)= pod_coef_old !pod_coef_all_obv(kk,:)pod_coef_all(s,:) !  
      call rbf_interp_nd ( m, nd, xd, r0, phi4, wm(:,1), ni, xi, fi ) 
     ! fi(1)=fi(1)*0.0000000000009869233/1000
      interpolate_pressure(u_nodesi,snapshotsi,numphasei)=fi(1)
   enddo !Do u_nodesi  
  enddo !Do snapshotsi 
  enddo !numphasei 
  
  
  
  deallocate ( fi )
  deallocate ( xi ) 
  deallocate(fd) 
  deallocate(wm)
  deallocate(xd)
  !deallocate(w)
  end subroutine interpolate_snapmatrix_pressure_rbf_8d

 subroutine interpolate_snapmatrix_volumefraction_rbf_8d( snapmatrix_volumefraction, interpolate_volumefraction)  
  
!
  implicit none
  !type(state_type), intent(in), dimension(:,:,:) :: state
  real, dimension(:,:,:,:), allocatable :: snapmatrix_volumefraction ! allocate(snapmatrix_pressure(no_smolyak_nodes,p_nodes,snapshots,numphase))
  real, dimension(:,:,:), intent(out), allocatable ::  interpolate_volumefraction 
   
  integer ( kind = 4 ) h 
  integer ( kind = 4 ) nd 
  integer ( kind = 4 ) ni 
   real, allocatable :: wm(:,:)
 
  integer ( kind = 4 ), allocatable :: w(:)
  real ( kind = 8 ), allocatable :: xd(:,:)
  real ( kind = 8 ), allocatable :: xi(:,:)
  real , allocatable ::  fd(:) 
   real ( kind = 8 ), allocatable :: fi(:) 
   real :: r0
  integer :: i, ii, total, total_timestep,num_pod,j,k,jj,kk, smolyak_total_points,u_nodesi,snapshotsi,dimi,numphasei
      ni = 1
   nd = no_smolyak_nodes
   allocate ( xi(1:m,1:ni) )
   allocate ( fi(1:ni) ) 
   allocate (fd(nd)) 
   allocate (wm(nd,1))
   allocate (xd(m,nd))   
     
  allocate(interpolate_volumefraction(p_nodes,snapshots,numphase))
  !smolynode, u_nodes,snapshots,dim,numphase
  Do u_nodesi=1, p_nodes  
   Do snapshotsi=1, snapshots  
   do numphasei=1, numphase
      
    xd(1,1)=0.3! 
    xd(2,1)=0.3! 
    xd(3,1)=0.3
    xd(4,1)=0.3
    xd(5,1)=0.3
    xd(6,1)=0.3
    xd(7,1)=0.3
    xd(8,1)=0.3
     fd(1)=snapmatrix_volumefraction(1,u_nodesi,snapshotsi,numphasei)
 
    xd(1,2)=0.1! 
    xd(2,2)=0.3! 
    xd(3,2)=0.3
    xd(4,2)=0.3
    xd(5,2)=0.3
    xd(6,2)=0.3
    xd(7,2)=0.3
    xd(8,2)=0.3
    
      fd(2)=snapmatrix_volumefraction(2,u_nodesi,snapshotsi,numphasei)

    xd(1,3)=0.5! 
    xd(2,3)=0.3! 
    xd(3,3)=0.3
    xd(4,3)=0.3
    xd(5,3)=0.3
    xd(6,3)=0.3
    xd(7,3)=0.3
    xd(8,3)=0.3
      fd(3)=snapmatrix_volumefraction(3,u_nodesi,snapshotsi,numphasei)

    xd(1,4)=0.3   
    xd(2,4)=0.1   
    xd(3,4)=0.3
    xd(4,4)=0.3
    xd(5,4)=0.3
    xd(6,4)=0.3
    xd(7,4)=0.3
    xd(8,4)=0.3 
      fd(4)=snapmatrix_volumefraction(4,u_nodesi,snapshotsi,numphasei)
     
    xd(1,5)=0.3! 
    xd(2,5)=0.5! 
    xd(3,5)=0.3
    xd(4,5)=0.3
    xd(5,5)=0.3
    xd(6,5)=0.3
    xd(7,5)=0.3
    xd(8,5)=0.3
      fd(5)=snapmatrix_volumefraction(5,u_nodesi,snapshotsi,numphasei)

    xd(1,6)=0.3! 
    xd(2,6)=0.3! 
    xd(3,6)=0.1
    xd(4,6)=0.3
    xd(5,6)=0.3
    xd(6,6)=0.3
    xd(7,6)=0.3
    xd(8,6)=0.3
      fd(6)=snapmatrix_volumefraction(6,u_nodesi,snapshotsi,numphasei)

    xd(1,7)=0.3 
    xd(2,7)=0.3  
    xd(3,7)=0.5
    xd(4,7)=0.3
    xd(5,7)=0.3
    xd(6,7)=0.3
    xd(7,7)=0.3
    xd(8,7)=0.3
      fd(7)=snapmatrix_volumefraction(7,u_nodesi,snapshotsi,numphasei)
    xd(1,8)=0.3! 
    xd(2,8)=0.3! 
    xd(3,8)=0.3
    xd(4,8)=0.1
    xd(5,8)=0.3
    xd(6,8)=0.3
    xd(7,8)=0.3
    xd(8,8)=0.3
        fd(8)=snapmatrix_volumefraction(8,u_nodesi,snapshotsi,numphasei)

    xd(1,9)=0.3! 
    xd(2,9)=0.3! 
    xd(3,9)=0.3
    xd(4,9)=0.5
    xd(5,9)=0.3
    xd(6,9)=0.3
    xd(7,9)=0.3
    xd(8,9)=0.3
      fd(9)=snapmatrix_volumefraction(9,u_nodesi,snapshotsi,numphasei)
    xd(1,10)=0.3 
    xd(2,10)=0.3! 
    xd(3,10)=0.3
    xd(4,10)=0.3
    xd(5,10)=0.1
    xd(6,10)=0.3
    xd(7,10)=0.3
    xd(8,10)=0.3
        fd(10)=snapmatrix_volumefraction(10,u_nodesi,snapshotsi,numphasei)
    xd(1,11)=0.3 
    xd(2,11)=0.3! 
    xd(3,11)=0.3
    xd(4,11)=0.3
    xd(5,11)=0.5
    xd(6,11)=0.3
    xd(7,11)=0.3
    xd(8,11)=0.3
      fd(11)=snapmatrix_volumefraction(11,u_nodesi,snapshotsi,numphasei)
    xd(1,12)=0.3 
    xd(2,12)=0.3! 
    xd(3,12)=0.3
    xd(4,12)=0.3
    xd(5,12)=0.3
    xd(6,12)=0.1
    xd(7,12)=0.3
    xd(8,12)=0.3
      fd(12)=snapmatrix_volumefraction(12,u_nodesi,snapshotsi,numphasei)
    xd(1,13)=0.1! 
    xd(2,13)=0.3! 
    xd(3,13)=0.3
    xd(4,13)=0.3
    xd(5,13)=0.3
    xd(6,13)=0.5
    xd(7,13)=0.3
    xd(8,13)=0.3
      fd(13)=snapmatrix_volumefraction(13,u_nodesi,snapshotsi,numphasei)
    xd(1,14)=0.3 
    xd(2,14)=0.3! 
    xd(3,14)=0.3
    xd(4,14)=0.3
    xd(5,14)=0.3
    xd(6,14)=0.3
    xd(7,14)=0.1
    xd(8,14)=0.3
        fd(14)=snapmatrix_volumefraction(14,u_nodesi,snapshotsi,numphasei)

    xd(1,15)=0.3! 
    xd(2,15)=0.3! 
    xd(3,15)=0.3
    xd(4,15)=0.3
    xd(5,15)=0.3
    xd(6,15)=0.3
    xd(7,15)=0.5
    xd(8,15)=0.3
        fd(15)=snapmatrix_volumefraction(15,u_nodesi,snapshotsi,numphasei)
    xd(1,16)=0.3 
    xd(2,16)=0.3! 
    xd(3,16)=0.3
    xd(4,16)=0.3
    xd(5,16)=0.3
    xd(6,16)=0.3
    xd(7,16)=0.3
    xd(8,16)=0.1
        fd(16)=snapmatrix_volumefraction(16,u_nodesi,snapshotsi,numphasei)
    xd(1,17)=0.3 
    xd(2,17)=0.3! 
    xd(3,17)=0.3
    xd(4,17)=0.3
    xd(5,17)=0.3
    xd(6,17)=0.3
    xd(7,17)=0.3
    xd(8,17)=0.5
      fd(17)=snapmatrix_volumefraction(17,u_nodesi,snapshotsi,numphasei)
   
     ! print * ,'xd',xd
     ! r0 = ( 9.869233E-013 - 3.2897443e-13 ) / real ( nd, kind = 8 )
     r0 = ( 1000- 50 ) / real ( nd, kind = 8 )
     !r0=1
      call rbf_weight ( m, nd, xd, r0, phi4, fd, wm(:,1) )
     ! enddo 
        xi(1,1)=0.1 
        xi(2,1)=0.3 
        xi(3,1)=0.3 
        xi(4,1)=0.3
        xi(5,1)=0.1 
        xi(6,1)=0.3
        xi(7,1)=0.3 
        xi(8,1)=0.3 
      call rbf_interp_nd ( m, nd, xd, r0, phi4, wm(:,1), ni, xi, fi ) 
      !fi(1)=fi(1)*0.0000000000009869233/1000
   interpolate_volumefraction(u_nodesi,snapshotsi,numphasei)=fi(1)
   enddo !Do u_nodesi  
  enddo !Do snapshotsi 
  enddo !numphasei 
  deallocate ( fi )
  deallocate ( xi ) 
  deallocate(fd) 
  deallocate(wm)
  deallocate(xd)
  !deallocate(w)
  end subroutine interpolate_snapmatrix_volumefraction_rbf_8d

 
 subroutine interpolate_snapmatrix_velocity_rbf_22d( snapmatrix_velocity, interpolate_velocity)     ! CMAME-material-property
  implicit none
  
  real, dimension(:,:,:,:,:), allocatable :: snapmatrix_velocity ! allocate(snapmatrix_velocity(no_smolyak_nodes,u_nodes,snapshots,dim,numphase))
  real, dimension(:,:,:,:), intent(out),allocatable ::  interpolate_velocity 
 
  integer ( kind = 4 ) nd
  
  integer ( kind = 4 ) ni 
  integer ( kind = 4 ), allocatable :: w(:) 
  
  integer :: i, ii, total, total_timestep,num_pod,j,k,jj,kk, smolyak_total_points,u_nodesi,snapshotsi,dimi,numphasei, mi, no_smolyak_nodesi
  real, allocatable :: wm(:,:)
  real :: r0
  real , allocatable :: xd(:,:) 
  real , allocatable :: xi(:,:)
  real , allocatable ::  fd(:) 
  real ( kind = 8 ), allocatable :: fi(:) 
  ni = 1
   
  !call r8mat_uniform_abvec ( m, ni, a, b, seed, xi )
 !  print * , 'xi', xi
  print *, 'm,ni', m, ni 
 ! stop 345
 
   ! call get_option(&
    !     '/reduced_model/variables_to_interpolate/dimension_one', xi(1,1))
  ! call get_option(&
   !      '/reduced_model/variables_to_interpolate/dimension_two', xi(2,1))
  ! call get_option(&
  !        '/reduced_model/variables_to_interpolate/dimension_three', xi(3,1))
  ! call get_option(&
  !       '/reduced_model/variables_to_interpolate/dimension_four', xi(4,1))

  allocate (interpolate_velocity(u_nodes,snapshots,dimen,numphase))  
   
   ni = 1
   nd = no_smolyak_nodes
   allocate ( xi(1:m,1:ni) )
   allocate ( fi(1:ni) ) 
   allocate (fd(nd)) 
   allocate (wm(nd,1))
   allocate (xd(m,nd))  
  !smolynode, u_nodes,snapshots,dim,numphase 
  Do u_nodesi=1, u_nodes 
     print *, 'u_nodesi=1, u_nodes ', u_nodesi, u_nodes 
   Do snapshotsi=1, snapshots  
   do dimi=1, size(snapmatrix_velocity,4)
   do numphasei=1, 1!size(snapmatrix_velocity,5)       
        do no_smolyak_nodesi=1, no_smolyak_nodes
            do mi=1, m
              xd(mi, no_smolyak_nodesi)=input_parameters(no_smolyak_nodesi,mi)
            enddo
        fd(no_smolyak_nodesi) = snapmatrix_velocity(no_smolyak_nodesi,u_nodesi,snapshotsi,dimi,numphasei)
    
    
      r0 = ( 10000- 5 ) / real ( nd, kind = 8 )
     ! call rbf_weight ( m, nd, xd, r0, phi4, fd, wm(:,1) )
     
       xi(:,1)=input_parameters(no_smolyak_nodesi,:)
          
      
      ! call rbf_interp_nd ( m, nd, xd, r0, phi4, wm(:,1), ni, xi, fi )
    
      ! fi(1)=fi(1)*0.0000000000009869233/1000
      ! PHI4 evaluates the gaussian radial basis function, best
      ! PHI3 evaluates the thin-plate spline radial basis function. terrible
      ! PHI2 evaluates the inverse multiquadric radial basis function.,! 
      ! PHI1 evaluates the multiquadric radial basis function. terrible.! 
      ! pod_coef_new(l)=fi(1) 
     interpolate_velocity(u_nodesi,snapshotsi,dimi,numphasei)=0!fi(1) !zi(ni)
     !enddo  !do l=1, m     
   
    enddo ! no_smolyak_nodesi
   enddo  !numphasei 
  enddo !dimi
  enddo !Do snapshotsi
  enddo !Do u_nodesi 
  
  deallocate ( fi )
  deallocate ( xi ) 
  deallocate(fd) 
  deallocate(wm)
  deallocate(xd)
  !deallocate(w)
  end subroutine interpolate_snapmatrix_velocity_rbf_22d

 subroutine interpolate_snapmatrix_pressure_rbf_22d( snapmatrix_pressure, interpolate_pressure)  
    
  implicit none
  !type(state_type), intent(in), dimension(:,:,:) :: state
  real, dimension(:,:,:,:), allocatable :: snapmatrix_pressure ! allocate(snapmatrix_pressure(no_smolyak_nodes,p_nodes,snapshots,numphase))
  real, dimension(:,:,:), intent(out), allocatable ::  interpolate_pressure 
   
  integer ( kind = 4 ) h 
  integer ( kind = 4 ) nd
  
  integer ( kind = 4 ) ni
  real :: r0
  integer ( kind = 4 ) t
  integer ( kind = 4 ), allocatable :: w(:)
  real ( kind = 8 ), allocatable :: xd(:,:)
  real ( kind = 8 ), allocatable :: xi(:,:) 
   real , allocatable ::  fd(:) 
   real ( kind = 8 ), allocatable :: fi(:) 
 real, allocatable :: wm(:,:)
  integer :: i, ii, total, total_timestep,num_pod,j,k,jj,kk, smolyak_total_points,u_nodesi,snapshotsi,dimi,numphasei, mi, no_smolyak_nodesi
    ni = 1
   nd = no_smolyak_nodes
   allocate ( xi(1:m,1) )
   allocate ( fi(1:ni) ) 
   allocate (fd(nd)) 
   allocate (wm(nd,1))
   allocate (xd(m,nd))   
   
  
  allocate(interpolate_pressure(p_nodes,snapshots,numphase)) 
   
  !smolynode, u_nodes,snapshots,dim,numphase
  Do u_nodesi=1, p_nodes  
   print * ,'pressure nodes' , u_nodesi, p_nodes 
   Do snapshotsi=1, snapshots  
   do numphasei=1, 1!numphase
     
     do no_smolyak_nodesi=1, no_smolyak_nodes
            do mi=1, m
             xd(mi, no_smolyak_nodesi)=input_parameters(no_smolyak_nodesi,mi)
            enddo
        fd(no_smolyak_nodesi) = snapmatrix_pressure(no_smolyak_nodesi,u_nodesi,snapshotsi,numphasei)
           
    
     
     ! r0 = ( 9.869233E-013 - 3.2897443e-13 ) / real ( nd, kind = 8 )
       r0 = ( 10000- 5 ) / real ( nd, kind = 8 ) 
  !    call rbf_weight ( m, nd, xd, r0, phi4, fd, wm(:,1) )
     
       xi(:,1)=input_parameters(no_smolyak_nodesi,:)
        
         
   !    call rbf_interp_nd ( m, nd, xd, r0, phi4, wm(:,1), ni, xi, fi ) 
     
      interpolate_pressure(u_nodesi,snapshotsi,numphasei)=fi(1)
   enddo !Do u_nodesi  
  enddo !Do snapshotsi 
  enddo !numphasei 
  enddo
  
  
  deallocate ( fi )
  deallocate ( xi ) 
  deallocate(fd) 
  deallocate(wm)
  deallocate(xd)
  !deallocate(w)
  end subroutine interpolate_snapmatrix_pressure_rbf_22d

 subroutine interpolate_snapmatrix_volumefraction_rbf_22d( snapmatrix_volumefraction, interpolate_volumefraction)  
  
!
  implicit none
  !type(state_type), intent(in), dimension(:,:,:) :: state
  real, dimension(:,:,:,:), allocatable :: snapmatrix_volumefraction ! allocate(snapmatrix_pressure(no_smolyak_nodes,p_nodes,snapshots,numphase))
  real, dimension(:,:,:), intent(out), allocatable ::  interpolate_volumefraction 
   
  integer ( kind = 4 ) h 
  integer ( kind = 4 ) nd 
  integer ( kind = 4 ) ni 
   real, allocatable :: wm(:,:)
 
  integer ( kind = 4 ), allocatable :: w(:)
  real ( kind = 8 ), allocatable :: xd(:,:)
  real ( kind = 8 ), allocatable :: xi(:,:)
  real , allocatable ::  fd(:) 
   real ( kind = 8 ), allocatable :: fi(:) 
   real :: r0
  integer :: i, ii, total, total_timestep,num_pod,j,k,jj,kk, smolyak_total_points,u_nodesi,snapshotsi,dimi,numphasei, mi, no_smolyak_nodesi
      ni = 1
   nd = no_smolyak_nodes+1
   allocate ( xi(1:m,1:ni) )
   allocate ( fi(1:ni) ) 
   allocate (fd(nd)) 
   allocate (wm(nd,1))
   allocate (xd(m,nd))   
     
  allocate(interpolate_volumefraction(p_nodes,snapshots,numphase))
  !smolynode, u_nodes,snapshots,dim,numphase
  Do u_nodesi=1, p_nodes  
   Do snapshotsi=1, snapshots  
   do numphasei=1, 1! numphase 
   do no_smolyak_nodesi=1, no_smolyak_nodes
            do mi=1, m
             xd(mi, no_smolyak_nodesi)=input_parameters(no_smolyak_nodesi,mi)
            enddo
        fd(no_smolyak_nodesi) = snapmatrix_volumefraction(no_smolyak_nodesi,u_nodesi,snapshotsi,numphasei)
     
   
    ! print * ,'volumefraction' 
     ! r0 = ( 9.869233E-013 - 3.2897443e-13 ) / real ( nd, kind = 8 )
     r0 = ( 10000- 5 ) / real ( nd, kind = 8 )
       call rbf_weight ( m, nd, xd, r0, phi4, fd, wm(:,1) )
      
       xi(:,1)=input_parameters(no_smolyak_nodesi,:)
       call rbf_interp_nd ( m, nd, xd, r0, phi4, wm(:,1), ni, xi, fi ) 
      
   interpolate_volumefraction(u_nodesi,snapshotsi,numphasei)=fi(1)
   enddo !Do u_nodesi  
  enddo !Do snapshotsi 
  enddo !numphasei
  enddo 
  deallocate ( fi )
  deallocate ( xi ) 
  deallocate(fd) 
  deallocate(wm)
  deallocate(xd)
  !deallocate(w)
  end subroutine interpolate_snapmatrix_volumefraction_rbf_22d



 subroutine interpolate_snapmatrix_velocity_4dim( snapmatrix_velocity, interpolate_velocity)  
    
! top-bottom between 1-2, pod_coef_all_obv is used for storing some podbasis to avoid to much modifications.
 
!  Parameters:
!
!    Local, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) SPARSE_MAX, the maximum sparse grid level to try.
!
!  Local Parameters:
!
!    Local, real ( kind = 8 ) A(M), B(M), the upper and lower variable limits 
!    in each dimension.
!
!    Local, real ( kind = 8 ) APP_ERROR, the averaged Euclidean norm of the 
!    difference between the sparse interpolant and the exact function at 
!    the interpolation points.
!
!    Local, integer ( kind = 4 ) C(0:L_MAX), the sparse grid coefficient vector.
!    Results at level L are weighted by C(L).
!
!    Local, integer ( kind = 4 ) IND(M), the 1D indices defining a Lagrange grid.
!    Each entry is a 1d "level" that specifies the order of a 
!    Clenshaw Curtis 1D grid.
!
!    Local, integer ( kind = 4 ) L, the current Lagrange grid level.
!
!    Local, integer ( kind = 4 ) L_MAX, the current sparse grid level.
!
!    Local, integer ( kind = 4 ) MORE, used to control the enumeration of all the
!    Lagrange grids at a current grid level.
!
!    Local, integer ( kind = 4 ) ND, the number of points used in a Lagrange grid.
!
!    Local, integer ( kind = 4 ) ND_TOTAL, the total number of points used in all the
!    Lagrange interpolants for a given sparse interpolant points that occur
!    in more than one Lagrange grid are counted multiple times.
!
!    Local, integer ( kind = 4 ) NI, the number of interpolant evaluation points.
!
!    Local, integer ( kind = 4 ) SPARSE_MIN, the minimum sparse grid level to try.
!
!    Local, real ( kind = 8 ) XD(M,ND), the data points for a Lagrange grid.
!
!    Local, real ( kind = 8 ) XI(M,NI), the interpolant evaluation points.
!
!    Local, real ( kind = 8 ) ZD(ND), the data values for a Lagrange grid.
!
!    Local, real ( kind = 8 ) ZE(NI), the exact function values at XI.
!
!    Local, real ( kind = 8 ) ZI(NI), the sparse interpolant values at XI.
!
!    Local, real ( kind = 8 ) ZPI(NI), one set of Lagrange interpolant values at XI.
!
  implicit none
  !type(state_type), intent(in), dimension(:,:,:) :: state
  real, dimension(:,:,:,:,:), allocatable :: snapmatrix_velocity ! allocate(snapmatrix_velocity(no_smolyak_nodes,u_nodes,snapshots,dim,numphase))
  real, dimension(:,:,:,:), intent(out),allocatable ::  interpolate_velocity 
  real ( kind = 8 ), allocatable :: a(:)
  real ( kind = 8 ) app_error
  real ( kind = 8 ), allocatable :: b(:)
  integer ( kind = 4 ), allocatable :: c(:)
  integer ( kind = 4 ) h
  integer ( kind = 4 ), allocatable :: ind(:)
  integer ( kind = 4 ) l
  integer ( kind = 4 ) l_max
  integer ( kind = 4 ) l_min
  real :: a1,a2,a3,a4,a5
  logical more
  integer ( kind = 4 ) nd
  integer ( kind = 4 ) nd_total
  integer ( kind = 4 ) ni
  real ( kind = 8 ) r8vec_norm_affine
  integer ( kind = 4 ) seed
 ! integer ( kind = 4 ) sparse_max
  integer ( kind = 4 ) sparse_min
  integer ( kind = 4 ) t
  integer ( kind = 4 ), allocatable :: w(:)
  real ( kind = 8 ), allocatable :: xd(:,:)
  real ( kind = 8 ), allocatable :: xi(:,:)
  real ( kind = 8 ), allocatable :: zd(:)
  real ( kind = 8 ), allocatable :: ze(:)
  real ( kind = 8 ), allocatable :: zi(:)
  real ( kind = 8 ), allocatable :: zpi(:)
  integer :: i, ii, total, total_timestep,num_pod,j,k,jj,kk, smolyak_total_points,u_nodesi,snapshotsi,dimi,numphasei
   
  allocate (ind(1:m))

  ! Define the region.
  allocate ( a(1:m) )
  allocate ( b(1:m) )
  
 
   
  !total_timestep=864
  !num_pod=12*36
  !smolyak_total_points=41!2**sparse_max+1
  !print *, '  smolyak_total_points=2**sparse_max+1=', smolyak_total_points 
  !  stop 2323
  !  Define the interpolation evaluation information.
  ni = 1
  seed = 123456789
  allocate ( xi(m,ni) )
  !call r8mat_uniform_abvec ( m, ni, a, b, seed, xi )
!  print * , 'xi', xi
  print *, 'm,ni', m, ni 
 ! stop 345
   

  if (.true.) then
  a(1:m) = 0.1
  b(1:m) = 0.5
  a1=0.1
  a2=0.158578
  a3=0.3
  a4=0.441421
  a5=0.5
 
  xi(1,1)= 0.1 
  xi(2,1)= 0.5
  xi(3,1)= 0.5
  xi(4,1)= 0.5! 0.10000000000000000       0.30000000000000000       0.50000000000000000       0.30000000000000000  
   call get_option(&
         '/reduced_model/variables_to_interpolate/dimension_one', xi(1,1))
   call get_option(&
         '/reduced_model/variables_to_interpolate/dimension_two', xi(2,1))
 call get_option(&
         '/reduced_model/variables_to_interpolate/dimension_three', xi(3,1))
 call get_option(&
         '/reduced_model/variables_to_interpolate/dimension_four', xi(4,1))   
  endif
  write ( *, '(a,i6)' ) '  Number of interpolation points is NI = ', ni
  allocate ( ze(1:ni) )
 ! call f_sinr ( m, ni, xi, ze )
  print *, 'actual function value is, will interpolate this value later on,then do comparison', xi,ze
!
!  Compute a sequence of sparse grid interpolants of increasing level.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '   L     Order    ApproxError'
  write ( *, '(a)' ) ''
  allocate (interpolate_velocity(u_nodes,snapshots,dimen,numphase))
  allocate ( zpi(1:ni) )
  allocate ( zi(1:ni) )
  !smolynode, u_nodes,snapshots,dim,numphase
  Do u_nodesi=1, u_nodes  
   Do snapshotsi=1, snapshots  
   do dimi=1, size(snapmatrix_velocity,4)
   do numphasei=1, size(snapmatrix_velocity,5)
  total = 0
  sparse_min = 0

  ! do l_max = sparse_min, sparse_max
    do l_max = sparse_max, sparse_max
    allocate ( c(0:l_max) )
    allocate ( w(0:l_max) )
    call smolyak_coefficients ( l_max, m, c, w )

    zi(1:ni) = 0.0D+00
    nd_total = 0
    
    l_min = max ( l_max + 1 - m, 0 )

    do l = l_min, l_max

      more = .false.
      
      do
     ! print *,'l_max,l,' ,l_max,l   
!  Define the next product grid.
!
        call comp_next ( l, m, ind, more, h, t )
!
!  Count the grid, find its points, evaluate the data there.
!
        call lagrange_interp_nd_size2 ( m, ind, nd )
        allocate ( xd(1:m,1:nd) )
        call lagrange_interp_nd_grid2 ( m, ind, a, b, nd, xd )        
         !  print *,'xd', xd, 点1 坐标为： xd(:,1)
        allocate ( zd(1:nd) )

       ! print *,'---------------------nd',nd  
        do ii=1, nd           
              if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3   .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3) then 
!( 1.5000000000000000    1.5000000000000000   1.5000000000000000  1.5000000000000000)
           zd(ii)=snapmatrix_velocity(1,u_nodesi,snapshotsi,dimi,numphasei)
         else if(xd(1,ii) .eq. a1 .and. xd(2,ii) .eq. a3   .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3) then 
!(1.0000000000000000   1.5000000000000000   1.5000000000000000   1.5000000000000000)
           zd(ii)=snapmatrix_velocity(2,u_nodesi,snapshotsi,dimi,numphasei)
         else if(xd(1,ii) .eq. a5 .and. xd(2,ii) .eq. a3   .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3) then
 !( 2.0000000000000000   1.5000000000000000   1.5000000000000000   1.5000000000000000)
           zd(ii)=snapmatrix_velocity(3,u_nodesi,snapshotsi,dimi,numphasei)
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a1   .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3) then 
!(1.5000000000000000   1.0000000000000000   1.5000000000000000   1.5000000000000000 )
           zd(ii)=snapmatrix_velocity(4,u_nodesi,snapshotsi,dimi,numphasei)
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a5 .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3) then 
!(1.5000000000000000   2.0000000000000000   1.5000000000000000   1.5000000000000000 )
           zd(ii)=snapmatrix_velocity(5,u_nodesi,snapshotsi,dimi,numphasei)
        else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3 .and. xd(3,ii) .eq. a1 .and. xd(4,ii) .eq. a3) then 
!(1.5000000000000000   1.5000000000000000   1.0000000000000000   1.5000000000000000  )
           zd(ii)=snapmatrix_velocity(6,u_nodesi,snapshotsi,dimi,numphasei)
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3 .and. xd(3,ii) .eq. a5 .and. xd(4,ii) .eq. a3) then 
!(1.5000000000000000   1.5000000000000000   2.0000000000000000   a3000000000000000  )
           zd(ii)=snapmatrix_velocity(7,u_nodesi,snapshotsi,dimi,numphasei)
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3 .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a1) then
 !(a3000000000000000   a3000000000000000   a3000000000000000   1.0000000000000000   )
           zd(ii)=snapmatrix_velocity(8,u_nodesi,snapshotsi,dimi,numphasei)
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3   .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a5) then 
!(  a3000000000000000   a3000000000000000   a3000000000000000   2.0000000000000000  )
           zd(ii)=snapmatrix_velocity(9,u_nodesi,snapshotsi,dimi,numphasei) 
         else if(int(xd(1,ii)*1000000)/1000000 .eq. a2 .and. xd(2,ii) .eq. a3 .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3) then 
               !(a24466094067263   a3000000000000000   a3000000000000000   a3000000000000000 .and. xd(1,ii) <= 1.15 ..and. xd(2:4,ii) .eq. a3  )
           zd(ii)=snapmatrix_velocity(10,u_nodesi,snapshotsi,dimi,numphasei)
         elseif(int(xd(1,ii)*1000000)/1000000 .eq. a4 .and. xd(2,ii) .eq. a3 .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3) then 
               !(a45533905932737   a3000000000000000   a3000000000000000   a3000000000000000.and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3 )
           zd(ii)=snapmatrix_velocity(11,u_nodesi,snapshotsi,dimi,numphasei)
         else if(xd(1,ii) .eq. a1 .and. xd(2,ii) .eq. a1 .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3) then !(  1.0000000000000000   1.0000000000000000   a3000000000000000   a3000000000000000 )
           zd(ii)=snapmatrix_velocity(12,u_nodesi,snapshotsi,dimi,numphasei)
         else if(xd(1,ii) .eq. a5 .and. xd(2,ii) .eq. a1 .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3) then !( 2.0000000000000000   1.0000000000000000   a3000000000000000   a3000000000000000)
           zd(ii)=snapmatrix_velocity(13,u_nodesi,snapshotsi,dimi,numphasei)
          else if(xd(1,ii) .eq. a1 .and. xd(2,ii) .eq.a5  .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3) then !(1.0000000000000000   2.0000000000000000   a3000000000000000   a3000000000000000)
           zd(ii)=snapmatrix_velocity(14,u_nodesi,snapshotsi,dimi,numphasei)
         else if(xd(1,ii) .eq. a5 .and. xd(2,ii) .eq.a5  .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3) then !(2.0000000000000000   2.0000000000000000   a3000000000000000   a3000000000000000)
           zd(ii)=snapmatrix_velocity(15,u_nodesi,snapshotsi,dimi,numphasei)
         else if(xd(1,ii) .eq. a3 .and. int(xd(2,ii)*1000000)/1000000 .eq. a2 .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3) then 
            !(a3000000000000000   a24466094067263   a3000000000000000   a3000000000000000 xd(2,ii) >=1.14  .and. xd(2,ii) <=1.15)
           zd(ii)=snapmatrix_velocity(16,u_nodesi,snapshotsi,dimi,numphasei)
         else if(xd(1,ii) .eq. a3 .and. int(xd(2,ii)*1000000)/1000000 .eq. a4 .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3) then
          !(a3000000000000000   a45533905932737   a3000000000000000   a3000000000000000)
           zd(ii)=snapmatrix_velocity(17,u_nodesi,snapshotsi,dimi,numphasei)
           else if(xd(1,ii) .eq. a1 .and. xd(2,ii) .eq. a3 .and. xd(3,ii) .eq. a1 .and. xd(4,ii) .eq. a3) then !( 1.0000000000000000   a3000000000000000   1.0000000000000000  a3000000000000000)
           zd(ii)=snapmatrix_velocity(18,u_nodesi,snapshotsi,dimi,numphasei)
         else if(xd(1,ii) .eq.a5.and. xd(2,ii) .eq. a3   .and. xd(3,ii) .eq. a1 .and. xd(4,ii) .eq. a3) then !(2.0000000000000000   a3000000000000000   1.0000000000000000   a3000000000000000)
           zd(ii)=snapmatrix_velocity(19,u_nodesi,snapshotsi,dimi,numphasei)
         else if(xd(1,ii) .eq. a1 .and. xd(2,ii) .eq. a3   .and. xd(3,ii) .eq. a5 .and. xd(4,ii) .eq. a3) then !(1.0000000000000000   a3000000000000000   2.0000000000000000   a3000000000000000)
           zd(ii)=snapmatrix_velocity(20,u_nodesi,snapshotsi,dimi,numphasei)
         else if(xd(1,ii) .eq. a5 .and. xd(2,ii) .eq. a3   .and. xd(3,ii) .eq. a5 .and. xd(4,ii) .eq. a3) then !(2.0000000000000000   a3000000000000000   2.0000000000000000   a3000000000000000)
           zd(ii)=snapmatrix_velocity(21,u_nodesi,snapshotsi,dimi,numphasei) 
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a1   .and. xd(3,ii) .eq. a1 .and. xd(4,ii) .eq. a3) then !(a3000000000000000   1.0000000000000000   1.0000000000000000   a3000000000000000)
           zd(ii)=snapmatrix_velocity(22,u_nodesi,snapshotsi,dimi,numphasei)
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a5   .and. xd(3,ii) .eq. a1 .and. xd(4,ii) .eq. a3) then !(a3000000000000000   2.0000000000000000   1.0000000000000000   a3000000000000000)
           zd(ii)=snapmatrix_velocity(23,u_nodesi,snapshotsi,dimi,numphasei)
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a1   .and. xd(3,ii) .eq. a5 .and. xd(4,ii) .eq. a3) then !(a3000000000000000   1.0000000000000000   2.0000000000000000   a3000000000000000)
           zd(ii)=snapmatrix_velocity(24,u_nodesi,snapshotsi,dimi,numphasei)
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a5   .and. xd(3,ii) .eq. a5 .and. xd(4,ii) .eq. a3) then !(a3000000000000000   2.0000000000000000   2.0000000000000000   a3000000000000000)
           zd(ii)=snapmatrix_velocity(25,u_nodesi,snapshotsi,dimi,numphasei)
          else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3 .and. int(xd(3,ii)*1000000)/1000000 .eq. a2 .and. xd(4,ii) .eq. a3) then 
       !(a3000000000000000   a3000000000000000  a24466094067263   a3000000000000000)
           zd(ii)=snapmatrix_velocity(26,u_nodesi,snapshotsi,dimi,numphasei)
         else if(xd(1,ii) .eq. a3  .and. xd(2,ii) .eq. a3  .and. int(xd(3,ii)*1000000)/1000000 .eq. a4 .and. xd(4,ii) .eq. a3) then 
       !( a3000000000000000   a3000000000000000   a45533905932737   a3000000000000000 )
           zd(ii)=snapmatrix_velocity(27,u_nodesi,snapshotsi,dimi,numphasei)
         else if(xd(1,ii) .eq. a1 .and. xd(2,ii) .eq. a3   .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a1) then !(1.0000000000000000   a3000000000000000   a3000000000000000   1.0000000000000000)
           zd(ii)=snapmatrix_velocity(28,u_nodesi,snapshotsi,dimi,numphasei)
         else if(xd(1,ii) .eq. a5 .and. xd(2,ii) .eq. a3   .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a1) then !(2.0000000000000000   a3000000000000000   a3000000000000000   1.0000000000000000)
           zd(ii)=snapmatrix_velocity(29,u_nodesi,snapshotsi,dimi,numphasei)
           else if(xd(1,ii) .eq. a1 .and. xd(2,ii) .eq. a3   .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a5) then !(1.0000000000000000   a3000000000000000   a3000000000000000   2.0000000000000000)
           zd(ii)=snapmatrix_velocity(30,u_nodesi,snapshotsi,dimi,numphasei)
         else if(xd(1,ii) .eq. a5 .and. xd(2,ii) .eq. a3   .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a5) then !(2.0000000000000000   a3000000000000000   a3000000000000000   2.0000000000000000)
           zd(ii)=snapmatrix_velocity(31,u_nodesi,snapshotsi,dimi,numphasei)
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a1   .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a1) then !(a3000000000000000   1.0000000000000000   a3000000000000000   1.0000000000000000)
           zd(ii)=snapmatrix_velocity(32,u_nodesi,snapshotsi,dimi,numphasei)
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a5  .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a1) then !(a3000000000000000   2.0000000000000000   a3000000000000000   1.0000000000000000)
           zd(ii)=snapmatrix_velocity(33,u_nodesi,snapshotsi,dimi,numphasei) 
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a1   .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a5) then !(a3000000000000000   1.0000000000000000   a3000000000000000   2.0000000000000000)
           zd(ii)=snapmatrix_velocity(34,u_nodesi,snapshotsi,dimi,numphasei)
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a5   .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a5) then !(a3000000000000000   2.0000000000000000   a3000000000000000   2.0000000000000000)
           zd(ii)=snapmatrix_velocity(35,u_nodesi,snapshotsi,dimi,numphasei)
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3   .and. xd(3,ii) .eq. a1 .and. xd(4,ii) .eq. a1) then !(a3000000000000000   a3000000000000000   1.0000000000000000   1.0000000000000000)
           zd(ii)=snapmatrix_velocity(36,u_nodesi,snapshotsi,dimi,numphasei)
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3   .and. xd(3,ii) .eq. a5 .and. xd(4,ii) .eq. a1) then !( a3000000000000000   a3000000000000000   2.0000000000000000   1.0000000000000000)
           zd(ii)=snapmatrix_velocity(37,u_nodesi,snapshotsi,dimi,numphasei)
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3   .and. xd(3,ii) .eq. a1 .and. xd(4,ii) .eq. a5) then !(a3000000000000000   a3000000000000000   1.0000000000000000   2.0000000000000000)
           zd(ii)=snapmatrix_velocity(38,u_nodesi,snapshotsi,dimi,numphasei)
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3   .and. xd(3,ii) .eq. a5 .and. xd(4,ii) .eq. a5) then !(a3000000000000000   a3000000000000000   2.0000000000000000   2.0000000000000000)
           zd(ii)=snapmatrix_velocity(39,u_nodesi,snapshotsi,dimi,numphasei)
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3 .and. xd(3,ii) .eq. a3 .and. int(xd(4,ii)*1000000)/1000000 .eq. a2) then 
      !(a3000000000000000   a3000000000000000   a3000000000000000   a24466094067263 )
           zd(ii)=snapmatrix_velocity(40,u_nodesi,snapshotsi,dimi,numphasei)
         else if(xd(1,ii) .eq. a3.and. xd(2,ii) .eq. a3 .and. xd(3,ii) .eq. a3.and. int(xd(4,ii)*1000000)/1000000 .eq. a4) then 
     !(a3000000000000000   a3000000000000000   a3000000000000000   1.8535533905932737     )
           zd(ii)=snapmatrix_velocity(41,u_nodesi,snapshotsi,dimi,numphasei)
         endif
      !  print *, 'zdzdzd', zd(ii)
        enddo
          ! call f_sinr ( m, nd, xd, zd )  
          ! print *, 'zd',  zd(:) 
          !  Use the grid to evaluate the interpolant.
        call lagrange_interp_nd_value2 ( m, ind, a, b, nd, zd, ni, xi, zpi ) 
          !  Weighted the interpolant values and add to the sparse grid interpolant. 
        nd_total = nd_total + nd
        zi(1:ni) = zi(1:ni) + c(l) * zpi(1:ni)
        !  print *, ' nd_total,level',  nd_total, l, c(l), nd
        !  print *, 'zizizizizizizic(l)zpi', zi(1)!,c(l),zpi(1)
        deallocate ( xd )
        deallocate ( zd )     
       
       if ( .not. more ) then
          exit
        end if    
 
      end do !DO

    end do !do l = l_min, l_max   
    deallocate ( c )
    deallocate ( w )
    
   end do ! do l_max = sparse_max, sparse_max
   interpolate_velocity(u_nodesi,snapshotsi,dimi,numphasei)=zi(ni)
   enddo !Do u_nodesi  
  enddo !Do snapshotsi
  enddo !dimi
  enddo !numphasei

 !  print *, ' total', total , nd_total
 ! open(unit=66,file='podbasis_interp') 
 !  do j=1, total_timestep 
 !     write(66,*)  pod_coef_all_obv_41i(j,1:num_pod)
 ! enddo
 ! close(66)
 
  deallocate ( a )
  deallocate ( b )
  deallocate ( ind )
  deallocate ( xi )
  deallocate ( ze )
  deallocate ( zi )
  deallocate ( zpi )
  end subroutine interpolate_snapmatrix_velocity_4dim

 subroutine interpolate_snapmatrix_pressure_4dim( snapmatrix_pressure, interpolate_pressure)  
    
! top-bottom between 1-2, pod_coef_all_obv is used for storing some podbasis to avoid to much modifications.
 
!  Parameters:
!
!    Local, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) SPARSE_MAX, the maximum sparse grid level to try.
!
!  Local Parameters:
!
!    Local, real ( kind = 8 ) A(M), B(M), the upper and lower variable limits 
!    in each dimension.
!
!    Local, real ( kind = 8 ) APP_ERROR, the averaged Euclidean norm of the 
!    difference between the sparse interpolant and the exact function at 
!    the interpolation points.
!
!    Local, integer ( kind = 4 ) C(0:L_MAX), the sparse grid coefficient vector.
!    Results at level L are weighted by C(L).
!
!    Local, integer ( kind = 4 ) IND(M), the 1D indices defining a Lagrange grid.
!    Each entry is a 1d "level" that specifies the order of a 
!    Clenshaw Curtis 1D grid.
!
!    Local, integer ( kind = 4 ) L, the current Lagrange grid level.
!
!    Local, integer ( kind = 4 ) L_MAX, the current sparse grid level.
!
!    Local, integer ( kind = 4 ) MORE, used to control the enumeration of all the
!    Lagrange grids at a current grid level.
!
!    Local, integer ( kind = 4 ) ND, the number of points used in a Lagrange grid.
!
!    Local, integer ( kind = 4 ) ND_TOTAL, the total number of points used in all the
!    Lagrange interpolants for a given sparse interpolant points that occur
!    in more than one Lagrange grid are counted multiple times.
!
!    Local, integer ( kind = 4 ) NI, the number of interpolant evaluation points.
!
!    Local, integer ( kind = 4 ) SPARSE_MIN, the minimum sparse grid level to try.
!
!    Local, real ( kind = 8 ) XD(M,ND), the data points for a Lagrange grid.
!
!    Local, real ( kind = 8 ) XI(M,NI), the interpolant evaluation points.
!
!    Local, real ( kind = 8 ) ZD(ND), the data values for a Lagrange grid.
!
!    Local, real ( kind = 8 ) ZE(NI), the exact function values at XI.
!
!    Local, real ( kind = 8 ) ZI(NI), the sparse interpolant values at XI.
!
!    Local, real ( kind = 8 ) ZPI(NI), one set of Lagrange interpolant values at XI.
!
  implicit none
  !type(state_type), intent(in), dimension(:,:,:) :: state
  real, dimension(:,:,:,:), allocatable :: snapmatrix_pressure ! allocate(snapmatrix_pressure(no_smolyak_nodes,p_nodes,snapshots,numphase))
  real, dimension(:,:,:), intent(out), allocatable ::  interpolate_pressure 
  real ( kind = 8 ), allocatable :: a(:)
  real ( kind = 8 ) app_error
  real ( kind = 8 ), allocatable :: b(:)
  integer ( kind = 4 ), allocatable :: c(:)
  integer ( kind = 4 ) h
  integer ( kind = 4 ), allocatable :: ind(:)
  integer ( kind = 4 ) l
  integer ( kind = 4 ) l_max
  integer ( kind = 4 ) l_min
  !integer ( kind = 4 ) m
  logical more
  integer ( kind = 4 ) nd
  integer ( kind = 4 ) nd_total
  integer ( kind = 4 ) ni
  real ( kind = 8 ) r8vec_norm_affine
  integer ( kind = 4 ) seed
 ! integer ( kind = 4 ) sparse_max
  integer ( kind = 4 ) sparse_min
  integer ( kind = 4 ) t
  integer ( kind = 4 ), allocatable :: w(:)
  real ( kind = 8 ), allocatable :: xd(:,:)
  real ( kind = 8 ), allocatable :: xi(:,:)
  real ( kind = 8 ), allocatable :: zd(:)
  real ( kind = 8 ), allocatable :: ze(:)
  real ( kind = 8 ), allocatable :: zi(:)
  real ( kind = 8 ), allocatable :: zpi(:)
  integer :: i, ii, total, total_timestep,num_pod,j,k,jj,kk, smolyak_total_points,u_nodesi,snapshotsi,dimi,numphasei
   real :: a1,a2,a3,a4,a5 
  allocate (ind(1:m))

  ! Define the region.
  allocate ( a(1:m) )
  allocate ( b(1:m) )

  
   
  !total_timestep=864
  !num_pod=12*36
  !smolyak_total_points=41!2**sparse_max+1
  !print *, '  smolyak_total_points=2**sparse_max+1=', smolyak_total_points 
  !  stop 2323
  !  Define the interpolation evaluation information.
  ni = 1
  seed = 123456789
  allocate ( xi(m,ni) )
  !call r8mat_uniform_abvec ( m, ni, a, b, seed, xi )
!  print * , 'xi', xi
 
   if (.false.) then
  a(1:m) = 1
  b(1:m) = 2
  a1=1
  a2=1.1464!5780000000000
  a3=1.5
  a4=1.8535!421000
  a5=2
  xi(1,1)= 2 
  xi(2,1)= 1.5
  xi(3,1)= 2
  xi(4,1)= 1.55 !1.7 !2.0000000000000000   1.5000000000000000   2.0000000000000000   1.5000000000000000,21
  endif

  if (.true.) then
  a(1:m) = 0.1
  b(1:m) = 0.5
  a1=0.1
  a2=0.158578
  a3=0.3
  a4=0.441421
  a5=0.5
  xi(1,1)= 0.1 
  xi(2,1)= 0.5
  xi(3,1)= 0.5
  xi(4,1)= 0.5 !1.7 !2.0000000000000000   1.5000000000000000   2.0000000000000000   1.5000000000000000,21
   call get_option(&
         '/reduced_model/variables_to_interpolate/dimension_one', xi(1,1))
   call get_option(&
         '/reduced_model/variables_to_interpolate/dimension_two', xi(2,1))
 call get_option(&
         '/reduced_model/variables_to_interpolate/dimension_three', xi(3,1))
 call get_option(&
         '/reduced_model/variables_to_interpolate/dimension_four', xi(4,1))
  endif
  write ( *, '(a,i6)' ) '  Number of interpolation points is NI = ', ni
  allocate ( ze(1:ni) )
 ! call f_sinr ( m, ni, xi, ze )
  print *, 'actual function value is, will interpolate this value later on,then do comparison', xi,ze
!
!  Compute a sequence of sparse grid interpolants of increasing level.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '   L     Order    ApproxError'
  write ( *, '(a)' ) ''
  allocate(interpolate_pressure(p_nodes,snapshots,numphase)) 
  allocate ( zpi(1:ni) )
  allocate ( zi(1:ni) )
  !smolynode, u_nodes,snapshots,dim,numphase
  Do u_nodesi=1, p_nodes  
   Do snapshotsi=1, snapshots  
   do numphasei=1, numphase
  total = 0
  sparse_min = 0

  ! do l_max = sparse_min, sparse_max
    do l_max = sparse_max, sparse_max
    allocate ( c(0:l_max) )
    allocate ( w(0:l_max) )
    call smolyak_coefficients ( l_max, m, c, w )

    zi(1:ni) = 0.0D+00
    nd_total = 0
    
    l_min = max ( l_max + 1 - m, 0 )

    do l = l_min, l_max

      more = .false.
      
      do
     ! print *,'l_max,l,' ,l_max,l   
!  Define the next product grid.
!
        call comp_next ( l, m, ind, more, h, t )
!
!  Count the grid, find its points, evaluate the data there.
!
        call lagrange_interp_nd_size2 ( m, ind, nd )
        allocate ( xd(1:m,1:nd) )
        call lagrange_interp_nd_grid2 ( m, ind, a, b, nd, xd )        
         !  print *,'xd', xd, 点1 坐标为： xd(:,1)
        allocate ( zd(1:nd) )

       ! print *,'---------------------nd',nd  
        do ii=1, nd           
                if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3   .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3) then 
!( 1.5000000000000000    1.5000000000000000   1.5000000000000000  1.5000000000000000)
           zd(ii)=snapmatrix_pressure(1,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a1 .and. xd(2,ii) .eq. a3   .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3) then 
!(1.0000000000000000   1.5000000000000000   1.5000000000000000   1.5000000000000000)
           zd(ii)=snapmatrix_pressure(2,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a5 .and. xd(2,ii) .eq. a3   .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3) then
 !( 2.0000000000000000   1.5000000000000000   1.5000000000000000   1.5000000000000000)
           zd(ii)=snapmatrix_pressure(3,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a1   .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3) then 
!(1.5000000000000000   1.0000000000000000   1.5000000000000000   1.5000000000000000 )
           zd(ii)=snapmatrix_pressure(4,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a5 .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3) then 
!(1.5000000000000000   2.0000000000000000   1.5000000000000000   1.5000000000000000 )
           zd(ii)=snapmatrix_pressure(5,u_nodesi,snapshotsi, numphasei)
        else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3 .and. xd(3,ii) .eq. a1 .and. xd(4,ii) .eq. a3) then 
!(1.5000000000000000   1.5000000000000000   1.0000000000000000   1.5000000000000000  )
           zd(ii)=snapmatrix_pressure(6,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3 .and. xd(3,ii) .eq. a5 .and. xd(4,ii) .eq. a3) then 
!(1.5000000000000000   1.5000000000000000   2.0000000000000000   a3000000000000000  )
           zd(ii)=snapmatrix_pressure(7,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3 .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a1) then
 !(a3000000000000000   a3000000000000000   a3000000000000000   1.0000000000000000   )
           zd(ii)=snapmatrix_pressure(8,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3   .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a5) then 
!(  a3000000000000000   a3000000000000000   a3000000000000000   2.0000000000000000  )
           zd(ii)=snapmatrix_pressure(9,u_nodesi,snapshotsi, numphasei) 
         else if(int(xd(1,ii) *1000000)/1000000 .eq. a2 .and. xd(2,ii) .eq. a3 .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3) then 
               !(a24466094067263   a3000000000000000   a3000000000000000   a3000000000000000 .and. xd(1,ii) <= 1.15 ..and. xd(2:4,ii) .eq. a3  )
           zd(ii)=snapmatrix_pressure(10,u_nodesi,snapshotsi, numphasei)
         elseif(int(xd(1,ii) *1000000)/1000000 .eq. a4 .and. xd(2,ii) .eq. a3 .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3) then 
               !(a45533905932737   a3000000000000000   a3000000000000000   a3000000000000000.and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3 )
           zd(ii)=snapmatrix_pressure(11,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a1 .and. xd(2,ii) .eq. a1 .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3) then !(  1.0000000000000000   1.0000000000000000   a3000000000000000   a3000000000000000 )
           zd(ii)=snapmatrix_pressure(12,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a5 .and. xd(2,ii) .eq. a1 .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3) then !( 2.0000000000000000   1.0000000000000000   a3000000000000000   a3000000000000000)
           zd(ii)=snapmatrix_pressure(13,u_nodesi,snapshotsi, numphasei)
          else if(xd(1,ii) .eq. a1 .and. xd(2,ii) .eq.a5  .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3) then !(1.0000000000000000   2.0000000000000000   a3000000000000000   a3000000000000000)
           zd(ii)=snapmatrix_pressure(14,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a5 .and. xd(2,ii) .eq.a5  .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3) then !(2.0000000000000000   2.0000000000000000   a3000000000000000   a3000000000000000)
           zd(ii)=snapmatrix_pressure(15,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a3 .and. int(xd(2,ii) *1000000)/1000000 .eq. a2 .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3) then 
            !(a3000000000000000   a24466094067263   a3000000000000000   a3000000000000000 xd(2,ii) >=1.14  .and. xd(2,ii) <=1.15)
           zd(ii)=snapmatrix_pressure(16,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a3 .and. int(xd(2,ii) *1000000)/1000000 .eq. a4 .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3) then
          !(a3000000000000000   a45533905932737   a3000000000000000   a3000000000000000)
           zd(ii)=snapmatrix_pressure(17,u_nodesi,snapshotsi, numphasei)
           else if(xd(1,ii) .eq. a1 .and. xd(2,ii) .eq. a3 .and. xd(3,ii) .eq. a1 .and. xd(4,ii) .eq. a3) then !( 1.0000000000000000   a3000000000000000   1.0000000000000000  a3000000000000000)
           zd(ii)=snapmatrix_pressure(18,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq.a5.and. xd(2,ii) .eq. a3   .and. xd(3,ii) .eq. a1 .and. xd(4,ii) .eq. a3) then !(2.0000000000000000   a3000000000000000   1.0000000000000000   a3000000000000000)
           zd(ii)=snapmatrix_pressure(19,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a1 .and. xd(2,ii) .eq. a3   .and. xd(3,ii) .eq. a5 .and. xd(4,ii) .eq. a3) then !(1.0000000000000000   a3000000000000000   2.0000000000000000   a3000000000000000)
           zd(ii)=snapmatrix_pressure(20,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a5 .and. xd(2,ii) .eq. a3   .and. xd(3,ii) .eq. a5 .and. xd(4,ii) .eq. a3) then !(2.0000000000000000   a3000000000000000   2.0000000000000000   a3000000000000000)
           zd(ii)=snapmatrix_pressure(21,u_nodesi,snapshotsi, numphasei) 
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a1   .and. xd(3,ii) .eq. a1 .and. xd(4,ii) .eq. a3) then !(a3000000000000000   1.0000000000000000   1.0000000000000000   a3000000000000000)
           zd(ii)=snapmatrix_pressure(22,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a5   .and. xd(3,ii) .eq. a1 .and. xd(4,ii) .eq. a3) then !(a3000000000000000   2.0000000000000000   1.0000000000000000   a3000000000000000)
           zd(ii)=snapmatrix_pressure(23,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a1   .and. xd(3,ii) .eq. a5 .and. xd(4,ii) .eq. a3) then !(a3000000000000000   1.0000000000000000   2.0000000000000000   a3000000000000000)
           zd(ii)=snapmatrix_pressure(24,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a5   .and. xd(3,ii) .eq. a5 .and. xd(4,ii) .eq. a3) then !(a3000000000000000   2.0000000000000000   2.0000000000000000   a3000000000000000)
           zd(ii)=snapmatrix_pressure(25,u_nodesi,snapshotsi, numphasei)
          else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3 .and. int(xd(3,ii) *1000000)/1000000 .eq. a2 .and. xd(4,ii) .eq. a3) then 
       !(a3000000000000000   a3000000000000000  a24466094067263   a3000000000000000)
           zd(ii)=snapmatrix_pressure(26,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a3  .and. xd(2,ii) .eq. a3  .and. int(xd(3,ii) *1000000)/1000000 .eq. a4 .and. xd(4,ii) .eq. a3) then 
       !( a3000000000000000   a3000000000000000   a45533905932737   a3000000000000000 )
           zd(ii)=snapmatrix_pressure(27,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a1 .and. xd(2,ii) .eq. a3   .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a1) then !(1.0000000000000000   a3000000000000000   a3000000000000000   1.0000000000000000)
           zd(ii)=snapmatrix_pressure(28,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a5 .and. xd(2,ii) .eq. a3   .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a1) then !(2.0000000000000000   a3000000000000000   a3000000000000000   1.0000000000000000)
           zd(ii)=snapmatrix_pressure(29,u_nodesi,snapshotsi, numphasei)
           else if(xd(1,ii) .eq. a1 .and. xd(2,ii) .eq. a3   .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a5) then !(1.0000000000000000   a3000000000000000   a3000000000000000   2.0000000000000000)
           zd(ii)=snapmatrix_pressure(30,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a5 .and. xd(2,ii) .eq. a3   .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a5) then !(2.0000000000000000   a3000000000000000   a3000000000000000   2.0000000000000000)
           zd(ii)=snapmatrix_pressure(31,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a1   .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a1) then !(a3000000000000000   1.0000000000000000   a3000000000000000   1.0000000000000000)
           zd(ii)=snapmatrix_pressure(32,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a5  .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a1) then !(a3000000000000000   2.0000000000000000   a3000000000000000   1.0000000000000000)
           zd(ii)=snapmatrix_pressure(33,u_nodesi,snapshotsi, numphasei) 
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a1   .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a5) then !(a3000000000000000   1.0000000000000000   a3000000000000000   2.0000000000000000)
           zd(ii)=snapmatrix_pressure(34,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a5   .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a5) then !(a3000000000000000   2.0000000000000000   a3000000000000000   2.0000000000000000)
           zd(ii)=snapmatrix_pressure(35,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3   .and. xd(3,ii) .eq. a1 .and. xd(4,ii) .eq. a1) then !(a3000000000000000   a3000000000000000   1.0000000000000000   1.0000000000000000)
           zd(ii)=snapmatrix_pressure(36,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3   .and. xd(3,ii) .eq. a5 .and. xd(4,ii) .eq. a1) then !( a3000000000000000   a3000000000000000   2.0000000000000000   1.0000000000000000)
           zd(ii)=snapmatrix_pressure(37,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3   .and. xd(3,ii) .eq. a1 .and. xd(4,ii) .eq. a5) then !(a3000000000000000   a3000000000000000   1.0000000000000000   2.0000000000000000)
           zd(ii)=snapmatrix_pressure(38,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3   .and. xd(3,ii) .eq. a5 .and. xd(4,ii) .eq. a5) then !(a3000000000000000   a3000000000000000   2.0000000000000000   2.0000000000000000)
           zd(ii)=snapmatrix_pressure(39,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3 .and. xd(3,ii) .eq. a3 .and. int(xd(4,ii) *1000000)/1000000 .eq. a2) then 
      !(a3000000000000000   a3000000000000000   a3000000000000000   a24466094067263 )
           zd(ii)=snapmatrix_pressure(40,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a3.and. xd(2,ii) .eq. a3 .and. xd(3,ii) .eq. a3.and. int(xd(4,ii) *1000000)/1000000 .eq. a4) then 
     !(a3000000000000000   a3000000000000000   a3000000000000000   1.8535533905932737     )
           zd(ii)=snapmatrix_pressure(41,u_nodesi,snapshotsi, numphasei)
         endif
        enddo
          !call f_sinr ( m, nd, xd, zd )       
           
          ! print *, 'zd',  zd(:)
            
!  Use the grid to evaluate the interpolant.
        call lagrange_interp_nd_value2 ( m, ind, a, b, nd, zd, ni, xi, zpi )
 
!  Weighted the interpolant values and add to the sparse grid interpolant.
 
        nd_total = nd_total + nd
        zi(1:ni) = zi(1:ni) + c(l) * zpi(1:ni)
        !  print *, ' nd_total,level',  nd_total, l, c(l), nd
        !  print *, 'zizizizizizizic(l)zpi', zi(1)!,c(l),zpi(1)
        deallocate ( xd )
        deallocate ( zd )     
       
       if ( .not. more ) then
          exit
        end if    
 
      end do !DO

    end do !do l = l_min, l_max   
    deallocate ( c )
    deallocate ( w )
    
   end do ! do l_max = sparse_max, sparse_max
   interpolate_pressure(u_nodesi,snapshotsi,numphasei)=zi(ni)
   enddo !Do u_nodesi  
  enddo !Do snapshotsi 
  enddo !numphasei

 !  print *, ' total', total , nd_total
 ! open(unit=66,file='podbasis_interp') 
 !  do j=1, total_timestep 
 !     write(66,*)  pod_coef_all_obv_41i(j,1:num_pod)
 ! enddo
 ! close(66)
 
  deallocate ( a )
  deallocate ( b )
  deallocate ( ind )
  deallocate ( xi )
  deallocate ( ze )
  deallocate ( zi )
  deallocate ( zpi )
  end subroutine interpolate_snapmatrix_pressure_4dim

 subroutine interpolate_snapmatrix_volumefraction_4dim( snapmatrix_volumefraction, interpolate_volumefraction)  
    
! top-bottom between 1-2, pod_coef_all_obv is used for storing some podbasis to avoid to much modifications.
 
!  Parameters:
!
!    Local, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) SPARSE_MAX, the maximum sparse grid level to try.
!
!  Local Parameters:
!
!    Local, real ( kind = 8 ) A(M), B(M), the upper and lower variable limits 
!    in each dimension.
!
!    Local, real ( kind = 8 ) APP_ERROR, the averaged Euclidean norm of the 
!    difference between the sparse interpolant and the exact function at 
!    the interpolation points.
!
!    Local, integer ( kind = 4 ) C(0:L_MAX), the sparse grid coefficient vector.
!    Results at level L are weighted by C(L).
!
!    Local, integer ( kind = 4 ) IND(M), the 1D indices defining a Lagrange grid.
!    Each entry is a 1d "level" that specifies the order of a 
!    Clenshaw Curtis 1D grid.
!
!    Local, integer ( kind = 4 ) L, the current Lagrange grid level.
!
!    Local, integer ( kind = 4 ) L_MAX, the current sparse grid level.
!
!    Local, integer ( kind = 4 ) MORE, used to control the enumeration of all the
!    Lagrange grids at a current grid level.
!
!    Local, integer ( kind = 4 ) ND, the number of points used in a Lagrange grid.
!
!    Local, integer ( kind = 4 ) ND_TOTAL, the total number of points used in all the
!    Lagrange interpolants for a given sparse interpolant points that occur
!    in more than one Lagrange grid are counted multiple times.
!
!    Local, integer ( kind = 4 ) NI, the number of interpolant evaluation points.
!
!    Local, integer ( kind = 4 ) SPARSE_MIN, the minimum sparse grid level to try.
!
!    Local, real ( kind = 8 ) XD(M,ND), the data points for a Lagrange grid.
!
!    Local, real ( kind = 8 ) XI(M,NI), the interpolant evaluation points.
!
!    Local, real ( kind = 8 ) ZD(ND), the data values for a Lagrange grid.
!
!    Local, real ( kind = 8 ) ZE(NI), the exact function values at XI.
!
!    Local, real ( kind = 8 ) ZI(NI), the sparse interpolant values at XI.
!
!    Local, real ( kind = 8 ) ZPI(NI), one set of Lagrange interpolant values at XI.
!
  implicit none
  !type(state_type), intent(in), dimension(:,:,:) :: state
  real, dimension(:,:,:,:), allocatable :: snapmatrix_volumefraction ! allocate(snapmatrix_pressure(no_smolyak_nodes,p_nodes,snapshots,numphase))
  real, dimension(:,:,:), intent(out), allocatable ::  interpolate_volumefraction 
  real ( kind = 8 ), allocatable :: a(:)
  real ( kind = 8 ) app_error
  real ( kind = 8 ), allocatable :: b(:)
  integer ( kind = 4 ), allocatable :: c(:)
  integer ( kind = 4 ) h
  integer ( kind = 4 ), allocatable :: ind(:)
  integer ( kind = 4 ) l
  integer ( kind = 4 ) l_max
  integer ( kind = 4 ) l_min
  !integer ( kind = 4 ) m
  logical more
  integer ( kind = 4 ) nd
  integer ( kind = 4 ) nd_total
  integer ( kind = 4 ) ni
  real ( kind = 8 ) r8vec_norm_affine
  integer ( kind = 4 ) seed
  !integer ( kind = 4 ) sparse_max
  integer ( kind = 4 ) sparse_min
  integer ( kind = 4 ) t
  integer ( kind = 4 ), allocatable :: w(:)
  real ( kind = 8 ), allocatable :: xd(:,:)
  real ( kind = 8 ), allocatable :: xi(:,:)
  real ( kind = 8 ), allocatable :: zd(:)
  real ( kind = 8 ), allocatable :: ze(:)
  real ( kind = 8 ), allocatable :: zi(:)
  real ( kind = 8 ), allocatable :: zpi(:)
  integer :: i, ii, total, total_timestep,num_pod,j,k,jj,kk, smolyak_total_points,u_nodesi,snapshotsi,dimi,numphasei
   real :: a1,a2,a3,a4,a5  
  allocate (ind(1:m))

  ! Define the region.
  allocate ( a(1:m) )
  allocate ( b(1:m) )

 
   
  !total_timestep=864
  !num_pod=12*36
  !smolyak_total_points=41!2**sparse_max+1
  !print *, '  smolyak_total_points=2**sparse_max+1=', smolyak_total_points 
  
  !  Define the interpolation evaluation information.
  ni = 1
  seed = 123456789
  allocate ( xi(m,ni) )
  !call r8mat_uniform_abvec ( m, ni, a, b, seed, xi )
!  print * , 'xi', xi
   if (.false.) then
  a(1:m) = 1
  b(1:m) = 2
  a1=1
  a2=1.1464!5780000000000
  a3=1.5
  a4=1.8535!421000
  a5=2
  xi(1,1)= 2 
  xi(2,1)= 1.5
  xi(3,1)= 2
  xi(4,1)= 1.55 !1.7 !2.0000000000000000   1.5000000000000000   2.0000000000000000   1.5000000000000000,21
  endif

  if (.true.) then
  a(1:m) = 0.1
  b(1:m) = 0.5
  a1=0.1
  a2=0.158578
  a3=0.3
  a4=0.441421
  a5=0.5
  xi(1,1)= 0.1 
  xi(2,1)= 0.5
  xi(3,1)= 0.5
  xi(4,1)= 0.5 !1.7 !2.0000000000000000   1.5000000000000000   2.0000000000000000   1.5000000000000000,21
   call get_option(&
         '/reduced_model/variables_to_interpolate/dimension_one', xi(1,1))
   call get_option(&
         '/reduced_model/variables_to_interpolate/dimension_two', xi(2,1))
 call get_option(&
         '/reduced_model/variables_to_interpolate/dimension_three', xi(3,1))
 call get_option(&
         '/reduced_model/variables_to_interpolate/dimension_four', xi(4,1))
  endif
  write ( *, '(a,i6)' ) '  Number of interpolation points is NI = ', ni
  allocate ( ze(1:ni) )
 ! call f_sinr ( m, ni, xi, ze )
  print *, 'actual function value is, will interpolate this value later on,then do comparison', xi,ze
!
!  Compute a sequence of sparse grid interpolants of increasing level.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '   L     Order    ApproxError'
  write ( *, '(a)' ) ''

  allocate ( zpi(1:ni) )
  allocate ( zi(1:ni) )
  allocate(interpolate_volumefraction(p_nodes,snapshots,numphase))
  !smolynode, u_nodes,snapshots,dim,numphase
  Do u_nodesi=1, p_nodes  
   Do snapshotsi=1, snapshots  
   do numphasei=1, numphase
  total = 0
  sparse_min = 0

  ! do l_max = sparse_min, sparse_max
    do l_max = sparse_max, sparse_max
    allocate ( c(0:l_max) )
    allocate ( w(0:l_max) )
    call smolyak_coefficients ( l_max, m, c, w )

    zi(1:ni) = 0.0D+00
    nd_total = 0
    
    l_min = max ( l_max + 1 - m, 0 )

    do l = l_min, l_max

      more = .false.
      
      do
     ! print *,'l_max,l,' ,l_max,l   
!  Define the next product grid.
!
        call comp_next ( l, m, ind, more, h, t )
!
!  Count the grid, find its points, evaluate the data there.
!
        call lagrange_interp_nd_size2 ( m, ind, nd )
        allocate ( xd(1:m,1:nd) )
        call lagrange_interp_nd_grid2 ( m, ind, a, b, nd, xd )        
         !  print *,'xd', xd, 点1 坐标为： xd(:,1)
        allocate ( zd(1:nd) )

       ! print *,'---------------------nd',nd  
        do ii=1, nd           
         if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3   .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3) then 
!( 1.5000000000000000    1.5000000000000000   1.5000000000000000  1.5000000000000000)
           zd(ii)=snapmatrix_volumefraction(1,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a1 .and. xd(2,ii) .eq. a3   .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3) then 
!(1.0000000000000000   1.5000000000000000   1.5000000000000000   1.5000000000000000)
           zd(ii)=snapmatrix_volumefraction(2,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a5 .and. xd(2,ii) .eq. a3   .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3) then
 !( 2.0000000000000000   1.5000000000000000   1.5000000000000000   1.5000000000000000)
           zd(ii)=snapmatrix_volumefraction(3,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a1   .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3) then 
!(1.5000000000000000   1.0000000000000000   1.5000000000000000   1.5000000000000000 )
           zd(ii)=snapmatrix_volumefraction(4,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a5 .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3) then 
!(1.5000000000000000   2.0000000000000000   1.5000000000000000   1.5000000000000000 )
           zd(ii)=snapmatrix_volumefraction(5,u_nodesi,snapshotsi, numphasei)
        else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3 .and. xd(3,ii) .eq. a1 .and. xd(4,ii) .eq. a3) then 
!(1.5000000000000000   1.5000000000000000   1.0000000000000000   1.5000000000000000  )
           zd(ii)=snapmatrix_volumefraction(6,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3 .and. xd(3,ii) .eq. a5 .and. xd(4,ii) .eq. a3) then 
!(1.5000000000000000   1.5000000000000000   2.0000000000000000   a3000000000000000  )
           zd(ii)=snapmatrix_volumefraction(7,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3 .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a1) then
 !(a3000000000000000   a3000000000000000   a3000000000000000   1.0000000000000000   )
           zd(ii)=snapmatrix_volumefraction(8,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3   .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a5) then 
!(  a3000000000000000   a3000000000000000   a3000000000000000   2.0000000000000000  )
           zd(ii)=snapmatrix_volumefraction(9,u_nodesi,snapshotsi, numphasei) 
         else if(int(xd(1,ii) *1000000)/1000000 .eq. a2 .and. xd(2,ii) .eq. a3 .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3) then 
               !(a24466094067263   a3000000000000000   a3000000000000000   a3000000000000000 .and. xd(1,ii) <= 1.15 ..and. xd(2:4,ii) .eq. a3  )
           zd(ii)=snapmatrix_volumefraction(10,u_nodesi,snapshotsi, numphasei)
         elseif(int(xd(1,ii) *1000000)/1000000 .eq. a4 .and. xd(2,ii) .eq. a3 .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3) then 
               !(a45533905932737   a3000000000000000   a3000000000000000   a3000000000000000.and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3 )
           zd(ii)=snapmatrix_volumefraction(11,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a1 .and. xd(2,ii) .eq. a1 .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3) then !(  1.0000000000000000   1.0000000000000000   a3000000000000000   a3000000000000000 )
           zd(ii)=snapmatrix_volumefraction(12,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a5 .and. xd(2,ii) .eq. a1 .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3) then !( 2.0000000000000000   1.0000000000000000   a3000000000000000   a3000000000000000)
           zd(ii)=snapmatrix_volumefraction(13,u_nodesi,snapshotsi, numphasei)
          else if(xd(1,ii) .eq. a1 .and. xd(2,ii) .eq.a5  .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3) then !(1.0000000000000000   2.0000000000000000   a3000000000000000   a3000000000000000)
           zd(ii)=snapmatrix_volumefraction(14,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a5 .and. xd(2,ii) .eq.a5  .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3) then !(2.0000000000000000   2.0000000000000000   a3000000000000000   a3000000000000000)
           zd(ii)=snapmatrix_volumefraction(15,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a3 .and. int(xd(2,ii) *1000000)/1000000 .eq. a2 .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3) then 
            !(a3000000000000000   a24466094067263   a3000000000000000   a3000000000000000 xd(2,ii) >=1.14  .and. xd(2,ii) <=1.15)
           zd(ii)=snapmatrix_volumefraction(16,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a3 .and. int(xd(2,ii) *1000000)/1000000 .eq. a4 .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3) then
          !(a3000000000000000   a45533905932737   a3000000000000000   a3000000000000000)
           zd(ii)=snapmatrix_volumefraction(17,u_nodesi,snapshotsi, numphasei)
           else if(xd(1,ii) .eq. a1 .and. xd(2,ii) .eq. a3 .and. xd(3,ii) .eq. a1 .and. xd(4,ii) .eq. a3) then !( 1.0000000000000000   a3000000000000000   1.0000000000000000  a3000000000000000)
           zd(ii)=snapmatrix_volumefraction(18,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq.a5.and. xd(2,ii) .eq. a3   .and. xd(3,ii) .eq. a1 .and. xd(4,ii) .eq. a3) then !(2.0000000000000000   a3000000000000000   1.0000000000000000   a3000000000000000)
           zd(ii)=snapmatrix_volumefraction(19,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a1 .and. xd(2,ii) .eq. a3   .and. xd(3,ii) .eq. a5 .and. xd(4,ii) .eq. a3) then !(1.0000000000000000   a3000000000000000   2.0000000000000000   a3000000000000000)
           zd(ii)=snapmatrix_volumefraction(20,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a5 .and. xd(2,ii) .eq. a3   .and. xd(3,ii) .eq. a5 .and. xd(4,ii) .eq. a3) then !(2.0000000000000000   a3000000000000000   2.0000000000000000   a3000000000000000)
           zd(ii)=snapmatrix_volumefraction(21,u_nodesi,snapshotsi, numphasei) 
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a1   .and. xd(3,ii) .eq. a1 .and. xd(4,ii) .eq. a3) then !(a3000000000000000   1.0000000000000000   1.0000000000000000   a3000000000000000)
           zd(ii)=snapmatrix_volumefraction(22,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a5   .and. xd(3,ii) .eq. a1 .and. xd(4,ii) .eq. a3) then !(a3000000000000000   2.0000000000000000   1.0000000000000000   a3000000000000000)
           zd(ii)=snapmatrix_volumefraction(23,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a1   .and. xd(3,ii) .eq. a5 .and. xd(4,ii) .eq. a3) then !(a3000000000000000   1.0000000000000000   2.0000000000000000   a3000000000000000)
           zd(ii)=snapmatrix_volumefraction(24,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a5   .and. xd(3,ii) .eq. a5 .and. xd(4,ii) .eq. a3) then !(a3000000000000000   2.0000000000000000   2.0000000000000000   a3000000000000000)
           zd(ii)=snapmatrix_volumefraction(25,u_nodesi,snapshotsi, numphasei)
          else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3 .and. int(xd(3,ii) *1000000)/1000000 .eq. a2 .and. xd(4,ii) .eq. a3) then 
       !(a3000000000000000   a3000000000000000  a24466094067263   a3000000000000000)
           zd(ii)=snapmatrix_volumefraction(26,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a3  .and. xd(2,ii) .eq. a3  .and. int(xd(3,ii) *1000000)/1000000 .eq. a4 .and. xd(4,ii) .eq. a3) then 
       !( a3000000000000000   a3000000000000000   a45533905932737   a3000000000000000 )
           zd(ii)=snapmatrix_volumefraction(27,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a1 .and. xd(2,ii) .eq. a3   .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a1) then !(1.0000000000000000   a3000000000000000   a3000000000000000   1.0000000000000000)
           zd(ii)=snapmatrix_volumefraction(28,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a5 .and. xd(2,ii) .eq. a3   .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a1) then !(2.0000000000000000   a3000000000000000   a3000000000000000   1.0000000000000000)
           zd(ii)=snapmatrix_volumefraction(29,u_nodesi,snapshotsi, numphasei)
           else if(xd(1,ii) .eq. a1 .and. xd(2,ii) .eq. a3   .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a5) then !(1.0000000000000000   a3000000000000000   a3000000000000000   2.0000000000000000)
           zd(ii)=snapmatrix_volumefraction(30,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a5 .and. xd(2,ii) .eq. a3   .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a5) then !(2.0000000000000000   a3000000000000000   a3000000000000000   2.0000000000000000)
           zd(ii)=snapmatrix_volumefraction(31,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a1   .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a1) then !(a3000000000000000   1.0000000000000000   a3000000000000000   1.0000000000000000)
           zd(ii)=snapmatrix_volumefraction(32,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a5  .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a1) then !(a3000000000000000   2.0000000000000000   a3000000000000000   1.0000000000000000)
           zd(ii)=snapmatrix_volumefraction(33,u_nodesi,snapshotsi, numphasei) 
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a1   .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a5) then !(a3000000000000000   1.0000000000000000   a3000000000000000   2.0000000000000000)
           zd(ii)=snapmatrix_volumefraction(34,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a5   .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a5) then !(a3000000000000000   2.0000000000000000   a3000000000000000   2.0000000000000000)
           zd(ii)=snapmatrix_volumefraction(35,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3   .and. xd(3,ii) .eq. a1 .and. xd(4,ii) .eq. a1) then !(a3000000000000000   a3000000000000000   1.0000000000000000   1.0000000000000000)
           zd(ii)=snapmatrix_volumefraction(36,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3   .and. xd(3,ii) .eq. a5 .and. xd(4,ii) .eq. a1) then !( a3000000000000000   a3000000000000000   2.0000000000000000   1.0000000000000000)
           zd(ii)=snapmatrix_volumefraction(37,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3   .and. xd(3,ii) .eq. a1 .and. xd(4,ii) .eq. a5) then !(a3000000000000000   a3000000000000000   1.0000000000000000   2.0000000000000000)
           zd(ii)=snapmatrix_volumefraction(38,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3   .and. xd(3,ii) .eq. a5 .and. xd(4,ii) .eq. a5) then !(a3000000000000000   a3000000000000000   2.0000000000000000   2.0000000000000000)
           zd(ii)=snapmatrix_volumefraction(39,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3 .and. xd(3,ii) .eq. a3 .and. int(xd(4,ii) *1000000)/1000000 .eq. a2) then 
      !(a3000000000000000   a3000000000000000   a3000000000000000   a24466094067263 )
           zd(ii)=snapmatrix_volumefraction(40,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a3.and. xd(2,ii) .eq. a3 .and. xd(3,ii) .eq. a3.and. int(xd(4,ii) *1000000)/1000000 .eq. a4) then 
     !(a3000000000000000   a3000000000000000   a3000000000000000   1.8535533905932737     )
           zd(ii)=snapmatrix_volumefraction(41,u_nodesi,snapshotsi, numphasei)
         endif
        enddo
          !call f_sinr ( m, nd, xd, zd )       
           
          ! print *, 'zd',  zd(:)
            
!  Use the grid to evaluate the interpolant.
        call lagrange_interp_nd_value2 ( m, ind, a, b, nd, zd, ni, xi, zpi )
 
!  Weighted the interpolant values and add to the sparse grid interpolant.
 
        nd_total = nd_total + nd
        zi(1:ni) = zi(1:ni) + c(l) * zpi(1:ni)
        !  print *, ' nd_total,level',  nd_total, l, c(l), nd
        !  print *, 'zizizizizizizic(l)zpi', zi(1)!,c(l),zpi(1)
        deallocate ( xd )
        deallocate ( zd )     
       
       if ( .not. more ) then
          exit
        end if    
 
      end do !DO

    end do !do l = l_min, l_max   
    deallocate ( c )
    deallocate ( w )
    
   end do ! do l_max = sparse_max, sparse_max
   interpolate_volumefraction(u_nodesi,snapshotsi,numphasei)=zi(ni)
   enddo !Do u_nodesi  
  enddo !Do snapshotsi 
  enddo !numphasei

 !  print *, ' total', total , nd_total
 ! open(unit=66,file='podbasis_interp') 
 !  do j=1, total_timestep 
 !     write(66,*)  pod_coef_all_obv_41i(j,1:num_pod)
 ! enddo
 ! close(66)
 
  deallocate ( a )
  deallocate ( b )
  deallocate ( ind )
  deallocate ( xi )
  deallocate ( ze )
  deallocate ( zi )
  deallocate ( zpi )
  end subroutine interpolate_snapmatrix_volumefraction_4dim

subroutine interpolate_snapmatrix_velocity_4dim_level1( snapmatrix_velocity, interpolate_velocity)  
    
! top-bottom between 1-2, pod_coef_all_obv is used for storing some podbasis to avoid to much modifications.
 
!  Parameters:
!
!    Local, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) SPARSE_MAX, the maximum sparse grid level to try.
!
!  Local Parameters:
!
!    Local, real ( kind = 8 ) A(M), B(M), the upper and lower variable limits 
!    in each dimension.
!
!    Local, real ( kind = 8 ) APP_ERROR, the averaged Euclidean norm of the 
!    difference between the sparse interpolant and the exact function at 
!    the interpolation points.
!
!    Local, integer ( kind = 4 ) C(0:L_MAX), the sparse grid coefficient vector.
!    Results at level L are weighted by C(L).
!
!    Local, integer ( kind = 4 ) IND(M), the 1D indices defining a Lagrange grid.
!    Each entry is a 1d "level" that specifies the order of a 
!    Clenshaw Curtis 1D grid.
!
!    Local, integer ( kind = 4 ) L, the current Lagrange grid level.
!
!    Local, integer ( kind = 4 ) L_MAX, the current sparse grid level.
!
!    Local, integer ( kind = 4 ) MORE, used to control the enumeration of all the
!    Lagrange grids at a current grid level.
!
!    Local, integer ( kind = 4 ) ND, the number of points used in a Lagrange grid.
!
!    Local, integer ( kind = 4 ) ND_TOTAL, the total number of points used in all the
!    Lagrange interpolants for a given sparse interpolant points that occur
!    in more than one Lagrange grid are counted multiple times.
!
!    Local, integer ( kind = 4 ) NI, the number of interpolant evaluation points.
!
!    Local, integer ( kind = 4 ) SPARSE_MIN, the minimum sparse grid level to try.
!
!    Local, real ( kind = 8 ) XD(M,ND), the data points for a Lagrange grid.
!
!    Local, real ( kind = 8 ) XI(M,NI), the interpolant evaluation points.
!
!    Local, real ( kind = 8 ) ZD(ND), the data values for a Lagrange grid.
!
!    Local, real ( kind = 8 ) ZE(NI), the exact function values at XI.
!
!    Local, real ( kind = 8 ) ZI(NI), the sparse interpolant values at XI.
!
!    Local, real ( kind = 8 ) ZPI(NI), one set of Lagrange interpolant values at XI.
!
  implicit none
  !type(state_type), intent(in), dimension(:,:,:) :: state
  real, dimension(:,:,:,:,:), allocatable :: snapmatrix_velocity ! allocate(snapmatrix_velocity(no_smolyak_nodes,u_nodes,snapshots,dim,numphase))
  real, dimension(:,:,:,:), intent(out),allocatable ::  interpolate_velocity 
  real ( kind = 8 ), allocatable :: a(:)
  real ( kind = 8 ) app_error
  real ( kind = 8 ), allocatable :: b(:)
  integer ( kind = 4 ), allocatable :: c(:)
  integer ( kind = 4 ) h
  integer ( kind = 4 ), allocatable :: ind(:)
  integer ( kind = 4 ) l
  integer ( kind = 4 ) l_max
  integer ( kind = 4 ) l_min
  real :: a1,a2,a3,a4,a5
  logical more
  integer ( kind = 4 ) nd
  integer ( kind = 4 ) nd_total
  integer ( kind = 4 ) ni
  real ( kind = 8 ) r8vec_norm_affine
  integer ( kind = 4 ) seed
 ! integer ( kind = 4 ) sparse_max
  integer ( kind = 4 ) sparse_min
  integer ( kind = 4 ) t
  integer ( kind = 4 ), allocatable :: w(:)
  real ( kind = 8 ), allocatable :: xd(:,:)
  real ( kind = 8 ), allocatable :: xi(:,:)
  real ( kind = 8 ), allocatable :: zd(:)
  real ( kind = 8 ), allocatable :: ze(:)
  real ( kind = 8 ), allocatable :: zi(:)
  real ( kind = 8 ), allocatable :: zpi(:)
  integer :: i, ii, total, total_timestep,num_pod,j,k,jj,kk, smolyak_total_points,u_nodesi,snapshotsi,dimi,numphasei
   
  allocate (ind(1:m))

  ! Define the region.
  allocate ( a(1:m) )
  allocate ( b(1:m) )
  
 
   
  !total_timestep=864
  !num_pod=12*36
  !smolyak_total_points=41!2**sparse_max+1
  !print *, '  smolyak_total_points=2**sparse_max+1=', smolyak_total_points 
  !  stop 2323
  !  Define the interpolation evaluation information.
  ni = 1
  seed = 123456789
  allocate ( xi(m,ni) )
  !call r8mat_uniform_abvec ( m, ni, a, b, seed, xi )
!  print * , 'xi', xi
  print *, 'm,ni', m, ni 
 ! stop 345
   if (.false.) then
  a(1:m) = 1
  b(1:m) = 2
  a1=1
  a2=1.1464578
  a3=1.5
  a4=1.8535421
  a5=2
  xi(1,1)= 2 
  xi(2,1)= 1.5
  xi(3,1)= 2
  xi(4,1)= 1.55 !1.7 !2.0000000000000000   1.5000000000000000   2.0000000000000000   1.5000000000000000,21
  endif

  if (.true.) then
  a(1:m) = 0.1
  b(1:m) = 0.5
  a1=0.1
  a2=0.158578
  a3=0.3
  a4=0.441421
  a5=0.5
  xi(1,1)= 0.1 
  xi(2,1)= 0.5
  xi(3,1)= 0.5
  xi(4,1)= 0.5! 0.10000000000000000       0.30000000000000000       0.50000000000000000       0.30000000000000000    
   call get_option(&
         '/reduced_model/variables_to_interpolate/dimension_one', xi(1,1))
   call get_option(&
         '/reduced_model/variables_to_interpolate/dimension_two', xi(2,1))
 call get_option(&
         '/reduced_model/variables_to_interpolate/dimension_three', xi(3,1))
 call get_option(&
         '/reduced_model/variables_to_interpolate/dimension_four', xi(4,1)) 
  endif
  write ( *, '(a,i6)' ) '  Number of interpolation points is NI = ', ni
  allocate ( ze(1:ni) )
 ! call f_sinr ( m, ni, xi, ze )
  print *, 'actual function value is, will interpolate this value later on,then do comparison', xi,ze
!
!  Compute a sequence of sparse grid interpolants of increasing level.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '   L     Order    ApproxError'
  write ( *, '(a)' ) ''
  allocate (interpolate_velocity(u_nodes,snapshots,dimen,numphase))
  allocate ( zpi(1:ni) )
  allocate ( zi(1:ni) )
  !smolynode, u_nodes,snapshots,dim,numphase
  Do u_nodesi=1, u_nodes  
   Do snapshotsi=1, snapshots  
   do dimi=1, size(snapmatrix_velocity,4)
   do numphasei=1, size(snapmatrix_velocity,5)
  total = 0
  sparse_min = 0

  ! do l_max = sparse_min, sparse_max
    do l_max = sparse_max, sparse_max
    allocate ( c(0:l_max) )
    allocate ( w(0:l_max) )
    call smolyak_coefficients ( l_max, m, c, w )

    zi(1:ni) = 0.0D+00
    nd_total = 0
    
    l_min = max ( l_max + 1 - m, 0 )

    do l = l_min, l_max

      more = .false.
      
      do
     ! print *,'l_max,l,' ,l_max,l   
!  Define the next product grid.
!
        call comp_next ( l, m, ind, more, h, t )
!
!  Count the grid, find its points, evaluate the data there.
!
        call lagrange_interp_nd_size2 ( m, ind, nd )
        allocate ( xd(1:m,1:nd) )
        call lagrange_interp_nd_grid2 ( m, ind, a, b, nd, xd )        
         !  print *,'xd', xd, 点1 坐标为： xd(:,1)
        allocate ( zd(1:nd) )

       ! print *,'---------------------nd',nd  
        do ii=1, nd           
              if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3   .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3) then 
!( 1.5000000000000000    1.5000000000000000   1.5000000000000000  1.5000000000000000)
           zd(ii)=snapmatrix_velocity(1,u_nodesi,snapshotsi,dimi,numphasei)
         else if(xd(1,ii) .eq. a1 .and. xd(2,ii) .eq. a3   .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3) then 
!(1.0000000000000000   1.5000000000000000   1.5000000000000000   1.5000000000000000)
           zd(ii)=snapmatrix_velocity(2,u_nodesi,snapshotsi,dimi,numphasei)
         else if(xd(1,ii) .eq. a5 .and. xd(2,ii) .eq. a3   .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3) then
 !( 2.0000000000000000   1.5000000000000000   1.5000000000000000   1.5000000000000000)
           zd(ii)=snapmatrix_velocity(3,u_nodesi,snapshotsi,dimi,numphasei)
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a1   .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3) then 
!(1.5000000000000000   1.0000000000000000   1.5000000000000000   1.5000000000000000 )
           zd(ii)=snapmatrix_velocity(4,u_nodesi,snapshotsi,dimi,numphasei)
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a5 .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3) then 
!(1.5000000000000000   2.0000000000000000   1.5000000000000000   1.5000000000000000 )
           zd(ii)=snapmatrix_velocity(5,u_nodesi,snapshotsi,dimi,numphasei)
        else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3 .and. xd(3,ii) .eq. a1 .and. xd(4,ii) .eq. a3) then 
!(1.5000000000000000   1.5000000000000000   1.0000000000000000   1.5000000000000000  )
           zd(ii)=snapmatrix_velocity(6,u_nodesi,snapshotsi,dimi,numphasei)
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3 .and. xd(3,ii) .eq. a5 .and. xd(4,ii) .eq. a3) then 
!(1.5000000000000000   1.5000000000000000   2.0000000000000000   a3000000000000000  )
           zd(ii)=snapmatrix_velocity(7,u_nodesi,snapshotsi,dimi,numphasei)
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3 .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a1) then
 !(a3000000000000000   a3000000000000000   a3000000000000000   1.0000000000000000   )
           zd(ii)=snapmatrix_velocity(8,u_nodesi,snapshotsi,dimi,numphasei)
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3   .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a5) then 
!(  a3000000000000000   a3000000000000000   a3000000000000000   2.0000000000000000  )
           zd(ii)=snapmatrix_velocity(9,u_nodesi,snapshotsi,dimi,numphasei) 
         else if(int(xd(1,ii)*1000000)/1000000 .eq. a2 .and. xd(2,ii) .eq. a3 .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3) then 
               !(a24466094067263   a3000000000000000   a3000000000000000   a3000000000000000 .and. xd(1,ii) <= 1.15 ..and. xd(2:4,ii) .eq. a3  )
           zd(ii)=snapmatrix_velocity(10,u_nodesi,snapshotsi,dimi,numphasei)
         elseif(int(xd(1,ii)*1000000)/1000000 .eq. a4 .and. xd(2,ii) .eq. a3 .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3) then 
               !(a45533905932737   a3000000000000000   a3000000000000000   a3000000000000000.and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3 )
           zd(ii)=snapmatrix_velocity(11,u_nodesi,snapshotsi,dimi,numphasei)
         else if(xd(1,ii) .eq. a1 .and. xd(2,ii) .eq. a1 .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3) then !(  1.0000000000000000   1.0000000000000000   a3000000000000000   a3000000000000000 )
           zd(ii)=snapmatrix_velocity(12,u_nodesi,snapshotsi,dimi,numphasei)
         else if(xd(1,ii) .eq. a5 .and. xd(2,ii) .eq. a1 .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3) then !( 2.0000000000000000   1.0000000000000000   a3000000000000000   a3000000000000000)
           zd(ii)=snapmatrix_velocity(13,u_nodesi,snapshotsi,dimi,numphasei)
          else if(xd(1,ii) .eq. a1 .and. xd(2,ii) .eq.a5  .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3) then !(1.0000000000000000   2.0000000000000000   a3000000000000000   a3000000000000000)
           zd(ii)=snapmatrix_velocity(14,u_nodesi,snapshotsi,dimi,numphasei)
         else if(xd(1,ii) .eq. a5 .and. xd(2,ii) .eq.a5  .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3) then !(2.0000000000000000   2.0000000000000000   a3000000000000000   a3000000000000000)
           zd(ii)=snapmatrix_velocity(15,u_nodesi,snapshotsi,dimi,numphasei)
         else if(xd(1,ii) .eq. a3 .and. int(xd(2,ii)*1000000)/1000000 .eq. a2 .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3) then 
            !(a3000000000000000   a24466094067263   a3000000000000000   a3000000000000000 xd(2,ii) >=1.14  .and. xd(2,ii) <=1.15)
           zd(ii)=snapmatrix_velocity(16,u_nodesi,snapshotsi,dimi,numphasei)
         else if(xd(1,ii) .eq. a3 .and. int(xd(2,ii)*1000000)/1000000 .eq. a4 .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3) then
          !(a3000000000000000   a45533905932737   a3000000000000000   a3000000000000000)
           zd(ii)=snapmatrix_velocity(17,u_nodesi,snapshotsi,dimi,numphasei)
           else if(xd(1,ii) .eq. a1 .and. xd(2,ii) .eq. a3 .and. xd(3,ii) .eq. a1 .and. xd(4,ii) .eq. a3) then !( 1.0000000000000000   a3000000000000000   1.0000000000000000  a3000000000000000)
           zd(ii)=snapmatrix_velocity(18,u_nodesi,snapshotsi,dimi,numphasei)
         else if(xd(1,ii) .eq.a5.and. xd(2,ii) .eq. a3   .and. xd(3,ii) .eq. a1 .and. xd(4,ii) .eq. a3) then !(2.0000000000000000   a3000000000000000   1.0000000000000000   a3000000000000000)
           zd(ii)=snapmatrix_velocity(19,u_nodesi,snapshotsi,dimi,numphasei)
         else if(xd(1,ii) .eq. a1 .and. xd(2,ii) .eq. a3   .and. xd(3,ii) .eq. a5 .and. xd(4,ii) .eq. a3) then !(1.0000000000000000   a3000000000000000   2.0000000000000000   a3000000000000000)
           zd(ii)=snapmatrix_velocity(20,u_nodesi,snapshotsi,dimi,numphasei)
         else if(xd(1,ii) .eq. a5 .and. xd(2,ii) .eq. a3   .and. xd(3,ii) .eq. a5 .and. xd(4,ii) .eq. a3) then !(2.0000000000000000   a3000000000000000   2.0000000000000000   a3000000000000000)
           zd(ii)=snapmatrix_velocity(21,u_nodesi,snapshotsi,dimi,numphasei) 
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a1   .and. xd(3,ii) .eq. a1 .and. xd(4,ii) .eq. a3) then !(a3000000000000000   1.0000000000000000   1.0000000000000000   a3000000000000000)
           zd(ii)=snapmatrix_velocity(22,u_nodesi,snapshotsi,dimi,numphasei)
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a5   .and. xd(3,ii) .eq. a1 .and. xd(4,ii) .eq. a3) then !(a3000000000000000   2.0000000000000000   1.0000000000000000   a3000000000000000)
           zd(ii)=snapmatrix_velocity(23,u_nodesi,snapshotsi,dimi,numphasei)
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a1   .and. xd(3,ii) .eq. a5 .and. xd(4,ii) .eq. a3) then !(a3000000000000000   1.0000000000000000   2.0000000000000000   a3000000000000000)
           zd(ii)=snapmatrix_velocity(24,u_nodesi,snapshotsi,dimi,numphasei)
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a5   .and. xd(3,ii) .eq. a5 .and. xd(4,ii) .eq. a3) then !(a3000000000000000   2.0000000000000000   2.0000000000000000   a3000000000000000)
           zd(ii)=snapmatrix_velocity(25,u_nodesi,snapshotsi,dimi,numphasei)
          else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3 .and. int(xd(3,ii)*1000000)/1000000 .eq. a2 .and. xd(4,ii) .eq. a3) then 
       !(a3000000000000000   a3000000000000000  a24466094067263   a3000000000000000)
           zd(ii)=snapmatrix_velocity(26,u_nodesi,snapshotsi,dimi,numphasei)
         else if(xd(1,ii) .eq. a3  .and. xd(2,ii) .eq. a3  .and. int(xd(3,ii)*1000000)/1000000 .eq. a4 .and. xd(4,ii) .eq. a3) then 
       !( a3000000000000000   a3000000000000000   a45533905932737   a3000000000000000 )
           zd(ii)=snapmatrix_velocity(27,u_nodesi,snapshotsi,dimi,numphasei)
         else if(xd(1,ii) .eq. a1 .and. xd(2,ii) .eq. a3   .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a1) then !(1.0000000000000000   a3000000000000000   a3000000000000000   1.0000000000000000)
           zd(ii)=snapmatrix_velocity(28,u_nodesi,snapshotsi,dimi,numphasei)
         else if(xd(1,ii) .eq. a5 .and. xd(2,ii) .eq. a3   .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a1) then !(2.0000000000000000   a3000000000000000   a3000000000000000   1.0000000000000000)
           zd(ii)=snapmatrix_velocity(29,u_nodesi,snapshotsi,dimi,numphasei)
           else if(xd(1,ii) .eq. a1 .and. xd(2,ii) .eq. a3   .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a5) then !(1.0000000000000000   a3000000000000000   a3000000000000000   2.0000000000000000)
           zd(ii)=snapmatrix_velocity(30,u_nodesi,snapshotsi,dimi,numphasei)
         else if(xd(1,ii) .eq. a5 .and. xd(2,ii) .eq. a3   .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a5) then !(2.0000000000000000   a3000000000000000   a3000000000000000   2.0000000000000000)
           zd(ii)=snapmatrix_velocity(31,u_nodesi,snapshotsi,dimi,numphasei)
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a1   .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a1) then !(a3000000000000000   1.0000000000000000   a3000000000000000   1.0000000000000000)
           zd(ii)=snapmatrix_velocity(32,u_nodesi,snapshotsi,dimi,numphasei)
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a5  .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a1) then !(a3000000000000000   2.0000000000000000   a3000000000000000   1.0000000000000000)
           zd(ii)=snapmatrix_velocity(33,u_nodesi,snapshotsi,dimi,numphasei) 
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a1   .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a5) then !(a3000000000000000   1.0000000000000000   a3000000000000000   2.0000000000000000)
           zd(ii)=snapmatrix_velocity(34,u_nodesi,snapshotsi,dimi,numphasei)
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a5   .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a5) then !(a3000000000000000   2.0000000000000000   a3000000000000000   2.0000000000000000)
           zd(ii)=snapmatrix_velocity(35,u_nodesi,snapshotsi,dimi,numphasei)
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3   .and. xd(3,ii) .eq. a1 .and. xd(4,ii) .eq. a1) then !(a3000000000000000   a3000000000000000   1.0000000000000000   1.0000000000000000)
           zd(ii)=snapmatrix_velocity(36,u_nodesi,snapshotsi,dimi,numphasei)
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3   .and. xd(3,ii) .eq. a5 .and. xd(4,ii) .eq. a1) then !( a3000000000000000   a3000000000000000   2.0000000000000000   1.0000000000000000)
           zd(ii)=snapmatrix_velocity(37,u_nodesi,snapshotsi,dimi,numphasei)
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3   .and. xd(3,ii) .eq. a1 .and. xd(4,ii) .eq. a5) then !(a3000000000000000   a3000000000000000   1.0000000000000000   2.0000000000000000)
           zd(ii)=snapmatrix_velocity(38,u_nodesi,snapshotsi,dimi,numphasei)
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3   .and. xd(3,ii) .eq. a5 .and. xd(4,ii) .eq. a5) then !(a3000000000000000   a3000000000000000   2.0000000000000000   2.0000000000000000)
           zd(ii)=snapmatrix_velocity(39,u_nodesi,snapshotsi,dimi,numphasei)
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3 .and. xd(3,ii) .eq. a3 .and. int(xd(4,ii)*1000000)/1000000 .eq. a2) then 
      !(a3000000000000000   a3000000000000000   a3000000000000000   a24466094067263 )
           zd(ii)=snapmatrix_velocity(40,u_nodesi,snapshotsi,dimi,numphasei)
         else if(xd(1,ii) .eq. a3.and. xd(2,ii) .eq. a3 .and. xd(3,ii) .eq. a3.and. int(xd(4,ii)*1000000)/1000000 .eq. a4) then 
     !(a3000000000000000   a3000000000000000   a3000000000000000   1.8535533905932737     )
           zd(ii)=snapmatrix_velocity(41,u_nodesi,snapshotsi,dimi,numphasei)
         endif
      !  print *, 'zdzdzd', zd(ii)
        enddo
          ! call f_sinr ( m, nd, xd, zd )  
          ! print *, 'zd',  zd(:) 
          !  Use the grid to evaluate the interpolant.
        call lagrange_interp_nd_value2 ( m, ind, a, b, nd, zd, ni, xi, zpi ) 
          !  Weighted the interpolant values and add to the sparse grid interpolant. 
        nd_total = nd_total + nd
        zi(1:ni) = zi(1:ni) + c(l) * zpi(1:ni)
        !  print *, ' nd_total,level',  nd_total, l, c(l), nd
        !  print *, 'zizizizizizizic(l)zpi', zi(1)!,c(l),zpi(1)
        deallocate ( xd )
        deallocate ( zd )     
       
       if ( .not. more ) then
          exit
        end if    
 
      end do !DO

    end do !do l = l_min, l_max   
    deallocate ( c )
    deallocate ( w )
    
   end do ! do l_max = sparse_max, sparse_max
   interpolate_velocity(u_nodesi,snapshotsi,dimi,numphasei)=zi(ni)
   enddo !Do u_nodesi  
  enddo !Do snapshotsi
  enddo !dimi
  enddo !numphasei

 !  print *, ' total', total , nd_total
 ! open(unit=66,file='podbasis_interp') 
 !  do j=1, total_timestep 
 !     write(66,*)  pod_coef_all_obv_41i(j,1:num_pod)
 ! enddo
 ! close(66)
 
  deallocate ( a )
  deallocate ( b )
  deallocate ( ind )
  deallocate ( xi )
  deallocate ( ze )
  deallocate ( zi )
  deallocate ( zpi )
  end subroutine interpolate_snapmatrix_velocity_4dim_level1

subroutine interpolate_snapmatrix_velocity_8dim_level1( snapmatrix_velocity, interpolate_velocity)  
    
! top-bottom between 1-2, pod_coef_all_obv is used for storing some podbasis to avoid to much modifications.
 
!  Parameters:
!
!    Local, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) SPARSE_MAX, the maximum sparse grid level to try.
!
!  Local Parameters:
!
!    Local, real ( kind = 8 ) A(M), B(M), the upper and lower variable limits 
!    in each dimension.
!
!    Local, real ( kind = 8 ) APP_ERROR, the averaged Euclidean norm of the 
!    difference between the sparse interpolant and the exact function at 
!    the interpolation points.
!
!    Local, integer ( kind = 4 ) C(0:L_MAX), the sparse grid coefficient vector.
!    Results at level L are weighted by C(L).
!
!    Local, integer ( kind = 4 ) IND(M), the 1D indices defining a Lagrange grid.
!    Each entry is a 1d "level" that specifies the order of a 
!    Clenshaw Curtis 1D grid.
!
!    Local, integer ( kind = 4 ) L, the current Lagrange grid level.
!
!    Local, integer ( kind = 4 ) L_MAX, the current sparse grid level.
!
!    Local, integer ( kind = 4 ) MORE, used to control the enumeration of all the
!    Lagrange grids at a current grid level.
!
!    Local, integer ( kind = 4 ) ND, the number of points used in a Lagrange grid.
!
!    Local, integer ( kind = 4 ) ND_TOTAL, the total number of points used in all the
!    Lagrange interpolants for a given sparse interpolant points that occur
!    in more than one Lagrange grid are counted multiple times.
!
!    Local, integer ( kind = 4 ) NI, the number of interpolant evaluation points.
!
!    Local, integer ( kind = 4 ) SPARSE_MIN, the minimum sparse grid level to try.
!
!    Local, real ( kind = 8 ) XD(M,ND), the data points for a Lagrange grid.
!
!    Local, real ( kind = 8 ) XI(M,NI), the interpolant evaluation points.
!
!    Local, real ( kind = 8 ) ZD(ND), the data values for a Lagrange grid.
!
!    Local, real ( kind = 8 ) ZE(NI), the exact function values at XI.
!
!    Local, real ( kind = 8 ) ZI(NI), the sparse interpolant values at XI.
!
!    Local, real ( kind = 8 ) ZPI(NI), one set of Lagrange interpolant values at XI.
!
  implicit none
  !type(state_type), intent(in), dimension(:,:,:) :: state
  real, dimension(:,:,:,:,:), allocatable :: snapmatrix_velocity ! allocate(snapmatrix_velocity(no_smolyak_nodes,u_nodes,snapshots,dim,numphase))
  real, dimension(:,:,:,:), intent(out),allocatable ::  interpolate_velocity 
  real ( kind = 8 ), allocatable :: a(:)
  real ( kind = 8 ) app_error
  real ( kind = 8 ), allocatable :: b(:)
  integer ( kind = 4 ), allocatable :: c(:)
  integer ( kind = 4 ) h
  integer ( kind = 4 ), allocatable :: ind(:)
  integer ( kind = 4 ) l
  integer ( kind = 4 ) l_max
  integer ( kind = 4 ) l_min
  real :: a1,a2,a3,a4,a5
  logical more
  integer ( kind = 4 ) nd
  integer ( kind = 4 ) nd_total
  integer ( kind = 4 ) ni
  real ( kind = 8 ) r8vec_norm_affine
  integer ( kind = 4 ) seed
 ! integer ( kind = 4 ) sparse_max
  integer ( kind = 4 ) sparse_min
  integer ( kind = 4 ) t
  integer ( kind = 4 ), allocatable :: w(:)
  real ( kind = 8 ), allocatable :: xd(:,:)
  real ( kind = 8 ), allocatable :: xi(:,:)
  real ( kind = 8 ), allocatable :: zd(:)
  real ( kind = 8 ), allocatable :: ze(:)
  real ( kind = 8 ), allocatable :: zi(:)
  real ( kind = 8 ), allocatable :: zpi(:)
  integer :: i, ii, total, total_timestep,num_pod,j,k,jj,kk, smolyak_total_points,u_nodesi,snapshotsi,dimi,numphasei
   
  allocate (ind(1:m))

  ! Define the region.
  allocate ( a(1:m) )
  allocate ( b(1:m) )
  
 
   
  !total_timestep=864
  !num_pod=12*36
  !smolyak_total_points=41!2**sparse_max+1
  !print *, '  smolyak_total_points=2**sparse_max+1=', smolyak_total_points 
  !  stop 2323
  !  Define the interpolation evaluation information.
  ni = 1
  seed = 123456789
  allocate ( xi(m,ni) )
  !call r8mat_uniform_abvec ( m, ni, a, b, seed, xi )
!  print * , 'xi', xi
  print *, 'm,ni', m, ni 
 ! stop 345
   if (.false.) then
  a(1:m) = 1
  b(1:m) = 2
  a1=1
  a2=1.1464578
  a3=1.5
  a4=1.8535421
  a5=2
  xi(1,1)= 2 
  xi(2,1)= 1.5
  xi(3,1)= 2
  xi(4,1)= 1.55 !1.7 !2.0000000000000000   1.5000000000000000   2.0000000000000000   1.5000000000000000,21
  endif

  if (.true.) then
  a(1:m) = 0.1
  b(1:m) = 0.5
  a1=0.1
  a2=0.158578
  a3=0.3
  a4=0.441421
  a5=0.5
  xi(1,1)= 0.1 
   xi(2,1)= 0.5
   xi(3,1)= 0.5
  xi(4,1)= 0.5! 0.10000000000000000       0.30000000000000000       0.50000000000000000       0.30000000000000000    
 
   call get_option(&
         '/reduced_model/variables_to_interpolate/dimension_one', xi(1,1))
   call get_option(&
         '/reduced_model/variables_to_interpolate/dimension_two', xi(2,1))
   call get_option(&
         '/reduced_model/variables_to_interpolate/dimension_three', xi(3,1))
   call get_option(&
         '/reduced_model/variables_to_interpolate/dimension_four', xi(4,1)) 
  endif
   
   xi(5,1)= 0.1 
   xi(6,1)= 0.5
   xi(7,1)= 0.5
   xi(8,1)= 0.5  ! 0.10000000000000000       0.30000000000000000       0.50000000000000000       0.30000000000000000   
 
  write ( *, '(a,i6)' ) '  Number of interpolation points is NI = ', ni
  allocate ( ze(1:ni) )
 ! call f_sinr ( m, ni, xi, ze )
  print *, 'actual function value is, will interpolate this value later on,then do comparison', xi,ze
!
!  Compute a sequence of sparse grid interpolants of increasing level.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '   L     Order    ApproxError'
  write ( *, '(a)' ) ''
  allocate (interpolate_velocity(u_nodes,snapshots,dimen,numphase))
  allocate ( zpi(1:ni) )
  allocate ( zi(1:ni) )
  !smolynode, u_nodes,snapshots,dim,numphase
  Do u_nodesi=1, u_nodes  
   Do snapshotsi=1, snapshots  
   do dimi=1, size(snapmatrix_velocity,4)
   do numphasei=1, size(snapmatrix_velocity,5)
  total = 0
  sparse_min = 0

  ! do l_max = sparse_min, sparse_max
    do l_max = sparse_max, sparse_max
    allocate ( c(0:l_max) )
    allocate ( w(0:l_max) )
    call smolyak_coefficients ( l_max, m, c, w )

    zi(1:ni) = 0.0D+00
    nd_total = 0
    
    l_min = max ( l_max + 1 - m, 0 )

    do l = l_min, l_max

      more = .false.
      
      do
     ! print *,'l_max,l,' ,l_max,l   
!  Define the next product grid.
!
        call comp_next ( l, m, ind, more, h, t )
!
!  Count the grid, find its points, evaluate the data there.
!
        call lagrange_interp_nd_size2 ( m, ind, nd )
        allocate ( xd(1:m,1:nd) )
        call lagrange_interp_nd_grid2 ( m, ind, a, b, nd, xd )        
         !  print *,'xd', xd, 点1 坐标为： xd(:,1)
        allocate ( zd(1:nd) )

       ! print *,'---------------------nd',nd  

!   0.30000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000     
 !  0.10000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000     
 !   0.50000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000     
!    0.30000000000000000  0.10000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000      
 !   0.30000000000000000  0.50000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000     
 !    0.30000000000000000  0.30000000000000000  0.10000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000      
!   0.30000000000000000  0.30000000000000000  0.50000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000     
 !  0.30000000000000000  0.30000000000000000  0.30000000000000000  0.10000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000      
 !  0.30000000000000000  0.30000000000000000  0.30000000000000000  0.50000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000     
  ! 0.30000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000  0.10000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000     
  ! 0.30000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000  0.50000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000     
  ! 0.30000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000  0.10000000000000000  0.30000000000000000  0.30000000000000000      
  ! 0.30000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000  0.50000000000000000  0.30000000000000000  0.30000000000000000     
  ! 0.30000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000  0.10000000000000000  0.30000000000000000     
  ! 0.30000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000  0.50000000000000000  0.30000000000000000     
  ! 0.30000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000  0.10000000000000000     
 !  0.30000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000  0.50000000000000000  
        do ii=1, nd           
               if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3 .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3 .and. xd(5,ii) .eq. a3 .and. xd(6,ii).eq. a3 .and. xd(7,ii) .eq. a3 .and. xd(8,ii) .eq. a3 ) then  
           zd(ii)=snapmatrix_velocity(1,u_nodesi,snapshotsi,dimi, numphasei)
   else if(xd(1,ii) .eq. a1 .and. xd(2,ii) .eq. a3 .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3 .and. xd(5,ii) .eq. a3 .and. xd(6,ii).eq. a3 .and. xd(7,ii) .eq. a3 .and. xd(8,ii) .eq. a3 ) then  
           zd(ii)=snapmatrix_velocity(2,u_nodesi,snapshotsi,dimi, numphasei)
   else if(xd(1,ii) .eq. a5 .and. xd(2,ii) .eq. a3 .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3 .and. xd(5,ii) .eq. a3 .and. xd(6,ii).eq. a3 .and. xd(7,ii) .eq. a3 .and. xd(8,ii) .eq. a3 ) then 
           zd(ii)=snapmatrix_velocity(3,u_nodesi,snapshotsi,dimi, numphasei)
   else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a1 .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3 .and. xd(5,ii) .eq. a3 .and. xd(6,ii).eq. a3 .and. xd(7,ii) .eq. a3 .and. xd(8,ii) .eq. a3 ) then 
           zd(ii)=snapmatrix_velocity(4,u_nodesi,snapshotsi,dimi, numphasei)
   else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a5 .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3 .and. xd(5,ii) .eq. a3 .and. xd(6,ii).eq. a3 .and. xd(7,ii) .eq. a3 .and. xd(8,ii) .eq. a3 ) then 
           zd(ii)=snapmatrix_velocity(5,u_nodesi,snapshotsi,dimi, numphasei)
   else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3 .and. xd(3,ii) .eq. a1 .and. xd(4,ii) .eq. a3 .and. xd(5,ii) .eq. a3 .and. xd(6,ii).eq. a3 .and. xd(7,ii) .eq. a3 .and. xd(8,ii) .eq. a3 ) then 
           zd(ii)=snapmatrix_velocity(6,u_nodesi,snapshotsi,dimi, numphasei)
   else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3 .and. xd(3,ii) .eq. a5 .and. xd(4,ii) .eq. a3 .and. xd(5,ii) .eq. a3 .and. xd(6,ii).eq. a3 .and. xd(7,ii) .eq. a3 .and. xd(8,ii) .eq. a3) then 
           zd(ii)=snapmatrix_velocity(7,u_nodesi,snapshotsi,dimi, numphasei)
   else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3 .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a1 .and. xd(5,ii) .eq. a3 .and. xd(6,ii).eq. a3 .and. xd(7,ii) .eq. a3 .and. xd(8,ii) .eq. a3) then 
           zd(ii)=snapmatrix_velocity(8,u_nodesi,snapshotsi,dimi, numphasei)
   else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3 .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a5 .and. xd(5,ii) .eq. a3 .and. xd(6,ii).eq. a3 .and. xd(7,ii) .eq. a3 .and. xd(8,ii) .eq. a3) then
           zd(ii)=snapmatrix_velocity(9,u_nodesi,snapshotsi,dimi, numphasei)
   else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3 .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3 .and. xd(5,ii) .eq. a1 .and. xd(6,ii).eq. a3 .and. xd(7,ii) .eq. a3 .and. xd(8,ii) .eq. a3) then              
           zd(ii)=snapmatrix_velocity(10,u_nodesi,snapshotsi,dimi, numphasei)
  else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3 .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3 .and. xd(5,ii) .eq. a5 .and. xd(6,ii).eq. a3 .and. xd(7,ii) .eq. a3 .and. xd(8,ii) .eq. a3) then           
           zd(ii)=snapmatrix_velocity(11,u_nodesi,snapshotsi,dimi, numphasei)
 else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3 .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3 .and. xd(5,ii) .eq. a3 .and. xd(6,ii).eq. a1 .and. xd(7,ii) .eq. a3 .and. xd(8,ii) .eq. a3) then
           zd(ii)=snapmatrix_velocity(12,u_nodesi,snapshotsi,dimi, numphasei)
 else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3 .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3 .and. xd(5,ii) .eq. a3 .and. xd(6,ii).eq. a5 .and. xd(7,ii) .eq. a3 .and. xd(8,ii) .eq. a3) then
           zd(ii)=snapmatrix_velocity(13,u_nodesi,snapshotsi,dimi, numphasei)
 else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3 .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3 .and. xd(5,ii) .eq. a3 .and. xd(6,ii).eq. a3 .and. xd(7,ii) .eq. a1 .and. xd(8,ii) .eq. a3) then
           zd(ii)=snapmatrix_velocity(14,u_nodesi,snapshotsi,dimi, numphasei)
 else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3 .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3 .and. xd(5,ii) .eq. a3 .and. xd(6,ii).eq. a3 .and. xd(7,ii) .eq. a5 .and. xd(8,ii) .eq. a3) then
           zd(ii)=snapmatrix_velocity(15,u_nodesi,snapshotsi,dimi, numphasei)
 else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3 .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3 .and. xd(5,ii) .eq. a3 .and. xd(6,ii).eq. a3 .and. xd(7,ii) .eq. a3 .and. xd(8,ii) .eq. a1) then
            zd(ii)=snapmatrix_velocity(16,u_nodesi,snapshotsi,dimi, numphasei)
 else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3 .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3 .and. xd(5,ii) .eq. a3 .and. xd(6,ii).eq. a3 .and. xd(7,ii) .eq. a3 .and. xd(8,ii) .eq. a5) then
           zd(ii)=snapmatrix_velocity(17,u_nodesi,snapshotsi,dimi, numphasei)
         endif
      !  print *, 'zdzdzd', zd(ii)
        enddo
          ! call f_sinr ( m, nd, xd, zd )  
          ! print *, 'zd',  zd(:) 
          !  Use the grid to evaluate the interpolant.
        call lagrange_interp_nd_value2 ( m, ind, a, b, nd, zd, ni, xi, zpi ) 
          !  Weighted the interpolant values and add to the sparse grid interpolant. 
        nd_total = nd_total + nd
        zi(1:ni) = zi(1:ni) + c(l) * zpi(1:ni)
        !  print *, ' nd_total,level',  nd_total, l, c(l), nd
        !  print *, 'zizizizizizizic(l)zpi', zi(1)!,c(l),zpi(1)
        deallocate ( xd )
        deallocate ( zd )     
       
       if ( .not. more ) then
          exit
        end if    
 
      end do !DO

    end do !do l = l_min, l_max   
    deallocate ( c )
    deallocate ( w )
    
   end do ! do l_max = sparse_max, sparse_max
   interpolate_velocity(u_nodesi,snapshotsi,dimi,numphasei)=zi(ni)
   enddo !Do u_nodesi  
  enddo !Do snapshotsi
  enddo !dimi
  enddo !numphasei

  
 
  deallocate ( a )
  deallocate ( b )
  deallocate ( ind )
  deallocate ( xi )
  deallocate ( ze )
  deallocate ( zi )
  deallocate ( zpi )
  end subroutine interpolate_snapmatrix_velocity_8dim_level1

subroutine interpolate_snapmatrix_pressure_8dim_level1( snapmatrix_pressure, interpolate_pressure)  
    
! top-bottom between 1-2, pod_coef_all_obv is used for storing some podbasis to avoid to much modifications.
 
!  Parameters:
!
!    Local, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) SPARSE_MAX, the maximum sparse grid level to try.
!
!  Local Parameters:
!
!    Local, real ( kind = 8 ) A(M), B(M), the upper and lower variable limits 
!    in each dimension.
!
!    Local, real ( kind = 8 ) APP_ERROR, the averaged Euclidean norm of the 
!    difference between the sparse interpolant and the exact function at 
!    the interpolation points.
!
!    Local, integer ( kind = 4 ) C(0:L_MAX), the sparse grid coefficient vector.
!    Results at level L are weighted by C(L).
!
!    Local, integer ( kind = 4 ) IND(M), the 1D indices defining a Lagrange grid.
!    Each entry is a 1d "level" that specifies the order of a 
!    Clenshaw Curtis 1D grid.
!
!    Local, integer ( kind = 4 ) L, the current Lagrange grid level.
!
!    Local, integer ( kind = 4 ) L_MAX, the current sparse grid level.
!
!    Local, integer ( kind = 4 ) MORE, used to control the enumeration of all the
!    Lagrange grids at a current grid level.
!
!    Local, integer ( kind = 4 ) ND, the number of points used in a Lagrange grid.
!
!    Local, integer ( kind = 4 ) ND_TOTAL, the total number of points used in all the
!    Lagrange interpolants for a given sparse interpolant points that occur
!    in more than one Lagrange grid are counted multiple times.
!
!    Local, integer ( kind = 4 ) NI, the number of interpolant evaluation points.
!
!    Local, integer ( kind = 4 ) SPARSE_MIN, the minimum sparse grid level to try.
!
!    Local, real ( kind = 8 ) XD(M,ND), the data points for a Lagrange grid.
!
!    Local, real ( kind = 8 ) XI(M,NI), the interpolant evaluation points.
!
!    Local, real ( kind = 8 ) ZD(ND), the data values for a Lagrange grid.
!
!    Local, real ( kind = 8 ) ZE(NI), the exact function values at XI.
!
!    Local, real ( kind = 8 ) ZI(NI), the sparse interpolant values at XI.
!
!    Local, real ( kind = 8 ) ZPI(NI), one set of Lagrange interpolant values at XI.
!
  implicit none
  !type(state_type), intent(in), dimension(:,:,:) :: state
  real, dimension(:,:,:,:), allocatable :: snapmatrix_pressure ! allocate(snapmatrix_pressure(no_smolyak_nodes,p_nodes,snapshots,numphase))
  real, dimension(:,:,:), intent(out), allocatable ::  interpolate_pressure 
  real ( kind = 8 ), allocatable :: a(:)
  real ( kind = 8 ) app_error
  real ( kind = 8 ), allocatable :: b(:)
  integer ( kind = 4 ), allocatable :: c(:)
  integer ( kind = 4 ) h
  integer ( kind = 4 ), allocatable :: ind(:)
  integer ( kind = 4 ) l
  integer ( kind = 4 ) l_max
  integer ( kind = 4 ) l_min
  !integer ( kind = 4 ) m
  logical more
  integer ( kind = 4 ) nd
  integer ( kind = 4 ) nd_total
  integer ( kind = 4 ) ni
  real ( kind = 8 ) r8vec_norm_affine
  integer ( kind = 4 ) seed
 ! integer ( kind = 4 ) sparse_max
  integer ( kind = 4 ) sparse_min
  integer ( kind = 4 ) t
  integer ( kind = 4 ), allocatable :: w(:)
  real ( kind = 8 ), allocatable :: xd(:,:)
  real ( kind = 8 ), allocatable :: xi(:,:)
  real ( kind = 8 ), allocatable :: zd(:)
  real ( kind = 8 ), allocatable :: ze(:)
  real ( kind = 8 ), allocatable :: zi(:)
  real ( kind = 8 ), allocatable :: zpi(:)
  integer :: i, ii, total, total_timestep,num_pod,j,k,jj,kk, smolyak_total_points,u_nodesi,snapshotsi,dimi,numphasei
   real :: a1,a2,a3,a4,a5 
  allocate (ind(1:m))

  ! Define the region.
  allocate ( a(1:m) )
  allocate ( b(1:m) )

  
   
  !total_timestep=864
  !num_pod=12*36
  !smolyak_total_points=41!2**sparse_max+1
  !print *, '  smolyak_total_points=2**sparse_max+1=', smolyak_total_points 
  !  stop 2323
  !  Define the interpolation evaluation information.
  ni = 1
  seed = 123456789
  allocate ( xi(m,ni) )
  !call r8mat_uniform_abvec ( m, ni, a, b, seed, xi )
!  print * , 'xi', xi
 
   if (.false.) then
  a(1:m) = 1
  b(1:m) = 2
  a1=1
  a2=1.1464!5780000000000
  a3=1.5
  a4=1.8535!421000
  a5=2
  xi(1,1)= 2 
  xi(2,1)= 1.5
  xi(3,1)= 2
  xi(4,1)= 1.55 !1.7 !2.0000000000000000   1.5000000000000000   2.0000000000000000   1.5000000000000000,21
  endif

  if (.true.) then
  a(1:m) = 0.1
  b(1:m) = 0.5
  a1=0.1
  a2=0.158578
  a3=0.3
  a4=0.441421
  a5=0.5
  xi(1,1)= 0.1 
   xi(2,1)= 0.5
   xi(3,1)= 0.5
  xi(4,1)= 0.5! 0.10000000000000000       0.30000000000000000       0.50000000000000000       0.30000000000000000    
 
   call get_option(&
         '/reduced_model/variables_to_interpolate/dimension_one', xi(1,1))
   call get_option(&
         '/reduced_model/variables_to_interpolate/dimension_two', xi(2,1))
   call get_option(&
         '/reduced_model/variables_to_interpolate/dimension_three', xi(3,1))
   call get_option(&
         '/reduced_model/variables_to_interpolate/dimension_four', xi(4,1)) 
  endif
   
   xi(5,1)= 0.1 
   xi(6,1)= 0.5
   xi(7,1)= 0.5
   xi(8,1)= 0.5  ! 0.10000000000000000       0.30000000000000000       0.50000000000000000       0.30000000000000000   
 
  write ( *, '(a,i6)' ) '  Number of interpolation points is NI = ', ni
  allocate ( ze(1:ni) )
 ! call f_sinr ( m, ni, xi, ze )
  print *, 'actual function value is, will interpolate this value later on,then do comparison', xi,ze
!
!  Compute a sequence of sparse grid interpolants of increasing level.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '   L     Order    ApproxError'
  write ( *, '(a)' ) ''
  allocate(interpolate_pressure(p_nodes,snapshots,numphase)) 
  allocate ( zpi(1:ni) )
  allocate ( zi(1:ni) )
  !smolynode, u_nodes,snapshots,dim,numphase
  Do u_nodesi=1, p_nodes  
   Do snapshotsi=1, snapshots  
   do numphasei=1, numphase
  total = 0
  sparse_min = 0

  ! do l_max = sparse_min, sparse_max
    do l_max = sparse_max, sparse_max
    allocate ( c(0:l_max) )
    allocate ( w(0:l_max) )
    call smolyak_coefficients ( l_max, m, c, w )

    zi(1:ni) = 0.0D+00
    nd_total = 0
    
    l_min = max ( l_max + 1 - m, 0 )

    do l = l_min, l_max

      more = .false.
      
      do
     ! print *,'l_max,l,' ,l_max,l   
!  Define the next product grid.
!
        call comp_next ( l, m, ind, more, h, t )
!
!  Count the grid, find its points, evaluate the data there.
!
        call lagrange_interp_nd_size2 ( m, ind, nd )
        allocate ( xd(1:m,1:nd) )
        call lagrange_interp_nd_grid2 ( m, ind, a, b, nd, xd )        
         !  print *,'xd', xd, 点1 坐标为： xd(:,1)
        allocate ( zd(1:nd) )

       ! print *,'---------------------nd',nd  


!   0.30000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000     
 !  0.10000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000     
 !   0.50000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000     
!    0.30000000000000000  0.10000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000      
 !   0.30000000000000000  0.50000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000     
 !    0.30000000000000000  0.30000000000000000  0.10000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000      
!   0.30000000000000000  0.30000000000000000  0.50000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000     
 !  0.30000000000000000  0.30000000000000000  0.30000000000000000  0.10000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000      
 !  0.30000000000000000  0.30000000000000000  0.30000000000000000  0.50000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000     
  ! 0.30000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000  0.10000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000     
  ! 0.30000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000  0.50000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000     
  ! 0.30000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000  0.10000000000000000  0.30000000000000000  0.30000000000000000      
  ! 0.30000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000  0.50000000000000000  0.30000000000000000  0.30000000000000000     
  ! 0.30000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000  0.10000000000000000  0.30000000000000000     
  ! 0.30000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000  0.50000000000000000  0.30000000000000000     
  ! 0.30000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000  0.10000000000000000     
 !  0.30000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000  0.30000000000000000  0.50000000000000000  
  do ii=1, nd           
        if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3 .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3 .and. xd(5,ii) .eq. a3 .and. xd(6,ii).eq. a3 .and. xd(7,ii) .eq. a3 .and. xd(8,ii) .eq. a3 ) then  
           zd(ii)=snapmatrix_pressure(1,u_nodesi,snapshotsi, numphasei)
   else if(xd(1,ii) .eq. a1 .and. xd(2,ii) .eq. a3 .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3 .and. xd(5,ii) .eq. a3 .and. xd(6,ii).eq. a3 .and. xd(7,ii) .eq. a3 .and. xd(8,ii) .eq. a3 ) then  
           zd(ii)=snapmatrix_pressure(2,u_nodesi,snapshotsi, numphasei)
   else if(xd(1,ii) .eq. a5 .and. xd(2,ii) .eq. a3 .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3 .and. xd(5,ii) .eq. a3 .and. xd(6,ii).eq. a3 .and. xd(7,ii) .eq. a3 .and. xd(8,ii) .eq. a3 ) then 
           zd(ii)=snapmatrix_pressure(3,u_nodesi,snapshotsi, numphasei)
   else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a1 .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3 .and. xd(5,ii) .eq. a3 .and. xd(6,ii).eq. a3 .and. xd(7,ii) .eq. a3 .and. xd(8,ii) .eq. a3 ) then 
           zd(ii)=snapmatrix_pressure(4,u_nodesi,snapshotsi, numphasei)
   else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a5 .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3 .and. xd(5,ii) .eq. a3 .and. xd(6,ii).eq. a3 .and. xd(7,ii) .eq. a3 .and. xd(8,ii) .eq. a3 ) then 
           zd(ii)=snapmatrix_pressure(5,u_nodesi,snapshotsi, numphasei)
   else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3 .and. xd(3,ii) .eq. a1 .and. xd(4,ii) .eq. a3 .and. xd(5,ii) .eq. a3 .and. xd(6,ii).eq. a3 .and. xd(7,ii) .eq. a3 .and. xd(8,ii) .eq. a3 ) then 
           zd(ii)=snapmatrix_pressure(6,u_nodesi,snapshotsi, numphasei)
   else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3 .and. xd(3,ii) .eq. a5 .and. xd(4,ii) .eq. a3 .and. xd(5,ii) .eq. a3 .and. xd(6,ii).eq. a3 .and. xd(7,ii) .eq. a3 .and. xd(8,ii) .eq. a3) then 
           zd(ii)=snapmatrix_pressure(7,u_nodesi,snapshotsi, numphasei)
   else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3 .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a1 .and. xd(5,ii) .eq. a3 .and. xd(6,ii).eq. a3 .and. xd(7,ii) .eq. a3 .and. xd(8,ii) .eq. a3) then 
           zd(ii)=snapmatrix_pressure(8,u_nodesi,snapshotsi, numphasei)
   else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3 .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a5 .and. xd(5,ii) .eq. a3 .and. xd(6,ii).eq. a3 .and. xd(7,ii) .eq. a3 .and. xd(8,ii) .eq. a3) then
           zd(ii)=snapmatrix_pressure(9,u_nodesi,snapshotsi, numphasei)
   else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3 .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3 .and. xd(5,ii) .eq. a1 .and. xd(6,ii).eq. a3 .and. xd(7,ii) .eq. a3 .and. xd(8,ii) .eq. a3) then              
           zd(ii)=snapmatrix_pressure(10,u_nodesi,snapshotsi, numphasei)
  else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3 .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3 .and. xd(5,ii) .eq. a5 .and. xd(6,ii).eq. a3 .and. xd(7,ii) .eq. a3 .and. xd(8,ii) .eq. a3) then           
           zd(ii)=snapmatrix_pressure(11,u_nodesi,snapshotsi, numphasei)
 else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3 .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3 .and. xd(5,ii) .eq. a3 .and. xd(6,ii).eq. a1 .and. xd(7,ii) .eq. a3 .and. xd(8,ii) .eq. a3) then
           zd(ii)=snapmatrix_pressure(12,u_nodesi,snapshotsi, numphasei)
 else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3 .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3 .and. xd(5,ii) .eq. a3 .and. xd(6,ii).eq. a5 .and. xd(7,ii) .eq. a3 .and. xd(8,ii) .eq. a3) then
           zd(ii)=snapmatrix_pressure(13,u_nodesi,snapshotsi, numphasei)
 else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3 .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3 .and. xd(5,ii) .eq. a3 .and. xd(6,ii).eq. a3 .and. xd(7,ii) .eq. a1 .and. xd(8,ii) .eq. a3) then
           zd(ii)=snapmatrix_pressure(14,u_nodesi,snapshotsi, numphasei)
 else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3 .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3 .and. xd(5,ii) .eq. a3 .and. xd(6,ii).eq. a3 .and. xd(7,ii) .eq. a5 .and. xd(8,ii) .eq. a3) then
           zd(ii)=snapmatrix_pressure(15,u_nodesi,snapshotsi, numphasei)
 else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3 .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3 .and. xd(5,ii) .eq. a3 .and. xd(6,ii).eq. a3 .and. xd(7,ii) .eq. a3 .and. xd(8,ii) .eq. a1) then
            zd(ii)=snapmatrix_pressure(16,u_nodesi,snapshotsi, numphasei)
 else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3 .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3 .and. xd(5,ii) .eq. a3 .and. xd(6,ii).eq. a3 .and. xd(7,ii) .eq. a3 .and. xd(8,ii) .eq. a5) then
           zd(ii)=snapmatrix_pressure(17,u_nodesi,snapshotsi, numphasei)
         endif
        enddo
          !call f_sinr ( m, nd, xd, zd )       
           
          ! print *, 'zd',  zd(:)
            
!  Use the grid to evaluate the interpolant.
        call lagrange_interp_nd_value2 ( m, ind, a, b, nd, zd, ni, xi, zpi )
 
!  Weighted the interpolant values and add to the sparse grid interpolant.
 
        nd_total = nd_total + nd
        zi(1:ni) = zi(1:ni) + c(l) * zpi(1:ni)
        !  print *, ' nd_total,level',  nd_total, l, c(l), nd
        !  print *, 'zizizizizizizic(l)zpi', zi(1)!,c(l),zpi(1)
        deallocate ( xd )
        deallocate ( zd )     
       
       if ( .not. more ) then
          exit
        end if    
 
      end do !DO

    end do !do l = l_min, l_max   
    deallocate ( c )
    deallocate ( w )
    
   end do ! do l_max = sparse_max, sparse_max
   interpolate_pressure(u_nodesi,snapshotsi,numphasei)=zi(ni)
   enddo !Do u_nodesi  
  enddo !Do snapshotsi 
  enddo !numphasei

 !  print *, ' total', total , nd_total
 ! open(unit=66,file='podbasis_interp') 
 !  do j=1, total_timestep 
 !     write(66,*)  pod_coef_all_obv_41i(j,1:num_pod)
 ! enddo
 ! close(66)
 
  deallocate ( a )
  deallocate ( b )
  deallocate ( ind )
  deallocate ( xi )
  deallocate ( ze )
  deallocate ( zi )
  deallocate ( zpi )
  end subroutine interpolate_snapmatrix_pressure_8dim_level1

subroutine interpolate_snapmatrix_volumefraction_8dim_level1( snapmatrix_volumefraction, interpolate_volumefraction)  
    
! top-bottom between 1-2, pod_coef_all_obv is used for storing some podbasis to avoid to much modifications.
 
!  Parameters:
!
!    Local, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) SPARSE_MAX, the maximum sparse grid level to try.
!
!  Local Parameters:
!
!    Local, real ( kind = 8 ) A(M), B(M), the upper and lower variable limits 
!    in each dimension.
!
!    Local, real ( kind = 8 ) APP_ERROR, the averaged Euclidean norm of the 
!    difference between the sparse interpolant and the exact function at 
!    the interpolation points.
!
!    Local, integer ( kind = 4 ) C(0:L_MAX), the sparse grid coefficient vector.
!    Results at level L are weighted by C(L).
!
!    Local, integer ( kind = 4 ) IND(M), the 1D indices defining a Lagrange grid.
!    Each entry is a 1d "level" that specifies the order of a 
!    Clenshaw Curtis 1D grid.
!
!    Local, integer ( kind = 4 ) L, the current Lagrange grid level.
!
!    Local, integer ( kind = 4 ) L_MAX, the current sparse grid level.
!
!    Local, integer ( kind = 4 ) MORE, used to control the enumeration of all the
!    Lagrange grids at a current grid level.
!
!    Local, integer ( kind = 4 ) ND, the number of points used in a Lagrange grid.
!
!    Local, integer ( kind = 4 ) ND_TOTAL, the total number of points used in all the
!    Lagrange interpolants for a given sparse interpolant points that occur
!    in more than one Lagrange grid are counted multiple times.
!
!    Local, integer ( kind = 4 ) NI, the number of interpolant evaluation points.
!
!    Local, integer ( kind = 4 ) SPARSE_MIN, the minimum sparse grid level to try.
!
!    Local, real ( kind = 8 ) XD(M,ND), the data points for a Lagrange grid.
!
!    Local, real ( kind = 8 ) XI(M,NI), the interpolant evaluation points.
!
!    Local, real ( kind = 8 ) ZD(ND), the data values for a Lagrange grid.
!
!    Local, real ( kind = 8 ) ZE(NI), the exact function values at XI.
!
!    Local, real ( kind = 8 ) ZI(NI), the sparse interpolant values at XI.
!
!    Local, real ( kind = 8 ) ZPI(NI), one set of Lagrange interpolant values at XI.
!
  implicit none
  !type(state_type), intent(in), dimension(:,:,:) :: state
  real, dimension(:,:,:,:), allocatable :: snapmatrix_volumefraction ! allocate(snapmatrix_pressure(no_smolyak_nodes,p_nodes,snapshots,numphase))
  real, dimension(:,:,:), intent(out), allocatable ::  interpolate_volumefraction 
  real ( kind = 8 ), allocatable :: a(:)
  real ( kind = 8 ) app_error
  real ( kind = 8 ), allocatable :: b(:)
  integer ( kind = 4 ), allocatable :: c(:)
  integer ( kind = 4 ) h
  integer ( kind = 4 ), allocatable :: ind(:)
  integer ( kind = 4 ) l
  integer ( kind = 4 ) l_max
  integer ( kind = 4 ) l_min
  !integer ( kind = 4 ) m
  logical more
  integer ( kind = 4 ) nd
  integer ( kind = 4 ) nd_total
  integer ( kind = 4 ) ni
  real ( kind = 8 ) r8vec_norm_affine
  integer ( kind = 4 ) seed
  !integer ( kind = 4 ) sparse_max
  integer ( kind = 4 ) sparse_min
  integer ( kind = 4 ) t
  integer ( kind = 4 ), allocatable :: w(:)
  real ( kind = 8 ), allocatable :: xd(:,:)
  real ( kind = 8 ), allocatable :: xi(:,:)
  real ( kind = 8 ), allocatable :: zd(:)
  real ( kind = 8 ), allocatable :: ze(:)
  real ( kind = 8 ), allocatable :: zi(:)
  real ( kind = 8 ), allocatable :: zpi(:)
  integer :: i, ii, total, total_timestep,num_pod,j,k,jj,kk, smolyak_total_points,u_nodesi,snapshotsi,dimi,numphasei
   real :: a1,a2,a3,a4,a5  
  allocate (ind(1:m))

  ! Define the region.
  allocate ( a(1:m) )
  allocate ( b(1:m) )

 
   
  !total_timestep=864
  !num_pod=12*36
  !smolyak_total_points=41!2**sparse_max+1
  !print *, '  smolyak_total_points=2**sparse_max+1=', smolyak_total_points 
  
  !  Define the interpolation evaluation information.
  ni = 1
  seed = 123456789
  allocate ( xi(m,ni) )
  !call r8mat_uniform_abvec ( m, ni, a, b, seed, xi )
!  print * , 'xi', xi
   if (.false.) then
  a(1:m) = 1
  b(1:m) = 2
  a1=1
  a2=1.1464!5780000000000
  a3=1.5
  a4=1.8535!421000
  a5=2
  xi(1,1)= 2 
  xi(2,1)= 1.5
  xi(3,1)= 2
  xi(4,1)= 1.55 !1.7 !2.0000000000000000   1.5000000000000000   2.0000000000000000   1.5000000000000000,21
  endif

  if (.true.) then
  a(1:m) = 0.1
  b(1:m) = 0.5
  a1=0.1
  a2=0.158578
  a3=0.3
  a4=0.441421
  a5=0.5
  xi(1,1)= 0.1 
  xi(2,1)= 0.5
  xi(3,1)= 0.5
  xi(4,1)= 0.5 !1.7 !2.0000000000000000   1.5000000000000000   2.0000000000000000   1.5000000000000000,21
   call get_option(&
         '/reduced_model/variables_to_interpolate/dimension_one', xi(1,1))
   call get_option(&
         '/reduced_model/variables_to_interpolate/dimension_two', xi(2,1))
 call get_option(&
         '/reduced_model/variables_to_interpolate/dimension_three', xi(3,1))
 call get_option(&
         '/reduced_model/variables_to_interpolate/dimension_four', xi(4,1))
    xi(5,1)= 0.1 
   xi(6,1)= 0.5
   xi(7,1)= 0.5
   xi(8,1)= 0.5  ! 0.10000000000000000       0.30000000000000000       0.50000000000000000       0.30000000000000000   
  endif
  write ( *, '(a,i6)' ) '  Number of interpolation points is NI = ', ni
  allocate ( ze(1:ni) )
 ! call f_sinr ( m, ni, xi, ze )
  print *, 'actual function value is, will interpolate this value later on,then do comparison', xi,ze
!
!  Compute a sequence of sparse grid interpolants of increasing level.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '   L     Order    ApproxError'
  write ( *, '(a)' ) ''

  allocate ( zpi(1:ni) )
  allocate ( zi(1:ni) )
  allocate(interpolate_volumefraction(p_nodes,snapshots,numphase))
  !smolynode, u_nodes,snapshots,dim,numphase
  Do u_nodesi=1, p_nodes  
   Do snapshotsi=1, snapshots  
   do numphasei=1, numphase
  total = 0
  sparse_min = 0

  ! do l_max = sparse_min, sparse_max
    do l_max = sparse_max, sparse_max
    allocate ( c(0:l_max) )
    allocate ( w(0:l_max) )
    call smolyak_coefficients ( l_max, m, c, w )

    zi(1:ni) = 0.0D+00
    nd_total = 0
    
    l_min = max ( l_max + 1 - m, 0 )

    do l = l_min, l_max

      more = .false.
      
      do
     ! print *,'l_max,l,' ,l_max,l   
!  Define the next product grid.
!
        call comp_next ( l, m, ind, more, h, t )
!
!  Count the grid, find its points, evaluate the data there.
!
        call lagrange_interp_nd_size2 ( m, ind, nd )
        allocate ( xd(1:m,1:nd) )
        call lagrange_interp_nd_grid2 ( m, ind, a, b, nd, xd )        
         !  print *,'xd', xd, 点1 坐标为： xd(:,1)
        allocate ( zd(1:nd) )

       ! print *,'---------------------nd',nd  
        do ii=1, nd           
        if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3 .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3 .and. xd(5,ii) .eq. a3 .and. xd(6,ii).eq. a3 .and. xd(7,ii) .eq. a3 .and. xd(8,ii) .eq. a3 ) then  
           zd(ii)=snapmatrix_volumefraction(1,u_nodesi,snapshotsi, numphasei)
   else if(xd(1,ii) .eq. a1 .and. xd(2,ii) .eq. a3 .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3 .and. xd(5,ii) .eq. a3 .and. xd(6,ii).eq. a3 .and. xd(7,ii) .eq. a3 .and. xd(8,ii) .eq. a3 ) then  
           zd(ii)=snapmatrix_volumefraction(2,u_nodesi,snapshotsi, numphasei)
   else if(xd(1,ii) .eq. a5 .and. xd(2,ii) .eq. a3 .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3 .and. xd(5,ii) .eq. a3 .and. xd(6,ii).eq. a3 .and. xd(7,ii) .eq. a3 .and. xd(8,ii) .eq. a3 ) then 
           zd(ii)=snapmatrix_volumefraction(3,u_nodesi,snapshotsi, numphasei)
   else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a1 .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3 .and. xd(5,ii) .eq. a3 .and. xd(6,ii).eq. a3 .and. xd(7,ii) .eq. a3 .and. xd(8,ii) .eq. a3 ) then 
           zd(ii)=snapmatrix_volumefraction(4,u_nodesi,snapshotsi, numphasei)
   else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a5 .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3 .and. xd(5,ii) .eq. a3 .and. xd(6,ii).eq. a3 .and. xd(7,ii) .eq. a3 .and. xd(8,ii) .eq. a3 ) then 
           zd(ii)=snapmatrix_volumefraction(5,u_nodesi,snapshotsi, numphasei)
   else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3 .and. xd(3,ii) .eq. a1 .and. xd(4,ii) .eq. a3 .and. xd(5,ii) .eq. a3 .and. xd(6,ii).eq. a3 .and. xd(7,ii) .eq. a3 .and. xd(8,ii) .eq. a3 ) then 
           zd(ii)=snapmatrix_volumefraction(6,u_nodesi,snapshotsi, numphasei)
   else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3 .and. xd(3,ii) .eq. a5 .and. xd(4,ii) .eq. a3 .and. xd(5,ii) .eq. a3 .and. xd(6,ii).eq. a3 .and. xd(7,ii) .eq. a3 .and. xd(8,ii) .eq. a3) then 
           zd(ii)=snapmatrix_volumefraction(7,u_nodesi,snapshotsi, numphasei)
   else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3 .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a1 .and. xd(5,ii) .eq. a3 .and. xd(6,ii).eq. a3 .and. xd(7,ii) .eq. a3 .and. xd(8,ii) .eq. a3) then 
           zd(ii)=snapmatrix_volumefraction(8,u_nodesi,snapshotsi, numphasei)
   else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3 .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a5 .and. xd(5,ii) .eq. a3 .and. xd(6,ii).eq. a3 .and. xd(7,ii) .eq. a3 .and. xd(8,ii) .eq. a3) then
           zd(ii)=snapmatrix_volumefraction(9,u_nodesi,snapshotsi, numphasei)
   else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3 .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3 .and. xd(5,ii) .eq. a1 .and. xd(6,ii).eq. a3 .and. xd(7,ii) .eq. a3 .and. xd(8,ii) .eq. a3) then              
           zd(ii)=snapmatrix_volumefraction(10,u_nodesi,snapshotsi, numphasei)
  else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3 .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3 .and. xd(5,ii) .eq. a5 .and. xd(6,ii).eq. a3 .and. xd(7,ii) .eq. a3 .and. xd(8,ii) .eq. a3) then           
           zd(ii)=snapmatrix_volumefraction(11,u_nodesi,snapshotsi, numphasei)
 else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3 .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3 .and. xd(5,ii) .eq. a3 .and. xd(6,ii).eq. a1 .and. xd(7,ii) .eq. a3 .and. xd(8,ii) .eq. a3) then
           zd(ii)=snapmatrix_volumefraction(12,u_nodesi,snapshotsi, numphasei)
 else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3 .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3 .and. xd(5,ii) .eq. a3 .and. xd(6,ii).eq. a5 .and. xd(7,ii) .eq. a3 .and. xd(8,ii) .eq. a3) then
           zd(ii)=snapmatrix_volumefraction(13,u_nodesi,snapshotsi, numphasei)
 else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3 .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3 .and. xd(5,ii) .eq. a3 .and. xd(6,ii).eq. a3 .and. xd(7,ii) .eq. a1 .and. xd(8,ii) .eq. a3) then
           zd(ii)=snapmatrix_volumefraction(14,u_nodesi,snapshotsi, numphasei)
 else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3 .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3 .and. xd(5,ii) .eq. a3 .and. xd(6,ii).eq. a3 .and. xd(7,ii) .eq. a5 .and. xd(8,ii) .eq. a3) then
           zd(ii)=snapmatrix_volumefraction(15,u_nodesi,snapshotsi, numphasei)
 else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3 .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3 .and. xd(5,ii) .eq. a3 .and. xd(6,ii).eq. a3 .and. xd(7,ii) .eq. a3 .and. xd(8,ii) .eq. a1) then
            zd(ii)=snapmatrix_volumefraction(16,u_nodesi,snapshotsi, numphasei)
 else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3 .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3 .and. xd(5,ii) .eq. a3 .and. xd(6,ii).eq. a3 .and. xd(7,ii) .eq. a3 .and. xd(8,ii) .eq. a5) then
           zd(ii)=snapmatrix_volumefraction(17,u_nodesi,snapshotsi, numphasei)
         endif
        enddo
          !call f_sinr ( m, nd, xd, zd )       
           
          ! print *, 'zd',  zd(:)
            
!  Use the grid to evaluate the interpolant.
        call lagrange_interp_nd_value2 ( m, ind, a, b, nd, zd, ni, xi, zpi )
 
!  Weighted the interpolant values and add to the sparse grid interpolant.
 
        nd_total = nd_total + nd
        zi(1:ni) = zi(1:ni) + c(l) * zpi(1:ni)
        !  print *, ' nd_total,level',  nd_total, l, c(l), nd
        !  print *, 'zizizizizizizic(l)zpi', zi(1)!,c(l),zpi(1)
        deallocate ( xd )
        deallocate ( zd )     
       
       if ( .not. more ) then
          exit
        end if    
 
      end do !DO

    end do !do l = l_min, l_max   
    deallocate ( c )
    deallocate ( w )
    
   end do ! do l_max = sparse_max, sparse_max
   interpolate_volumefraction(u_nodesi,snapshotsi,numphasei)=zi(ni)
   enddo !Do u_nodesi  
  enddo !Do snapshotsi 
  enddo !numphasei

 !  print *, ' total', total , nd_total
 ! open(unit=66,file='podbasis_interp') 
 !  do j=1, total_timestep 
 !     write(66,*)  pod_coef_all_obv_41i(j,1:num_pod)
 ! enddo
 ! close(66)
 
  deallocate ( a )
  deallocate ( b )
  deallocate ( ind )
  deallocate ( xi )
  deallocate ( ze )
  deallocate ( zi )
  deallocate ( zpi )
  end subroutine interpolate_snapmatrix_volumefraction_8dim_level1

 subroutine interpolate_snapmatrix_pressure_4dim_level1( snapmatrix_pressure, interpolate_pressure)  
    
! top-bottom between 1-2, pod_coef_all_obv is used for storing some podbasis to avoid to much modifications.
 
!  Parameters:
!
!    Local, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) SPARSE_MAX, the maximum sparse grid level to try.
!
!  Local Parameters:
!
!    Local, real ( kind = 8 ) A(M), B(M), the upper and lower variable limits 
!    in each dimension.
!
!    Local, real ( kind = 8 ) APP_ERROR, the averaged Euclidean norm of the 
!    difference between the sparse interpolant and the exact function at 
!    the interpolation points.
!
!    Local, integer ( kind = 4 ) C(0:L_MAX), the sparse grid coefficient vector.
!    Results at level L are weighted by C(L).
!
!    Local, integer ( kind = 4 ) IND(M), the 1D indices defining a Lagrange grid.
!    Each entry is a 1d "level" that specifies the order of a 
!    Clenshaw Curtis 1D grid.
!
!    Local, integer ( kind = 4 ) L, the current Lagrange grid level.
!
!    Local, integer ( kind = 4 ) L_MAX, the current sparse grid level.
!
!    Local, integer ( kind = 4 ) MORE, used to control the enumeration of all the
!    Lagrange grids at a current grid level.
!
!    Local, integer ( kind = 4 ) ND, the number of points used in a Lagrange grid.
!
!    Local, integer ( kind = 4 ) ND_TOTAL, the total number of points used in all the
!    Lagrange interpolants for a given sparse interpolant points that occur
!    in more than one Lagrange grid are counted multiple times.
!
!    Local, integer ( kind = 4 ) NI, the number of interpolant evaluation points.
!
!    Local, integer ( kind = 4 ) SPARSE_MIN, the minimum sparse grid level to try.
!
!    Local, real ( kind = 8 ) XD(M,ND), the data points for a Lagrange grid.
!
!    Local, real ( kind = 8 ) XI(M,NI), the interpolant evaluation points.
!
!    Local, real ( kind = 8 ) ZD(ND), the data values for a Lagrange grid.
!
!    Local, real ( kind = 8 ) ZE(NI), the exact function values at XI.
!
!    Local, real ( kind = 8 ) ZI(NI), the sparse interpolant values at XI.
!
!    Local, real ( kind = 8 ) ZPI(NI), one set of Lagrange interpolant values at XI.
!
  implicit none
  !type(state_type), intent(in), dimension(:,:,:) :: state
  real, dimension(:,:,:,:), allocatable :: snapmatrix_pressure ! allocate(snapmatrix_pressure(no_smolyak_nodes,p_nodes,snapshots,numphase))
  real, dimension(:,:,:), intent(out), allocatable ::  interpolate_pressure 
  real ( kind = 8 ), allocatable :: a(:)
  real ( kind = 8 ) app_error
  real ( kind = 8 ), allocatable :: b(:)
  integer ( kind = 4 ), allocatable :: c(:)
  integer ( kind = 4 ) h
  integer ( kind = 4 ), allocatable :: ind(:)
  integer ( kind = 4 ) l
  integer ( kind = 4 ) l_max
  integer ( kind = 4 ) l_min
  !integer ( kind = 4 ) m
  logical more
  integer ( kind = 4 ) nd
  integer ( kind = 4 ) nd_total
  integer ( kind = 4 ) ni
  real ( kind = 8 ) r8vec_norm_affine
  integer ( kind = 4 ) seed
 ! integer ( kind = 4 ) sparse_max
  integer ( kind = 4 ) sparse_min
  integer ( kind = 4 ) t
  integer ( kind = 4 ), allocatable :: w(:)
  real ( kind = 8 ), allocatable :: xd(:,:)
  real ( kind = 8 ), allocatable :: xi(:,:)
  real ( kind = 8 ), allocatable :: zd(:)
  real ( kind = 8 ), allocatable :: ze(:)
  real ( kind = 8 ), allocatable :: zi(:)
  real ( kind = 8 ), allocatable :: zpi(:)
  integer :: i, ii, total, total_timestep,num_pod,j,k,jj,kk, smolyak_total_points,u_nodesi,snapshotsi,dimi,numphasei
   real :: a1,a2,a3,a4,a5 
  allocate (ind(1:m))

  ! Define the region.
  allocate ( a(1:m) )
  allocate ( b(1:m) )

  
   
  !total_timestep=864
  !num_pod=12*36
  !smolyak_total_points=41!2**sparse_max+1
  !print *, '  smolyak_total_points=2**sparse_max+1=', smolyak_total_points 
  !  stop 2323
  !  Define the interpolation evaluation information.
  ni = 1
  seed = 123456789
  allocate ( xi(m,ni) )
  !call r8mat_uniform_abvec ( m, ni, a, b, seed, xi )
!  print * , 'xi', xi
 
   if (.false.) then
  a(1:m) = 1
  b(1:m) = 2
  a1=1
  a2=1.1464!5780000000000
  a3=1.5
  a4=1.8535!421000
  a5=2
  xi(1,1)= 2 
  xi(2,1)= 1.5
  xi(3,1)= 2
  xi(4,1)= 1.55 !1.7 !2.0000000000000000   1.5000000000000000   2.0000000000000000   1.5000000000000000,21
  endif

  if (.true.) then
  a(1:m) = 0.1
  b(1:m) = 0.5
  a1=0.1
  a2=0.158578
  a3=0.3
  a4=0.441421
  a5=0.5
  xi(1,1)= 0.1 
  xi(2,1)= 0.5
  xi(3,1)= 0.5
  xi(4,1)= 0.5 !1.7 !2.0000000000000000   1.5000000000000000   2.0000000000000000   1.5000000000000000,21
   call get_option(&
         '/reduced_model/variables_to_interpolate/dimension_one', xi(1,1))
   call get_option(&
         '/reduced_model/variables_to_interpolate/dimension_two', xi(2,1))
 call get_option(&
         '/reduced_model/variables_to_interpolate/dimension_three', xi(3,1))
 call get_option(&
         '/reduced_model/variables_to_interpolate/dimension_four', xi(4,1))
  endif
  write ( *, '(a,i6)' ) '  Number of interpolation points is NI = ', ni
  allocate ( ze(1:ni) )
 ! call f_sinr ( m, ni, xi, ze )
  print *, 'actual function value is, will interpolate this value later on,then do comparison', xi,ze
!
!  Compute a sequence of sparse grid interpolants of increasing level.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '   L     Order    ApproxError'
  write ( *, '(a)' ) ''
  allocate(interpolate_pressure(p_nodes,snapshots,numphase)) 
  allocate ( zpi(1:ni) )
  allocate ( zi(1:ni) )
  !smolynode, u_nodes,snapshots,dim,numphase
  Do u_nodesi=1, p_nodes  
   Do snapshotsi=1, snapshots  
   do numphasei=1, numphase
  total = 0
  sparse_min = 0

  ! do l_max = sparse_min, sparse_max
    do l_max = sparse_max, sparse_max
    allocate ( c(0:l_max) )
    allocate ( w(0:l_max) )
    call smolyak_coefficients ( l_max, m, c, w )

    zi(1:ni) = 0.0D+00
    nd_total = 0
    
    l_min = max ( l_max + 1 - m, 0 )

    do l = l_min, l_max

      more = .false.
      
      do
     ! print *,'l_max,l,' ,l_max,l   
!  Define the next product grid.
!
        call comp_next ( l, m, ind, more, h, t )
!
!  Count the grid, find its points, evaluate the data there.
!
        call lagrange_interp_nd_size2 ( m, ind, nd )
        allocate ( xd(1:m,1:nd) )
        call lagrange_interp_nd_grid2 ( m, ind, a, b, nd, xd )        
         !  print *,'xd', xd, 点1 坐标为： xd(:,1)
        allocate ( zd(1:nd) )

       ! print *,'---------------------nd',nd  
        do ii=1, nd           
                if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3   .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3) then 
!( 1.5000000000000000    1.5000000000000000   1.5000000000000000  1.5000000000000000)
           zd(ii)=snapmatrix_pressure(1,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a1 .and. xd(2,ii) .eq. a3   .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3) then 
!(1.0000000000000000   1.5000000000000000   1.5000000000000000   1.5000000000000000)
           zd(ii)=snapmatrix_pressure(2,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a5 .and. xd(2,ii) .eq. a3   .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3) then
 !( 2.0000000000000000   1.5000000000000000   1.5000000000000000   1.5000000000000000)
           zd(ii)=snapmatrix_pressure(3,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a1   .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3) then 
!(1.5000000000000000   1.0000000000000000   1.5000000000000000   1.5000000000000000 )
           zd(ii)=snapmatrix_pressure(4,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a5 .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3) then 
!(1.5000000000000000   2.0000000000000000   1.5000000000000000   1.5000000000000000 )
           zd(ii)=snapmatrix_pressure(5,u_nodesi,snapshotsi, numphasei)
        else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3 .and. xd(3,ii) .eq. a1 .and. xd(4,ii) .eq. a3) then 
!(1.5000000000000000   1.5000000000000000   1.0000000000000000   1.5000000000000000  )
           zd(ii)=snapmatrix_pressure(6,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3 .and. xd(3,ii) .eq. a5 .and. xd(4,ii) .eq. a3) then 
!(1.5000000000000000   1.5000000000000000   2.0000000000000000   a3000000000000000  )
           zd(ii)=snapmatrix_pressure(7,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3 .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a1) then
 !(a3000000000000000   a3000000000000000   a3000000000000000   1.0000000000000000   )
           zd(ii)=snapmatrix_pressure(8,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3   .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a5) then 
!(  a3000000000000000   a3000000000000000   a3000000000000000   2.0000000000000000  )
           zd(ii)=snapmatrix_pressure(9,u_nodesi,snapshotsi, numphasei) 
         else if(int(xd(1,ii) *1000000)/1000000 .eq. a2 .and. xd(2,ii) .eq. a3 .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3) then 
               !(a24466094067263   a3000000000000000   a3000000000000000   a3000000000000000 .and. xd(1,ii) <= 1.15 ..and. xd(2:4,ii) .eq. a3  )
           zd(ii)=snapmatrix_pressure(10,u_nodesi,snapshotsi, numphasei)
         elseif(int(xd(1,ii) *1000000)/1000000 .eq. a4 .and. xd(2,ii) .eq. a3 .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3) then 
               !(a45533905932737   a3000000000000000   a3000000000000000   a3000000000000000.and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3 )
           zd(ii)=snapmatrix_pressure(11,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a1 .and. xd(2,ii) .eq. a1 .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3) then !(  1.0000000000000000   1.0000000000000000   a3000000000000000   a3000000000000000 )
           zd(ii)=snapmatrix_pressure(12,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a5 .and. xd(2,ii) .eq. a1 .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3) then !( 2.0000000000000000   1.0000000000000000   a3000000000000000   a3000000000000000)
           zd(ii)=snapmatrix_pressure(13,u_nodesi,snapshotsi, numphasei)
          else if(xd(1,ii) .eq. a1 .and. xd(2,ii) .eq.a5  .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3) then !(1.0000000000000000   2.0000000000000000   a3000000000000000   a3000000000000000)
           zd(ii)=snapmatrix_pressure(14,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a5 .and. xd(2,ii) .eq.a5  .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3) then !(2.0000000000000000   2.0000000000000000   a3000000000000000   a3000000000000000)
           zd(ii)=snapmatrix_pressure(15,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a3 .and. int(xd(2,ii) *1000000)/1000000 .eq. a2 .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3) then 
            !(a3000000000000000   a24466094067263   a3000000000000000   a3000000000000000 xd(2,ii) >=1.14  .and. xd(2,ii) <=1.15)
           zd(ii)=snapmatrix_pressure(16,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a3 .and. int(xd(2,ii) *1000000)/1000000 .eq. a4 .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3) then
          !(a3000000000000000   a45533905932737   a3000000000000000   a3000000000000000)
           zd(ii)=snapmatrix_pressure(17,u_nodesi,snapshotsi, numphasei)
           else if(xd(1,ii) .eq. a1 .and. xd(2,ii) .eq. a3 .and. xd(3,ii) .eq. a1 .and. xd(4,ii) .eq. a3) then !( 1.0000000000000000   a3000000000000000   1.0000000000000000  a3000000000000000)
           zd(ii)=snapmatrix_pressure(18,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq.a5.and. xd(2,ii) .eq. a3   .and. xd(3,ii) .eq. a1 .and. xd(4,ii) .eq. a3) then !(2.0000000000000000   a3000000000000000   1.0000000000000000   a3000000000000000)
           zd(ii)=snapmatrix_pressure(19,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a1 .and. xd(2,ii) .eq. a3   .and. xd(3,ii) .eq. a5 .and. xd(4,ii) .eq. a3) then !(1.0000000000000000   a3000000000000000   2.0000000000000000   a3000000000000000)
           zd(ii)=snapmatrix_pressure(20,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a5 .and. xd(2,ii) .eq. a3   .and. xd(3,ii) .eq. a5 .and. xd(4,ii) .eq. a3) then !(2.0000000000000000   a3000000000000000   2.0000000000000000   a3000000000000000)
           zd(ii)=snapmatrix_pressure(21,u_nodesi,snapshotsi, numphasei) 
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a1   .and. xd(3,ii) .eq. a1 .and. xd(4,ii) .eq. a3) then !(a3000000000000000   1.0000000000000000   1.0000000000000000   a3000000000000000)
           zd(ii)=snapmatrix_pressure(22,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a5   .and. xd(3,ii) .eq. a1 .and. xd(4,ii) .eq. a3) then !(a3000000000000000   2.0000000000000000   1.0000000000000000   a3000000000000000)
           zd(ii)=snapmatrix_pressure(23,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a1   .and. xd(3,ii) .eq. a5 .and. xd(4,ii) .eq. a3) then !(a3000000000000000   1.0000000000000000   2.0000000000000000   a3000000000000000)
           zd(ii)=snapmatrix_pressure(24,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a5   .and. xd(3,ii) .eq. a5 .and. xd(4,ii) .eq. a3) then !(a3000000000000000   2.0000000000000000   2.0000000000000000   a3000000000000000)
           zd(ii)=snapmatrix_pressure(25,u_nodesi,snapshotsi, numphasei)
          else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3 .and. int(xd(3,ii) *1000000)/1000000 .eq. a2 .and. xd(4,ii) .eq. a3) then 
       !(a3000000000000000   a3000000000000000  a24466094067263   a3000000000000000)
           zd(ii)=snapmatrix_pressure(26,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a3  .and. xd(2,ii) .eq. a3  .and. int(xd(3,ii) *1000000)/1000000 .eq. a4 .and. xd(4,ii) .eq. a3) then 
       !( a3000000000000000   a3000000000000000   a45533905932737   a3000000000000000 )
           zd(ii)=snapmatrix_pressure(27,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a1 .and. xd(2,ii) .eq. a3   .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a1) then !(1.0000000000000000   a3000000000000000   a3000000000000000   1.0000000000000000)
           zd(ii)=snapmatrix_pressure(28,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a5 .and. xd(2,ii) .eq. a3   .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a1) then !(2.0000000000000000   a3000000000000000   a3000000000000000   1.0000000000000000)
           zd(ii)=snapmatrix_pressure(29,u_nodesi,snapshotsi, numphasei)
           else if(xd(1,ii) .eq. a1 .and. xd(2,ii) .eq. a3   .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a5) then !(1.0000000000000000   a3000000000000000   a3000000000000000   2.0000000000000000)
           zd(ii)=snapmatrix_pressure(30,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a5 .and. xd(2,ii) .eq. a3   .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a5) then !(2.0000000000000000   a3000000000000000   a3000000000000000   2.0000000000000000)
           zd(ii)=snapmatrix_pressure(31,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a1   .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a1) then !(a3000000000000000   1.0000000000000000   a3000000000000000   1.0000000000000000)
           zd(ii)=snapmatrix_pressure(32,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a5  .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a1) then !(a3000000000000000   2.0000000000000000   a3000000000000000   1.0000000000000000)
           zd(ii)=snapmatrix_pressure(33,u_nodesi,snapshotsi, numphasei) 
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a1   .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a5) then !(a3000000000000000   1.0000000000000000   a3000000000000000   2.0000000000000000)
           zd(ii)=snapmatrix_pressure(34,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a5   .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a5) then !(a3000000000000000   2.0000000000000000   a3000000000000000   2.0000000000000000)
           zd(ii)=snapmatrix_pressure(35,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3   .and. xd(3,ii) .eq. a1 .and. xd(4,ii) .eq. a1) then !(a3000000000000000   a3000000000000000   1.0000000000000000   1.0000000000000000)
           zd(ii)=snapmatrix_pressure(36,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3   .and. xd(3,ii) .eq. a5 .and. xd(4,ii) .eq. a1) then !( a3000000000000000   a3000000000000000   2.0000000000000000   1.0000000000000000)
           zd(ii)=snapmatrix_pressure(37,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3   .and. xd(3,ii) .eq. a1 .and. xd(4,ii) .eq. a5) then !(a3000000000000000   a3000000000000000   1.0000000000000000   2.0000000000000000)
           zd(ii)=snapmatrix_pressure(38,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3   .and. xd(3,ii) .eq. a5 .and. xd(4,ii) .eq. a5) then !(a3000000000000000   a3000000000000000   2.0000000000000000   2.0000000000000000)
           zd(ii)=snapmatrix_pressure(39,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3 .and. xd(3,ii) .eq. a3 .and. int(xd(4,ii) *1000000)/1000000 .eq. a2) then 
      !(a3000000000000000   a3000000000000000   a3000000000000000   a24466094067263 )
           zd(ii)=snapmatrix_pressure(40,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a3.and. xd(2,ii) .eq. a3 .and. xd(3,ii) .eq. a3.and. int(xd(4,ii) *1000000)/1000000 .eq. a4) then 
     !(a3000000000000000   a3000000000000000   a3000000000000000   1.8535533905932737     )
           zd(ii)=snapmatrix_pressure(41,u_nodesi,snapshotsi, numphasei)
         endif
        enddo
          !call f_sinr ( m, nd, xd, zd )       
           
          ! print *, 'zd',  zd(:)
            
!  Use the grid to evaluate the interpolant.
        call lagrange_interp_nd_value2 ( m, ind, a, b, nd, zd, ni, xi, zpi )
 
!  Weighted the interpolant values and add to the sparse grid interpolant.
 
        nd_total = nd_total + nd
        zi(1:ni) = zi(1:ni) + c(l) * zpi(1:ni)
        !  print *, ' nd_total,level',  nd_total, l, c(l), nd
        !  print *, 'zizizizizizizic(l)zpi', zi(1)!,c(l),zpi(1)
        deallocate ( xd )
        deallocate ( zd )     
       
       if ( .not. more ) then
          exit
        end if    
 
      end do !DO

    end do !do l = l_min, l_max   
    deallocate ( c )
    deallocate ( w )
    
   end do ! do l_max = sparse_max, sparse_max
   interpolate_pressure(u_nodesi,snapshotsi,numphasei)=zi(ni)
   enddo !Do u_nodesi  
  enddo !Do snapshotsi 
  enddo !numphasei

 !  print *, ' total', total , nd_total
 ! open(unit=66,file='podbasis_interp') 
 !  do j=1, total_timestep 
 !     write(66,*)  pod_coef_all_obv_41i(j,1:num_pod)
 ! enddo
 ! close(66)
 
  deallocate ( a )
  deallocate ( b )
  deallocate ( ind )
  deallocate ( xi )
  deallocate ( ze )
  deallocate ( zi )
  deallocate ( zpi )
  end subroutine interpolate_snapmatrix_pressure_4dim_level1

 subroutine interpolate_snapmatrix_volumefraction_4dim_level1( snapmatrix_volumefraction, interpolate_volumefraction)  
    
! top-bottom between 1-2, pod_coef_all_obv is used for storing some podbasis to avoid to much modifications.
 
!  Parameters:
!
!    Local, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) SPARSE_MAX, the maximum sparse grid level to try.
!
!  Local Parameters:
!
!    Local, real ( kind = 8 ) A(M), B(M), the upper and lower variable limits 
!    in each dimension.
!
!    Local, real ( kind = 8 ) APP_ERROR, the averaged Euclidean norm of the 
!    difference between the sparse interpolant and the exact function at 
!    the interpolation points.
!
!    Local, integer ( kind = 4 ) C(0:L_MAX), the sparse grid coefficient vector.
!    Results at level L are weighted by C(L).
!
!    Local, integer ( kind = 4 ) IND(M), the 1D indices defining a Lagrange grid.
!    Each entry is a 1d "level" that specifies the order of a 
!    Clenshaw Curtis 1D grid.
!
!    Local, integer ( kind = 4 ) L, the current Lagrange grid level.
!
!    Local, integer ( kind = 4 ) L_MAX, the current sparse grid level.
!
!    Local, integer ( kind = 4 ) MORE, used to control the enumeration of all the
!    Lagrange grids at a current grid level.
!
!    Local, integer ( kind = 4 ) ND, the number of points used in a Lagrange grid.
!
!    Local, integer ( kind = 4 ) ND_TOTAL, the total number of points used in all the
!    Lagrange interpolants for a given sparse interpolant points that occur
!    in more than one Lagrange grid are counted multiple times.
!
!    Local, integer ( kind = 4 ) NI, the number of interpolant evaluation points.
!
!    Local, integer ( kind = 4 ) SPARSE_MIN, the minimum sparse grid level to try.
!
!    Local, real ( kind = 8 ) XD(M,ND), the data points for a Lagrange grid.
!
!    Local, real ( kind = 8 ) XI(M,NI), the interpolant evaluation points.
!
!    Local, real ( kind = 8 ) ZD(ND), the data values for a Lagrange grid.
!
!    Local, real ( kind = 8 ) ZE(NI), the exact function values at XI.
!
!    Local, real ( kind = 8 ) ZI(NI), the sparse interpolant values at XI.
!
!    Local, real ( kind = 8 ) ZPI(NI), one set of Lagrange interpolant values at XI.
!
  implicit none
  !type(state_type), intent(in), dimension(:,:,:) :: state
  real, dimension(:,:,:,:), allocatable :: snapmatrix_volumefraction ! allocate(snapmatrix_pressure(no_smolyak_nodes,p_nodes,snapshots,numphase))
  real, dimension(:,:,:), intent(out), allocatable ::  interpolate_volumefraction 
  real ( kind = 8 ), allocatable :: a(:)
  real ( kind = 8 ) app_error
  real ( kind = 8 ), allocatable :: b(:)
  integer ( kind = 4 ), allocatable :: c(:)
  integer ( kind = 4 ) h
  integer ( kind = 4 ), allocatable :: ind(:)
  integer ( kind = 4 ) l
  integer ( kind = 4 ) l_max
  integer ( kind = 4 ) l_min
  !integer ( kind = 4 ) m
  logical more
  integer ( kind = 4 ) nd
  integer ( kind = 4 ) nd_total
  integer ( kind = 4 ) ni
  real ( kind = 8 ) r8vec_norm_affine
  integer ( kind = 4 ) seed
  !integer ( kind = 4 ) sparse_max
  integer ( kind = 4 ) sparse_min
  integer ( kind = 4 ) t
  integer ( kind = 4 ), allocatable :: w(:)
  real ( kind = 8 ), allocatable :: xd(:,:)
  real ( kind = 8 ), allocatable :: xi(:,:)
  real ( kind = 8 ), allocatable :: zd(:)
  real ( kind = 8 ), allocatable :: ze(:)
  real ( kind = 8 ), allocatable :: zi(:)
  real ( kind = 8 ), allocatable :: zpi(:)
  integer :: i, ii, total, total_timestep,num_pod,j,k,jj,kk, smolyak_total_points,u_nodesi,snapshotsi,dimi,numphasei
   real :: a1,a2,a3,a4,a5  
  allocate (ind(1:m))

  ! Define the region.
  allocate ( a(1:m) )
  allocate ( b(1:m) )

 
   
  !total_timestep=864
  !num_pod=12*36
  !smolyak_total_points=41!2**sparse_max+1
  !print *, '  smolyak_total_points=2**sparse_max+1=', smolyak_total_points 
  
  !  Define the interpolation evaluation information.
  ni = 1
  seed = 123456789
  allocate ( xi(m,ni) )
  !call r8mat_uniform_abvec ( m, ni, a, b, seed, xi )
!  print * , 'xi', xi
   if (.false.) then
  a(1:m) = 1
  b(1:m) = 2
  a1=1
  a2=1.1464!5780000000000
  a3=1.5
  a4=1.8535!421000
  a5=2
  xi(1,1)= 2 
  xi(2,1)= 1.5
  xi(3,1)= 2
  xi(4,1)= 1.55 !1.7 !2.0000000000000000   1.5000000000000000   2.0000000000000000   1.5000000000000000,21
  endif

  if (.true.) then
  a(1:m) = 0.1
  b(1:m) = 0.5
  a1=0.1
  a2=0.158578
  a3=0.3
  a4=0.441421
  a5=0.5
  xi(1,1)= 0.1 
  xi(2,1)= 0.5
  xi(3,1)= 0.5
  xi(4,1)= 0.5 !1.7 !2.0000000000000000   1.5000000000000000   2.0000000000000000   1.5000000000000000,21
   call get_option(&
         '/reduced_model/variables_to_interpolate/dimension_one', xi(1,1))
   call get_option(&
         '/reduced_model/variables_to_interpolate/dimension_two', xi(2,1))
 call get_option(&
         '/reduced_model/variables_to_interpolate/dimension_three', xi(3,1))
 call get_option(&
         '/reduced_model/variables_to_interpolate/dimension_four', xi(4,1))
  endif
  write ( *, '(a,i6)' ) '  Number of interpolation points is NI = ', ni
  allocate ( ze(1:ni) )
 ! call f_sinr ( m, ni, xi, ze )
  print *, 'actual function value is, will interpolate this value later on,then do comparison', xi,ze
!
!  Compute a sequence of sparse grid interpolants of increasing level.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '   L     Order    ApproxError'
  write ( *, '(a)' ) ''

  allocate ( zpi(1:ni) )
  allocate ( zi(1:ni) )
  allocate(interpolate_volumefraction(p_nodes,snapshots,numphase))
  !smolynode, u_nodes,snapshots,dim,numphase
  Do u_nodesi=1, p_nodes  
   Do snapshotsi=1, snapshots  
   do numphasei=1, numphase
  total = 0
  sparse_min = 0

  ! do l_max = sparse_min, sparse_max
    do l_max = sparse_max, sparse_max
    allocate ( c(0:l_max) )
    allocate ( w(0:l_max) )
    call smolyak_coefficients ( l_max, m, c, w )

    zi(1:ni) = 0.0D+00
    nd_total = 0
    
    l_min = max ( l_max + 1 - m, 0 )

    do l = l_min, l_max

      more = .false.
      
      do
     ! print *,'l_max,l,' ,l_max,l   
!  Define the next product grid.
!
        call comp_next ( l, m, ind, more, h, t )
!
!  Count the grid, find its points, evaluate the data there.
!
        call lagrange_interp_nd_size2 ( m, ind, nd )
        allocate ( xd(1:m,1:nd) )
        call lagrange_interp_nd_grid2 ( m, ind, a, b, nd, xd )        
         !  print *,'xd', xd, 点1 坐标为： xd(:,1)
        allocate ( zd(1:nd) )

       ! print *,'---------------------nd',nd  
        do ii=1, nd           
         if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3   .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3) then 
!( 1.5000000000000000    1.5000000000000000   1.5000000000000000  1.5000000000000000)
           zd(ii)=snapmatrix_volumefraction(1,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a1 .and. xd(2,ii) .eq. a3   .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3) then 
!(1.0000000000000000   1.5000000000000000   1.5000000000000000   1.5000000000000000)
           zd(ii)=snapmatrix_volumefraction(2,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a5 .and. xd(2,ii) .eq. a3   .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3) then
 !( 2.0000000000000000   1.5000000000000000   1.5000000000000000   1.5000000000000000)
           zd(ii)=snapmatrix_volumefraction(3,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a1   .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3) then 
!(1.5000000000000000   1.0000000000000000   1.5000000000000000   1.5000000000000000 )
           zd(ii)=snapmatrix_volumefraction(4,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a5 .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3) then 
!(1.5000000000000000   2.0000000000000000   1.5000000000000000   1.5000000000000000 )
           zd(ii)=snapmatrix_volumefraction(5,u_nodesi,snapshotsi, numphasei)
        else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3 .and. xd(3,ii) .eq. a1 .and. xd(4,ii) .eq. a3) then 
!(1.5000000000000000   1.5000000000000000   1.0000000000000000   1.5000000000000000  )
           zd(ii)=snapmatrix_volumefraction(6,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3 .and. xd(3,ii) .eq. a5 .and. xd(4,ii) .eq. a3) then 
!(1.5000000000000000   1.5000000000000000   2.0000000000000000   a3000000000000000  )
           zd(ii)=snapmatrix_volumefraction(7,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3 .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a1) then
 !(a3000000000000000   a3000000000000000   a3000000000000000   1.0000000000000000   )
           zd(ii)=snapmatrix_volumefraction(8,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3   .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a5) then 
!(  a3000000000000000   a3000000000000000   a3000000000000000   2.0000000000000000  )
           zd(ii)=snapmatrix_volumefraction(9,u_nodesi,snapshotsi, numphasei) 
         else if(int(xd(1,ii) *1000000)/1000000 .eq. a2 .and. xd(2,ii) .eq. a3 .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3) then 
               !(a24466094067263   a3000000000000000   a3000000000000000   a3000000000000000 .and. xd(1,ii) <= 1.15 ..and. xd(2:4,ii) .eq. a3  )
           zd(ii)=snapmatrix_volumefraction(10,u_nodesi,snapshotsi, numphasei)
         elseif(int(xd(1,ii) *1000000)/1000000 .eq. a4 .and. xd(2,ii) .eq. a3 .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3) then 
               !(a45533905932737   a3000000000000000   a3000000000000000   a3000000000000000.and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3 )
           zd(ii)=snapmatrix_volumefraction(11,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a1 .and. xd(2,ii) .eq. a1 .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3) then !(  1.0000000000000000   1.0000000000000000   a3000000000000000   a3000000000000000 )
           zd(ii)=snapmatrix_volumefraction(12,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a5 .and. xd(2,ii) .eq. a1 .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3) then !( 2.0000000000000000   1.0000000000000000   a3000000000000000   a3000000000000000)
           zd(ii)=snapmatrix_volumefraction(13,u_nodesi,snapshotsi, numphasei)
          else if(xd(1,ii) .eq. a1 .and. xd(2,ii) .eq.a5  .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3) then !(1.0000000000000000   2.0000000000000000   a3000000000000000   a3000000000000000)
           zd(ii)=snapmatrix_volumefraction(14,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a5 .and. xd(2,ii) .eq.a5  .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3) then !(2.0000000000000000   2.0000000000000000   a3000000000000000   a3000000000000000)
           zd(ii)=snapmatrix_volumefraction(15,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a3 .and. int(xd(2,ii) *1000000)/1000000 .eq. a2 .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3) then 
            !(a3000000000000000   a24466094067263   a3000000000000000   a3000000000000000 xd(2,ii) >=1.14  .and. xd(2,ii) <=1.15)
           zd(ii)=snapmatrix_volumefraction(16,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a3 .and. int(xd(2,ii) *1000000)/1000000 .eq. a4 .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a3) then
          !(a3000000000000000   a45533905932737   a3000000000000000   a3000000000000000)
           zd(ii)=snapmatrix_volumefraction(17,u_nodesi,snapshotsi, numphasei)
           else if(xd(1,ii) .eq. a1 .and. xd(2,ii) .eq. a3 .and. xd(3,ii) .eq. a1 .and. xd(4,ii) .eq. a3) then !( 1.0000000000000000   a3000000000000000   1.0000000000000000  a3000000000000000)
           zd(ii)=snapmatrix_volumefraction(18,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq.a5.and. xd(2,ii) .eq. a3   .and. xd(3,ii) .eq. a1 .and. xd(4,ii) .eq. a3) then !(2.0000000000000000   a3000000000000000   1.0000000000000000   a3000000000000000)
           zd(ii)=snapmatrix_volumefraction(19,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a1 .and. xd(2,ii) .eq. a3   .and. xd(3,ii) .eq. a5 .and. xd(4,ii) .eq. a3) then !(1.0000000000000000   a3000000000000000   2.0000000000000000   a3000000000000000)
           zd(ii)=snapmatrix_volumefraction(20,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a5 .and. xd(2,ii) .eq. a3   .and. xd(3,ii) .eq. a5 .and. xd(4,ii) .eq. a3) then !(2.0000000000000000   a3000000000000000   2.0000000000000000   a3000000000000000)
           zd(ii)=snapmatrix_volumefraction(21,u_nodesi,snapshotsi, numphasei) 
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a1   .and. xd(3,ii) .eq. a1 .and. xd(4,ii) .eq. a3) then !(a3000000000000000   1.0000000000000000   1.0000000000000000   a3000000000000000)
           zd(ii)=snapmatrix_volumefraction(22,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a5   .and. xd(3,ii) .eq. a1 .and. xd(4,ii) .eq. a3) then !(a3000000000000000   2.0000000000000000   1.0000000000000000   a3000000000000000)
           zd(ii)=snapmatrix_volumefraction(23,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a1   .and. xd(3,ii) .eq. a5 .and. xd(4,ii) .eq. a3) then !(a3000000000000000   1.0000000000000000   2.0000000000000000   a3000000000000000)
           zd(ii)=snapmatrix_volumefraction(24,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a5   .and. xd(3,ii) .eq. a5 .and. xd(4,ii) .eq. a3) then !(a3000000000000000   2.0000000000000000   2.0000000000000000   a3000000000000000)
           zd(ii)=snapmatrix_volumefraction(25,u_nodesi,snapshotsi, numphasei)
          else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3 .and. int(xd(3,ii) *1000000)/1000000 .eq. a2 .and. xd(4,ii) .eq. a3) then 
       !(a3000000000000000   a3000000000000000  a24466094067263   a3000000000000000)
           zd(ii)=snapmatrix_volumefraction(26,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a3  .and. xd(2,ii) .eq. a3  .and. int(xd(3,ii) *1000000)/1000000 .eq. a4 .and. xd(4,ii) .eq. a3) then 
       !( a3000000000000000   a3000000000000000   a45533905932737   a3000000000000000 )
           zd(ii)=snapmatrix_volumefraction(27,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a1 .and. xd(2,ii) .eq. a3   .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a1) then !(1.0000000000000000   a3000000000000000   a3000000000000000   1.0000000000000000)
           zd(ii)=snapmatrix_volumefraction(28,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a5 .and. xd(2,ii) .eq. a3   .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a1) then !(2.0000000000000000   a3000000000000000   a3000000000000000   1.0000000000000000)
           zd(ii)=snapmatrix_volumefraction(29,u_nodesi,snapshotsi, numphasei)
           else if(xd(1,ii) .eq. a1 .and. xd(2,ii) .eq. a3   .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a5) then !(1.0000000000000000   a3000000000000000   a3000000000000000   2.0000000000000000)
           zd(ii)=snapmatrix_volumefraction(30,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a5 .and. xd(2,ii) .eq. a3   .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a5) then !(2.0000000000000000   a3000000000000000   a3000000000000000   2.0000000000000000)
           zd(ii)=snapmatrix_volumefraction(31,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a1   .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a1) then !(a3000000000000000   1.0000000000000000   a3000000000000000   1.0000000000000000)
           zd(ii)=snapmatrix_volumefraction(32,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a5  .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a1) then !(a3000000000000000   2.0000000000000000   a3000000000000000   1.0000000000000000)
           zd(ii)=snapmatrix_volumefraction(33,u_nodesi,snapshotsi, numphasei) 
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a1   .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a5) then !(a3000000000000000   1.0000000000000000   a3000000000000000   2.0000000000000000)
           zd(ii)=snapmatrix_volumefraction(34,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a5   .and. xd(3,ii) .eq. a3 .and. xd(4,ii) .eq. a5) then !(a3000000000000000   2.0000000000000000   a3000000000000000   2.0000000000000000)
           zd(ii)=snapmatrix_volumefraction(35,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3   .and. xd(3,ii) .eq. a1 .and. xd(4,ii) .eq. a1) then !(a3000000000000000   a3000000000000000   1.0000000000000000   1.0000000000000000)
           zd(ii)=snapmatrix_volumefraction(36,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3   .and. xd(3,ii) .eq. a5 .and. xd(4,ii) .eq. a1) then !( a3000000000000000   a3000000000000000   2.0000000000000000   1.0000000000000000)
           zd(ii)=snapmatrix_volumefraction(37,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3   .and. xd(3,ii) .eq. a1 .and. xd(4,ii) .eq. a5) then !(a3000000000000000   a3000000000000000   1.0000000000000000   2.0000000000000000)
           zd(ii)=snapmatrix_volumefraction(38,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3   .and. xd(3,ii) .eq. a5 .and. xd(4,ii) .eq. a5) then !(a3000000000000000   a3000000000000000   2.0000000000000000   2.0000000000000000)
           zd(ii)=snapmatrix_volumefraction(39,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a3 .and. xd(2,ii) .eq. a3 .and. xd(3,ii) .eq. a3 .and. int(xd(4,ii) *1000000)/1000000 .eq. a2) then 
      !(a3000000000000000   a3000000000000000   a3000000000000000   a24466094067263 )
           zd(ii)=snapmatrix_volumefraction(40,u_nodesi,snapshotsi, numphasei)
         else if(xd(1,ii) .eq. a3.and. xd(2,ii) .eq. a3 .and. xd(3,ii) .eq. a3.and. int(xd(4,ii) *1000000)/1000000 .eq. a4) then 
     !(a3000000000000000   a3000000000000000   a3000000000000000   1.8535533905932737     )
           zd(ii)=snapmatrix_volumefraction(41,u_nodesi,snapshotsi, numphasei)
         endif
        enddo
          !call f_sinr ( m, nd, xd, zd )       
           
          ! print *, 'zd',  zd(:)
            
!  Use the grid to evaluate the interpolant.
        call lagrange_interp_nd_value2 ( m, ind, a, b, nd, zd, ni, xi, zpi )
 
!  Weighted the interpolant values and add to the sparse grid interpolant.
 
        nd_total = nd_total + nd
        zi(1:ni) = zi(1:ni) + c(l) * zpi(1:ni)
        !  print *, ' nd_total,level',  nd_total, l, c(l), nd
        !  print *, 'zizizizizizizic(l)zpi', zi(1)!,c(l),zpi(1)
        deallocate ( xd )
        deallocate ( zd )     
       
       if ( .not. more ) then
          exit
        end if    
 
      end do !DO

    end do !do l = l_min, l_max   
    deallocate ( c )
    deallocate ( w )
    
   end do ! do l_max = sparse_max, sparse_max
   interpolate_volumefraction(u_nodesi,snapshotsi,numphasei)=zi(ni)
   enddo !Do u_nodesi  
  enddo !Do snapshotsi 
  enddo !numphasei

 !  print *, ' total', total , nd_total
 ! open(unit=66,file='podbasis_interp') 
 !  do j=1, total_timestep 
 !     write(66,*)  pod_coef_all_obv_41i(j,1:num_pod)
 ! enddo
 ! close(66)
 
  deallocate ( a )
  deallocate ( b )
  deallocate ( ind )
  deallocate ( xi )
  deallocate ( ze )
  deallocate ( zi )
  deallocate ( zpi )
  end subroutine interpolate_snapmatrix_volumefraction_4dim_level1


  function count_dumps(dump_sampling_period) result (count)
    !! Work out how many dumps we're going to read in.
    integer :: count,dump_sampling_period

    logical :: exists
    !      character(len=FILE_NAME_LEN) :: filename
    character(len=1024) :: filename, phase_name
    numphase=option_count("/material_phase")

    count=1
    call get_option("/material_phase["//int2str(0)//"]/name", phase_name)
    do 
       !! Note that this won't work in parallel. Have to look for the pvtu in that case.
       !! not ' but "
       if(numphase .eq. 1) then
          write(filename, '(a, i0, a,a)')  trim(simulation_name)//'_VelocityMesh_', (count)*dump_sampling_period,'_checkpoint',".vtu" 
       else
          write(filename, '(a, i0, a,a)')  trim(simulation_name)//'_'//trim(phase_name)//'_VelocityMesh_', (count)*dump_sampling_period,'_checkpoint',".vtu" 
          ! write(filename, '(a, i0, a)') trim(simulation_name)//'_', (count-1)*dump_sampling_period,".vtu" 
       endif
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
       FLExit("Sorry! There are no vtu files found in function count_dump!")
    end if

  end function count_dumps

 function count_vtu(dump_sampling_period) result (count)
    !! Work out how many dumps we're going to read in.
    integer :: count,dump_sampling_period

    logical :: exists 
    character(len=1024) :: filename, phase_name
    numphase=option_count("/material_phase")

    count=1
    call get_option("/material_phase["//int2str(0)//"]/name", phase_name)
    do 
       !! Note that this won't work in parallel. Have to look for the pvtu in that case.
       !! not ' but "
        
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
       FLExit("Sorry! There are no vtu files found in function count_dump!")
    end if

  end function count_vtu
 

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
   
subroutine comp_next ( n, k, a, more, h, t )

!*****************************************************************************
!
!! COMP_NEXT computes the compositions of the integer N into K parts.
!
!  Discussion:
!
!    A composition of the integer N into K parts is an ordered sequence
!    of K nonnegative integers which sum to N.  The compositions (1,2,1)
!    and (1,1,2) are considered to be distinct.
!
!    The routine computes one composition on each call until there are no more.
!    For instance, one composition of 6 into 3 parts is
!    3+2+1, another would be 6+0+0.
!
!    On the first call to this routine, set MORE = FALSE.  The routine
!    will compute the first element in the sequence of compositions, and
!    return it, as well as setting MORE = TRUE.  If more compositions
!    are desired, call again, and again.  Each time, the routine will
!    return with a new composition.
!
!    However, when the LAST composition in the sequence is computed
!    and returned, the routine will reset MORE to FALSE, signaling that
!    the end of the sequence has been reached.
!
!    This routine originally used a SAVE statement to maintain the
!    variables H and T.  I have decided that it is safer
!    to pass these variables as arguments, even though the user should
!    never alter them.  This allows this routine to safely shuffle
!    between several ongoing calculations.
!
!
!    There are 28 compositions of 6 into three parts.  This routine will
!    produce those compositions in the following order:
!
!     I         A
!     -     ---------
!     1     6   0   0
!     2     5   1   0
!     3     4   2   0
!     4     3   3   0
!     5     2   4   0
!     6     1   5   0
!     7     0   6   0
!     8     5   0   1
!     9     4   1   1
!    10     3   2   1
!    11     2   3   1
!    12     1   4   1
!    13     0   5   1
!    14     4   0   2
!    15     3   1   2
!    16     2   2   2
!    17     1   3   2
!    18     0   4   2
!    19     3   0   3
!    20     2   1   3
!    21     1   2   3
!    22     0   3   3
!    23     2   0   4
!    24     1   1   4
!    25     0   2   4
!    26     1   0   5
!    27     0   1   5
!    28     0   0   6
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 July 2008
!
!  Author:
!
!    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms for Computers and Calculators,
!    Second Edition,
!    Academic Press, 1978,
!    ISBN: 0-12-519260-6,
!    LC: QA164.N54.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the integer whose compositions are desired.
!
!    Input, integer ( kind = 4 ) K, the number of parts in the composition.
!
!    Input/output, integer ( kind = 4 ) A(K), the parts of the composition.
!
!    Input/output, logical MORE, set by the user to start the computation,
!    and by the routine to terminate it.
!
!    Input/output, integer ( kind = 4 ) H, T, two internal parameters needed
!    for the computation.  The user should allocate space for these in the
!    calling program, include them in the calling sequence, but never alter
!    them!
!
  implicit none

  integer ( kind = 4 ) k

  integer ( kind = 4 ) a(k)
  integer ( kind = 4 ) h
  logical more
  integer ( kind = 4 ) n
  integer ( kind = 4 ) t
!
!  The first computation.
!
  if ( .not. more ) then

    t = n
    h = 0
    a(1) = n
    a(2:k) = 0
!
!  The next computation.
!
  else
!
!  If the first entry A(1) is positive, then set H to zero,
!  so that when we increment H, it points to A(1); we will decrement A(1) by 1
!  and increment A(2).
!
    if ( 1 < t ) then
      h = 0
    end if
!
!  Otherwise, A(1) is 0.  Then by H + 1 is the entry we incremented last time.
!  Set H = H + 1, zero A(H), adding all but one of its value to A(1),
!  and incrementing A(H+1) by 1.
!
    h = h + 1
    t = a(h)
    a(h) = 0
    a(1) = t - 1
    a(h+1) = a(h+1) + 1

  end if
!
!  This is the last element of the sequence if all the
!  items are in the last slot.
!
  more = ( a(k) /= n )

  return
end

function i4_mop ( i ) result (i4_mop_r)

!*****************************************************************************80
!
!! I4_MOP returns the I-th power of -1 as an I4.
!
!  Discussion:
!
!    An I4 is an integer ( kind = 4 ) value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 November 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the power of -1.
!
!    Output, integer ( kind = 4 ) I4_MOP, the I-th power of -1.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_mop_r

  if ( mod ( i, 2 ) == 0 ) then
    i4_mop_r = 1
  else
    i4_mop_r = -1
  end if

 ! return
end function i4_mop

function i4_choose ( n, k ) ! return(i4_choose_r)

!*****************************************************************************
!
!! I4_CHOOSE computes the binomial coefficient C(N,K) as an I4.
!
!  Discussion:
!
!    The value is calculated in such a way as to avoid overflow and
!    roundoff.  The calculation is done in integer arithmetic.
!
!    The formula used is:
!
!      C(N,K) = N! / ( K! * (N-K)! )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    ML Wolfson, HV Wright,
!    Algorithm 160:
!    Combinatorial of M Things Taken N at a Time,
!    Communications of the ACM,
!    Volume 6, Number 4, April 1963, page 161.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, K, are the values of N and K.
!
!    Output, integer ( kind = 4 ) I4_CHOOSE, the number of combinations of N
!    things taken K at a time.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_choose
  integer ( kind = 4 ) k
  integer ( kind = 4 ) mn
  integer ( kind = 4 ) mx
  integer ( kind = 4 ) n
  integer ( kind = 4 ) value

  mn = min ( k, n - k )

  if ( mn < 0 ) then

    value = 0

  else if ( mn == 0 ) then

    value = 1

  else

    mx = max ( k, n - k )
    value = mx + 1

    do i = 2, mn
      value = ( value * ( mx + i ) ) / i
    end do

  end if

  i4_choose = value

  return
end !function i4_choose

subroutine smolyak_coefficients ( l_max, m, c, w )

!*****************************************************************************80
!
!! SMOLYAK_COEFFICIENTS returns the Smolyak coefficients and counts.
!
!  Discussion:
!
!    The Smolyak sparse interpolant can be written as:
!
!      A(L,M)(X) = sum ( L-M+1 <= |L| <= L_max ) 
!        C(|L|) * g(l1)(x1) * g(l2)(x2) * ... * g(lm)(xm).
!
!    where:
!
!    * L=(l1,l2,...,lm) is a vector of M nonnegative integers;
!    * |L| is the sum of the entries of L;
!    * X=(x1,x2,...,xm) is an M-dimensional point in a product space;
!    * g(i)(xj) is the i-th 1-d interpolation function in dimension j;
!
!    Note that:
!
!    * W(|L|) will represent the number of distinct interpolants for which
!      the sublevel, or sum of the L vector entries, is |L|;
!
!    * the coefficients C and counts W will be zero for sublevels 
!      0 through L_MAX - M (and MATLAB indices 1 through L_MAX-M+1).
!
!    * it will be the case that W' * C = 1, essentially because the interpolant
!      to the identity function must be the identity function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 July 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) L_MAX, the (maximum) level.
!    0 <= L_MAX.
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!    1 <= M.
!
!    Output, integer ( kind = 4 ) C(0:L_MAX), the coefficients for objects 
!    at sublevels 0 through L_MAX.
!
!    Output, integer ( kind = 4 ) W(0:L_MAX), the number of objects at 
!    sublevels 0 through L_MAX.
!
  implicit none

  integer ( kind = 4 ) l_max

  integer ( kind = 4 ) c(0:l_max)
  !integer ( kind = 4 ) i4_choose
  !integer ( kind = 4 ) i4_mop
  integer ( kind = 4 ) l
  integer ( kind = 4 ) l_min
  integer ( kind = 4 ) m
  integer ( kind = 4 ) w(0:l_max)

  l_min = max ( l_max - m + 1, 0 )

  c(0:l_min-1) = 0
  do l = l_min, l_max
    c(l) = i4_mop ( l_max - l ) * i4_choose ( m - 1, l_max - l )
  end do

  w(0:l_min-1) = 0
  do l = l_min, l_max
    w(l) = i4_choose ( l + m - 1, m - 1 )
  end do

  return
end

subroutine lagrange_interp_nd_size2 ( m, ind, nd )

!*****************************************************************************80
!
!! LAGRANGE_INTERP_ND_SIZE2 sizes an M-dimensional Lagrange interpolant.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 September 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) IND(M), the index or level of the 1D rule 
!    to be used in each dimension.
!
!    Output, integer ( kind = 4 ) ND, the number of points in the product grid.
!
  implicit none

  integer ( kind = 4 ) m

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ind(m)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nd
!
!  Determine the number of data points.
!
  nd = 1
  do i = 1, m
    call order_from_level_135 ( ind(i), n )
    nd = nd * n
  end do

  return
end

subroutine lagrange_interp_nd_grid2 ( m, ind, a, b, nd, xd )

!*****************************************************************************80
!
!! LAGRANGE_INTERP_ND_GRID2 sets an M-dimensional Lagrange interpolant grid.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 September 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) IND(M), the index or level of the 1D rule 
!    to be used in each dimension.
!
!    Input, real ( kind = 8 ) A(M), B(M), the lower and upper limits.
!
!    Input, integer ( kind = 4 ) ND, the number of points in the product grid.
!
!    Output, real ( kind = 8 ) XD(M,ND), the points at which data was sampled.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) nd

  real ( kind = 8 ) a(m)
  real ( kind = 8 ) b(m)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ind(m)
  integer ( kind = 4 ) n
  real ( kind = 8 ), allocatable :: x_1d(:)
  real ( kind = 8 ) xd(m,nd)
!
!  Compute the data points.
!
  xd(1:m,1:nd) = 0.0D+00
  do i = 1, m
    call order_from_level_135 ( ind(i), n ) !get n,for each dimension
    allocate ( x_1d(1:n) )
    call cc_compute_points ( n, x_1d ) !get coordinates save it to x_ld
    x_1d(1:n) = 0.5D+00 * ( ( 1.0D+00 - x_1d(1:n) ) * a(i) &
                          + ( 1.0D+00 + x_1d(1:n) ) * b(i) )
    call r8vec_direct_product5 ( i, n, x_1d, m, nd, xd )
    deallocate ( x_1d )
  end do

  return
end

subroutine lagrange_interp_nd_value2 ( m, ind, a, b, nd, zd, ni, xi, zi )

!*****************************************************************************80
!
!! LAGRANGE_INTERP_ND_VALUE2 evaluates an ND Lagrange interpolant.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 September 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) IND(M), the index or level of the 1D rule 
!    to be used in each dimension.
!
!    Input, real ( kind = 8 ) A(M), B(M), the lower and upper limits.
!
!    Input, integer ( kind = 4 ) ND, the number of points in the product grid.
!
!    Input, real ( kind = 8 ) ZD(ND), the function evaluated at the points XD.
!
!    Input, integer ( kind = 4 ) NI, the number of points at which the 
!    interpolant is to be evaluated.
!
!    Input, real ( kind = 8 ) XI(M,NI), the points at which the interpolant 
!    is to be evaluated.
!
!    Output, real ( kind = 8 ) ZI(NI), the interpolant evaluated at the 
!    points XI.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) nd
  integer ( kind = 4 ) ni

  real ( kind = 8 ) a(m)
  real ( kind = 8 ) b(m)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ind(m)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n
  real ( kind = 8 ), allocatable :: value(:)
  real ( kind = 8 ) w(nd)
  real ( kind = 8 ), allocatable :: x_1d(:)
  real ( kind = 8 ) xi(m,ni)
  real ( kind = 8 ) zd(nd)
  real ( kind = 8 ) zi(ni)

  do j = 1, ni

    w(1:nd) = 1.0D+00

    do i = 1, m
      call order_from_level_135 ( ind(i), n )
      allocate ( x_1d(1:n) )
      allocate ( value(1:n) )
      call cc_compute_points ( n, x_1d )
      x_1d(1:n) = 0.5D+00 * ( ( 1.0D+00 - x_1d(1:n) ) * a(i) &
                            + ( 1.0D+00 + x_1d(1:n) ) * b(i) )
      call lagrange_basis_1d ( n, x_1d, 1, xi(i,j), value )
      call r8vec_direct_product3 ( i, n, value, m, nd, w )
      deallocate ( value )
      deallocate ( x_1d )
    end do

    zi(j) = dot_product ( w, zd )

  end do

  return
end

subroutine order_from_level_135 ( l, n )

!*****************************************************************************80
!
!! ORDER_FROM_LEVEL_135 evaluates the 135 level-to-order relationship.
!
!  Discussion:
!
!    Clenshaw Curtis rules, and some others, often use the following
!    scheme:
!
!    L: 0  1  2  3   4   5
!    N: 1  3  5  9  17  33 ... 2^L+1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 July 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) L, the level, which should be 0 or greater.
!
!    Output, integer ( kind = 4 ) N, the order.
!
  implicit none

  integer ( kind = 4 ) l
  integer ( kind = 4 ) n

  if ( l < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ORDER_FROM_LEVEL_135 - Fatal error!'
    write ( *, '(a)' ) '  Illegal input value of L!'
    stop
  else if ( l == 0 ) then
    n = 1
  else
    n = ( 2 ** l ) + 1
  end if

  return
end





subroutine cc_compute_points ( n, points )

!*****************************************************************************
!
!! CC_COMPUTE_POINTS: abscissas of a Clenshaw Curtis rule.
!
!  Discussion:
!
!    Our convention is that the abscissas are numbered from left to right.
!
!    The rule is defined on [-1,1].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!    1 <= N.
!
!    Output, real ( kind = 8 ) POINTS(N), the abscissas.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) points(n)

  if ( n < 1 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CC_COMPUTE_POINTS - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of N = ', n
    stop

  else if ( n == 1 ) then

    points(1) = 0.0D+00

  else

    do i = 1, n
      points(i) = cos ( real ( n - i, kind = 8 ) * pi &
                      / real ( n - 1, kind = 8 ) )
    end do

    points(1) = -1.0D+00
    if ( mod ( n, 2 ) == 1 ) then
      points((n+1)/2) = 0.0D+00
    end if
    points(n) = +1.0D+00

  end if

  return
end






subroutine lagrange_basis_1d ( nd, xd, ni, xi, lb ) 

!*****************************************************************************80
!
!! LAGRANGE_BASIS_1D evaluates a 1D Lagrange basis.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 October 2012
!
!  Author:
!
!    John Burkardt

!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ND, the number of data points.
!
!    Input, real ( kind = 8 ) XD(ND), the interpolation nodes.
!
!    Input, integer ( kind = 4 ) NI, the number of evaluation points.
!
!    Input, real ( kind = 8 ) XI(NI), the evaluation points.
!
!    Output, real ( kind = 8 ) LB(NI,ND), the value, at the I-th point XI, 
!    of the Jth basis function.
!
  implicit none

  integer ( kind = 4 ) nd
  integer ( kind = 4 ) ni

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) lb(ni,nd)
  real ( kind = 8 ) xd(nd)
  real ( kind = 8 ) xi(ni)
  
  do i = 1, ni
    do j = 1, nd
      lb(i,j) = product ( ( xi(i) - xd(1:j-1)  ) / ( xd(j) - xd(1:j-1)  ) ) &
              * product ( ( xi(i) - xd(j+1:nd) ) / ( xd(j) - xd(j+1:nd) ) )
    end do
  end do

  return
end






subroutine r8vec_direct_product3 ( factor_index, factor_order, factor_value, &
  factor_num, point_num, w )

!*****************************************************************************80
!
!! R8VEC_DIRECT_PRODUCT2 creates a direct product of R8VEC's.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    To explain what is going on here, suppose we had to construct
!    a multidimensional quadrature rule as the product of K rules
!    for 1D quadrature.
!
!    The product rule will be represented as a list of points and weights.
!
!    The J-th item in the product rule will be associated with
!      item J1 of 1D rule 1,
!      item J2 of 1D rule 2,
!      ...,
!      item JK of 1D rule K.
!
!    In particular,
!      X(J) = ( X(1,J1), X(2,J2), ..., X(K,JK))
!    and
!      W(J) = W(1,J1) * W(2,J2) * ... * W(K,JK)
!
!    So we can construct the quadrature rule if we can properly
!    distribute the information in the 1D quadrature rules.
!
!    This routine carries out the task involving the weights W.
!
!    Another way to do this would be to compute, one by one, the
!    set of all possible indices (J1,J2,...,JK), and then index
!    the appropriate information.  An advantage of the method shown
!    here is that you can process the K-th set of information and
!    then discard it.
!
!  Example:
!
!    Rule 1:
!      Order = 4
!      W(1:4) = ( 2, 3, 5, 7 )
!
!    Rule 2:
!      Order = 3
!      W(1:3) = ( 11, 13, 17 )
!
!    Rule 3:
!      Order = 2
!      W(1:2) = ( 19, 23 )
!
!    Product Rule:
!      Order = 24
!      W(1:24) =
!        ( 2 * 11 * 19 )
!        ( 3 * 11 * 19 )
!        ( 4 * 11 * 19 )
!        ( 7 * 11 * 19 )
!        ( 2 * 13 * 19 )
!        ( 3 * 13 * 19 )
!        ( 5 * 13 * 19 )
!        ( 7 * 13 * 19 )
!        ( 2 * 17 * 19 )
!        ( 3 * 17 * 19 )
!        ( 5 * 17 * 19 )
!        ( 7 * 17 * 19 )
!        ( 2 * 11 * 23 )
!        ( 3 * 11 * 23 )
!        ( 5 * 11 * 23 )
!        ( 7 * 11 * 23 )
!        ( 2 * 13 * 23 )
!        ( 3 * 13 * 23 )
!        ( 5 * 13 * 23 )
!        ( 7 * 13 * 23 )
!        ( 2 * 17 * 23 )
!        ( 3 * 17 * 23 )
!        ( 5 * 17 * 23 )
!        ( 7 * 17 * 23 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 April 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) FACTOR_INDEX, the index of the factor being
!    processed.  The first factor processed must be factor 1!
!
!    Input, integer ( kind = 4 ) FACTOR_ORDER, the order of the factor.
!
!    Input, real ( kind = 8 ) FACTOR_VALUE(FACTOR_ORDER), the factor values
!    for factor FACTOR_INDEX.
!
!    Input, integer ( kind = 4 ) FACTOR_NUM, the number of factors.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of elements in the
!    direct product.
!
!    Input/output, real ( kind = 8 ) W(POINT_NUM), the elements of the
!    direct product, which are built up gradually.
!
!  Local Parameters:
!
!    Local, integer ( kind = 4 ) START, the first location of a block of values
!    to set.
!
!    Local, integer ( kind = 4 ) CONTIG, the number of consecutive values
!    to set.
!
!    Local, integer ( kind = 4 ) SKIP, the distance from the current value 
!    of START to the next location of a block of values to set.
!
!    Local, integer ( kind = 4 ) REP, the number of blocks of values to set.
!
  implicit none

  integer ( kind = 4 ) factor_num
  integer ( kind = 4 ) factor_order
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ), save :: contig
  integer ( kind = 4 ) factor_index
  real ( kind = 8 ) factor_value(factor_order)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ), save :: rep
  integer ( kind = 4 ), save :: skip
  integer ( kind = 4 ) start
  real ( kind = 8 ) w(point_num)

  if ( factor_index == 1 ) then
    contig = 1
    skip = 1
    rep = point_num
    w(1:point_num) = 1.0D+00
  end if

  rep = rep / factor_order
  skip = skip * factor_order

  do j = 1, factor_order

    start = 1 + ( j - 1 ) * contig

    do k = 1, rep
      w(start:start+contig-1) = w(start:start+contig-1) * factor_value(j)
      start = start + skip
    end do

  end do

  contig = contig * factor_order

  return
end
 

subroutine r8vec_direct_product5 ( factor_index, factor_order, factor_value, &
  factor_num, point_num, x )

!*****************************************************************************80
!
!! R8VEC_DIRECT_PRODUCT creates a direct product of R8VEC's.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    To explain what is going on here, suppose we had to construct
!    a multidimensional quadrature rule as the product of K rules
!    for 1D quadrature.
!
!    The product rule will be represented as a list of points and weights.
!
!    The J-th item in the product rule will be associated with
!      item J1 of 1D rule 1,
!      item J2 of 1D rule 2,
!      ...,
!      item JK of 1D rule K.
!
!    In particular,
!      X(J) = ( X(1,J1), X(2,J2), ..., X(K,JK))
!    and
!      W(J) = W(1,J1) * W(2,J2) * ... * W(K,JK)
!
!    So we can construct the quadrature rule if we can properly
!    distribute the information in the 1D quadrature rules.
!
!    This routine carries out that task for the abscissas X.
!
!    Another way to do this would be to compute, one by one, the
!    set of all possible indices (J1,J2,...,JK), and then index
!    the appropriate information.  An advantage of the method shown
!    here is that you can process the K-th set of information and
!    then discard it.
!
!  Example:
!
!    Rule 1:
!      Order = 4
!      X(1:4) = ( 1, 2, 3, 4 )
!
!    Rule 2:
!      Order = 3
!      X(1:3) = ( 10, 20, 30 )
!
!    Rule 3:
!      Order = 2
!      X(1:2) = ( 100, 200 )
!
!    Product Rule:
!      Order = 24
!      X(1:24) =
!        ( 1, 10, 100 )
!        ( 2, 10, 100 )
!        ( 3, 10, 100 )
!        ( 4, 10, 100 )
!        ( 1, 20, 100 )
!        ( 2, 20, 100 )
!        ( 3, 20, 100 )
!        ( 4, 20, 100 )
!        ( 1, 30, 100 )
!        ( 2, 30, 100 )
!        ( 3, 30, 100 )
!        ( 4, 30, 100 )
!        ( 1, 10, 200 )
!        ( 2, 10, 200 )
!        ( 3, 10, 200 )
!        ( 4, 10, 200 )
!        ( 1, 20, 200 )
!        ( 2, 20, 200 )
!        ( 3, 20, 200 )
!        ( 4, 20, 200 )
!        ( 1, 30, 200 )
!        ( 2, 30, 200 )
!        ( 3, 30, 200 )
!        ( 4, 30, 200 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 April 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) FACTOR_INDEX, the index of the factor being
!    processed.  The first factor processed must be factor 1!
!
!    Input, integer ( kind = 4 ) FACTOR_ORDER, the order of the factor.
!
!    Input, real ( kind = 8 ) FACTOR_VALUE(FACTOR_ORDER), the factor values
!    for factor FACTOR_INDEX.
!
!    Input, integer ( kind = 4 ) FACTOR_NUM, the number of factors.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of elements in the
!    direct product.
!
!    Input/output, real ( kind = 8 ) X(FACTOR_NUM,POINT_NUM), the elements of
!    the direct product, which are built up gradually.
!
!  Local Parameters:
!
!    Local, integer ( kind = 4 ) START, the first location of a block of 
!    values to set.
!
!    Local, integer ( kind = 4 ) CONTIG, the number of consecutive values 
!    to set.
!
!    Local, integer ( kind = 4 ) SKIP, the distance from the current value 
!    of START to the next location of a block of values to set.
!
!    Local, integer ( kind = 4 ) REP, the number of blocks of values to set.
!
  implicit none

  integer ( kind = 4 ) factor_num
  integer ( kind = 4 ) factor_order
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ), save :: contig
  integer ( kind = 4 ) factor_index
  real ( kind = 8 ) factor_value(factor_order)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ), save :: rep
  integer ( kind = 4 ), save :: skip
  integer ( kind = 4 ) start
  real ( kind = 8 ) x(factor_num,point_num)

  if ( factor_index == 1 ) then
    contig = 1
    skip = 1
    rep = point_num
    x(1:factor_num,1:point_num) = 0.0D+00
  end if

  rep = rep / factor_order
  skip = skip * factor_order

  do j = 1, factor_order

    start = 1 + ( j - 1 ) * contig

    do k = 1, rep
      x(factor_index,start:start+contig-1) = factor_value(j)
      start = start + skip
    end do

  end do

  contig = contig * factor_order

  return
end

subroutine rbf_weight ( m, nd, xd, r0, phi, fd, w )

!*****************************************************************************80
!
!! RBF_WEIGHT computes weights for radial basis function interpolation.
!
!  Discussion:
!
!    We assume that there are N (nonsingular) equations in N unknowns.
!
!    However, it should be clear that, if we are willing to do some kind
!    of least squares calculation, we could allow for singularity,
!    inconsistency, or underdetermine systems.  This could be associated
!    with data points that are very close or repeated, a smaller number
!    of data points than function values, or some other ill-conditioning
!    of the system arising from a peculiarity in the point spacing.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!    
!
!    30 June 2012
!
!      
!
!     
!
!  Reference:
!
!    William Press, Brian Flannery, Saul Teukolsky, William Vetterling,
!    Numerical Recipes in FORTRAN: The Art of Scientific Computing,
!    Third Edition,
!    Cambridge University Press, 2007,
!    ISBN13: 978-0-521-88068-8,
!    LC: QA297.N866.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) ND, the number of data points.
!
!    Input, real ( kind = 8 ) XD(M,ND), the data points.
!
!    Input, real ( kind = 8 ) R0, a scale factor.  R0 should be larger than 
!    the typical separation between points, but smaller than the maximum 
!    separation.  The value of R0 has a significant effect on the resulting 
!    interpolant.
!
!    Input, subroutine PHI ( N, R, R0, V ), a subroutine to evaluate the radial
!    basis functions.
!
!    Input, real ( kind = 8 ) FD(ND), the function values at the data points.
!
!    Output, real ( kind = 8 ) W(ND), the weights.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) nd

  real ( kind = 8 ) a(nd,nd)
  real ( kind = 8 ) fd(nd)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  external phi
  real ( kind = 8 ) r(nd)
  real ( kind = 8 ) r0
  real ( kind = 8 ) v(nd)
  real ( kind = 8 ) w(nd)
  real ( kind = 8 ) xd(m,nd)

  do i = 1, nd

    do j = 1, nd
      r(j) = sqrt ( sum ( ( xd(1:m,i) - xd(1:m,j) )**2 ) )
    end do

    call phi ( nd, r, r0, v )

    a(i,1:nd) = v(1:nd)

  end do
!
!  Solve for the weights.
!
  call r8mat_solve_svd ( nd, nd, a, fd, w )

  return
end
subroutine phi1 ( n, r, r0, v )

!*****************************************************************************80
!
!! PHI1 evaluates the multiquadric radial basis function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!    
!
!    30 June 2012
!
!      
!
!     
!
!  Reference:
!
!    William Press, Brian Flannery, Saul Teukolsky, William Vetterling,
!    Numerical Recipes in FORTRAN: The Art of Scientific Computing,
!    Third Edition,
!    Cambridge University Press, 2007,
!    ISBN13: 978-0-521-88068-8,
!    LC: QA297.N866.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) R(N), the radial separation.
!    0 < R.
!
!    Input, real ( kind = 8 ) R0, a scale factor.
!
!    Output, real ( kind = 8 ) V(N), the value of the radial basis function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) r(n)
  real ( kind = 8 ) r0
  real ( kind = 8 ) v(n)

  v(1:n) = sqrt ( r(1:n)**2 + r0**2 )

  return
end
subroutine phi2 ( n, r, r0, v )

!*****************************************************************************80
!
!! PHI2 evaluates the inverse multiquadric radial basis function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!    
!
!    30 June 2012
!
!      
!
!     
!
!  Reference:
!
!    William Press, Brian Flannery, Saul Teukolsky, William Vetterling,
!    Numerical Recipes in FORTRAN: The Art of Scientific Computing,
!    Third Edition,
!    Cambridge University Press, 2007,
!    ISBN13: 978-0-521-88068-8,
!    LC: QA297.N866.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) R(N), the radial separation.
!    0 < R.
!
!    Input, real ( kind = 8 ) R0, a scale factor.
!
!    Output, real ( kind = 8 ) V(N), the value of the radial basis function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) r(n)
  real ( kind = 8 ) r0
  real ( kind = 8 ) v(n)

  v = 1.0D+00 / sqrt ( r**2 + r0**2 )

  return
end
subroutine phi3 ( n, r, r0, v )

!*****************************************************************************80
!
!! PHI3 evaluates the thin-plate spline radial basis function.
!
!  Discussion:
!
!    Note that PHI3(R,R0) is negative if R < R0.  Thus, for this basis function,
!    it may be desirable to choose a value of R0 smaller than any possible R.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!    
!
!    30 June 2012
!
!      
!
!     
!
!  Reference:
!
!    William Press, Brian Flannery, Saul Teukolsky, William Vetterling,
!    Numerical Recipes in FORTRAN: The Art of Scientific Computing,
!    Third Edition,
!    Cambridge University Press, 2007,
!    ISBN13: 978-0-521-88068-8,
!    LC: QA297.N866.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) R(N), the radial separation.
!    0 < R.
!
!    Input, real ( kind = 8 ) R0, a scale factor.
!
!    Output, real ( kind = 8 ) V(N), the value of the radial basis function.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  real ( kind = 8 ) r(n)
  real ( kind = 8 ) r0
  real ( kind = 8 ) v(n)

  do i = 1, n
    if ( r(i) .le. 0.0D+00 ) then
      v(i) = 0.0D+00
    else
      v(i) = r(i)**2 * log ( r(i) / r0 )
    end if
  end do

  return
end
subroutine phi4 ( n, r, r0, v )

!*****************************************************************************80
!
!! PHI4 evaluates the gaussian radial basis function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!    
!
!    30 June 2012
!
!      
!
!     
!
!  Reference:
!
!    William Press, Brian Flannery, Saul Teukolsky, William Vetterling,
!    Numerical Recipes in FORTRAN: The Art of Scientific Computing,
!    Third Edition,
!    Cambridge University Press, 2007,
!    ISBN13: 978-0-521-88068-8,
!    LC: QA297.N866.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) R(N), the radial separation.
!    0 < R.
!
!    Input, real ( kind = 8 ) R0, a scale factor.
!
!    Output, real ( kind = 8 ) V(N), the value of the radial basis function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) r(n)
  real ( kind = 8 ) r0
  real ( kind = 8 ) v(n)

  v(1:n) = exp ( - 0.5D+00 * r(1:n)**2 / r0**2 )

  return
end
subroutine r8mat_solve_svd ( m, n, a, b, x )

!*****************************************************************************80
!
!! R8MAT_SOLVE_SVD solves a linear system A*x=b using the SVD.
!
!  Discussion:
!
!    When the system is determined, the solution is the solution in the
!    ordinary sense, and A*x = b.
!
!    When the system is overdetermined, the solution minimizes the
!    L2 norm of the residual ||A*x-b||.
!
!    When the system is underdetermined, ||A*x-b|| should be zero, and
!    the solution is the solution of minimum L2 norm, ||x||.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!    30 June 2012 
!    Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns
!    in the matrix A.
!
!    Input, real ( kind = 8 ) A(M,N), the matrix.
!
!    Input, real ( kind = 8 ) B(M), the right hand side.
!
!    Output, real ( kind = 8 ) X(N), the solution.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) a_copy(m,n)
  real ( kind = 8 ) a_pseudo(n,m)
  real ( kind = 8 ) b(m)
  real ( kind = 8 ) e(max(m+1,n))
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) lda
  integer ( kind = 4 ) ldu
  integer ( kind = 4 ) ldv
  integer ( kind = 4 ) job
  integer ( kind = 4 ) lwork
  real ( kind = 8 ) s(m,n)
  real ( kind = 8 ) sp(n,m)
  real ( kind = 8 ) sdiag(max(m+1,n))
  real ( kind = 8 ) u(m,m)
  real ( kind = 8 ) v(n,n)
  real ( kind = 8 ) work(m)
  real ( kind = 8 ) x(n)
!
!  Compute the SVD decomposition.
!
  job = 11
  lda = m
  ldu = m
  ldv = n

  a_copy(1:m,1:n) = a(1:m,1:n)

  call dsvdc ( a_copy, lda, m, n, sdiag, e, u, ldu, v, ldv, work, job, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8MAT_SOLVE_SVD - Failure!'
    write ( *, '(a)' ) '  The SVD could not be calculated.'
    write ( *, '(a)' ) '  LINPACK routine DSVDC returned a nonzero'
    write ( *, '(a,i8)' ) '  value of the error flag, INFO = ', info
    return
  end if

  s(1:m,1:n) = 0.0D+00
  do i = 1, min ( m, n )
    s(i,i) = sdiag(i)
  end do
!
!  Compute the pseudo inverse.
!
  sp(1:n,1:m) = 0.0D+00
  do i = 1, min ( m, n )
    if ( s(i,i) /= 0.0D+00 ) then
      sp(i,i) = 1.0D+00 / s(i,i)
    end if
  end do

  a_pseudo(1:n,1:m) = matmul ( v(1:n,1:n), &
    matmul ( sp(1:n,1:m), transpose ( u(1:m,1:m) ) ) )
!
!  Compute x = A_pseudo * b
!
  x(1:n) = matmul ( a_pseudo(1:n,1:m), b(1:m) )

  return
end
subroutine rbf_interp_nd ( m, nd, xd, r0, phi, w, ni, xi, fi )

!*****************************************************************************80
!
!! RBF_INTERP_ND evaluates a radial basis function interpolant.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!    
!
!    30 June 2012
!  Reference:
!
!    William Press, Brian Flannery, Saul Teukolsky, William Vetterling,
!    Numerical Recipes in FORTRAN: The Art of Scientific Computing,
!    Third Edition,
!    Cambridge University Press, 2007,
!    ISBN13: 978-0-521-88068-8,
!    LC: QA297.N866.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) ND, the number of data points.
!
!    Input, real ( kind = 8 ) XD(M,ND), the data points.
!
!    Input, real ( kind = 8 ) R0, a scale factor.  R0 should be larger than 
!    the typical separation between points, but smaller than the maximum 
!    separation.  The value of R0 has a significant effect on the resulting 
!    interpolant.
!
!    Input, subroutine PHI ( N, R, R0, V ), a subroutine to evaluate the radial
!    basis functions.
!
!    Input, real ( kind = 8 ) W(ND), the weights, as computed by RBF_WEIGHTS.
!
!    Input, integer ( kind = 4 ) NI, the number of interpolation points.
!
!    Input, real ( kind = 8 ) XI(M,NI), the interpolation points.
!
!    Output, real ( kind = 8 ) FI(NI), the interpolated values.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) nd
  integer ( kind = 4 ) ni

  real ( kind = 8 ) fi(ni)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  external phi
  real ( kind = 8 ) r(nd)
  real ( kind = 8 ) r0
  real ( kind = 8 ) v(nd)
  real ( kind = 8 ) w(nd)
  real ( kind = 8 ) xd(m,nd)
  real ( kind = 8 ) xi(m,ni)

  do i = 1, ni

    do j = 1, nd
      r(j) = sqrt ( sum ( ( xi(1:m,i) - xd(1:m,j) )**2 ) )
    end do

    call phi ( nd, r, r0, v )

    fi(i) = dot_product ( v, w )

  end do

  return
end
subroutine r8mat_transpose_print ( m, n, a, title )

!*****************************************************************************80
!
!! R8MAT_TRANSPOSE_PRINT prints an R8MAT, transposed.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!    
!
!    14 June 2004
!
!      
!
!     
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, real ( kind = 8 ) A(M,N), an M by N matrix to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  character ( len = * ) title

  call r8mat_transpose_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine r8mat_transpose_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8MAT_TRANSPOSE_PRINT_SOME prints some of an R8MAT, transposed.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!    
!
!    10 September 2009
!
!      
!
!     
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, real ( kind = 8 ) A(M,N), an M by N matrix to be printed.
!
!    Input, integer ( kind = 4 ) ILO, JLO, the first row and column to print.
!
!    Input, integer ( kind = 4 ) IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 5
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  character ( len = 14 ) ctemp(incx)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  if ( m <= 0 .or. n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  (None)'
    return
  end if

  do i2lo = max ( ilo, 1 ), min ( ihi, m ), incx

    i2hi = i2lo + incx - 1
    i2hi = min ( i2hi, m )
    i2hi = min ( i2hi, ihi )

    inc = i2hi + 1 - i2lo

    write ( *, '(a)' ) ' '

    do i = i2lo, i2hi
      i2 = i + 1 - i2lo
      write ( ctemp(i2), '(i8,6x)' ) i
    end do

    write ( *, '(''  Row   '',5a14)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Col'
    write ( *, '(a)' ) ' '

    j2lo = max ( jlo, 1 )
    j2hi = min ( jhi, n )

    do j = j2lo, j2hi

      do i2 = 1, inc
        i = i2lo - 1 + i2
        write ( ctemp(i2), '(g14.6)' ) a(i,j)
      end do

      write ( *, '(i5,a,5a14)' ) j, ':', ( ctemp(i), i = 1, inc )

    end do

  end do

  return
end
subroutine dsvdc ( a, lda, m, n, s, e, u, ldu, v, ldv, work, job, info )

!*****************************************************************************80
!
!! DSVDC computes the singular value decomposition of a real rectangular matrix.
!
!  Discussion:
!
!    This routine reduces an M by N matrix A to diagonal form by orthogonal
!    transformations U and V.  The diagonal elements S(I) are the singular
!    values of A.  The columns of U are the corresponding left singular
!    vectors, and the columns of V the right singular vectors.
!
!    The form of the singular value decomposition is then
!
!      A(MxN) = U(MxM) * S(MxN) * V(NxN)'
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!    
!
!    16 September 2006
!
!      
!
!    Original FORTRAN77 version by Jack Dongarra, Jim Bunch, Cleve Moler, 
!    Pete Stewart.
!    FORTRAN90 version by  .
!
!  Reference:
!
!    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979,
!    ISBN13: 978-0-898711-72-1,
!    LC: QA214.L56.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) A(LDA,N).  On input, the M by N
!    matrix whose singular value decomposition is to be computed.
!    On output, the matrix has been destroyed.  Depending on the user's
!    requests, the matrix may contain other useful information.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of the array A.
!    LDA must be at least N.
!
!    Input, integer ( kind = 4 ) M, the number of rows of the matrix.
!
!    Input, integer ( kind = 4 ) N, the number of columns of the matrix A.
!
!    Output, real ( kind = 8 ) S(MM), where MM = max(M+1,N).  The first
!    min(M,N) entries of S contain the singular values of A arranged in
!    descending order of magnitude.
!
!    Output, real ( kind = 8 ) E(MM), where MM = max(M+1,N).  Ordinarily
!    contains zeros.  However see the discussion of INFO for exceptions.
!
!    Output, real ( kind = 8 ) U(LDU,K).  If JOBA = 1 then K = M;
!    if 2 <= JOBA, then K = min(M,N).  U contains the M by M matrix of
!    left singular vectors.  U is not referenced if JOBA = 0.  If M <= N
!    or if JOBA = 2, then U may be identified with A in the subroutine call.
!
!    Input, integer ( kind = 4 ) LDU, the leading dimension of the array U.
!    LDU must be at least M.
!
!    Output, real ( kind = 8 ) V(LDV,N), the N by N matrix of right singular
!    vectors.  V is not referenced if JOB is 0.  If N <= M, then V may be
!    identified with A in the subroutine call.
!
!    Input, integer ( kind = 4 ) LDV, the leading dimension of the array V.
!    LDV must be at least N.
!
!    Workspace, real ( kind = 8 ) WORK(M).
!
!    Input, integer ( kind = 4 ) JOB, controls the computation of the singular
!    vectors.  It has the decimal expansion AB with the following meaning:
!      A =  0, do not compute the left singular vectors.
!      A =  1, return the M left singular vectors in U.
!      A >= 2, return the first min(M,N) singular vectors in U.
!      B =  0, do not compute the right singular vectors.
!      B =  1, return the right singular vectors in V.
!
!    Output, integer ( kind = 4 ) INFO, status indicator.
!    The singular values (and their corresponding singular vectors)
!    S(INFO+1), S(INFO+2),...,S(MN) are correct.  Here MN = min ( M, N ).
!    Thus if INFO is 0, all the singular values and their vectors are
!    correct.  In any event, the matrix B = U' * A * V is the bidiagonal
!    matrix with the elements of S on its diagonal and the elements of E on
!    its superdiagonal.  Thus the singular values of A and B are the same.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) ldu
  integer ( kind = 4 ) ldv
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(lda,n)
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) cs
  real ( kind = 8 ) e(*)
  real ( kind = 8 ) el
  real ( kind = 8 ) emm1
  real ( kind = 8 ) f
  real ( kind = 8 ) g
  integer ( kind = 4 ) info
  integer ( kind = 4 ) iter
  integer ( kind = 4 ) j
  integer ( kind = 4 ) job
  integer ( kind = 4 ) jobu
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kase
  integer ( kind = 4 ) kk
  integer ( kind = 4 ) l
  integer ( kind = 4 ) ll
  integer ( kind = 4 ) lls
  integer ( kind = 4 ) ls
  integer ( kind = 4 ) lu
  integer ( kind = 4 ), parameter :: maxit = 30
  integer ( kind = 4 ) mm
  integer ( kind = 4 ) mm1
  integer ( kind = 4 ) mn
  integer ( kind = 4 ) nct
  integer ( kind = 4 ) nctp1
  integer ( kind = 4 ) ncu
  integer ( kind = 4 ) nrt
  integer ( kind = 4 ) nrtp1
  real ( kind = 8 ) s(*)
  real ( kind = 8 ) scale
  real ( kind = 8 ) ddot
  real ( kind = 8 ) shift
  real ( kind = 8 ) sl
  real ( kind = 8 ) sm
  real ( kind = 8 ) smm1
  real ( kind = 8 ) sn
  real ( kind = 8 ) dnrm2
  real ( kind = 8 ) t
  real ( kind = 8 ) t1
  real ( kind = 8 ) test
  real ( kind = 8 ) u(ldu,m)
  real ( kind = 8 ) v(ldv,n)
  logical wantu
  logical wantv
  real ( kind = 8 ) work(m)
  real ( kind = 8 ) ztest
!
!  Determine what is to be computed.
!
  wantu = .false.
  wantv = .false.
  jobu = mod ( job, 100 ) / 10

  if ( 1 < jobu ) then
    ncu = min ( m, n )
  else
    ncu = m
  end if

  if ( jobu /= 0 ) then
    wantu = .true.
  end if

  if ( mod ( job, 10 ) /= 0 ) then
    wantv = .true.
  end if
!
!  Reduce A to bidiagonal form, storing the diagonal elements
!  in S and the super-diagonal elements in E.
!
  info = 0
  nct = min ( m-1, n )
  nrt = max ( 0, min ( m, n-2 ) )
  lu = max ( nct, nrt )

  do l = 1, lu
!
!  Compute the transformation for the L-th column and
!  place the L-th diagonal in S(L).
!
    if ( l <= nct ) then

      s(l) = dnrm2 ( m-l+1, a(l,l), 1 )

      if ( s(l) /= 0.0D+00 ) then
        if ( a(l,l) /= 0.0D+00 ) then
          s(l) = sign ( s(l), a(l,l) )
        end if
        call dscal ( m-l+1, 1.0D+00 / s(l), a(l,l), 1 )
        a(l,l) = 1.0D+00 + a(l,l)
      end if

      s(l) = -s(l)

    end if

    do j = l+1, n
!
!  Apply the transformation.
!
      if ( l <= nct .and. s(l) /= 0.0D+00 ) then
        t = -ddot ( m-l+1, a(l,l), 1, a(l,j), 1 ) / a(l,l)
        call daxpy ( m-l+1, t, a(l,l), 1, a(l,j), 1 )
      end if
!
!  Place the L-th row of A into E for the
!  subsequent calculation of the row transformation.
!
      e(j) = a(l,j)

    end do
!
!  Place the transformation in U for subsequent back multiplication.
!
    if ( wantu .and. l <= nct ) then
      u(l:m,l) = a(l:m,l)
    end if
!
!  Compute the L-th row transformation and place the
!  L-th superdiagonal in E(L).
!
    if ( l <= nrt ) then

      e(l) = dnrm2 ( n-l, e(l+1), 1 )

      if ( e(l) /= 0.0D+00 ) then
        if ( e(l+1) /= 0.0D+00 ) then
          e(l) = sign ( e(l), e(l+1) )
        end if
        call dscal ( n-l, 1.0D+00 / e(l), e(l+1), 1 )
        e(l+1) = 1.0D+00 + e(l+1)
      end if

      e(l) = -e(l)
!
!  Apply the transformation.
!
      if ( l + 1 <= m .and. e(l) /= 0.0D+00 ) then

        work(l+1:m) = 0.0D+00

        do j = l+1, n
          call daxpy ( m-l, e(j), a(l+1,j), 1, work(l+1), 1 )
        end do

        do j = l+1, n
          call daxpy ( m-l, -e(j)/e(l+1), work(l+1), 1, a(l+1,j), 1 )
        end do

      end if
!
!  Place the transformation in V for subsequent back multiplication.
!
      if ( wantv ) then
        v(l+1:n,l) = e(l+1:n)
      end if

    end if

  end do
!
!  Set up the final bidiagonal matrix of order MN.
!
  mn = min ( m + 1, n )
  nctp1 = nct + 1
  nrtp1 = nrt + 1

  if ( nct < n ) then
    s(nctp1) = a(nctp1,nctp1)
  end if

  if ( m < mn ) then
    s(mn) = 0.0D+00
  end if

  if ( nrtp1 < mn ) then
    e(nrtp1) = a(nrtp1,mn)
  end if

  e(mn) = 0.0D+00
!
!  If required, generate U.
!
  if ( wantu ) then

    u(1:m,nctp1:ncu) = 0.0D+00

    do j = nctp1, ncu
      u(j,j) = 1.0D+00
    end do

    do ll = 1, nct

      l = nct - ll + 1

      if ( s(l) /= 0.0D+00 ) then

        do j = l+1, ncu
          t = -ddot ( m-l+1, u(l,l), 1, u(l,j), 1 ) / u(l,l)
          call daxpy ( m-l+1, t, u(l,l), 1, u(l,j), 1 )
        end do

        u(l:m,l) = -u(l:m,l)
        u(l,l) = 1.0D+00 + u(l,l)
        u(1:l-1,l) = 0.0D+00

      else

        u(1:m,l) = 0.0D+00
        u(l,l) = 1.0D+00

      end if

    end do

  end if
!
!  If it is required, generate V.
!
  if ( wantv ) then

    do ll = 1, n

      l = n - ll + 1

      if ( l <= nrt .and. e(l) /= 0.0D+00 ) then

        do j = l + 1, n
          t = -ddot ( n-l, v(l+1,l), 1, v(l+1,j), 1 ) / v(l+1,l)
          call daxpy ( n-l, t, v(l+1,l), 1, v(l+1,j), 1 )
        end do

      end if

      v(1:n,l) = 0.0D+00
      v(l,l) = 1.0D+00

    end do

  end if
!
!  Main iteration loop for the singular values.
!
  mm = mn
  iter = 0

  do while ( 0 < mn )
!
!  If too many iterations have been performed, set flag and return.
!
    if ( maxit <= iter ) then
      info = mn
      return
    end if
!
!  This section of the program inspects for
!  negligible elements in the S and E arrays.
!
!  On completion the variables KASE and L are set as follows:
!
!  KASE = 1     if S(MN) and E(L-1) are negligible and L < MN
!  KASE = 2     if S(L) is negligible and L < MN
!  KASE = 3     if E(L-1) is negligible, L < MN, and
!               S(L), ..., S(MN) are not negligible (QR step).
!  KASE = 4     if E(MN-1) is negligible (convergence).
!
    do ll = 1, mn

      l = mn - ll

      if ( l == 0 ) then
        exit
      end if

      test = abs ( s(l) ) + abs ( s(l+1) )
      ztest = test + abs ( e(l) )

      if ( ztest == test ) then
        e(l) = 0.0D+00
        exit
      end if

    end do

    if ( l == mn - 1 ) then

      kase = 4

    else

      do lls = l + 1, mn + 1

        ls = mn - lls + l + 1

        if ( ls == l ) then
          exit
        end if

        test = 0.0D+00
        if ( ls /= mn ) then
          test = test + abs ( e(ls) )
        end if

        if ( ls /= l + 1 ) then
          test = test + abs ( e(ls-1) )
        end if

        ztest = test + abs ( s(ls) )

        if ( ztest == test ) then
          s(ls) = 0.0D+00
          exit
        end if

      end do

      if ( ls == l ) then
        kase = 3
      else if ( ls == mn ) then
        kase = 1
      else
        kase = 2
        l = ls
      end if

    end if

    l = l + 1
!
!  Deflate negligible S(MN).
!
    if ( kase == 1 ) then

      mm1 = mn - 1
      f = e(mn-1)
      e(mn-1) = 0.0D+00

      do kk = l, mm1

        k = mm1 - kk + l
        t1 = s(k)
        call drotg ( t1, f, cs, sn )
        s(k) = t1

        if ( k /= l ) then
          f = -sn * e(k-1)
          e(k-1) = cs * e(k-1)
        end if

        if ( wantv ) then
          call drot ( n, v(1,k), 1, v(1,mn), 1, cs, sn )
        end if

      end do
!
!  Split at negligible S(L).
!
    else if ( kase == 2 ) then

      f = e(l-1)
      e(l-1) = 0.0D+00

      do k = l, mn

        t1 = s(k)
        call drotg ( t1, f, cs, sn )
        s(k) = t1
        f = -sn * e(k)
        e(k) = cs * e(k)
        if ( wantu ) then
          call drot ( m, u(1,k), 1, u(1,l-1), 1, cs, sn )
        end if

      end do
!
!  Perform one QR step.
!
    else if ( kase == 3 ) then
!
!  Calculate the shift.
!
      scale = max ( abs ( s(mn) ), abs ( s(mn-1) ), abs ( e(mn-1) ), &
                    abs ( s(l) ), abs ( e(l) ) )

      sm = s(mn) / scale
      smm1 = s(mn-1) / scale
      emm1 = e(mn-1) / scale
      sl = s(l) / scale
      el = e(l) / scale
      b = ( ( smm1 + sm ) * ( smm1 - sm ) + emm1 * emm1 ) / 2.0D+00
      c = sm  * sm * emm1 * emm1
      shift = 0.0D+00

      if ( b /= 0.0D+00 .or. c /= 0.0D+00 ) then
        shift = sqrt ( b * b + c )
        if ( b < 0.0D+00 ) then
          shift = -shift
        end if
        shift = c / ( b + shift )
      end if

      f = ( sl + sm ) * ( sl - sm ) + shift
      g = sl * el
!
!  Chase zeros.
!
      mm1 = mn - 1

      do k = l, mm1

        call drotg ( f, g, cs, sn )

        if ( k /= l ) then
          e(k-1) = f
        end if

        f = cs * s(k) + sn * e(k)
        e(k) = cs * e(k) - sn * s(k)
        g = sn * s(k+1)
        s(k+1) = cs * s(k+1)

        if ( wantv ) then
          call drot ( n, v(1,k), 1, v(1,k+1), 1, cs, sn )
        end if

        call drotg ( f, g, cs, sn )
        s(k) = f
        f = cs * e(k) + sn * s(k+1)
        s(k+1) = -sn * e(k) + cs * s(k+1)
        g = sn * e(k+1)
        e(k+1) = cs * e(k+1)

        if ( wantu .and. k < m ) then
          call drot ( m, u(1,k), 1, u(1,k+1), 1, cs, sn )
        end if

      end do

      e(mn-1) = f
      iter = iter + 1
!
!  Convergence.
!
    else if ( kase == 4 ) then
!
!  Make the singular value nonnegative.
!
      if ( s(l) < 0.0D+00 ) then
        s(l) = -s(l)
        if ( wantv ) then
          v(1:n,l) = -v(1:n,l)
        end if
      end if
!
!  Order the singular value.
!
      do

        if ( l == mm ) then
          exit
        end if

        if ( s(l+1) <= s(l) ) then
          exit
        end if

        t = s(l)
        s(l) = s(l+1)
        s(l+1) = t

        if ( wantv .and. l < n ) then
          call dswap ( n, v(1,l), 1, v(1,l+1), 1 )
        end if

        if ( wantu .and. l < m ) then
          call dswap ( m, u(1,l), 1, u(1,l+1), 1 )
        end if

        l = l + 1

      end do

      iter = 0
      mn = mn - 1

    end if

  end do

  return
end
end program non_intrusive_high_materials
