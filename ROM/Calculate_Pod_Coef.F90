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
program calculate_pod_coef
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
  !use reduced_model_runtime
  !use momentum_equation_reduced
  !USE IFPORT
  !Use DFLib
  !use checkpoint

  implicit none
#ifdef HAVE_PETSC
#include "finclude/petsc.h"
#endif
  integer, dimension(:), allocatable :: indices,indices_tmp
  type(state_type), dimension(:,:,:), allocatable :: state_p  
  integer :: timestep, dimen
  integer :: ierr
  integer :: numphase,total_dumps
  character(len = OPTION_PATH_LEN) :: simulation_name
  integer :: no_smolyak_nodes, m, sparse_max,snapshots, u_nodes, statei,p_nodes,v_nodes
  type(state_type), dimension(:), pointer :: state_show=> null()
  type(state_type), dimension(:,:,:), allocatable :: pod_state 
  type(state_type), dimension(:), allocatable :: state

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
  integer :: nsvd
  real, dimension(:,:,:), allocatable :: pod_coef_all_obv
#ifdef HAVE_MPI
  call mpi_init(ierr)
#endif

#ifdef HAVE_PETSC
  call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
#endif

  call python_init()
  call read_command_line()
  no_smolyak_nodes=41
  m = 4
  sparse_max =2
  call calculate_pod_coef_main() 
  do statei=1, size(state,1)
  !call deallocate(state(statei,:,:))
  call deallocate(state_p(statei,:,:))
  enddo
  call deallocate(state)
  deallocate(pod_state)

  deallocate(pod_coef_all_obv)
  call deallocate_transform_cache()

 ! call print_references(0)
#ifdef HAVE_MEMORY_STATS
  call print_current_memory_stats(0)
#endif

#ifdef HAVE_MPI
  call mpi_finalize(ierr)
#endif

contains 

  subroutine calculate_pod_coef_main() 
     type(vector_field), pointer :: snapmean_velocity 
      type(scalar_field), pointer :: snapmean_pressure , snapmean_volumefraction 
      type(vector_field), pointer :: PODBasis_velocity 
      type(scalar_field), pointer :: PODBasis_pressure, podbasis_volumefraction
    type(vector_field), pointer :: velocity_field
    type(scalar_field), pointer :: pressure_field, volumefraction_field
   integer :: i,j,k,itime
   call get_option('/simulation_name',simulation_name)
   call initialise_write_state
    ! Read state from .flml file
  ! call populate_state(state)
   call read_input_states(state)
   call read_pod_basis_interp(POD_state, state)
    call print_state(state(1))
    numphase=option_count("/material_phase")
   call get_option(&
         '/reduced_model/pod_basis_formation/pod_basis_count', nsvd)
   allocate(pod_coef_all_obv(total_dumps,2*nsvd,numphase)) 
   do itime=1, total_dumps
      if (have_option("/reduced_model/Velocity_Interpolation")) then
                 print *,'need coding'
              
                else 
                              
                 do k=1, numphase 
                    volumefraction_field=>extract_scalar_field (state(itime), "phase"//trim(int2str(k))//"::PhaseVolumeFraction") ! 1 donates phase 1.                              
                     pressure_field=>extract_scalar_field (state(itime), "phase1::Pressure")                     
                     velocity_field=>extract_vector_field (state(itime), "phase"//trim(int2str(k))//"::Velocity")   
                   do i=1,size(POD_state,1)                                               
                             snapmean_pressure=>extract_Scalar_field(POD_state(i,2,k),"SnapmeanPressure")
                             snapmean_volumefraction=>extract_Scalar_field(POD_state(i,3,k),"SnapmeanVolumefraction")                           
                             PODBasis_pressure=>extract_scalar_field(POD_state(i,2,k), "PODPressure")
                             PODBasis_volumefraction=>extract_scalar_field(POD_state(i,3,k), "PODVolumefraction")
                             pod_coef_all_obv(itime,i,k)= &
                                        dot_product(PODBasis_pressure%val(:),(pressure_field%val(:)-snapmean_pressure%val(:)))                             
                             pod_coef_all_obv(itime,i+nsvd,k)= &
                                        dot_product(PODBasis_volumefraction%val(:),(volumefraction_field%val(:)-snapmean_volumefraction%val(:)))
                     enddo  
                    
                  ! write(40,*)(pod_coef(i),i=1,(u%dim+1)*size(POD_state,1))
                    if((itime .eq. 1).and.(k.eq.1)) then                       
                    	open(11,file="coef_pod_all_obv"//"_"//trim(simulation_name))
                          write(11,*)('')
                         
                   	write(11,*)(pod_coef_all_obv(itime,:,k))
                        print *, 'total_timestep----timestep===coef_pod_all_obv', total_dumps, itime
                    
                   	close(11) 

                   
                    else 
                      open(11,file="coef_pod_all_obv"//"_"//trim(simulation_name), position='append',ACTION='WRITE')   
                                	
                      write(11,*)(pod_coef_all_obv(itime,:,k))
                      print *, 'total_timestep----timestep===coef_pod_all_obv', total_dumps, itime
                    
                      close(11)
                    
                   endif
                enddo !numphase

                ! endif     ! if(its==nonlinear_iterations) then  
               endif !if (have_option("/reduced_model/Velocity_Interpolation") then
               
           enddo 

           
      
  end subroutine calculate_pod_coef_main

  


 subroutine read_input_states(state)
 
    type(state_type), intent(out), dimension(:), allocatable :: state
    character(len=1024) :: filename

    integer :: dump_sampling_period, quadrature_degree
    integer :: i,j,k,stable_dumps

    call get_option('/reduced_model/pod_basis_formation/dump_sampling_period',dump_sampling_period)
    call get_option('/geometry/quadrature/degree', quadrature_degree)

    ewrite(3,*) 'dump_sampling_period',dump_sampling_period
 
    total_dumps=count_dumps(dump_sampling_period)-1
    print *, 'total_snapshots=,', total_dumps
    allocate(state(total_dumps)) 

    do i=1, total_dumps  
      
       write(filename, '(a, i0, a)') trim(simulation_name)//'_', (i)*dump_sampling_period,".vtu"   
       call vtk_read_state(filename, state(i), quadrature_degree)  
        
     end do
   
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

 subroutine read_pod_basis_interp(POD_state, state)
    
    !! Read the podbasis from the set of vtu files.
   
    character(len=1024) :: simulation_name, filename, phase_name 
    integer :: dump_period, quadrature_degree, numphase
    integer :: i,j,POD_num,k,nfield
    ! type(state_type), dimension(:,:,:), allocatable :: POD_state
    type(state_type), dimension(:,:,:), intent(out), allocatable :: POD_state !POD_state(podnum,:,numphase)
    type(state_type), dimension(:,:,:), allocatable ::  POD_state_p,POD_state_v!POD_state(podnum,:,numphase)
    type(state_type), dimension(:) :: state
    type(vector_field) :: podVelocity, newpodVelocity 
    type(scalar_field) :: podPressure, newpodPressure, podVolumefraction, newpodVolumefraction, podTemperature, newpodTemperature ,snapmean_pressure, snapmean_volumefraction
    type(mesh_type) :: VelocityMesh, PressureMesh,TemperatureMesh
    
    type(scalar_field), pointer :: pres
    integer :: pod_pnodes,pod_unodes,p_nodes,u_nodes
    !type(mesh_type) ,pointer ::pmesh
    call get_option('/simulation_name', simulation_name)
    print *, 'simulation_name=', simulation_name
    ! if (have_option("/reduced_model/execute_reduced_model")) then
   ! If(have_option("/reduced_model/Smolyak")) then
    if (have_option("/reduced_model/Smolyak") .or. have_option("/reduced_model/RBF_interpolation")) then
        simulation_name(len_trim(simulation_name)-3:len_trim(simulation_name))="" !delete _POD in the name
    endif
    print *, 'simulation_name==', simulation_name
    call get_option('/geometry/quadrature/degree', quadrature_degree)
    call get_option(&
         "/reduced_model/pod_basis_formation/pod_basis_count", POD_num) 
    numphase=option_count("/material_phase")
    !total_dumps=POD_num!count_dumps(simulation_name)
    nfield = vector_field_count( state(1) )+scalar_field_count( state(1))+tensor_field_count(state(1))
    
   ! allocate(pod_state(POD_num,nfield,numphase))
    allocate(pod_state_p(POD_num,nfield,numphase))
    allocate(pod_state_v(POD_num,nfield,numphase))
    allocate(pod_state(POD_num,nfield,numphase))
    VelocityMesh=extract_velocitymesh_from_state(state(1))
    PressureMesh=extract_pressuremesh_from_state(state(1))  
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
 
 function extract_pressuremesh_from_state(state, stat)
    type(state_type), intent(in):: state
    type(mesh_type), pointer:: extract_pressuremesh_from_state
    integer, optional, intent(out):: stat
      
    type(scalar_field), pointer:: p
      
    p => extract_scalar_field(state, "phase1::Pressure", stat=stat)
    if (associated(p)) then
      extract_pressuremesh_from_state => p%mesh
    else
      nullify(extract_pressuremesh_from_state)
    end if
  
  end function extract_pressuremesh_from_state

  
  
  function extract_velocitymesh_from_state(state, stat)
    type(state_type), intent(in):: state
    type(mesh_type), pointer:: extract_velocitymesh_from_state
    integer, optional, intent(out):: stat
      
    type(vector_field), pointer:: u
      
    u => extract_vector_field(state, "phase1::Velocity", stat=stat)
    if (associated(u)) then
      extract_velocitymesh_from_state => u%mesh
    else
      nullify(extract_velocitymesh_from_state)
    end if
  
  end function extract_velocitymesh_from_state

end program calculate_pod_coef
