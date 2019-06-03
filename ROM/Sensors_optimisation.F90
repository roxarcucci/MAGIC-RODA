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
  !    This code do non-intrusive ROM directly through Vtu files. 
  !    Read Vtu files into states, then form podbasis through those snapshots, 
  !    then non-intrusive, 28,May,2015 by Dunhui Xiao
  !    You should have received a copy of the GNU Lesser General Public
  !    License along with this library; if not, write to the Free Software
  !    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
  !    USA
  !

!  
#include "fdebug.h"


program sensors_optimisation
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
 


  implicit none
!#ifdef HAVE_PETSC
!#include "finclude/petsc.h"
!#endif
  integer, dimension(:), allocatable :: indices,indices_tmp
  type(state_type), dimension(:,:,:), allocatable :: state ,state_p  
  integer :: timestep, dimen
  integer :: ierr
  integer :: numphase,total_dumps
  real :: ADD_EPSILON2
  character(len = OPTION_PATH_LEN) :: simulation_name
  integer :: no_smolyak_nodes, sparse_max, snapshots, u_nodes, statei,p_nodes,v_nodes
  type(state_type), dimension(:), pointer :: state_show=> null()
  type(state_type), dimension(:,:,:), allocatable :: pod_state 
  type(state_type), dimension(:,:,:), allocatable :: pod_state_p, pod_state_v
  real, dimension(:,:,:,:), allocatable ::  interpolate_velocity
  real, dimension(:,:,:), allocatable ::  C_in
  real, dimension(:,:,:,:), allocatable :: snapmatrix_pressure, snapmatrix_volumefraction
  real, dimension(:,:,:), allocatable ::  interpolate_pressure,interpolate_volumefraction
  real, dimension(:,:,:,:), allocatable :: leftsvd_velocity 
  real, dimension(:,:,:), allocatable :: leftsvd_pressure,  leftsvd_volumefraction
  real, dimension(:,:,:), allocatable :: svdval_velocity 
  real, dimension(:,:), allocatable :: svdval_pressure, svdval_pressure_deim,svdval_volumefraction
  real, dimension(:,:,:), allocatable :: snapmean_velocity  
  real, dimension(:,:), allocatable :: snapmean_pressure, snapmean_volumefraction 
  real, dimension(:,:,:), allocatable :: pod_coef_all_obv
  real, dimension(:,:), allocatable :: G_in
  real, dimension(:), allocatable :: G_SUM_in
  real, dimension(:), allocatable :: WEIGHT_NOD
  type(vector_field), pointer :: nodes_positions
  integer :: nsvd,NCON, NWINDOW
  integer :: i, dump_no, d, j,k ,m,i2,i3,ph,p,s
  real :: minx,miny, maxx, maxy, minz, maxz
  character(len=1024) :: TracerFront, TracerBack,TracerObelisk,TracerWestminster,TracerWaterloo,TracerBlackfriars,TracerLambeth,TracerLondon,TracerGeorge,tracer_name
  type(vector_field), pointer :: velocity_field
  type(scalar_field), pointer :: pressure_field, volumefraction_field
    

#ifdef HAVE_MPI
  call mpi_init(ierr)
#endif

#ifdef HAVE_PETSC
  !call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
#endif

  call python_init()
  call read_command_line()
  no_smolyak_nodes=1 !no use in this program, but keep it for consitency with interpolate_snapshots_nonintruve_material_poeperty.
    numphase=option_count("/material_phase")
    numphase=1
    NCON=9
    NWINDOW=10
  ! m = 4
  ! sparse_max =12
  ! call get_option('/reduced_model/smolyak_order', sparse_max) 
  !call importance_map_main() 
  call initialise_write_state
   call get_option('/simulation_name',simulation_name) 
    ! Read state from .flml file
    call populate_state(state_show)
    call print_state(state_show(1))   
    !call put_vtu_to_states(state,state_p)
    call put_pvtu_to_states(state,state_p)
     nodes_positions=>extract_vector_field(state(1,1,1), "Coordinate")
     minx=nodes_positions%val(1,1)
     miny=nodes_positions%val(1,1)
     minz=nodes_positions%val(1,1)
     maxx=nodes_positions%val(1,1)
     maxy=nodes_positions%val(1,1)
     maxz=nodes_positions%val(1,1)
     
    
     Do i=1, size(nodes_positions%val,2) !number of nodes
           if(minx>nodes_positions%val(1,i)) then
              minx=nodes_positions%val(1,i)
           endif
           if(maxx<nodes_positions%val(1,i)) then
              maxx=nodes_positions%val(1,i)
           endif
      enddo

     do i=1, size(nodes_positions%val,2)
           if(miny>nodes_positions%val(2,i)) then
              miny=nodes_positions%val(2,i)
           endif
           if(maxy<nodes_positions%val(2,i)) then
              maxy=nodes_positions%val(2,i)
           endif
      enddo
    
      do i=1, size(nodes_positions%val,2)
           if(minz>nodes_positions%val(3,i)) then
              minz=nodes_positions%val(3,i)
           endif
           if(maxz<nodes_positions%val(3,i)) then
              maxz=nodes_positions%val(3,i)
           endif
      enddo
     print *, 'size(nodes_positions%val,2)', size(nodes_positions%val,2), p_nodes, minx,maxx, miny, maxy,minz, maxz
    !stop 345
      allocate(G_in(p_nodes, NWINDOW))
    allocate(G_SUM_in(p_nodes))
    allocate(C_in(p_nodes, NCON, snapshots))
   allocate(interpolate_pressure(p_nodes,snapshots,NCON))  
  write(tracer_name, '(a)')   "TracerFront"
  call importance_map_main_loop(tracer_name,1)
  write(tracer_name, '(a)')   "TracerBack"
   call importance_map_main_loop(tracer_name,2)  
  write(tracer_name, '(a)')   "TracerObelisk"
   call importance_map_main_loop(tracer_name,3)
   write(tracer_name, '(a)')   "TracerWestminster"
  call importance_map_main_loop(tracer_name,4)
  write(tracer_name, '(a)')   "TracerWaterloo"
   call importance_map_main_loop(tracer_name,5)  
  write(tracer_name, '(a)')   "TracerBlackfriars"
   call importance_map_main_loop(tracer_name,6)
   write(tracer_name, '(a)')   "TracerLambeth"
  ! write(tracer_name, '(a)')   "TracerFront"
  call importance_map_main_loop(tracer_name,7)
  write(tracer_name, '(a)')   "TracerLondon"
   call importance_map_main_loop(tracer_name,8)  
  write(tracer_name, '(a)')   "TracerGeorge"
   call importance_map_main_loop(tracer_name,9)
    
     allocate(WEIGHT_NOD(p_nodes))
    WEIGHT_NOD(:)=0 
    ADD_EPSILON2=1.E-12
     do d=1, NCON 
     do i=1, snapshots
       C_in(:, d, i)=interpolate_pressure(:,i,d)
       do j=1, size(nodes_positions%val,2)  
            if ((nodes_positions%val(1,j) .GT. 0.55*(maxx-minx)+minx) .and. (nodes_positions%val(1,j) .LE. 0.58*(maxx-minx)+minx) .and. (nodes_positions%val(2,j) .GT. 0.40*(maxy-miny)+miny) .and. (nodes_positions%val(2,j) .LE. 0.45*(maxy-miny)+miny) .and. (nodes_positions%val(3,j) .GT. 0.2) .and. (nodes_positions%val(3,j) .LE. 0.5)) then !0.55
          ! WEIGHT_NOD(j)=1
           print *, 'id', j
           endif     
          if (nodes_positions%val(1,j) .LE. 0.45*(maxx-minx)+minx) then
           C_in(j, d, i)=0 !print *, 'id', j
           endif
           if (nodes_positions%val(1,j) .GT. 0.655*(maxx-minx)+minx ) then !0.55
           C_in(j, d, i)=0 !print *, 'id', j 
           endif 
           if (nodes_positions%val(2,j) .GT. 0.6*(maxy-miny)+miny ) then !0.55
           C_in(j, d, i)=0 
           !print *, 'id', j
           endif  
       enddo 
     ! stop 34567
      enddo
     enddo 
  ! print *, 'size location of sensors', size(nodes_positions%val,1), size(nodes_positions%val,2) 
  ! print *, 'location of sensors',nodes_positions%val(:, 25841), 'location of sensors', nodes_positions%val(:,25812),   nodes_positions%val(:,20106), nodes_positions%val(:,13598)
   !stop 11
   !C_in=0
   
    G_SUM_in=0    
   
    WEIGHT_NOD(91046)=0
   WEIGHT_NOD(91047)=1
   print *, 'location of sensors,91047',nodes_positions%val(:, 91047), 'location of sensors,91046', nodes_positions%val(:,91046)
  ! location of sensors,91047   158.96729688999994       -188.94147613000032       0.20000000298000001      location of sensors,91046   149.39312487694247       -174.07346281543175       0.20000000298002982
  !stop 235
    G_in=0
    call G_CALC(C_in, WEIGHT_NOD, p_nodes, NCON, total_dumps, NWINDOW, G_in, G_SUM_in, .true., ADD_EPSILON2 )
    print *, 'G calculating correct'
   !  print *, 'G_in', G_in
  
    print *, 'start G calculating'
   !visulisiation!
     if (.true.) then
      DO s=1, NWINDOW 
       pressure_field => extract_scalar_field(state(1,s,1), tracer_name ) 
       pressure_field%val(:)=0 
       pressure_field%val(:)= G_in(:, s)   
       dump_no=NWINDOW-s+1 
       call vtk_write_state(filename=trim(tracer_name)//'_importance_map_G', index=dump_no, model="Mesh", state=state(:,s,1))  
     ENDDO ! Do s=1,  NWINDOW 
     endif  
       s=1
       pressure_field => extract_scalar_field(state(1,s,1), tracer_name ) 
       pressure_field%val(:)=0 
       pressure_field%val(:)= G_SUM_in(:)
       ! do j=1, size(nodes_positions%val,2) 
        !   print *, 'nodes_positions%val(1,j)',j, nodes_positions%val(1,j), nodes_positions%val(2,j), nodes_positions%val(3,j)!,minx,maxx,miny,maxy
          ! if (nodes_positions%val(1,j) .LE. 0.45*(maxx-minx)+minx) then
           ! pressure_field%val(j)=0 
          ! endif
          ! if (nodes_positions%val(1,j) .GT. 0.55*(maxx-minx)+minx ) then
          ! pressure_field%val(j)=0 
          ! endif
        !enddo        
       dump_no=s   
       call vtk_write_state(filename='Importance_map_G_sum', index=dump_no, model="Mesh", state=state(:,s,1))  
   
  do statei=1, size(state,1)
  call deallocate(state(statei,:,:)) 
  enddo
   ! deallocate(interpolate_velocity) ! need more
    deallocate(interpolate_pressure)
    deallocate(G_in)
    deallocate(G_SUM_in)
    deallocate(C_in) 
    
  call deallocate_transform_cache()

  !call print_references(0)
#ifdef HAVE_MEMORY_STATS
  call print_current_memory_stats(0)
#endif

#ifdef HAVE_MPI
  call mpi_finalize(ierr)
#endif

contains 
    subroutine importance_map_main_loop(tracer_name, no_tracer)
    character(len=1024), intent(in) :: tracer_name
    real, dimension(:,:,:,:,:), allocatable :: snapmatrix_velocity  
    
    integer :: stat       
    real :: cpu_basis_start, cpu_basis_end
    integer :: udim,deim_number_tmp  
    type(vector_field), pointer :: velocityudim 
    type(vector_field), pointer :: vfield
    character(len=1024) :: phase_name
    character(len = OPTION_PATH_LEN) :: simulation_name1 
    integer, intent(in):: no_tracer
    
    !stop 2121
    call retrieve_tracer_nodesvalue_from_states(state, snapmatrix_pressure, tracer_name) 
     
    !allocate(interpolate_velocity(u_nodes,snapshots,dimen,numphase)) ! need more
     
     !interpolate_velocity=snapmatrix_velocity(1,:,:,:,:) 
     interpolate_pressure(:,:,no_tracer)=snapmatrix_pressure(1,:,:,1) !-1000000
     
    !print *, C_in!,WEIGHT_NOD, p_nodes, NCON, total_dumps, NWINDOW, G_in, G_SUM_in , ADD_EPSILON2 
   ! print *,  WEIGHT_NOD, p_nodes, NCON, total_dumps, NWINDOW, G_in, G_SUM_in , ADD_EPSILON2
    !stop 111 
   ! deallocate(snapmatrix_velocity)
    deallocate(snapmatrix_pressure) 
   
  end subroutine importance_map_main_loop
 

  SUBROUTINE G_CALC(C, WEIGHT_NOD, NONODS, NCON, NTIME, NWINDOW, G,G_SUM, SUBTRACT_MEAN_C, ADD_EPSILON2 )
! This sub calculates the value of G at every time level of 
! the time window AND also G_SUM. It does this by assuming the flow is chaotic and 
! thus there are natual fluctuations that can be used to calculate 
! sensitivities with. 
! If SUBTRACT_MEAN_C then subtract out the mean from the concentration field 
! as this then follows the perturbation theory a bit more may be - try both.
! ADD_EPSILON = the no added onto diagonal diagonal for conditioning of the Mel Petros sudo inverse...
! if -ve then use the defaut value. 
         IMPLICIT NONE
         INTEGER, intent(in) :: NONODS, NCON, NTIME, NWINDOW
         LOGICAL, intent(in) :: SUBTRACT_MEAN_C
         real, intent(in) :: ADD_EPSILON2
         REAL, dimension(NONODS, NCON, NTIME), intent(in) :: C 
         REAL, dimension(NONODS), intent(in) :: WEIGHT_NOD
! C is the concentration fields
! weight_nod contains the weights of the nodes used to calculate 
! the functional F=sum(WEIGHT_NOD,C(:,ICON,ITIME))
! NONODS is the no of nodes/CVs in the mesh
! NCON is the no of concentration fields solved for each time level
! NTIM is the no of time levels. 
! NWINDOW is the no of time levels in the time window used to calculate 
! the senitivies (importance map) so there will be NWINDOW sensitivties e.g. 40
         REAL, dimension(NONODS,NWINDOW), intent(inout) :: G
         REAL, dimension(NONODS), intent(inout) :: G_SUM
! local variables...
         REAL, PARAMETER :: ADD_EPSILON_DEFAULT = 1.E-12
         REAL, dimension( :, : ), allocatable :: C_WIN, C_MEAN
         REAL, dimension( : ), allocatable :: DF
         integer N_ENSEMBLE, IWIN, I_ENSEMB, ITIME, ITIME2, ICON
         real ADD_EPSILON

         ADD_EPSILON=ADD_EPSILON2
         IF(ADD_EPSILON2<0.0) ADD_EPSILON=ADD_EPSILON_DEFAULT

         N_ENSEMBLE=NCON*(NTIME-NWINDOW+1)
         allocate(C_WIN(NONODS,N_ENSEMBLE))
         allocate(DF(N_ENSEMBLE))
         allocate(C_MEAN(NONODS,NCON) )
! F is the functional evaluated for each ensemble member. 
         C_MEAN = 0.0
         IF(SUBTRACT_MEAN_C) THEN
            DO ICON=1,NCON
               DO ITIME=1, NTIME
                  C_MEAN(:,ICON) = C_MEAN(:,ICON) + C(:,ICON,ITIME)
               END DO
                  C_MEAN(:,ICON) = C_MEAN(:,ICON)/REAL(NTIME)
            END DO
         ENDIF ! ENDOF IF(SUBTRACT_MEAN_C) THEN

         G_SUM=0.0
         DO IWIN=1,NWINDOW
! Calaculate G for this time window...
! form C for time window IWIN
            print *, 'NWINDOW, NTIME=', IWIN, NTIME
         ! stop 1212
            I_ENSEMB=0
            DO ITIME2 = NWINDOW, NTIME ! ITIME2 time level is over the last time level of the time window.
              ITIME=ITIME2 - NWINDOW + IWIN ! real time level
              DO ICON=1,NCON
                 I_ENSEMB=I_ENSEMB+1
                 C_WIN(:,I_ENSEMB) = C(:,ICON,ITIME) - C_MEAN(:,ICON)
                    !  print *, C_WIN 
                 DF(I_ENSEMB)=sum( WEIGHT_NOD(:)*(C(:,ICON,ITIME2)-C_MEAN(:,ICON)) ) 
                  !  print *, 'DF is not zero, have tested', DF
              END DO
            END DO ! ENDOF DO ITIME2 = NWINDOW, NTIME
              ! stop 1212
! calculate G for time window...
            CALL G_WINDOW(C_WIN,G(:,IWIN),DF,NONODS,N_ENSEMBLE, ADD_EPSILON)
                print *, 'before sum' 
            if(IWIN.ne.NWINDOW) then
              G_SUM(:)=G_SUM(:)+G(:,IWIN) ! Miss out last time window as no info here.
            endif
           ! print *, 'after sum' ,IWIN
         END DO ! DO IWIN=1,NWINDOW
         print *, 'finish calculation'  
        ! DEALLOCATE(C_WIN, DF,C_MEAN)
         deallocate(C_WIN)
          deallocate(DF)
          DEALLOCATE(C_MEAN)
         RETURN
   END SUBROUTINE G_CALC
!
 
!
         SUBROUTINE G_WINDOW(C_WIN_in,G,DF_in,NONODS,N_ENSEMBLE, ADD_EPSILON)
! This sub calculates the value of g at every time level of 
! the time window.
         IMPLICIT NONE
         INTEGER, intent(in) :: NONODS, N_ENSEMBLE
         real, intent(in) :: ADD_EPSILON
         REAL, dimension(NONODS,N_ENSEMBLE), intent(in) :: C_WIN_in 
         REAL, dimension(N_ENSEMBLE), intent(in) :: DF_in
         REAL, dimension(NONODS), intent(inout) :: G
! Local variables...
         REAL MTM(N_ENSEMBLE,N_ENSEMBLE),MTM_INV(N_ENSEMBLE,N_ENSEMBLE),G_SHORT(N_ENSEMBLE)
         INTEGER I_ENSEM,J_ENSEM
         print *, 'N_ENSEMBLE=', N_ENSEMBLE
         !stop 11111
         DO I_ENSEM=1,N_ENSEMBLE
            DO J_ENSEM=1,N_ENSEMBLE
              MTM(I_ENSEM,J_ENSEM)=SUM( C_WIN_in(:,I_ENSEM)*C_WIN_in(:,J_ENSEM) )
            END DO
         END DO
       ! print *,'MT has tested, not zero', MTM
        !  stop 1212
! Add a small no onto diagonal for conditioning of the Mel Petros sudo inverse...
         DO I_ENSEM=1,N_ENSEMBLE
           MTM(I_ENSEM,I_ENSEM)=MTM(I_ENSEM,I_ENSEM) + ADD_EPSILON
         END DO
         MTM_INV=MTM
         call invert(MTM_INV)
       ! print *, 'mtm_inv has tested, not zero', MTM_INV
        !print *, 'DF_IN', DF_in
       ! stop 1212
! Calculate g: G=K*DF with K=C_WIN^T * MTM_INV
         G_SHORT(:)=MATMUL(MTM_INV,DF_in)
         G=0.0
         print *, 'calculating G in G_WIN'
        !  print *, 'G_SHORT(I_ENSEM), tested is zero, why?', G_SHORT(:)
      ! stop 233
         DO I_ENSEM=1,N_ENSEMBLE
             print *, 'in G_WIN loop', I_ENSEM
            !   print *, 'C_WIN(:,I_ENSEM)', C_WIN(:,I_ENSEM)
             
             G(:) = G(:) + C_WIN_in(:,I_ENSEM)*G_SHORT(I_ENSEM) 
            
            ! print *, 'G=', G(:)
         END DO
        ! print *, 'finish calculating G in G_WIN'
         RETURN
       END SUBROUTINE G_WINDOW
  
 
  subroutine retrieve_nodesvalue_from_states(state, snapmatrix_velocity, snapmatrix_pressure, snapmatrix_volumefraction)             
   
    type(state_type), intent(in), dimension(:,:,:) :: state
    real, dimension(:,:,:,:,:), allocatable :: snapmatrix_velocity
    real, dimension(:,:,:,:), allocatable :: snapmatrix_pressure,snapmatrix_volumefraction
     
    ! integer :: numphase
    integer ::  i, d, p,k
    type(vector_field), pointer :: velocity
    type(scalar_field), pointer :: pressure, volumefraction  
   ! print *, 'sizeofstate', size(state,1), size(state,2), size(state,3)
  !  velocity => extract_vector_field(state(1,1,1), "Velocity")
    pressure => extract_scalar_field(state(1,1,1), "TracerLondon") !(state_p(1,1,1), "Pressure") 
    volumefraction=> extract_scalar_field(state(1,1,1), "TracerLondon")!(state_p(1,1,1), "PhaseVolumeFraction")
    dimen=velocity%dim
    p_nodes=node_count(pressure)
    u_nodes=node_count(velocity)
    v_nodes=node_count(volumefraction)
    snapshots=total_dumps!size(state(1,:,1))
    allocate(snapmatrix_velocity(no_smolyak_nodes,u_nodes,snapshots,dimen,numphase))
    allocate(snapmatrix_pressure(no_smolyak_nodes,p_nodes,snapshots,numphase)) 
    allocate(snapmatrix_volumefraction(no_smolyak_nodes,v_nodes,snapshots,numphase))
    print *, 'p_nodes=node_count(pressure)', p_nodes
    snapmatrix_velocity=0.0
    snapmatrix_pressure=0.0  
    snapmatrix_volumefraction=0.0 
   do k=1, no_smolyak_nodes
   do p = 1, numphase
    do i = 1, snapshots
     !print *, 'iiiiiiiii', i
     ! call print_state(state(k,i,p))
      ! velocity => extract_vector_field(state(k,i,p), "Velocity")
           !pressure => extract_scalar_field(state_p(k,i,1), "Pressure")
           !volumefraction=> extract_scalar_field(state_p(k,i,p), "PhaseVolumeFraction")
       pressure => extract_scalar_field(state(k,i,1), "TracerLondon")
       volumefraction=> extract_scalar_field(state(k,i,p), "TracerLondon")
       do d = 1, dimen 
           snapmatrix_velocity(k,:,i,d,p)=velocity%val(d,:)
       end do 
           snapmatrix_pressure(k,:,i,p)=pressure%val
           snapmatrix_volumefraction(k,:,i,p)=volumefraction%val !(p,:)
    end do
   enddo 
   enddo !end numphase      
  end subroutine retrieve_nodesvalue_from_states 

  subroutine retrieve_tracer_nodesvalue_from_states(state,   snapmatrix_pressure, tracer_name)             
    character(len=1024),intent(in) :: tracer_name
    type(state_type), intent(in), dimension(:,:,:) :: state
    real, dimension(:,:,:,:,:), allocatable :: snapmatrix_velocity
    real, dimension(:,:,:,:), allocatable :: snapmatrix_pressure,snapmatrix_volumefraction
     
    ! integer :: numphase
    integer ::  i, d, p,k
    !type(vector_field), pointer :: velocity
    type(scalar_field), pointer :: pressure  
  !  velocity => extract_vector_field(state(1,1,1), "Velocity")
    pressure => extract_scalar_field(state(1,1,1), tracer_name) !(state_p(1,1,1), "Pressure")  
    !dimen=velocity%dim
    p_nodes=node_count(pressure)
    !u_nodes=node_count(velocity) 
    snapshots=total_dumps!size(state(1,:,1))
  !  allocate(snapmatrix_velocity(no_smolyak_nodes,u_nodes,snapshots,dimen,numphase))
    allocate(snapmatrix_pressure(no_smolyak_nodes,p_nodes,snapshots,numphase)) 
     
    
  !  snapmatrix_velocity=0.0
    snapmatrix_pressure=0.0  
     
   do k=1, no_smolyak_nodes
   do p = 1, numphase
    do i = 1, snapshots 
      ! velocity => extract_vector_field(state(k,i,p), "Velocity") 
       pressure => extract_scalar_field(state(k,i,1), tracer_name) 
      ! do d = 1, dimen 
        !   snapmatrix_velocity(k,:,i,d,p)=velocity%val(d,:)
      ! end do 
           snapmatrix_pressure(k,:,i,p)=pressure%val 
    end do
   enddo 
   enddo !end numphase      
  end subroutine retrieve_tracer_nodesvalue_from_states

  subroutine read_input_states(state)
 !!!  for case has the same mesh like P1_P1
     type(state_type), intent(inout), dimension(:,:,:), allocatable :: state 
    character(len=1024) :: filename

    integer :: dump_sampling_period, quadrature_degree
    integer :: i,j,k,stable_dumps 
    character(len=1024) :: phase_name

    call get_option('/reduced_model/pod_basis_formation/dump_sampling_period',dump_sampling_period)
    call get_option('/geometry/quadrature/degree', quadrature_degree)
    ewrite(3,*) 'dump_sampling_period',dump_sampling_period
 
    total_dumps=count_dumps(dump_sampling_period)-1
    allocate(state(no_smolyak_nodes,total_dumps,numphase)) 
  
     do k=1, no_smolyak_nodes
     do j=1, numphase
        ! call get_option("/material_phase["//int2str(j-1)//"]/name", phase_name) ! not'but" 
         ! call get_option("/material_phase["//int2str(m)//"]/name", mat_name)
       !  print *,'phase_name=',  int2str(j), phase_name
    do i=1, total_dumps       
       write(filename, '(a, i0, a)') trim(simulation_name)//'_', (i)*dump_sampling_period,".vtu"   
       call vtk_read_state(filename, state(k,i,j), quadrature_degree) 
       ! write(filename, '(a, i0, a)') trim(simulation_name)//'_POD_', (i)*dump_sampling_period,".vtu"   
       ! call vtk_read_state(filename, state_pod(i), quadrature_degree) 
     end do
    enddo 
   enddo
  end subroutine read_input_states

  subroutine put_vtu_to_states(state,state_p)
    !!< Read the input states from the vtu dumps of the forward run. ! different mesh such as P1_dg_P2
    type(state_type), intent(out), dimension(:,:,:), allocatable :: state
    type(state_type), intent(out), dimension(:,:,:), allocatable :: state_p
    character(len=1024) :: filename
     type(scalar_field) ::pressure,  volumefraction
    integer :: dump_sampling_period, quadrature_degree
    integer :: i,j,k,stable_dumps
    type(mesh_type) ::pressuremesh
   ! integer :: numphase
    character(len=1024) :: phase_name
    
    call get_option('/reduced_model/pod_basis_formation/dump_sampling_period',dump_sampling_period)
    call get_option('/geometry/quadrature/degree', quadrature_degree)

    ewrite(3,*) 'dump_sampling_period',dump_sampling_period
    
    !substract gyre_0.vtu
    total_dumps=count_dumps(dump_sampling_period)
    ! numphase=option_count("/material_phase")
    allocate(state(no_smolyak_nodes,total_dumps,numphase))
    allocate(state_p(no_smolyak_nodes,total_dumps,numphase))
    
    !    vtu--->state
    do k=1, no_smolyak_nodes
      do j=1, numphase
       !  call get_option("/material_phase["//int2str(j-1)//"]/name", phase_name) ! not'but" 
         ! call get_option("/material_phase["//int2str(m)//"]/name", mat_name)
        ! print *,'phase_name=',  int2str(j), phase_name
     do i=1, total_dumps
         ! if(numphase .eq. 1) then
             write(filename, '(a, i0, a,a)')  trim(simulation_name)//'_'//'VelocityMesh_', (i)*dump_sampling_period,'_checkpoint',".vtu" 
             call vtk_read_state(filename, state(k,i,j), quadrature_degree)
             print *,'filename=', filename  
            
             ! write(filename, '(a, i0, a,a)')  trim(simulation_name)//'_'//'PressureMesh_', (i)*dump_sampling_period,'_checkpoint',".vtu" 
           !  call vtk_read_state(filename, state_p(k,i,j), quadrature_degree) 
            !  print *,'filename=', filename
               write(filename, '(a, i0, a)')  trim(simulation_name)//'_', (i)*dump_sampling_period,".vtu" 
              call vtk_read_state(filename, state_p(k,i,j), quadrature_degree)  
              print *,'filename=', filename
              call print_state(state_p(1,1,1))
        !  else
          !   write(filename, '(a, i0, a,a)')  trim(k)//'_'//trim(phase_name)//'_VelocityMesh_', (i)*dump_sampling_period,'_checkpoint',".vtu"
           !  call vtk_read_state(filename, state(k,i,j), quadrature_degree)
           !  write(filename, '(a, i0, a,a)')  trim(k)//'_'//trim(phase_name)//'_PressureMesh_', (i)*dump_sampling_period,'_checkpoint',".vtu" 
           !  call vtk_read_state(filename, state_p(k,i,j), quadrature_degree)
        !  endif           
          pressuremesh =extract_mesh(state_p(k,i,j), "Mesh") 
        call insert(state(k,i,j), pressuremesh, name="PressureMesh") 
      
        pressure=extract_scalar_field(state_p(k,i,j),"Pressure")
       call insert(state(k,i,j),  pressure, name="Pressure")
     
       volumefraction=extract_scalar_field(state_p(k,i,j),"Pressure")
        call insert(state(k,i,j),  volumefraction, name="Pressure")
       ! call insert(state(k,i,j), pressure_field, "Pressure")   change 5 places if changing one field
    enddo
    enddo  
    enddo
   ! call print_state(state(1,1,1))
    call print_state(state_p(1,1,1))
  !  stop 12
   end subroutine put_vtu_to_states
 
 

  subroutine put_pvtu_to_states(state,state_p)
    !!< Read the input states from the pvtu dumps of the forward run. ! different mesh such as P1_dg_P2
    type(state_type), intent(out), dimension(:,:,:), allocatable :: state
    type(state_type), intent(out), dimension(:,:,:), allocatable :: state_p
    character(len=1024) :: filename
     type(scalar_field) ::  volumefraction
    integer :: dump_sampling_period, quadrature_degree
    integer :: i,j,k,stable_dumps
    type(mesh_type) ::pressuremesh
   ! integer :: numphase
    character(len=1024) :: phase_name
    type(scalar_field), pointer :: pressure
    
    call get_option('/reduced_model/pod_basis_formation/dump_sampling_period',dump_sampling_period)
    call get_option('/geometry/quadrature/degree', quadrature_degree)

    ewrite(3,*) 'dump_sampling_period',dump_sampling_period
    
    !substract gyre_0.vtu
    total_dumps=490!count_dumps_pvtu(dump_sampling_period)
    snapshots=total_dumps
   ! numphase=option_count("/material_phase")
    allocate(state(no_smolyak_nodes,total_dumps,numphase))
    !allocate(state_p(no_smolyak_nodes,total_dumps,numphase))
    
    !    vtu--->state
    do k=1, no_smolyak_nodes 
     do i=1, total_dumps 
            ! write(filename, '(a, i0, a,a)')  trim(simulation_name)//'_'//'VelocityMesh_', (i*5)*dump_sampling_period,'_checkpoint',".pvtu" 
              write(filename, '(a, i0, a,a)')  trim(simulation_name)//'_', (i)*dump_sampling_period,".pvtu" 
             call vtk_read_state(filename, state(k,i,1), quadrature_degree)
             print *,'filename=', filename  
            ! write(filename, '(a, i0, a,a)')  trim(simulation_name)//'_'//'PressureMesh_', (i*5)*dump_sampling_period,'_checkpoint',".pvtu" 
            ! call vtk_read_state(filename, state_p(k,i,1), quadrature_degree) 
           !  print *,'filename=', filename  
    enddo  
    enddo
    pressure => extract_scalar_field(state(1,1,1), "TracerFront") !(state_p(1,1,1), "Pressure")   
    p_nodes=node_count(pressure)
   end subroutine put_pvtu_to_states
 
 
  
  function count_dumps(dump_sampling_period) result (count)
    !! Work out how many dumps we're going to read in.
    integer :: count,dump_sampling_period

    logical :: exists
    !      character(len=FILE_NAME_LEN) :: filename
    character(len=1024) :: filename, phase_name
   ! numphase=option_count("/material_phase")

    count=1
    call get_option("/material_phase["//int2str(0)//"]/name", phase_name)
    do 
       !! Note that this won't work in parallel. Have to look for the pvtu in that case.
       !! not ' but "
       if(numphase .eq. 1) then
           write(filename, '(a, i0, a)') trim(simulation_name)//'_', (count)*dump_sampling_period,".vtu"   
         ! write(filename, '(a, i0, a,a)')  trim(simulation_name)//'_VelocityMesh_', (count)*dump_sampling_period,'_checkpoint',".vtu" 
       else
         ! write(filename, '(a, i0, a,a)')  trim(simulation_name)//'_'//trim(phase_name)//'_VelocityMesh_', (count)*dump_sampling_period,'_checkpoint',".vtu" 
            write(filename, '(a, i0, a)') trim(simulation_name)//'_', (count)*dump_sampling_period,".vtu" 
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

 function count_dumps_pvtu(dump_sampling_period) result (count)
    !! Work out how many dumps we're going to read in.
    integer :: count,dump_sampling_period

    logical :: exists
    !      character(len=FILE_NAME_LEN) :: filename
    character(len=1024) :: filename, phase_name
   ! numphase=option_count("/material_phase")

    count=1
    call get_option("/material_phase["//int2str(0)//"]/name", phase_name)
    do 
       !! Note that this won't work in parallel. Have to look for the pvtu in that case.
       !! not ' but "
       if(numphase .eq. 1) then
           write(filename, '(a, i0, a)') trim(simulation_name)//'_', (count)*dump_sampling_period,".pvtu"   
         ! write(filename, '(a, i0, a,a)')  trim(simulation_name)//'_VelocityMesh_', (count)*dump_sampling_period,'_checkpoint',".vtu" 
       else
         ! write(filename, '(a, i0, a,a)')  trim(simulation_name)//'_'//trim(phase_name)//'_VelocityMesh_', (count)*dump_sampling_period,'_checkpoint',".vtu" 
            write(filename, '(a, i0, a)') trim(simulation_name)//'_', (count)*dump_sampling_period,".pvtu" 
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

  end function count_dumps_pvtu

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
 
end program sensors_optimisation
