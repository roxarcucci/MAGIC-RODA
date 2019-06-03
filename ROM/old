 






Scalar fields: 
 +phase1::Pressure (phase1::Pressure)  on Mesh
 +phase1::Density (phase1::Density)  on Mesh
 +phase1::PhaseVolumeFraction (phase1::PhaseVolumeFraction)  on Mesh
 +phase2::Density (phase2::Density)  on Mesh
 +phase2::PhaseVolumeFraction (phase2::PhaseVolumeFraction)  on Mesh
 +phase1::Porosity (phase1::Porosity)  on P0Mesh
 +phase1::Permeability (phase1::Permeability)  on P0Mesh
Vector fields: 
 +Coordinate (Coordinate)  on Mesh
 +phase1::Velocity (phase1::Velocity)  on Mesh
 +phase1::VelocityAbsorption (phase1::VelocityAbsorption)  on Mesh
 +phase2::Velocity (phase2::Velocity)  on Mesh
 +phase2::VelocityAbsorption (phase2::VelocityAbsorption)  on Mesh
Tensor fields: 
 +phase1::Viscosity (phase1::Viscosity)  on Mesh
 +phase2::Viscosity (phase2::Viscosity)  on Mesh





     subroutine non_intrusive_rbf_main() 
        implicit none
!
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
 
    type(state_type), dimension(:), pointer :: state!=> null()
    !type( state_type ), dimension( : ), intent( inout ) :: state
    integer :: i, dump_no = 0
    integer :: ntsol, nonlinear_iterations
    character(len = option_path_len) :: filename
    character(len = option_path_len) :: simulation_name, dump_format
    real ::   nonlinear_iteration_tolerance
    
    ! type(state_type), dimension(:), allocatable :: state_adj
    type(state_type), dimension(:,:,:), allocatable :: POD_state
    integer :: ii, kj, numphase,nsvd,romdim,ph
    integer :: ierr,stat,j,k ,ij 
    integer :: total_timestep,total_dump_no
    type(vector_field), pointer :: snapmean_velocity
    type(scalar_field), pointer :: snapmean_pressure
    type(vector_field), pointer :: u ,POD_velocity
    type(scalar_field), pointer :: p,POD_pressure
    real, dimension(:), allocatable :: pod_coef_optimised
    logical :: if_optimal,if_timetwo
    ! type(state_type), dimension(:,:,:) :: POD_state 
    ! Change in velocity
    !type(vector_field) :: delta_u
    REAL, DIMENSION(:,:,:), allocatable :: delta_u
    real, dimension(:,:), allocatable :: pod_sol_velocity
    real, dimension(:), allocatable :: pod_sol_pressure,pod_sol_volumefraction
    type(vector_field), pointer :: velocity_field!,volumefraction_field
    ! non_intrusive_ROM
    type(state_type), dimension(:,:,:), allocatable :: Non_Intrusive_state
    real, dimension(:), allocatable ::  pod_coef_orig,pod_coef_new,pod_coef_old,test
    
    integer istate, ntime
    integer dump_no_tmp, int_dump_period
    real current_time,finish_time,dt
    !smolyak
    real:: cpu_start1,cpu_end1,cpu_start2,cpu_end2, project_start, project_end
    type(scalar_field) :: delta_p, delta_vf
  real ( kind = 8 ), allocatable :: a(:)
  real ( kind = 8 ) app_error,app_error1
  real ( kind = 8 ), allocatable :: b(:)
  integer ( kind = 4 ), allocatable :: c(:)
  integer ( kind = 4 ) h
  integer ( kind = 4 ), allocatable :: ind(:)
  integer ( kind = 4 ) l
  integer ( kind = 4 ) l_max
  integer ( kind = 4 ) l_min
  real ( kind = 8 ), allocatable :: fi(:) 
  logical more
  integer ( kind = 4 ) nd,nd_lq
  integer ( kind = 4 ) ni
  ! real ( kind = 8 ) r8vec_norm_affine
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) sparse_max
  integer ( kind = 4 ) sparse_min
  integer ( kind = 4 ) t
  real, allocatable :: w(:)
  real, allocatable :: wm(:,:)
  real :: r0
  real , allocatable :: xd(:,:),xd_lq(:,:), pod_coef_all_obv(:,:), smoylak_all(:,:), pod_coef_all_obv2(:,:), pod_coef_all_obv3(:,:),pod_coef_obv_vpv(:,:)
  real , allocatable :: xi(:,:)
  real , allocatable :: zd(:),fd(:)
  real , allocatable :: ze(:),temp_coef(:)
  real , allocatable :: zi(:),minin(:),maxin(:)
  real , allocatable :: zpi(:)
  real , allocatable :: fe(:,:) !exact values of function produced from fluidity
  integer :: m,ll,kk,u_nonods
  !real ( kind = 8 ) w_lq(total_timestep)
  !real ( kind = 8 ) fd(nd_lq)
  integer ::  podnum
  real :: mean, positive,negtive,num_pos,num_neg
  type(scalar_field), pointer :: pressure_field, volumefraction_field 
  type(scalar_field), pointer :: snapmean_volumefraction
  
  real :: maxall,minall
  call get_option('/simulation_name', simulation_name)
  print *, 'simulation_name' , simulation_name,  len_trim(simulation_name)
  ! if (have_option("/reduced_model/execute_reduced_model")) then
    if (have_option("/reduced_model/Smolyak") .or. have_option("/reduced_model/RBF_interpolation")) then
       simulation_name(len_trim(simulation_name)-3:len_trim(simulation_name))="" !delete _POD in the name
    endif
  print *, 'simulation_name=' , simulation_name
  ! real minin(podnum),maxin(podnum) 
  call get_option(&
            '/reduced_model/pod_basis_formation/pod_basis_count', nsvd) 
  call delete_option("/reduced_model/execute_reduced_model") ! switch to full model 
  call get_option("/timestepping/current_time", current_time)
  call get_option("/timestepping/finish_time", finish_time)       
  call get_option("/timestepping/timestep", dt)
  call get_option('/geometry/dimension/', romdim) 
   u => extract_vector_field(state(1), "Velocity", stat)
   p => extract_scalar_field(state(1), "Pressure") 
   numphase=option_count("/material_phase")
    ! numphase=1
      if(have_option("/io/dump_period_in_timesteps/constant")) then
       call get_option("/io/dump_period_in_timesteps/constant", int_dump_period)
    else
       int_dump_period = 1
    endif
   ! numphase=1
   ntime=int((finish_time-current_time)/dt)
   total_timestep=int((finish_time-current_time)/dt)
   total_timestep=total_timestep/int_dump_period-1
   if (have_option("/reduced_model/Velocity_Interpolation")) then
      ! m = (romdim+2)*nsvd*numphase
       m = (romdim+1)*nsvd*numphase
       allocate(pod_coef_all_obv3(total_timestep,m*2)) 
    else 
       m = 2*nsvd*numphase
   endif
   podnum =m    
   allocate(pod_coef_all_obv(total_timestep,podnum))   
   allocate(pod_coef_all_obv2(total_timestep,2*nsvd*numphase)) 
   allocate(pod_coef_obv_vpv(total_timestep,(romdim+2)*nsvd*numphase)) 
   allocate(test(m)) 
    open(unit=61,file="test")
    !   read(61,*) test(:) !one phase
    !   print *, 'test', test
    close(61)

   if (have_option("/reduced_model/Velocity_Interpolation")) then 
      open(unit=61,file="coef_pod_all_obv"//"_"//trim(simulation_name))
    
       read(61,*)((pod_coef_all_obv(j,k),k=1,podnum),j=1,total_timestep)
     close(61)
     !  pod_coef_all_obv(:,:)=pod_coef_all_obv3(:,1:podnum) !one phase
   else
    open(unit=61,file="coef_pod_all_obv"//"_"//trim(simulation_name))
    
      read(61,*)((pod_coef_all_obv(j,k),k=1,2*nsvd*numphase),j=1,total_timestep) !read(61,*)((pod_coef_all_obv(j,k),k=1,2*nsvd),j=1,total_timestep)
     
    close(61)
    
    endif
  
     Do i=1, size(pod_coef_all_obv,1)-1
     pod_coef_all_obv(i,:)=pod_coef_all_obv(i+1,:) ! because some propblem produced for the first timestep's coef in multiphase. 
     enddo
   ! 
     Do i=1,m   ! cannot use minval because some of the vector is null.          
          minall=pod_coef_all_obv(1,1)
          maxall=pod_coef_all_obv(1,1)      
     do j=1,total_timestep
           if(minall>pod_coef_all_obv(j,i)) then
              minall=pod_coef_all_obv(j,i)          
           endif
           if(maxall<pod_coef_all_obv(j,i)) then
              maxall=pod_coef_all_obv(j,i)
           endif
      enddo
   ENDDO 
  
   
   ni = 1
   nd = total_timestep-1
   allocate ( xi(1:m,1:ni) )
   allocate ( fi(1:ni) ) 
   allocate (fd(nd))
   allocate (w(nd))
   allocate (wm(nd,m))
   allocate (xd(podnum,nd)) 
   u_nonods = node_count(u) 
   if (have_option("/reduced_model/multi_cases")) then
       call read_pod_basis_mp_allshots(POD_state, state) 
    else
       call read_pod_basis_mp(POD_state, state)
   endif
    allocate(pod_coef_orig(m))
    allocate(pod_coef_new(m)) 
    allocate(pod_coef_old(m)) 
    ! Allocate the change in pressure field
    call allocate(delta_p, p%mesh, "DeltaP")
    delta_p%option_path = trim(p%option_path)
    call zero(delta_p)
    
     call allocate(delta_vf, p%mesh, "DeltaVF")
     !delta_vf%option_path = trim(vf%option_path)
     call zero(delta_vf)
    ! get the initial ROM solution coef.
    !----------------------------------
    istate = 1
    if (have_option("/reduced_model/Velocity_Interpolation")) then
      do ph=1,numphase 
       ! call project_from_full_to_pod_mp(ph,pod_state, state, pod_coef_orig((ph-1)*(romdim+2)*nsvd+1:ph*(romdim+2)*nsvd))!pod_coef_orig) 
        ! call project_from_full_to_pod_mp_vp(ph,pod_state, state, pod_coef_orig)!pod_coef_orig)   
        ! pod_coef_old((ph-1)*2*nsvd+1:ph*2*nsvd)=pod_coef_orig
        pod_coef_old=pod_coef_orig
      enddo  
    else
     do ph=1,numphase 
      call project_from_full_to_pod_mp_pv(ph,pod_state, state, pod_coef_orig)        
      ! stop 111   
      pod_coef_old((ph-1)*2*nsvd+1:ph*2*nsvd)=pod_coef_orig(1:2*nsvd)
      ! print *,'initial pod_coef', pod_coef_old  
     enddo  
    endif 

       do j=1,m
         do k=1,nd
           xd(j,k)=pod_coef_all_obv(k,j)
         enddo
       enddo
     r0 = ( maxall - minall ) / real ( nd, kind = 8 )
    
    pod_coef_old=pod_coef_all_obv(1,:) !test the initial condition
    call cpu_time(cpu_start1)
     do l=1, m
        do k=1,nd
          fd(k)=pod_coef_all_obv(k+1,l)  ! target function value is the l^th next timestep's pod_coefficient  
          
        enddo 
      call rbf_weight ( m, nd, xd, r0, phi1, fd, wm(:,l) )
     enddo
   call cpu_time(cpu_end1)
    print *,'weighting',cpu_end1-cpu_start1
    do kk =1, total_timestep !ntime!,int_dump_period 
       
         xi(:,1)=  pod_coef_old !pod_coef_all_obv(kk,:)
     
     do l=1, m
      
       call cpu_time(cpu_start1) 
      call rbf_interp_nd ( m, nd, xd, r0, phi1, wm(:,l), ni, xi, fi )
      pod_coef_new(l)=fi(1)  
     enddo  !do l=1, m 
  
        ALLOCATE(delta_u(romdim, numphase, u_nonods))
        delta_u=0
      do ph=1, numphase 
        if(have_option("/reduced_model/Velocity_Interpolation")) then  !interpolate velocity and pressure
        pod_coef_old=pod_coef_new  
        call project_full_mp_vp(delta_u(:,ph,:), delta_p,  pod_sol_velocity, pod_sol_pressure,  POD_state(:,:,ph), pod_coef_new) 
        snapmean_velocity=>extract_vector_field(POD_state(1,1,ph),"SnapmeanVelocity")
        snapmean_pressure=>extract_Scalar_field(POD_state(1,2,ph),"SnapmeanPressure")
        snapmean_volumefraction=>extract_Scalar_field(POD_state(1,3,ph),"SnapmeanVolumefraction")
        pressure_field => extract_scalar_field( state(ph), "Pressure" )
             
        velocity_field=>extract_vector_field(state(ph),"Velocity")
   
        velocity_field%val=snapmean_velocity%val
        velocity_field%val=velocity_field%val+delta_u(:,ph,:)
    
        print *, 'size: 2, 1, 1152, 1', size(snapmean_pressure%val), size(pressure_field%val )        
        pressure_field%val=snapmean_pressure%val
        call addto(pressure_field, delta_p)
    
         
        call insert(state(ph), pressure_field, "pressure")
   
        call insert(state(ph), velocity_field, "Velocity")
  

         else          

        snapmean_velocity=>extract_vector_field(POD_state(1,1,ph),"SnapmeanVelocity")
        snapmean_pressure=>extract_Scalar_field(POD_state(1,2,ph),"SnapmeanPressure")
        snapmean_volumefraction=>extract_Scalar_field(POD_state(1,3,ph),"SnapmeanVolumefraction")
         call cpu_time(project_start)
         if (have_option("/reduced_model/interpolate_podbasis")) then
          call project_full_mp_pv_indep_pod( delta_p, delta_vf,pod_sol_velocity, pod_sol_pressure, pod_sol_volumefraction, POD_state(:,:,ph), pod_coef_new)
            open(19,file='snapmean_interp') 
           
               do j=1,romdim              
                 read(19,*)(snapmean_velocity%val(j,i),i=1,size(snapmean_velocity%val,2)) 
               enddo                       
               read(19,*)(snapmean_pressure%val(i),i=1,size(snapmean_pressure%val))
               read(19,*)(snapmean_volumefraction%val(i),i=1,size(snapmean_volumefraction%val))
           
        close(19)
      
        else 
          call project_full_mp_pv( delta_p, delta_vf,pod_sol_velocity, pod_sol_pressure, pod_sol_volumefraction, POD_state(:,:,ph), pod_coef_new)
        endif
       call cpu_time(project_end)
       print *, 'projection_time', project_end-project_start
       pressure_field => extract_scalar_field( state(ph), "Pressure" )
       volumefraction_field=> extract_scalar_field(state(ph), "PhaseVolumeFraction" )
       print *, 'size: 2, 1, 1152, 1', size(snapmean_pressure%val), size(pressure_field%val ) 
       
       pressure_field%val=snapmean_pressure%val
       call addto(pressure_field, delta_p)
       volumefraction_field%val(:)=snapmean_volumefraction%val(:)
       volumefraction_field%val(:)=volumefraction_field%val(:)+ delta_vf%val(:)
       call insert(state(ph), pressure_field, "pressure")
        !volumefraction_field%val(1,2,:)=1-volumefraction_field%val(1,1,:)
       call insert(state(ph), volumefraction_field, "PhaseVolumeFraction")

       endif     !if(have_option("/reduced_model/Velocity_Interpolation")) then  
       !call print_state(state(1))
       dump_no=kk+2!/int_dump_period 
      enddo ! do ph=1, numphase

       call cpu_time(cpu_end1)
       print *,'cpurbfint',cpu_end1-cpu_start1
       call calculate_diagnostic_variables( state, exclude_nonrecalculated = .true. )
       call calculate_diagnostic_variables_new( state, exclude_nonrecalculated = .true. )
       call write_state(dump_no, state)
       deallocate(delta_u) 
     
    enddo !kk,timestep loop
     
    !close(20)
    !----------------------------------
    ! deallocate variables
    !----------------------------------
   if (allocated(pod_state)) then
       do i=1, size(pod_state,1)
          do j=1,size(pod_state,2)
             do k=1,size(pod_state,3)
                call deallocate(pod_state(i,j,k))
             enddo
          enddo
       end do
    end if
     deallocate ( fi )
     deallocate ( xi ) 
     deallocate(fd)
     deallocate(w)
     deallocate(wm)
     deallocate(xd)
    deallocate(pod_coef_orig)
    deallocate(pod_coef_new)
    deallocate(test)
    deallocate(pod_coef_old)
    deallocate(pod_coef_all_obv)
    deallocate(pod_coef_all_obv2) 
    deallocate(pod_coef_obv_vpv)
  
 end subroutine non_intrusive_rbf_main


  subroutine project_back_interp_pv( delta_p, delta_vf,pod_sol_velocity, pod_sol_pressure, pod_sol_volumefraction, POD_state, pod_coef)

        !type(vector_field), intent(inout) :: delta_u
       ! REAL, DIMENSION(:,:) :: delta_u
        type(scalar_field), intent(inout) :: delta_p, delta_vf
        real, dimension(:,:), intent(out), allocatable :: pod_sol_velocity
        real, dimension(:), intent(out), allocatable :: pod_sol_pressure, pod_sol_volumefraction
        type(state_type), dimension(:,:), intent(in) :: POD_state
        real, dimension(:), intent(in) :: pod_coef

        type(vector_field), pointer :: POD_velocity
        type(scalar_field), pointer :: POD_pressure, POD_volumefraction

        type(vector_field), pointer :: snapmean_velocity
        type(scalar_field), pointer :: snapmean_pressure,snapmean_volumefraction

        integer :: d, POD_num, i, u_nodes, p_nodes, stat,ph,k

        type(mesh_type), pointer :: pod_umesh, pod_pmesh, pod_vmesh
        


        
        !POD_pressure=>extract_scalar_field(POD_state(1,2), "PODPressure")
        !POD_volumefraction=>extract_scalar_field(POD_state(1,3), "PODVolumefraction") 
        !snapmean_pressure=>extract_scalar_field(POD_state(1,2), "SnapmeanPressure")
        ! snapmean_volumefraction=>extract_scalar_field(POD_state(1,3), "SnapmeanVolumefraction")

        u_nodes=node_count(POD_velocity)
        p_nodes=node_count(POD_pressure)

        !allocate(pod_sol_velocity(u_nodes,POD_velocity%dim))
        allocate(pod_sol_pressure(p_nodes))
        allocate(pod_sol_volumefraction(p_nodes))
        pod_sol_velocity=0.0
        pod_sol_pressure=0.0
        pod_sol_volumefraction=0.0
        POD_num=size(POD_state,1)

         do k=1, numphase          
           do i=1,POD_num 
         snapmean_pressure=>extract_Scalar_field(POD_state(i,2,k),"SnapmeanPressure")
         snapmean_volumefraction=>extract_Scalar_field(POD_state(i,3,k),"SnapmeanVolumefraction")                           
         POD_pressure=>extract_scalar_field(POD_state(i,2,k), "PODPressure")
         POD_volumefraction=>extract_scalar_field(POD_state(i,3,k), "PODVolumefraction")
       
         !  POD_pressure=>extract_scalar_field(POD_state(i,2), "PODPressure")
          ! POD_volumefraction=>extract_scalar_field(POD_state(i,3), "PODVolumefraction")
        
           
           pod_sol_pressure(:)=pod_sol_pressure(:)+pod_coef(i)*POD_pressure%val(:)
           pod_sol_volumefraction(:)=pod_sol_volumefraction(:)+pod_coef(i+POD_num)*POD_volumefraction%val(:)
            
          ! pod_umesh => extract_mesh(pod_state(1,1,1), "Mesh", stat)
           pod_pmesh => extract_mesh(pod_state(1,2,1), "Mesh", stat)
           pod_vmesh=> extract_mesh(pod_state(1,3,1), "Mesh", stat)
 
           !  delta_u=0 ! no use here
           call set_all(delta_p, pod_sol_pressure(:))
           call set_all(delta_vf, pod_sol_volumefraction(:))
           !print *,'podbasis_velocity,pressure,volumn==', POD_velocity%val(d,1), POD_pressure%val(1),POD_volumefraction%val(1)
          ! print *,'delta_u(d,:),delta_p,delta_vf==',  delta_p%val(1),delta_vf%val(1)
        enddo
        enddo
          
      end subroutine project_back_interp_pv

