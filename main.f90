!Exercise 3b: FORTRAN program to solve the unsteady, 2D diffusion problem

PROGRAM diffusion

USE mod_diff ! contains all the declarations
USE mod_alloc! contains allocation subroutine
use mod_mpi_initialize

implicit none
!include 'mpif.h'

! define interface
INTERFACE

    SUBROUTINE file_out(Nx, Ny, dx, dy, T_new, tstep_count)
        USE mod_diff, ONLY:MK! contains allocation subroutine
        IMPLICIT none
        !integer :: MK
        INTEGER, INTENT(IN) :: Nx, Ny
        INTEGER, OPTIONAL :: tstep_count
        REAL, INTENT(IN) :: dx, dy
        REAL(MK), DIMENSION(:, :) :: T_new
        CHARACTER(LEN=20) :: filename
    END SUBROUTINE file_out
    
    ! subroutine diagnostic
    SUBROUTINE diagnostic(k, dt, nstep, T_new)
         USE mod_diff, ONLY:MK! contains allocation subroutine
        IMPLICIT none
        !integer :: MK
        INTEGER, INTENT(IN) :: k, nstep
        REAL, INTENT(IN) :: dt
        REAL:: time, min_T  
        REAL(MK), DIMENSION(:, :), INTENT(IN) :: T_new !assumed shape array
        LOGICAL :: first = .TRUE. ! saves the value for opening file
    END SUBROUTINE diagnostic

    ! subroutine update field
    ELEMENTAL SUBROUTINE elem_update_field(T_old, T_new)
        USE mod_diff, ONLY:MK! contains allocation subroutine
        IMPLICIT none
        !integer :: MK
        REAL(MK), INTENT(IN) :: T_new
        REAL(MK), INTENT(OUT) :: T_old
    END SUBROUTINE elem_update_field
    
    ! subroutine initialize
    SUBROUTINE  initialize(Lx, Ly, nstep, T_old, T_new, L, inp_file, hotstart_file, Nx, Ny, D, sim_time, nstep_start, dt, info)
        
        USE mod_diff, ONLY:MK! contains allocation subroutine
        implicit none
        !integer :: MK
        integer, intent(inout) :: nstep, nstep_start, Nx, Ny, D
        real, intent(inout) :: Lx, Ly, sim_time, dt
        integer :: Nx_tmp, Ny_tmp, info ! Nx, Ny from the hotstat file 
        real(MK), dimension(:, :), allocatable :: T_old, L, T_new
        real(MK), dimension(:,:), allocatable :: tmp_field
        character(len=*) :: inp_file, hotstart_file
        logical :: file_exists
    END SUBROUTINE initialize
    !subroutine to dave restart file
    subroutine save_restart(hotstart_file, Nx, Ny, D, sim_time, dt, itstep, T_old)
     USE mod_diff, ONLY:MK! contains allocation subroutine
    implicit none
    !integer :: MK
    character(len=*) :: hotstart_file
    integer :: Nx, Ny, D, itstep
    real :: sim_time, dt
    real(MK), dimension (:,:) :: T_old
end subroutine save_restart
   
END INTERFACE


call mpi_initialize(rank, size, name, status, ierror)


!Initialize
call initialize(Lx, Ly, nstep, T_old, T_new, L, inp_file, hotstart_file, &
                Nx, Ny, D, sim_time, nstep_start, dt, info)

! reallocate T_old with changed size (here it remains same)
!call alloc(L, T_new, T_old, Nx, Ny, info)

! set the dt, dx, dy

dx = Lx/REAL(Nx - 1) ! discrete length in x
dy = Ly/REAL(Ny - 1) ! discrete length in y

! Fourier limit check
! calculate dt_limit
dt_limit = MIN(dx,dy)**2/REAL(4*D)
!print*,'dt_limit= ', dt_limit !DEBUG
IF (dt.GE.dt_limit) THEN
    !print*, 'dt-dt_limit=',(dt-dt_limit) !DEBUG
    print*, 'WARNING! Fourier limit violated. Ensuring compliance by reducing time-step....'
    dt = dt_limit - 0.001*dt_limit !reduce by 0.1% from the dt limit
    nstep = int(sim_time/dt)

ENDIF        

! force total time steps for benchmark
!nstep = 100
!print*, 'Forcing total time steps = ',nstep

print*, 'Using the input values:' 
print*, 'sim_time=',sim_time,'[s], Nx=',Nx,&
        ', Ny=',Ny,', dt=',dt,'[s], No. of time steps=', nstep

!set the dirichlet boundary condition
!T_old(1:Nx,1) = 1.0
!T_old(1:Nx,Ny) = 1.0
!T_old(1,1:Ny) = 1.0
!T_old(Ny,1:Ny) = 1.0

!T_new(1:Nx,1) = 1.0
!T_new(1:Nx,Ny) = 1.0
!T_new(1,1:Ny) = 1.0
!T_new(Ny,1:Ny) = 1.0

! square the discrete lengths
sq_dx = dx**2
sq_dy = dy**2

!call cpu_time(cpu_t1)

!call system_clock(count=timer_start, count_rate=timer_rate)


! distribute the grid in the i-direction to different proc
! Ny
Ny_local = Ny
! Nx

if (size .ne. 1) then
    
    if (rank .ne. size-1) then
    
        Nx_local = int(Nx/size)
    else
        Nx_local = int(Nx - int(Nx/size)*(size-1))
    endif


    ! perform a ghost layer update according to the ranks
    if (rank .eq. 0 .or. rank .eq. size-1) then

        Nx_local = Nx_local + 1 ! 1 ghost layer either to the right or left
        
    else 
        Nx_local = Nx_local + 2
    
    endif

else

    Nx_local = Nx    
endif    

print*,'rank: ', rank, 'Nx: ', Nx_local

! allocate local T_old and T_new!!!!!!!!!!!!
call alloc(L, T_new_loc, T_old_loc, Nx_local, Ny_local, info)
!call alloc(L, sendmsg, recvmsg, Nx_local, Ny_local)


! initialize T_new_loc and T_old_loc!!!!!!!!!!
T_old_loc(1:Nx_local, 1:Ny_local) = 0.0
T_new_loc(1:Nx_local, 1:Ny_local) = 0.0

! enter the dirichlet boundary conditions!!!!!!!!!!!!!

if (size .ne. 1) then
    if (rank .eq. 0) then
        ! left, top and bottom
        T_old_loc(1, 1:Ny_local) = 1.0 !left wall
        T_old_loc(1:Nx_local, 1) = 1.0 !bottom wall
        T_old_loc(1:Nx_local, Ny_local) = 1.0 ! top wall
    
        ! same for the T_new_loc
        ! left, top and bottom
        T_new_loc(1, 1:Ny_local) = 1.0 !left wall
        T_new_loc(1:Nx_local, 1) = 1.0 !bottom wall
        T_new_loc(1:Nx_local, Ny_local) = 1.0 ! top wall
    
    elseif (rank .eq. size-1) then
        ! right, top and bottom
        T_old_loc(Nx_local, 1:Ny_local) = 1.0 !right wall
        T_old_loc(1:Nx_local, 1) = 1.0 !bottom wall
        T_old_loc(1:Nx_local, Ny_local) = 1.0 ! top wall
    
        ! same for the T_new_loc
        ! left, top and bottom
        T_new_loc(Nx_local, 1:Ny_local) = 1.0 !right wall
        T_new_loc(1:Nx_local, 1) = 1.0 !bottom wall
        T_new_loc(1:Nx_local, Ny_local) = 1.0 ! top wall
    
    else
        ! top and bottom
        T_old_loc(1:Nx_local, 1) = 1.0 !bottom wall
        T_old_loc(1:Nx_local, Ny_local) = 1.0 ! top wall
    
        ! same for the T_new_loc
        ! left, top and bottom
        T_new_loc(1:Nx_local, 1) = 1.0 !bottom wall
        T_new_loc(1:Nx_local, Ny_local) = 1.0 ! top wall
    
    endif

else
        ! left, top and bottom
        T_old_loc(1, 1:Ny_local) = 1.0 !left wall
        T_old_loc(1:Nx_local, 1) = 1.0 !bottom wall
        T_old_loc(Nx_local, 1:Ny_local) = 1.0 ! right wall
        T_old_loc(1:Nx_local, Ny_local) = 1.0 ! top wall
    
        ! same for the T_new_loc
        ! left, top and bottom
        T_new_loc(1, 1:Ny_local) = 1.0 !left wall
        T_new_loc(1:Nx_local, 1) = 1.0 !bottom wall
        T_new_loc(Nx_local, 1:Ny_local) = 1.0 ! right wall
        T_new_loc(1:Nx_local, Ny_local) = 1.0 ! top wall

endif    
! end Dirichlet boundary condition assignment  !!!!!!!!!!!!!!!  

!euler time integration
DO k=nstep_start,nstep
    
    ! Saving a hotstart file for T_old at each time-step
    !call save_restart(hotstart_file, Nx, Ny, D, sim_time, dt, k, T_old)

    ! condition to exit the do loop if tot time exceeds sim_time
    tot_time = k*dt 
    if (tot_time > sim_time) then 
        print*, 'Stopping simulation....Total simulation time exceeded!'
        exit    
    endif


!! NOTES: A] The local field arrays need to be updated in the following way:
!!           1) T_new value from Nx-1 of left domain ==> T_new value of 0 of right domain
!!           2) T_new value from 1 of right domain ==> T_new value of NX of left domain
!!        B] The messgae passing can happen only after T_new has been calculated as
!!           the T_old = T_new update will take care of the values for T_old
!!        C] T_old update isn't required when it is the first time-step as everything starts from 0.0
    
   ! send and receive T_old
    if (size .ne. 1) then

        if (rank .eq. 0) then
            count = Ny_local

            dest = 1
            ! send the T(Nx_local -1, :) to rank 1
            call MPI_Send(T_old_loc(Nx_local-1, 1:Ny_local), count, MPI_DOUBLE_PRECISION,&
                           dest, tag, MPI_COMM_WORLD, ierror)

            ! receive into T(Nx_local,:) from rank 1
            src = 1
            call MPI_Recv(T_old_loc(Nx_local, 1:Ny_local), count, MPI_DOUBLE_PRECISION,&
                           src, tag, MPI_COMM_WORLD, status, ierror)

        elseif(rank .eq. size-1) then

            count = Ny_local

            dest = rank - 1
            ! send the T(2, :) to rank 1
            call MPI_Send(T_old_loc(2, 1:Ny_local), count, MPI_DOUBLE_PRECISION,&
                           dest, tag, MPI_COMM_WORLD, ierror)

            ! receive into T(1,:) from rank 1
            src = rank - 1
            call MPI_Recv(T_old_loc(1, 1:Ny_local), count, MPI_DOUBLE_PRECISION,&
                           src, tag, MPI_COMM_WORLD, status, ierror)
            

        else

            count = Ny_local
            dest = rank - 1
            src = rank + 1
            ! send T(2,:) and receive T(Nx_local,:)
            call MPI_SENDRECV(T_old_loc(2,1:Ny_local), count, MPI_DOUBLE_PRECISION, dest, sendtag,&
                              T_old_loc(Nx_local, 1:Ny_local), count, MPI_DOUBLE_PRECISION,&
                              src, recvtag, MPI_COMM_WORLD, status, ierror)


            dest = rank + 1
            src = rank - 1
            ! send T(Nx_local-1,:) and receive T(1,:)
            call MPI_SENDRECV(T_old_loc(Nx_local-1,1:Ny_local), count, MPI_DOUBLE_PRECISION, dest, sendtag,&
                              T_old_loc(1, 1:Ny_local), count, MPI_DOUBLE_PRECISION,&
                              src, recvtag, MPI_COMM_WORLD, status, ierror)

        endif
    endif





    DO j=2, Ny_local-1
        DO i=2,Nx_local-1
            !laplacian
            laplacian = (T_old_loc(i+1,j) - 2*T_old_loc(i,j) + T_old_loc(i-1,j))/sq_dx + &
                        (T_old_loc(i,j+1) - 2*T_old_loc(i,j) + T_old_loc(i,j-1))/sq_dy
            
            !update T_new_loc
            T_new_loc(i,j) = D*laplacian*dt + T_old_loc(i,j)      

            ! write the ghost layer update sequence !!!!!!!!!!!!!!!!!!!!!!!
            
            !if (size .ne. 1) then
            !
            !    if (rank.eq.0) then
            !    
            !    
            !        if ( i.eq.Nx_local-1) then ! first rank
            !            ! send Nx-1 to rank+1
            !            !msg = T_new_loc(i, j)
            !            count = 1
            !            dest = rank + 1
            !            call MPI_Send( T_new_loc(i, j), count, MPI_DOUBLE_PRECISION, dest,&
            !                          tag, MPI_COMM_WORLD, ierror)
            !    
            !        !elseif ( i.eq.Nx_local) then
            !            ! receive from rank+1
            !            src = rank + 1
            !            count = 1
            !    
            !            call MPI_Recv(T_new_loc(i+1,j), count, MPI_DOUBLE_PRECISION, src, tag,&
            !                          MPI_COMM_WORLD, status, ierror)
            !    
            !            !T_new_loc(i+1,j) = msg
            !    
            !        endif                                                  
            !    
            !    elseif (rank .eq. size-1) then ! last rank
            !    
            !        if ( i.eq.2) then
            !            ! send Nx-1 to rank+1
            !            !msg = T_new_loc(i, j)
            !            count = 1
            !            dest = rank - 1
            !            call MPI_Send(T_new_loc(i, j), count, MPI_DOUBLE_PRECISION,&
            !                           dest, tag, MPI_COMM_WORLD, ierror)
            !    
            !        !elseif ( i.eq.0) then
            !            ! receive from rank+1
            !            src = rank - 1
            !            count = 1
            !    
            !            call MPI_Recv(T_new_loc(i-1,j), count, MPI_DOUBLE_PRECISION, &
            !                          src, tag, MPI_COMM_WORLD, status, ierror)
            !    
            !            !T_new_loc(i-1,j) = msg
            !    
            !        endif                                                  
            !    
            !    
            !    else ! middle ranks use sendrecv
            !        
            !        if (i.eq.2) then  
            !            !msg_send = T_new_loc(i,j)
            !            dest = rank - 1
            !            src = rank + 1
            !            count_send = 1
            !            count_recv = 1
            !    
            !            call MPI_SENDRECV(T_new_loc(i,j), count_send, MPI_DOUBLE_PRECISION, dest, sendtag,&
            !                              T_new_loc(Nx_local, j), count_recv, MPI_DOUBLE_PRECISION,&
            !                              src, recvtag, MPI_COMM_WORLD, status, ierror)
            !    
            !            !T_new_loc(Nx_local, j) = msg_recv
            !    
            !        elseif (i.eq.Nx_local-1) then
            !    
            !            !msg_send = T_new_loc(i,j)
            !            dest = rank + 1
            !            src = rank - 1
            !            count_send = 1
            !            count_recv = 1
            !    
            !            call MPI_SENDRECV(T_new_loc(i,j), count_send, MPI_DOUBLE_PRECISION, dest, sendtag, &
            !                              T_new_loc(1, j), count_recv, MPI_DOUBLE_PRECISION,&
            !                              src, recvtag, MPI_COMM_WORLD, status, ierror)
            !    
            !            !T_new_loc(1, j) = msg_recv
            !    
            !        endif                                                  
            !    
            !    endif    
            !    !! END THE GHOST lAYER UPDATE SEQUENCE!!!!!!!!!!!!!!!
            ! endif

        ENDDO
        ! end the do loop over i
    ENDDO
    ! end the do loop over j

    ! print diagnostic
    !if (mod(k,10)==0) then
    !    call diagnostic(k, dt, nstep, T_new)
    !endif
    ! call optional argument and write field at each step
    !call file_out(Nx, Ny, dx, dy, T_new, k)

    !update T_old_loc
    call elem_update_field(T_old_loc, T_new_loc)
    
    ! place a barrier before moving the iteration to the next time-step
    call MPI_Barrier(MPI_COMM_WORLD,ierror) 

ENDDO
print*,'loop done'

!!!GATHER ALL THE T_local into T_global!!!!!!!!!!!!!!

! print the field for each rank

print*,'T for rank ', rank
do j = 1, Ny
    print*, T_new_loc(:, j)
enddo

! creates the sub-array type               
if (rank .eq. 0) then

   !create the new sub-array type
   ndim = 2 ! dimensions of the array
   sizes = (/Nx,Ny/) ! size of the global array
   ! defining sub sizes for the incoming arrays sub-array
   size_norm = int(Nx/size)
   subsizes_norm = (/size_norm, Ny_local/) ! size of the local arrays   

   
   starts = (/0,0/) ! address where the array assingment starts in the
               ! global array.

    ! create a type for the rest of the ranks
    call MPI_Type_create_subarray(ndim, sizes, subsizes_norm, starts, &
                                  MPI_ORDER_FORTRAN, MPI_Double_Precision, &
                                 newtype_norm, ierror)

    call MPI_Type_commit(newtype_norm, ierror)


   if (size .ne. 1) then    

       size_last = int(Nx - int(Nx/size)*(size-1))   
       subsizes_last = (/size_last, Ny_local/) ! size of the local arrays

        ! create a type for the last rank
       call MPI_Type_create_subarray(ndim, sizes, subsizes_last, starts, &
                                     MPI_ORDER_FORTRAN, MPI_Double_Precision, &
                                    newtype_last, ierror)
       call MPI_Type_commit(newtype_last, ierror)
   endif

endif

!call file_out(Nx_local, Ny_local, dx, dy, T_new_loc)
! barrier to ensure that rank=0 creates the data types
call MPI_Barrier(MPI_COMM_WORLD, ierror)

! start sending the arrays from various ranks to root
if (size .ne. 1) then

    if (rank .ne. 0) then
    ! send all the arrays to root
        dest = 0
        
        if (rank .ne. size-1) then
    
            count = (Nx_local-2)*Ny_local
            
            call MPI_Send(T_new_loc(2:Nx_local-1, 1:Ny), count, &
                          MPI_Double_Precision, dest, &
                          tag, MPI_COMM_WORLD, ierror)
        else
            
            count = (Nx_local-1)*Ny_local
            
            call MPI_Send(T_new_loc(2:Nx_local, 1:Ny), count, &
                          MPI_Double_Precision, dest, &
                          tag, MPI_COMM_WORLD, ierror)
        endif        
    
    else
        ! receive all the arrays in the root
        do i=1, size-1
            
            indi = int(Nx/size)*i + 1
            indj = 1
    
            if (i .ne. size-1) then
    
                call MPI_Recv(T_new(indi, indj),1,newtype_norm,i,tag,&
                              MPI_COMM_WORLD, status, ierror)
    
            else
                call MPI_Recv(T_new(indi, indj),1,newtype_last,i,tag,&
                              MPI_COMM_WORLD, status, ierror)
    
            endif
        enddo
    
        ! call cpu clock
        !call cpu_time(cpu_t2)
        !call system_clock(count=timer_stop)
        !cputime_timestep = (cpu_t2 - cpu_t1)/real(k-1)
        !e_time = real(timer_stop - timer_start)/timer_rate
    
        !print*, 'CPU TIME= ', cputime_timestep*100,'[s]'
        !print*, 'CPU TIME per time-step = ', cputime_timestep
        !print*, 'Wall time = ', e_time, '[s]'
    
        !print*, 'total time steps = ',(k-1)
        !print*, 'total time = ',tot_time,'[s]'
    
    ! write the final field
         !call file_out(Nx, Ny, dx, dy, T_new)
    
    endif    
endif
    ! fill in the rank=0 matrix
! wait for the root to write the file
call MPI_Barrier(MPI_COMM_WORLD, ierror)

if (rank .eq. 0) then

    indi = int(Nx/size)
    T_new(1:indi, 1:Ny) = T_new_loc(1:Nx_local-1, 1:Ny)

    ! write the final field
    call file_out(Nx, Ny, dx, dy, T_new)

endif    



! wait for the root to write the file
call MPI_Barrier(MPI_COMM_WORLD, ierror)
! terminate MPI
call MPI_Finalize(ierror)

END PROGRAM diffusion










