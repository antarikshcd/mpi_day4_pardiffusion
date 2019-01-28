!module containing the global data for exercise 2b diffusion problem
module mod_diff
    implicit none
    INTEGER, PARAMETER :: MK = kind(1.0D0)
    INTEGER :: Nx, Ny, D, Nx_local, Ny_local
    REAL :: sim_time ! total simulation time in [s]
    REAL(MK), DIMENSION(:, :), allocatable :: T_new, T_new_loc
    REAL(MK), DIMENSION (:,:), allocatable :: T_old, T_old_loc
    real(mk) :: laplacian, rsq_dx, rsq_dy ! laplacian vector
    REAL :: dx,dy,Lx,Ly,dt,dt_limit, tot_time
    INTEGER ::i,j,k,nstep, info, nstep_start
    logical :: file_exists
    character(len=256) :: inp_file, hotstart_file
    character(8) :: date ! variable for storing date
    character(10) :: time ! variable for storing time
    integer :: timer_start, timer_stop ! system clock counters
    integer :: timer_rate ! system clock count_rate: depends on the system
    real :: e_time ! elapsed time measuring time for the double do loops
    real :: cpu_t1, cpu_t2, cputime_timestep ! cpu times
    
    ! MPI variables  
    integer :: rank, size, ierror, ierror, tag, &
               src, dest, count, count_send, count_recv, &
               sendtag, recvtag
    integer, dimension(:), allocatable :: status              
    !real(mk) :: msg, msg_send, msg_recv               
    character(len=256) :: name
    !real(mk), DIMENSION(:,:),allocatable :: sendmsg, recvmsg
    !sub-array
    real(mk) :: newtype_norm, newtype_last
    integer, dimension(2) :: sizes, subsizes_norm, subsizes_last, starts
    integer :: ndim, indi, indj, size_norm, size_last

end module mod_diff

