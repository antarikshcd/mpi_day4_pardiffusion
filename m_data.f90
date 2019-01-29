!module containing the global data for exercise 2b diffusion problem
module mod_diff
    implicit none
    INTEGER, PARAMETER :: MK = kind(1.0D0)
    INTEGER :: Nx, Ny, D, Nx_local, Ny_local, Nx_serial, Ny_serial
    REAL :: sim_time ! total simulation time in [s]
    REAL(MK), DIMENSION(:, :), allocatable :: T_new, T_new_loc
    REAL(MK), DIMENSION (:,:), allocatable :: T_old, T_old_loc
    
    ! load the saved serial counterpart 
    real(mk), dimension(:, :), allocatable :: T_serial
    
    real(mk) :: laplacian, rsq_dx, rsq_dy, diff, sum, rms ! laplacian vector
    REAL :: dx,dy,Lx,Ly,dt,dt_limit, tot_time
    INTEGER ::i,j,k,nstep, info, nstep_start, nstep_serial
    logical :: file_exists,flag_save
    character(len=256) :: inp_file, hotstart_file
    
    ! time variables
    real(mk) :: t1, t2, wall_time  

    ! MPI variables  
    integer :: rank, size, ierror, ierror, tag, &
               src, dest, count, count_send, count_recv, &
               sendtag, recvtag
    integer, dimension(:), allocatable :: status              
    !real(mk) :: msg, msg_send, msg_recv               
    character(len=256) :: name, savename, load_file
    !real(mk), DIMENSION(:,:),allocatable :: sendmsg, recvmsg
    !sub-array
    real(mk) :: newtype_norm, newtype_last
    integer, dimension(2) :: sizes, subsizes_norm, subsizes_last, starts
    integer :: ndim, indi, indj, size_norm, size_last

end module mod_diff

