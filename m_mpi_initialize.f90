module mod_mpi_initialize
    use mod_diff, only:mk
    include 'mpif.h'

    contains
        subroutine mpi_initialize(rank, size, name, status, ierror)
            implicit none
            integer :: rank, size, ilen, ierror
            character(len=256) :: name
            integer :: info
            integer, dimension(:), allocatable :: status


            ! initialize MPI
            call MPI_Init(ierror)
            !print*, 'MPI initializion error signal: ',ierror
            
            ! print the rnak of MPI_COMM_WORLD
            call MPI_Comm_RANK(MPI_COMM_WORLD, rank, ierror)
            !print*, 'MPI rank: ', rank, 'error status: ',ierror
            
            ! print the total number of sizes in the communicator
            call MPI_Comm_Size(MPI_COMM_WORLD, size, ierror)
            !print*, 'MPI size: ', size, 'error status: ',ierror
            
            ! print the processor name
            call MPI_Get_Processor_Name(name, ilen, ierror) 
            !print*,'MPI processor name', name(1:ilen), 'error status: ',ierror

            allocate(status(MPI_STATUS_SIZE), stat=info)

        end subroutine mpi_initialize

end module mod_mpi_initialize
