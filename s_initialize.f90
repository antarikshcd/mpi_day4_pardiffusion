!intitalization subroutine
subroutine  initialize(Lx, Ly, nstep, T_old, T_new, L, inp_file, hotstart_file, Nx, Ny, D, sim_time, nstep_start, dt, info)
    use mod_alloc
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
    
    ! inquire if file exists. If it does call the read_input subroutine
    ! NOTE: it is possibel that all the inputs are not given in the 
    ! input file and one or more are missing. If that's the case the default
    ! values that are intialized are carried forward.
    !default data
    inp_file = 'input.in'
    hotstart_file = 'hotstart.bck'
    Nx=21
    Ny=21
    D=1
    ! initialize constants
    Lx = 1.0 ! length in x
    Ly = 1.0 ! length in y
    nstep=200 ! number of time steps   
    sim_time = 0.125  ! total simulation time 
    dt = sim_time/real(nstep) ! time step

    inquire(FILE=inp_file, EXIST=file_exists)
    if (file_exists) then
        print*, 'Input file exists...Getting input from it..'
        print*, 'WARNING! Input variables not defined in the input file will take default values.'
        call read_input(inp_file, Nx, Ny, sim_time, D, dt)
        ! if dt is changed then new nstep is calculated else the  
        ! nstep is based on old dt ie nstep=200
        nstep = int(sim_time/dt) 
    else
        print*, 'Input file does not exist. Continuing with default values..'

    endif

    INQUIRE(FILE=hotstart_file, EXIST=file_exists)
    if (file_exists) then
        OPEN(17, FILE=hotstart_file, FORM='unformatted')
        READ(17) Nx_tmp
        READ(17) Ny_tmp
        READ(17) D
        READ(17) sim_time
        READ(17) dt    
        READ(17) nstep_start !store the time step
        allocate(tmp_field(Nx_tmp, Ny_tmp), stat=info)
        READ(17) tmp_field ! store the array
        CLOSE(17)
        ! check for mismatch between input file and hotstart_file
        ! NOTE: Values from he hotstart file take precedence
        
        Nx = Nx_tmp
        Ny = Ny_tmp
        ! re-calculate the total steps
        nstep = int(sim_time/dt)

        ! reallocate T_old, T_new  and L
        call alloc(L, T_new, T_old, Nx, Ny, info)            
        T_old = tmp_field

        deallocate(tmp_field, stat=info)
        !endif 
    else
        call alloc(L, T_new, T_old, Nx, Ny, info)
            !initital condition
        T_old(:,:) = 0.0 ! temperature field at time step n
        L(:,:) = 0.0 !laplacian array
        nstep_start = 1
    endif

end subroutine initialize