subroutine save_restart(hotstart_file, Nx, Ny, D, sim_time, dt, itstep, T_old)
    USE mod_diff, ONLY:MK! contains allocation subroutine

    implicit none
    !integer :: MK
    character(len=*) :: hotstart_file
    integer :: Nx, Ny, D, itstep
    real :: sim_time, dt
    real(MK), dimension (:,:) :: T_old


    ! Saving a hotstart file for T_old at each time-step
    OPEN(15, FILE=hotstart_file, FORM='unformatted', STATUS='replace')
    WRITE(15) Nx
    WRITE(15) Ny
    WRITE(15) D
    WRITE(15) sim_time
    WRITE(15) dt    
    WRITE(15) itstep !store the time step
    WRITE(15) T_old ! store the array
    CLOSE(15)

end subroutine save_restart
