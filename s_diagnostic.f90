subroutine diagnostic(k, dt, nstep, T_new)

    USE mod_diff, ONLY:MK! contains allocation subroutine
    implicit none
    !integer :: MK
    integer, intent(in) :: k, nstep
    real, intent(in) :: dt
    real:: time, min_T  
    real(MK), dimension(:, :), intent(in) :: T_new !assumed shape array
    logical :: first = .TRUE. ! saves the value for opening file
    
    ! calculate time
    time = k*dt
    ! minimum value of array
    min_T = minval(T_new)

    if (first) then
        first = .FALSE.
        open(20, file='diag.dat')
    endif

    write(20, '(3E12.4)') time, min_T

    if (k .EQ. nstep) then
        close(20)
    endif    
end subroutine diagnostic