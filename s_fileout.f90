! file output subroutine
subroutine file_out(Nx, Ny, dx, dy, T_new, tstep_count)
    USE mod_diff, ONLY:MK! contains allocation subroutine
    
    implicit none
    !integer :: MK
    integer, intent(in) :: Nx, Ny
    integer, optional :: tstep_count
    real, intent(in) :: dx, dy
    real(MK), dimension(:, :) :: T_new
    character(len=20) :: filename
    integer :: i, j, real_len

    if (present(tstep_count)) then
        ! write field at present time step 
        write(filename,'(A,I6.6,A)') 'diff_',tstep_count,'.dat'
        real_len = len_trim(filename)
        filename = filename(1:real_len)
        open(unit=10, file=filename)
    else
        ! write final field
        OPEN(UNIT=10, FILE='diff_final.dat')
    endif 

    !WRITE to file
    DO j=1,Ny
        DO i=1,Nx
            WRITE(10,'(3E12.4)') REAL(i-1)*dx,REAL(j-1)*dy,T_new(i,j)
        ENDDO
        WRITE(10,'(A)') !Will produce a new empty line
    ENDDO
    CLOSE(10)

end subroutine file_out