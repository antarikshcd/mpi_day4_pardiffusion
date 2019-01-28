subroutine update_field(Nx, Ny, T_old, T_new)
    USE mod_diff, ONLY:MK! contains allocation subroutine
    implicit none
    integer, intent(in) :: Nx, Ny
    real(MK), dimension(Nx, Ny), intent(in) :: T_new
    real(MK), dimension(Nx, Ny), intent(out) :: T_old

    T_old = T_new
    

end subroutine update_field