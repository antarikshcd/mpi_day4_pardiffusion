elemental subroutine elem_update_field(T_old, T_new)
     USE mod_diff, ONLY:MK! contains allocation subroutine
    implicit none
    !integer :: MK
    real(MK), intent(in) :: T_new
    real(MK), intent(out) :: T_old

    T_old = T_new
    

end subroutine elem_update_field