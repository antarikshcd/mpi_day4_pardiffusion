!module containing the global data for exercise 2b diffusion problem
module mod_alloc
    
    interface alloc
        module procedure salloc, dalloc
    end interface alloc    	
    contains
        subroutine salloc(L, T_new, T_old, Nx, Ny, info)
            implicit none
            integer, parameter :: MK = KIND(1.0E0)
            integer :: Nx,Ny,info
            real(MK), dimension(:, :), allocatable :: T_old, T_new, L
            real(MK), dimension(:, :), allocatable :: work !local temp array
            ! if already allocated do not reallocate
            if (.not.allocated(T_old)) then
                allocate(T_old(Nx,Ny), stat=info)
            
            else
                allocate(work(Nx, Ny), stat=info)
                !print*,'1 status:', info
                work = T_old
                deallocate(T_old, stat=info)
                !print*,'2 status:', info
                ! allocate T_old with new size
                allocate(T_old(Nx,Ny), stat=info)
                !print*,'3 status:', info
                T_old = work
                ! deallocate work
                deallocate(work, stat=info)
                !print*,'4 status:', info 
                                
            endif
            ! allocation of T_new
            if (.not.allocated(T_new)) then
                allocate(T_new(Nx,Ny), stat=info)
                !print*,'5 status:', info
            else
                allocate(work(Nx, Ny), stat=info)
                work = T_new
                deallocate(T_new, stat = info)

                allocate(T_new(Nx, Ny), stat=info)
                T_new = work

                !deallocate work
                deallocate(work, stat=info)                    
            endif
            ! allocation of L
            if (.not.allocated(L)) then
                allocate(L(Nx,Ny), stat=info)
                !print*,'6 status:', info
            else
                allocate(work(Nx, Ny), stat=info)
                work = L
                deallocate(L, stat=info)

                allocate(L(Nx, Ny), stat=info)
                L = work

                ! deallocate work
                deallocate(work, stat=info)                    
            endif

        end subroutine salloc
        
        subroutine dalloc(L, T_new, T_old, Nx, Ny, info)
            implicit none
            integer, parameter :: MK = KIND(1.0D0)
            integer :: Nx,Ny,info
            real(MK), dimension(:, :), allocatable :: T_old, T_new, L
            real(MK), dimension(:,:), allocatable :: work !local temp array
            ! if already allocated do not reallocate
            if (.not. allocated(T_old)) then
                allocate(T_old(Nx,Ny), stat=info)
            
            else
                allocate(work(Nx, Ny), stat=info)
                !print*,'1 status:', info
                work = T_old
                deallocate(T_old, stat=info)
                !print*,'2 status:', info
                ! allocate T_old with new size
                allocate(T_old(Nx,Ny), stat=info)
                !print*,'3 status:', info
                T_old = work
                ! deallocate work
                deallocate(work, stat=info)
                !print*,'4 status:', info 
            
            endif
            ! allocation of T_new
            if (.not.allocated(T_new)) then
                allocate(T_new(Nx,Ny), stat=info)
                !print*,'5 status:', info
            else
                allocate(work(Nx, Ny), stat=info)
                work = T_new
                deallocate(T_new, stat = info)

                allocate(T_new(Nx, Ny), stat=info)
                T_new = work

                !deallocate work
                deallocate(work, stat=info)                    
            endif
            ! allocation of L
            if (.not.allocated(L)) then
                allocate(L(Nx,Ny), stat=info)
                !print*,'6 status:', info
            else
                allocate(work(Nx, Ny), stat=info)
                work = L
                deallocate(L, stat=info)

                allocate(L(Nx, Ny), stat=info)
                L = work

                ! deallocate work
                deallocate(work, stat=info)                    
            endif    
        
        end subroutine dalloc
end module mod_alloc

