! read input file for the diffusion problem: Nx,Ny,D,nstep,dt,sim_time

subroutine read_input(filename, Nx, Ny, sim_time, D, dt)
    USE mod_diff, ONLY:MK! contains allocation subroutine
   implicit none
   !integer :: MK
   character(len=256) :: buffer, keyword, cvalue, filename
   integer :: i, freq, flen, slen, vlen, pos_eq, ilen, Nx, Ny, D, ios
   real :: sim_time, dt
       
   print*,'filename:', filename ! debug
   OPEN(10, FILE=filename, STATUS='OLD', ACTION='READ')
   i = 1
   do
       read(10, '(A)', iostat=ios, err=100, end=200) buffer
          if (buffer(1:1).ne.'#') then ! for comments
              ilen = len_trim(buffer)
              if (ilen .ne. 0) then
                  pos_eq = index(buffer, '=')
                  if (pos_eq .ne. 0) then
                      ! split the string into the variable and the value
                      keyword = trim(buffer(1:pos_eq-1))
                      cvalue = trim(buffer(pos_eq+1:ilen))
                      !print*,'keyword:', keyword !debug
                      !print*,'cvalue:', cvalue !debug
                      !print*, 'lenght of keyword:', len(keyword) !debug
                      slen = len_trim(keyword)
                      !print*,'slen:', slen !debug
                      if (keyword(1:slen) == 'Nx') then
                          read(cvalue,*) Nx
                          !print*, 'Nx= ', Nx ! debug
                      elseif (keyword(1:slen) == 'Ny') then
                          read(cvalue,*) Ny
                          !print*, 'Ny= ', Ny ! debug
                      elseif (keyword(1:slen) == 'sim_time') then
                          read(cvalue,*) sim_time
                          !print*, 'sim_time= ', sim_time !debug
                      elseif (keyword(1:slen) == 'D') then
                          read(cvalue,*) D
                          !print*, 'D= ', D !debug
                      elseif (keyword(1:slen) == 'dt') then
                          read(cvalue,*) dt
                          !print*, 'dt= ', dt !debug
                        
                      endif    

                  endif
              endif
          endif    
       i = i+1
   enddo
   
   100 CONTINUE ! error label
       write(*,'(2(A,I5))') 'Error reading file in line:', i, 'IOSTAT=', ios
   200 CONTINUE ! end of line label
       write(*,'(A,I6,A)') 'Read', i-1,'lines.'
   close(10)          

end subroutine read_input