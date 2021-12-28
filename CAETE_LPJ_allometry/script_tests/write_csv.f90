program p1
    integer ::  j    
    integer, parameter:: npls = 20 
    real, dimension(npls) :: dwood !wood density (g/cm-3) *Fearnside, 1997 - aleatory choices
    real, dimension(npls) :: dwood_aux
    real, dimension (100,10):: a

    dwood=(/0.24,0.53,0.39,0.32,0.31,0.44,0.66,0.42,0.74,0.39,0.82,0.40,0.26,0.79,0.39,0.52,0.41,0.44,0.86,0.42/)!,&

    open(unit=1,file='test.txt',status='unknown')
    do j = 1,npls
      
       dwood_aux(j)= dwood(j)
       write(1,*) dwood_aux(j) ! newline
    enddo    
    
    close(1)


    ! open(unit=1,file='test2.csv',status='unknown')
    ! do j = 1,npls
      
    !     dwood_aux(j) = dwood(j)
    !  enddo 
    ! write(1,'(20(g14.7,x))') dwood_aux
    ! close(1)

    
    ! ! precompute full table
    ! call random_number(a)
    !  ! write table
    ! open(unit=1,file='test3.csv',status='unknown')
    ! write(1,'(100(g14.7,x))') a
    ! close(1)

  end program