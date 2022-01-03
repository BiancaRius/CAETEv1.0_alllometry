program test_random_number
  integer:: j,k
  integer, parameter:: npls=20
  integer, parameter :: time =100
 
  real:: x(npls,time)
 

  


        
    
    do k=1,time
        do j=1,npls
            x(j,k) = xmin + (xmax-xmin)*x(j,k)
            print*,x(j,k),j,k
        enddo
    enddo   

       
end program