program light_competition

		integer,parameter::npls=14
		real, dimension(npls) :: height
		real :: max_height
		real :: num_layer
		integer::num_layer_round
		real :: layer_size
		real :: max_layer
		real :: heigher_layer
		real :: ipar


		real,allocatable :: layer_height(:)
		real,allocatable :: light_availability(:)
		real,allocatable :: light_used (:)
		real,allocatable :: light_incidence (:)
		real,allocatable :: height_aux(:,:)
		real,allocatable :: id(:)
! do i = 1, npls
! 	height(i) = 35 -(i*2)+3
! 	if (i.eq.1) then
! 		height(i)=3.
! 	endif
! 	 print*, 'height', height(i)
! enddo
		height=(/2.0,3.0,3.7,10.,10.5,12.,13.,20.,22.,27.,27.5,28.,29.,34./)
! do i=1, npls
! 	print*, 'height',height(i)
! enddo

		max_height = maxval(height)
!print*, 'max_height',max_height

		num_layer = max_height/5
!print*, 'num_layer',num_layer

		num_layer_round = nint(num_layer)
		n=num_layer_round

print*,'num_layer_round',num_layer_round,n

		layer_size = max_height/num_layer_round
!print*, 'layer_size', layer_size

		allocate(layer_height(1:n))
		allocate(light_availability(1:n))
		allocate(light_used(1:n))
		allocate(light_incidence(1:n))
		allocate(height_aux(1:n,1:npls)) !it identi
		allocate(id(1:n))


		do j=1,n
			layer_height(j)=0
		enddo

!JOÂO: algoritmos mais simples
		do j=1,n
			if (j.eq.1) then
				layer_height(j)=layer_size
			else
				layer_height(j)=layer_height(j-1)+layer_size
			endif
			print*, 'layer_height',layer_height(j)
		enddo

!do for attribute the layer id to each pls 
		do i=1,npls
			do j=1,n
				if ((height(i).le.layer_height(j)).and.(height(i).gt.(layer_height(j-1)))) then
					height_aux(j,i)=height(i)
					id(j)=j
					print*,height_aux(j,i),i,j,id(j)
				endif
			enddo
		enddo




		ipar=100.


  ! 		do i=1,npls
		! 	do j=n,1,-1
		! 		if(id(j).eq.n)then
		! 			light_incidence(j) = ipar
		! 			light_used(j) = light_incidence(j)*0.2
		! 			light_availability(j) = light_incidence(j)-light_used(j)
		! 		endif
		! 		if ((id(j).gt.0).and.(id(j).ne.n)) then
		! 			light_incidence(j) = light_availability(j+1)
		! 			light_used(j) = light_incidence(j)*0.2
		! 			light_availability(j) = light_incidence(j)-light_used(j)
		! 		else if ((light_availability(j).eq.0).and.id(j).ne.n)then
		! 			light_incidence(j-1)=light_availability(j+1)
 	! 				light_used(j-1) = light_incidence(j-1)*0.2
		! 			light_availability(j-1) = light_incidence(j-1)-light_used(j-1)
		! 		endif
		! 	enddo
		! enddo
   

do i=1,npls
	do j=n,1,-1
		if (j.eq.n) then
			light_incidence(j) = ipar
			light_used(j) = light_incidence(j)*0.2
			light_availability(j) = light_incidence(j)-light_used(j)
		else if(j.eq.1) then
				light_incidence(j) = light_availability(2)
				light_used(j) = light_incidence(j)*0.2
				light_availability(j) = light_incidence(j)-light_used(j)
		else if ((height(i).le.layer_height(j)).and.(height(i).gt.(layer_height(j-1)))) then
			 light_incidence(j) = light_availability(j+1)
			 light_used(j) = light_incidence(j)*0.2
			 light_availability(j) = light_incidence(j)-light_used(j)
		
		endif
	enddo
enddo


 do j=n,1,-1
 	if (light_availability(j).eq.0) then
 		
 	
 	 	light_incidence(j-1)=light_availability(j+1)
 	 	light_used(j-1) = light_incidence(j-1)*0.2
		light_availability(j-1) = light_incidence(j-1)-light_used(j-1)
 	endif
 enddo

!!!!ATENTION IF A LAYER DOES'NT HAVE A PLS
	do j=1,n
		print*, j,'avai',light_availability(j),'used',light_used(j),'inc',light_incidence(j)
	enddo
end program light_competition