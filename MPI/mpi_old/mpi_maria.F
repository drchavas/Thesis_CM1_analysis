	program mpi
	
	character*3 aaa, bbb
	character*4 cccc, dddd, eeee
	
	parameter(NA=17)
	real pres(NA), tver(NA), shum(NA)
	
	open (1,file='input_sounding',status='old')
	read(1,*) year, month, day
	read(1,*) aaa, sst
	read(1,*) bbb, psl
	write(*,*) psl
	read(1,*) cccc, npres
	read(1,*) pres
	read(1,*) dddd, ntver
	read(1,*) (tver(i),i=1,ntver)
	read(1,*) eeee, nshum
	read(1,*) (shum(i),i=1,nshum)
	close(1)

c	psl=psl*0.01 !per passar de Pa a mb

	do i=1, npres
c	pres(i)=pres(i)*0.01 !per passar de Pa a mb
	enddo

	
c	do i=1,nshum
c	shum(i)=shum(i)*1000.0 !per passar de kg/kg a g/kg
c	enddo
	if (nshum.lt.NA) then 
	do i=nshum+1,NA
	shum(i)=0.0
	enddo
	endif
	
	idim=min(na, npres, ntver, nshum)
	
	call pcmin(sst, psl, pres, tver, shum, NA, idim, pmin, vmax, ifl)
		
c	write(*,*) 'pmin ',pmin, ' vmax ', vmax, ' prova_ifl ',ifl


	open (2,file='sortidasubrutina.txt')
	write(2,*) pmin
	write(2,*) vmax
	write(2,*) ifl
	close(2)
		
	end
	
	include 'pcmin_revised.f'
