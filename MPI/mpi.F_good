        program mpi

	parameter(NL=200)        !max number of possible vertical levels
	real junk1, junk2, junk3
        real sst, psfc, thsfc, qvsfc
        real psnd(NL), tsnd(NL), qvsnd(NL), thsnd(NL)
        integer nlev

	!read in number of vertical levels in sounding
        open(3,file='num_lev.txt',status='old')
        read(3,*) nlev
        close(3)

	!extract sfc pressure
        open (1,file='input_sounding',status='old')
        read(1,*) psfc, junk1, junk2	!only need surface pressure [hPa]

        do k=1,nlev
          read(1,*) junk1,thsnd(k),qvsnd(k),junk2,junk3
        enddo
        close(1)

	!extract corresponding pressure data for sounding
        open (2,file='input_sounding_pres.txt',status='old')
        read (2,*) (psnd(k),k=1,nlev)
        do k=1,nlev
           psnd(k)=psnd(k)*.01	![Pa] --> [hPa]
        enddo

	!extract SST
        open (3,file='sstk.txt',status='old')
        read (3,*) sst	![K]
	sst=sst-273.15	![K]-->[C]

	!calculate actual temperature from potential temperature
        Rd=287.04	![J/kg/K]
        Cpd=1005.7 ![J/kg/K]
        do k=1,nlev
          tsnd(k)=thsnd(k)*(psnd(k)/psfc)**(Rd/Cpd)
          tsnd(k)=tsnd(k)-273.15
        enddo

        write(6,*) psfc, sst
        write(6,*)
        do k=1,nlev
           write(6,*)psnd(k), tsnd(k), qvsnd(k)
        enddo

	call pcmin(sst, psfc, psnd, tsnd, qvsnd, nlev, nlev, pmin, vmax, ifl)
		
c	write(*,*) 'pmin ',pmin, ' vmax ', vmax, ' prova_ifl ',ifl


	open (3,file='mpi_out.txt')
	write(3,*) pmin
	write(3,*) vmax
	write(3,*) ifl
	close(3)
	
       write(6,*) 'vmax = ', vmax, 'pmin = ', pmin
	
	end
	
	include 'pcmin_revised.f'
