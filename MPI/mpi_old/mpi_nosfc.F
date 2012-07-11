        program mpi

        parameter(NL=17)        !number of levels above surface to be read in
        real junk1, junk2, junk3
        real sst, psfc, thsfc, qvsfc
        real psnd(NL), tsnd(NL), qvsnd(NL), thsnd(NL)

	!extract sounding data for sfc pressure, pot temp, and wv mixing ratio
        open (1,file='input_sounding',status='old')
        read(1,*) psfc, thsfc, qvsfc	!need [hPa], [K], [g/kg]
        do k=1,NL
          read(1,*) junk1,thsnd(k),qvsnd(k),junk2,junk3
        enddo
        close(1)

	!extract corresponding pressure data for sounding
        open (2,file='input_sounding_pres.txt',status='old')
        read (2,*) (psnd(k),k=1,NL)
        do k=1,NL
           psnd(k)=psnd(k)*.01	![Pa] --> [hPa]
        enddo

	!SST
	sst=299.28-273.15	![K]-->[C]

	!calculate actual temperature from potential temperature
        Rd=287	![J/kg/K]
        Cpd=1004 ![J/kg/K]
        do k=1,NL
          tsnd(k)=thsnd(k)*(psnd(k)/psfc)**(Rd/Cpd)
          tsnd(k)=tsnd(k)-273.15
          write(6,*) tsnd(k)
        enddo

        write(6,*) psfc, thsfc, qvsfc, sst
        do k=1,NL
           write(6,*)psnd(k), tsnd(k), qvsnd(k)
        enddo

	call pcmin(sst, psfc, psnd, tsnd, qvsnd, NL, NL, pmin, vmax, ifl)
		
c	write(*,*) 'pmin ',pmin, ' vmax ', vmax, ' prova_ifl ',ifl


	open (3,file='mpi_out.txt')
	write(3,*) pmin
	write(3,*) vmax
	write(3,*) ifl
	close(3)
		
	end
	
	include 'pcmin_revised.f'
