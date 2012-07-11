#!/bin/csh

## $1 = SST [ex. 300.00] K; $2 = T_tpp [ex. 150] K; $3 = usfc [ex. 5] m/s; $4 = radcool [ex. 2] K/day
## input_sounding_3dRCE_nx48_SST300.00K_Tthresh150K_usfc5

set SSTs = ('275.00' '285.00' '290.00' '295.00' '297.50' '300.00' '302.50' '305.00' '310.00')
set TPPs = ('150' '175' '200' '225' '250')
set us = ('.5' '1' '2' '3' '4' '5' '10')
set radhalfs = ('0.125' '0.25' '0.5' '1.0' '2.0')

#set SSTs = ('300.00')
#set TPPs = ('200' '201')
#set us = ('3')
#set radhalfs = ('0.5')

rm mpi_list.txt

## Compile MPI code
gfortran -o mpi.exe mpi.F -ffixed-line-length-none

foreach SST (${SSTs})
foreach TPP (${TPPs})
foreach u (${us})
foreach radhalf (${radhalfs})

echo $SST $TPP $u $radhalf
./mpi.csh $SST $TPP $u $radhalf

end
end
end
end

