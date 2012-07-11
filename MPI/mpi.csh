#!/bin/csh

## $1 = SST [ex. 300.00] K; $2 = T_tpp [ex. 150] K; $3 = usfc [ex. 5] m/s; $4 = radcool [ex. 2] K/halfday
## input_sounding_3dRCE_nx48_SST300.00K_Tthresh150K_usfc5

if (-e input_sounding) then
  rm input_sounding input_sounding_pres
endif

## Copy over desired input_sounding and corresponding pressure file
if ($4 == 0.5) then
  if (-e ../input_soundings/input_sounding_3dRCE_nx48_SST${1}K_Tthresh${2}K_usfc${3}) then
    cp ../input_soundings/input_sounding_3dRCE_nx48_SST${1}K_Tthresh${2}K_usfc${3} input_sounding
    cp ../input_soundings/input_sounding_3dRCE_nx48_SST${1}K_Tthresh${2}K_usfc${3}_pressures input_sounding_pres.txt
  endif
else
  if (-e ../input_soundings/input_sounding_3dRCE_nx48_SST${1}K_Tthresh${2}K_usfc${3}_rad${4}K) then
    cp ../input_soundings/input_sounding_3dRCE_nx48_SST${1}K_Tthresh${2}K_usfc${3}_rad${4}K input_sounding
    cp ../input_soundings/input_sounding_3dRCE_nx48_SST${1}K_Tthresh${2}K_usfc${3}_rad${4}K_pressures input_sounding_pres.txt
  endif
endif

## Check if file exists -- if so, then run mpi code
if (-e input_sounding) then

## Update sstk.txt
echo $1 > sstk.txt

## Run code: requires input_sounding, input_sounding_pres.txt
grep -c "." input_sounding_pres.txt > num_lev.txt	#number of lines in sounding
./mpi.exe > dan.temp
tail -1 dan.temp
tail -1 dan.temp > dan2.temp

echo $1 $2 $3 $4 >> mpi_list.txt
tail -1 dan.temp >> mpi_list.txt
rm dan.temp

if ($4 == 0.5) then
  mv dan2.temp ../input_soundings/input_sounding_3dRCE_nx48_SST${1}K_Tthresh${2}K_usfc${3}_mpi
else
  mv dan2.temp ../input_soundings/input_sounding_3dRCE_nx48_SST${1}K_Tthresh${2}K_usfc${3}_rad${4}K_mpi
endif

rm input_sounding input_sounding_pres.txt sstk.txt

endif

