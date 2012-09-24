#!/bin/csh

#Syntax: ./mpi_runall.csh [SST]; ex: ./mpi_runall.csh 300.00

cd ../input_soundings
mkdir sndg_pres_mpi_temp
mv *_pressures *_mpi sndg_pres_mpi_temp/.

ls input_sounding*SST$1K* > file_list              
## find out how many files I have 
#@ nf = `cat file_list | wc -l`
#@ n = 1        # define a looping variable
#echo $nf

cd sndg_pres_mpi_temp/
mv *_pressures *_mpi ../.
cd ..
rm -r sndg_pres_mpi_temp

## loop through files
foreach sndfil ( "`cat file_list`" )
  cd ../mpi
  echo 'Calculating mpi for' $sndfil
  ./mpi.csh_usr $sndfil $1
  cd ../input_soundings
end

rm file_list

cd ../mpi

