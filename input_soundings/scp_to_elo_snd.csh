#!/bin/csh

mkdir sndg_pres_mpi_temp
mv *_pressures *_mpi sndg_pres_mpi_temp/.

mkdir sndg_scp_temp
cp input_sounding* sndg_scp_temp/.
rsync -rav sndg_scp_temp/* drchavas@elorenz.mit.edu:"/home/drchavas/scripts_CM1_3d_2/sounding_files_TC/."
rm -r sndg_scp_temp

mv sndg_pres_mpi_temp/* .
rm -r sndg_pres_mpi_temp

