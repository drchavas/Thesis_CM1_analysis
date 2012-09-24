#!/bin/csh

./mpi_runall.csh 275.00
./mpi_runall.csh 280.00
./mpi_runall.csh 285.00
./mpi_runall.csh 287.50
./mpi_runall.csh 290.00
./mpi_runall.csh 292.50
./mpi_runall.csh 295.00
./mpi_runall.csh 297.50
./mpi_runall.csh 300.00
./mpi_runall.csh 302.50
./mpi_runall.csh 305.00
./mpi_runall.csh 310.00

./mpi.csh_usr_varyCkCd input_sounding_3dRCE_nx48_SST300.00K_Tthresh200K_usfc3_drag 300.00 0.125 
./mpi.csh_usr_varyCkCd input_sounding_3dRCE_nx48_SST300.00K_Tthresh200K_usfc3_drag 300.00 0.25
./mpi.csh_usr_varyCkCd input_sounding_3dRCE_nx48_SST300.00K_Tthresh200K_usfc3_drag 300.00 0.5
./mpi.csh_usr_varyCkCd input_sounding_3dRCE_nx48_SST300.00K_Tthresh200K_usfc3_drag 300.00 2.0
./mpi.csh_usr_varyCkCd input_sounding_3dRCE_nx48_SST300.00K_Tthresh200K_usfc3_drag 300.00 4.0
./mpi.csh_usr_varyCkCd input_sounding_3dRCE_nx48_SST300.00K_Tthresh200K_usfc3_drag 300.00 8.0   
