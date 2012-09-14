#!/bin/csh

mkdir sndg_scp_temp
cp input_sounding*_drag sndg_scp_temp/.
rsync -rav sndg_scp_temp/* drchavas@hopper.nersc.gov:"/global/u1/d/drchavas/scripts_CM1_ax_hop/sounding_files/."
rsync -rav sndg_scp_temp/* drchavas@hopper.nersc.gov:"/global/u1/d/drchavas/scripts_CM1_3d_hop/sounding_files/."
rm -r sndg_scp_temp

