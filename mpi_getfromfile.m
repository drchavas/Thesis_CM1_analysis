%mpi_getfromfile.m

function [mpi] = mpi_getfromfile(mpi_file);

cd input_soundings/
assert(exist(mpi_file,'file')==2,'CANNOT FIND MPI, FILE DOES NOT EXIST')

fid=fopen(mpi_file);
temp=textscan(fid,'%s%s%s%s');
mpi=str2num(char(temp{3}(1)));
fclose(fid);

cd ..

end