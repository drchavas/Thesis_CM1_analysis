%sounding_to_3d.m

%Created: 17 Oct 2012, Dan Chavas

%Purpose: this function takes as input a VERTICAL profile of any quantity
%(e.g. from a sounding) and expands it into a three-dimensional array with
%the desired number of grid-points in the x and y directions

function [prof_3d] = sounding_to_3d(prof_z,nx,ny)

    temp = repmat(prof_z,nx,1);   %copy vertical profile in x-dimension
    temp = repmat(temp,[1 1 ny]);   %copy xz matrix in y-dimension
    prof_3d = permute(temp,[1 3 2]);   %permute so that z is final dimension, x is first, y is second

end