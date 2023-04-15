function [buildpath,Nx,Ny,Nz,dx,dy,dz]=genpath(dt, Lx, Ly, Lz)
% function to generate a path
% output is a Nx3 matrix with coordinates of each element starting from row 1


% in this test case, the head moves along x direction

v = 50e-3; % extrusion speed 50 mm/s taken from Garon's notes.

dz = 0.2e-3;
dx = v*dt  % parameter needs to be adjust, set to dz for now. Need to match with the build speed?
dy = 0.4e-3; %dz;

dim = [Lx/dx, Ly/dy, Lz/dz];
if any(dim ~= round(dim))
    fprintf('warning: Nx/Ny/Nz is not an integer');
    dim = round(dim);
end
fprintf('Nx = %d Ny = %d Nz = %d\n', dim(1), dim(2), dim(3))
buildpath = zeros(prod(dim),3); % initialize buildpath matrix
rownow = 0;

for iz = 1:dim(3)
    for iy = 1:dim(2)
        buildpath(rownow+1:rownow+dim(1),2) = iy;
        buildpath(rownow+1:rownow+dim(1),3) = iz;
        if mod(iy,2) ~= 0
            buildpath(rownow+1:rownow+dim(1),1) = 1:dim(1);

        else
            buildpath(rownow+1:rownow+dim(1),1) = dim(1):-1:1;

        end

        rownow = rownow+dim(1);
    end
end

Nx = dim(1); Ny = dim(2); Nz = dim(3);


