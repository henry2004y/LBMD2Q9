function [used_N, U, V, cbody] = set_geometry(nx,ny,U,V)
%set_geometry Set the geometry for the simulation domain.
%
%INPUTS:
% nx: number of grid in x 
% ny: number of grid in y
% U:  velocity in x
% V:  velocity in y
%OUTPUTS:
% used_N: logical for active nodes
% U:  velocity in x
% V:  velocity in y
% cbody: solid body centroid coordinate


% used_N(:,1)  = 0; % south wall/ down wall
% used_N(:,ny) = 0; % north wall/ upper wall

cbody = [floor(nx/5+1) floor(ny/2+3)];

[Gridx,Gridy] = ndgrid(1:nx,1:ny);

used_N = (Gridx - cbody(1)).^2 + (Gridy - cbody(2)).^2 >= Parameters.R^2;

U(~used_N) = 0;
V(~used_N) = 0;

end