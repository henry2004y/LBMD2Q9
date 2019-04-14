function [Fx,Fy,Fx_2,Fy_2] = calc_drag(f,used_N,iStep,Ures,istep)
%calc_drag Calculate the drag coefficient C_D.
%

nx = Parameters.nx;
ny = Parameters.ny;
ex = Parameters.ex;
ey = Parameters.ey;
matrix_f = Parameters.matrix_f;
Uinit = Parameters.Uinit;
R = Parameters.R; 

% C_D calculation
% Calculating fluid force on the cylinder: momentum exchange approach
Fx = 0;
Fy = 0;
Fx_2 = 0;
Fy_2 = 0;

for j=1:ny; for i=1:nx
      if ~used_N(i,j)
         f_x = 0;
         f_y = 0;
         f_x_2 = 0;
         f_y_2 = 0;
         for jP=-1:1; for iP=-1:1
               % for b and a = 1 -1 0 different cases a+2 and b+3 return
               % same results the desired ex and ey
               if used_N(i+iP,j+jP) == 1
                  % first method suggested by J. Gotz et al,
                  % Large scale simulation of fluid structure using LBM and
                  % the 'physics engine'
                  f_x = -ex(iP+2)*( 2*f(i+iP,j+jP,matrix_f(jP+2,iP+2)) ) +...
                     f_x; %*dx/dt
                  f_y = -ey(jP+3)*( 2*f(i+iP,j+jP,matrix_f(jP+2,iP+2)) ) +...
                     f_y; %*dx/dt same as above
                  f_x_2 = -ex(iP+2)*(f(i,j,matrix_f(-jP+2,-iP+2)) +...
                     f(i+iP,j+jP,matrix_f(jP+2,iP+2)) ) + f_x_2;%*dx/dt
                  f_y_2 = -ey(jP+3)*(f(i,j,matrix_f(-jP+2,-iP+2)) +...
                     f(i+iP,j+jP,matrix_f(jP+2,iP+2)) ) + f_y_2;%*dx/dt
               end
            end; end
         
         Fx = Fx + f_x;
         Fy = Fy + f_y;
         Fx_2 = Fx_2 + f_x_2;
         Fy_2 = Fy_2 + f_y_2;
      end
   end; end

% Based on second momentum method
Cd_lattice_2 = abs(Fx_2)/(Uinit)^2/R;
% When converging, average
% C_d_keeper(loop_counter) = Cd_lattice_2;
% cd can be calculated from values obtained by cd_keeper

if mod(iStep,Parameters.nLog) == 0
   fprintf('%i Force_x= %f, cd second lattice=%f, Uresidual = %f\n',...
      iStep,  Fx, Cd_lattice_2, Ures)
end

end