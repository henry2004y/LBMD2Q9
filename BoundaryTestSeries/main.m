% 2D Lattice Boltzmann (BGK) model of a fluid.                                  
%  c6  c2  c5  D2Q9 model. At each timestep, particle densities propagate      
%    \  |  /   outwards in the directions indicated in the figure. An          
%  c3 -c9 -c1  equivalent 'equilibrium' density is found, and the densities   
%    /  | \    relax towards that state, in a proportion governed by omega.
%  c7  c4  c8  
%
% Original author: Amirsaman Farrokhpanah
% Modified by Hongyang Zhou, hyzhou@umich.edu 03/31/2019

clear; clc
%% Parameters

nx = Parameters.nx;
ny = Parameters.ny;

ex = Parameters.ex;
ey = Parameters.ey;
w  = Parameters.w;
cs = Parameters.cs;
matrix_f = Parameters.matrix_f;

Gx = Parameters.Gx;
Gy = Parameters.Gy;

Uinit = Parameters.Uinit;
Rho0Out = Parameters.Rho0Out;
R = Parameters.R;

method = Parameters.method;

Tau = Parameters.Tau;

dt = Parameters.dt;
dx = Parameters.dx;

extemp = reshape(ex,1,1,[]);
eytemp = reshape(ey,1,1,[]);
Nxp1 = circshift(1:nx,1);
Nxm1 = circshift(1:nx,-1);
Nyp1 = circshift(1:ny,1);
Nym1 = circshift(1:ny,-1);

Utol = Parameters.Utol;	  % Tolerance of U
Vtol = Parameters.Vtol;	  % Tolerance of V
Ures = Parameters.Ures;   % Residual of U
Vres = Parameters.Vres;   % Residual of V 

% Define the parameters to be used
%ftemp = zeros(nx,ny,9);
feq   = zeros(nx,ny,9);
U     = zeros(nx,ny);
V     = zeros(nx,ny);
Rho   = zeros(nx,ny);

% Initialization
U(:,:)   = Uinit;
V(:,:)   = Parameters.Vinit;
Rho(:,:) = Parameters.Rhoinit;

Uold = U;
Vold = V;

ttt = 1000; % For Tecplot outputs
nStep = Parameters.nStep;


% Wall rotation available only in two methods
if method == 2 || method == 3 || method == 5 
    
    wall_rotation = 0;
    
    fprintf('wall rotation = ? (default: %d)\n', wall_rotation);
    replywall_rotation = input('');
    
    if isempty(replywall_rotation)
        replywall_rotation = wall_rotation ;
    end
    
    wall_rotation = replywall_rotation;
end

fprintf('2D Lattice Boltzmann Method Simulation\n')
fprintf('*****************************************\n')
fprintf('Re of flow  = %f              \n',Parameters.Re)
fprintf('Re Cylinder = %f              \n',Parameters.Re_cylinder)
fprintf('FD          = %f              \n',Parameters.FD)
fprintf('Tau         = %f              \n',Tau)
fprintf('time lattice= %f              \n',Parameters.dt)
fprintf('*****************************************\n')

%% Geometry

[used_N, U, V, cbody] = set_geometry(nx,ny,U,V);

%% Initial conditions

U(1,:) = Uinit;
V(1,:) = 0;
% V(nx,:) = 0;

% U(:,1) = Uinitial;
% V(:,1) = 0;
% U(:,ny)= Uinitial;
% V(:,ny)= 0;

% Initialize distribution function
% Calculate equlibrium distribution based on each node from initial
% macroscopic values
for ie=1:9
   feq(:,:,ie) = Rho.*w(ie).*(1 + 3*(ex(ie)*U + ey(ie)*V) + ...
      9/2*(ex(ie)*U + ey(ie)*V).^2 - ...
      3/2*(U.^2 + V.^2));
end

% Initialize f
f = feq;

%% Main

iStep = 1;

while iStep < nStep && (Ures > Utol) %|| Vres > Vtol)
   
   % Calculate g (hyzhou: is this the force term?)
   g = w.*(Gx.*ex + Gy.*ey) ./ cs^2;
   
   % Collision
   for j=1:ny; for i=1:nx
      if used_N(i,j)
         for ie=1:9
            f(i,j,ie) = f(i,j,ie) - (f(i,j,ie)-feq(i,j,ie))/Tau + g(ie);
         end
      end
   end; end

   % Post-collision values for curved walls, for boundary nodes x_b
   switch method
      case 2
         f = filippova(f,used_N,wall_rotation,cbody,U,V);
      case 3
         f = mei(f,used_N,wall_rotation,cbody,U,V);
      case 4
         f = bozidi(f,used_N,cbody);
      case 5
         f = yu(f,Rho,used_N,cbody,wall_rotation);
   end
   
   [Fx,Fy,Fx_2,Fy_2] = calc_drag(f,used_N,iStep,Ures);
   
   % Streaming   
   f(:,:,1) = f(Nxp1,:,   1);
   f(:,:,2) = f(:,   Nyp1,2);
   f(:,:,3) = f(Nxm1,:,   3);
   f(:,:,4) = f(:,   Nym1,4);
   f(:,:,5) = f(Nxp1,Nyp1,5);
   f(:,:,6) = f(Nxm1,Nyp1,6);
   f(:,:,7) = f(Nxm1,Nym1,7);
   f(:,:,8) = f(Nxp1,Nym1,8);   
  
   % Velocity calibration
   U(1,:) = Uinit;
   V(1,:) = 0;

   U(~used_N) = 0;
   V(~used_N) = 0;
   
   % Inlet & Outlet BC: Velocity and pressure boundary, He & Zou  
   for j=1:ny     
      % Inlet
      Rho0In = (f(1,j,9) + f(1,j,2) + f(1,j,4) + ...
         2*(f(1,j,3) + f(1,j,6) + f(1,j,7))) / (1 - U(1,j));
      ru = Rho0In*U(1,j);
      rv = Rho0In*V(1,j);
      f(1,j,1) = f(1,j,3) + (2/3)*ru;
      f(1,j,5) = f(1,j,7) + 1/6*ru + 0.5*rv + 0.5*(f(1,j,4) - f(1,j,2));
      f(1,j,8) = f(1,j,6) + 1/6*ru - 0.5*rv + 0.5*(f(1,j,2) - f(1,j,4));
 
      % Outlet
      U(nx,j) = -1 + (f(nx,j,9) + f(nx,j,2) + f(nx,j,4) + ...
         2*(f(nx,j,1) + f(nx,j,5) + f(nx,j,8))) / Rho0Out;
      ru = Rho0Out*U(nx,j);
      V(nx,j) = 0;
      rv = Rho0Out*V(nx,j);
      f(nx,j,3) = f(nx,j,1) - 2/3*ru;
      f(nx,j,6) = f(nx,j,8) + 0.5*(f(nx,j,4) - f(nx,j,2)) - 1/6*ru;
      f(nx,j,7) = f(nx,j,5) + 0.5*(f(nx,j,2) - f(nx,j,4)) - 1/6*ru;     
   end
   
   % Wall BC: bounce back
   if method == 1
      for j=1:ny; for i=1:nx
         if ~used_N(i,j) % lower and upper boundaries                           
            f(i,j,[3 1]) = f(i,j,[1 3]);
            f(i,j,[4 2]) = f(i,j,[2 4]);
            f(i,j,[7 5]) = f(i,j,[5 7]);
            f(i,j,[8 6]) = f(i,j,[6 8]);
         end
      end; end
   end
      
   % Calculate macroscopic quantities
   % Note: it doesn't matter to include the false nodes, because they are
   % not used anyway.
   Rho = sum(f,3);
   U = sum(f.*extemp,3) ./ Rho;
   V = sum(f.*eytemp,3) ./ Rho;
   
   % Velocity calibration
   U(1,:) = Uinit;
   V(1,:) = 0;
   
   U(~used_N) = 0;
   V(~used_N) = 0;
   
   % pressure difference calculations
   %calc_pressure(Rho,U);
   
   % Recalculate feq based on each node from initial macroscopic values
   for ie=1:9
      feq(:,:,ie) = Rho.*w(ie).*(1 + 3*(ex(ie)*U + ey(ie)*V) + ...
         4.5*(ex(ie)*U + ey(ie)*V).^2 - 1.5*(U.^2 + V.^2));
   end

   % Check the convergence
   % hyzhou: I think Vres is somewhere wrong...
   Ures = 0;
   Vres = 0;

   % This is somewhere in the middle?
   for j=1:ny
      for i=1:2*cbody(1)
         if used_N(i,j)
            if U(i,j) ~= 0
               Ures = Ures + abs((U(i,j) - Uold(i,j)) / U(i,j));
            end
            if V(i,j) ~= 0 
               Vres = Vres + abs((V(i,j) - Vold(i,j)) / V(i,j));
            end
         end
        
      end
   end
   
   Uold = U;
   Vold = V;
   
   % loop counter
   iStep = iStep + 1;
   
   % Animate
   if mod(iStep,Parameters.nPlot) == 0
      uMag = sqrt(U.^2+V.^2);

      imagesc(uMag');
      axis equal off
      drawnow
      fprintf('iStep = %i, Uresidual = %f, Vresidual = %f \n',...
         iStep, Ures, Vres)
   end
   
   % Tecplot format file writing 
   if mod(iStep,Parameters.nPlotSave) == 0
      s = num2str(ttt);
      A = fopen(s,'w');
      fprintf(A,[' TITLE= "STEP Output"   VARIABLES= "I","J","U", "V"'...
         'ZONE T="Flow Field", I= %d , J= %d , F=POINT \n '],nx,ny);
      for j=1:ny; for i=1:nx
         fprintf(A,'\n %f %f %f %f \n',i,j,U(i,j),V(i,j));
      end; end
      fclose('all');
      ttt = ttt + 1;
   end
   
end % end main loop

%% post processing

% Contour plot
y = (1:ny);
x = (1:nx);
[X,Y] = ndgrid(x,y);
starty = linspace(1,ny,10);
startx = ones(size(starty));
figure
contourf(X,Y,U,'LineColor','none'); colorbar; hold on
streamline(X',Y',U',V',startx,starty)
axis equal
title('contour of U velocity')
xlabel('Channel Length (X)'); ylabel('Channel Height (Y)')
set(gca,'linewidth',1.2,'fontsize',14)

% Output file
if Parameters.DoSaveTecplot
   A = fopen('fileee.plt','w');
   fprintf(A,[' TITLE= "STEP Output" m  VARIABLES= "I","J","U", "V" m  '...
      'ZONE T="Flow Field", I= %d , J= %d , F=POINT m'],nx,ny);
   for i=1:nx; for j=1:ny
      fprintf(A,'%f %f %f %fm',i,j,U(i,j),V(i,j));
   end; end
   fclose('all');
end
