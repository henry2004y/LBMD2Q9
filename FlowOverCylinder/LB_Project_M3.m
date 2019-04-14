% Lattice Bolzmann Simulation of Flow Around a Cylinder
%
% Milestone 3:
% Bounce back around cylinder, no slip at the walls
% Grad at outlet
% Inlet is a steady flow of constant profile

clear; clc
%% Lattice Parameters
nIter = 1000;
% # of iterations with f_eq constant in order to stabilise initial
% populations
pre_initialisation = 0;
Nx = 200;
Ny = 100;
x  = 1:Nx;  % Cells in x-direction
y  = 1:Ny;  % Cells in y-direction
[x,y] = meshgrid(x,y);
x = x'; y = y';
W  = [4/9, 1/9,1/9,1/9,1/9, 1/36,1/36,1/36,1/36];
% i     1,   2,  3,  4,  5,    6,   7,   8,   9
cx = [  0,   1,  0, -1,  0,    1,  -1,  -1,   1];
cy = [  0,   0,  1,  0, -1,    1,   1,  -1,  -1];
% opp contains the indices of opposite velocities
% i.e: cx(opp(i)) = -cx(i),    cy(opp(i)) = -cy(i)
opp= [  1,   4,  5,  2,  3,    8,   9,   6,   7];
% ymr contains the indices of y-mirrored velocities
ymr= [  1,   2,  5,  4,  3,    9,   8,   7,   6];

%% Storage Parameters
nInterval = 50;
% Amount of timepoints at which values are recorded
nOut = floor(nIter/nInterval); 

%% Obstacle & Wall Parameters
R  = Ny/30;
Ox = Nx/4+2;
Oy = Ny/2+2;
o  = ( (x-Ox).^2 + (y-Oy).^2 ) <= R.^2; % 2D obstacle as 1s and 0s
bodyNode = find(o); % Linear indexes of obstacle nodes
walls = (y == 1) | (y == Ny);
wallNode = find(walls); % linear indexes of wall nodes

%% FLow Parameters
cs   = 1/sqrt(3);
U    = 0.1;
Re   = 50;
visc = 2.*U*R/Re;
tau  = visc/cs^2;
beta = 1/(2*tau+1);
alpha = 2; % initial value for the entropic over-relaxation parameter

%% Initialization
rho = ones(Nx,Ny);
u   = U*ones(Nx,Ny);
v   = 0*ones(Nx,Ny);
feq = zeros(9,Nx,Ny);

% Initialize f_eq, f
for i=1:9
   feq(i,:,:) = rho * W(i) ...
      .* (2 - sqrt(1 + 3*u.^2)) ...
      .* (2 - sqrt(1 + 3*v.^2)) ...
      .* ((2*u+sqrt(1+3*u.^2))./(1-u)).^cx(i) ...
      .* ((2*v+sqrt(1+3*v.^2))./(1-v)).^cy(i);
end
f = feq;

% Outputs
rhoOut = zeros(Nx,Ny,nOut);
uOut   = zeros(Nx,Ny,nOut);
vOut   = zeros(Nx,Ny,nOut);
tOut   = zeros(nOut,1);
iOut   = 0;

%% Main
w = waitbar(0, 'Simulating...');
for it = 1:nIter
   % Output
   if mod(it-1,nInterval) == 0
      iOut = iOut + 1;
      rhoOut(:,:,iOut) = rho;
      uOut  (:,:,iOut) = u;
      vOut  (:,:,iOut) = v;
      tOut  (iOut)     = it;
      
      tmp = rho; tmp(bodyNode) = nan;
      rhoOut(:,:,iOut) = tmp;
      tmp = u; tmp(bodyNode) = nan;
      uOut(:,:,iOut) = tmp;
      tmp = v; tmp(bodyNode) = nan;
      vOut(:,:,iOut) = tmp;
      tOut(iOut) = it;
   end
   
   % Advection/Free-flight
   for i=1:9
      f(i,:,:) = circshift(f(i,:,:), [0,cx(i),cy(i)]);
   end
   
   % Boundary Conditions
   for i=[2 3 6 7]
      f([i opp(i)],bodyNode) = f([opp(i) i],bodyNode);
      f([i opp(i)],wallNode) = f([opp(i) i],wallNode);
   end
   
   % Relaxation/Collision
   rho = reshape(sum(f),Nx,Ny);
   u = reshape((cx * reshape(f,9,Nx*Ny)),Nx,Ny) ./ rho;
   v = reshape((cy * reshape(f,9,Nx*Ny)),Nx,Ny) ./ rho;
   % Inlet
   u(1,2:end-1) = U;
   v(1,2:end-1) = 0;
   % rho?
   % Outlet
   % Domain
   for i = 1:9
      if it >= pre_initialisation
         feq(i,:,:) = rho * W(i) ...
            .* (2 - sqrt(1 + 3*u.^2)) ...
            .* (2 - sqrt(1 + 3*v.^2)) ...
            .* ((2*u+sqrt(1+3*u.^2))./(1-u)).^cx(i) ...
            .* ((2*v+sqrt(1+3*v.^2))./(1-v)).^cy(i);
      end
      tmp = f(i,:,:);
      f(i,:,:) = f(i,:,:) - alpha*beta*( f(i,:,:) - feq(i,:,:) ) ;
      f(i,bodyNode)  = tmp(bodyNode);
      f(i,wallNode)  = tmp(wallNode);
   end
     
   % Update Entropic Relaxation Parameter
   %alpha = alpha;
   
   % Update Waitbar
   waitbar(it/nIter,w,['Simulating. - ', num2str(100*it/nIter), '% done']);
end
close(w)

%% Display
maxrho = max(rhoOut(:));
minrho = min(rhoOut(:));
if minrho == maxrho, maxrho = minrho*1.1; end
uMagOut = sqrt(uOut.^2+vOut.^2);

scrsz = get(0,'ScreenSize');
figure(1)
set(1, 'Position',[1 1 scrsz(3)/2 scrsz(4)/2])

for i = 1:size(rhoOut,3)   
   set(1,'Name', ['t = ', num2str(tOut(i)), ...
      ' - Re = ', num2str(Re),...
      ' - Press SPACE to advance'],...
      'NumberTitle', 'off')
   
   if i==1
      subplot(2,1,1)
      rho_plot = surf(x,y,rhoOut(:,:,i),'edgecolor','none');
      colorbar
      title('Density')
      zlim([minrho maxrho])
      view(2)
      
      subplot(2,1,2)
      u2_plot  = surf(x,y,uMagOut(:,:,i),'edgecolor', 'none');
      colorbar
      axis([0 Nx -0.001 Ny])
      view(0,90)
      hold on
      n = Ny/5;
      uv_plot = quiver3(x(1:n:end,1:n:end),...
         y(1:n:end,1:n:end),...
         ones(size(x(1:n:end,1:n:end))),...
         uOut(1:n:end,1:n:end,i),...
         vOut(1:n:end,1:n:end,i),...
         zeros(size(x(1:n:end,1:n:end))));
      title('Velocity')
      hold off
      
      rho_plot.ZDataSource = 'rhoOut(:,:,i)';
      u2_plot.ZDataSource = 'uMagOut(:,:,i)';
      uv_plot.UDataSource = 'uOut(1:n:end,1:n:end,i)';
      uv_plot.VDataSource = 'vOut(1:n:end,1:n:end,i)';
   end
   refreshdata(rho_plot)
   refreshdata(u2_plot)
   refreshdata(uv_plot)
   
   pause(0.1);
end
