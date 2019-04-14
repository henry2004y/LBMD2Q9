% Lattice Bolzmann Simulation of Flow Around a Cylinder
%
% Milestone VIII
% Bounce back around cylinder, no slip at the walls
% Grad at outlet
% Inlet is a steady flow of constant profile
% Moving Cylinder
% Grad at BB and Filling new fluid cells

clear; clc
%% Lattice Parameters
nIter  = 500;
% # of iterations with f_eq constant in order to stabilize initial
% populations
nPreInitStep = 0;
Nx = 400;
Ny = 200;
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
nOut = 100; % Amount of timepoints at which values are recorded
nInterval = floor(nIter/nOut);

%% Obstacle & Wall Parameters
R  = Ny/40;
Ox = Nx/4+2; Ax = 0;
Oy = Ny/2+2; Ay = 0.25;
F = 1.2; T = 570/F; % period of 1 movement cycle
o  = ( (x-Ox).^2 + (y-Oy).^2 ) <= R.^2; % 2D obstacle as 1s and 0s
bodyNode = find(o); % Linear indexes of obstacle nodes
walls = (y == 1) | (y == Ny);
wallNode = find(walls); % linear indexes of wall nodes
oi = o;
for i=1:9
   oi  = oi | circshift(o,[cx(i) cy(i)]);
end
border = find(oi-o); % linear indexes of border nodes

%% FLow Parameters
cs   = 1/sqrt(3);
U    = 0.1;
Re   = 200;
visc = 2.*U*R/Re;
tau  = visc/cs^2;
beta = 1/(2*tau+1);
alpha = 2; % initial value for the entropic over-relaxation parameter

%% Initialization
rho = ones(Nx,Ny);
u  = U*ones(Nx,Ny);
v  = 0*ones(Nx,Ny);
feq = zeros(9,Nx,Ny);

% Initialize feq, f ,fgr
for i=1:9
   feq(i,:,:) = rho * W(i) ...
      .* (2 - sqrt(1 + 3*u.^2)) ...
      .* (2 - sqrt(1 + 3*v.^2)) ...
      .* ((2*u+sqrt(1+3*u.^2))./(1-u)).^cx(i) ...
      .* ((2*v+sqrt(1+3*v.^2))./(1-v)).^cy(i);
end
f = feq;
fgr_next = feq;

% Outputs
rhoOut = zeros(Nx,Ny,nOut);
uOut   = zeros(Nx,Ny,nOut);
vOut   = zeros(Nx,Ny,nOut);
tOut   = zeros(nOut,1);
iOut   = 0;
FxOut  = zeros(nOut,1);
FyOut  = zeros(nOut,1);
vprobeOut = zeros(nOut,1);
OxOut  = zeros(nOut,1);
OyOut  = zeros(nOut,1);

% Initial Force Probe
Fx = 0; Fy = 0;
Ou = 0; Ov = 0; % TODO: organise

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

      tmp = rho; tmp([bodyNode; wallNode]) = nan;
      rhoOut(:,:,iOut) = tmp;
      tmp = u; tmp([bodyNode; wallNode]) = nan;
      uOut(:,:,iOut) = tmp;
      tmp = v; tmp([bodyNode; wallNode]) = nan;
      vOut(:,:,iOut) = tmp;
      tOut(iOut) = it;
      
      vprobeOut(iOut) = v(100,50);
      FxOut(iOut) = Fx; FyOut(iOut) = Fy;
      OxOut(iOut) = Ox; OyOut(iOut) = Oy;
   end
   
   % Advection/Free-flight
   for i=1:9
      f(i,:,:) = circshift(f(i,:,:), [0,cx(i),cy(i)]);
   end
   
   % Boundary Conditions
%    for i=1:9
%       f(i,wallNode) = f(ymr(i),wallNode);
%    end
%    for n=border
%       for i=1:9
%          if o(n+cx(i)+cy(i)*Nx) % if the ith velocity leads to a solid node
%             j = opp(i);
%             ns = n+cx(i)+cy(i)*Nx; % that solid node's index
%             f(j,n) = f(i,ns);
%          end
%       end
%    end
   for i=[2 3 6 7]
      f([i opp(i)],bodyNode) = f([opp(i) i],bodyNode);
      f([i opp(i)],wallNode) = f([opp(i) i],wallNode);
   end

   % TODO: Decide : Recalculate uu, vv, rho here?
   rhon = reshape(sum(f),Nx,Ny);
   un = reshape((cx * reshape(f,9,Nx*Ny)),Nx,Ny) ./ rhon;
   vn = reshape((cy * reshape(f,9,Nx*Ny)),Nx,Ny) ./ rhon;
   for n=border
      % interpolate velocity derivatives
      % TODO: This doesn't work if a boundary node is trapped
      %       between two solid nodes!
      for i=2:5
         j = opp(i);
         nf = n+cx(i)+cy(i)*Nx; %guess
         if o(nf)
            nf = n+cx(j)+cy(j)*Nx; %correction
            direction = -1;
         else
            direction = 1;
         end
         if cy(i) == 0
            dudx = cx(i) * (u(n) - u(nf)) * direction;
            dvdx = cx(i) * (v(n) - v(nf)) * direction;
         end
         if cx(i) == 0
            dudy = cy(i) * (u(n) - u(nf)) * direction;
            dvdy = cy(i) * (v(n) - v(nf)) * direction;
         end
      end
      
      nDbar = 0; uutgt = 0; vvtgt = 0; rhobb = 0; rhos = 0;
      for i=1:9
         j = opp(i);
         ns = n+cx(i)+cy(i)*Nx;
         nf = n+cx(j)+cy(j)*Nx;
         if o(ns)
            nDbar = nDbar+1;
            q(j) = 1; %TODO: interpolate?
            uutgt = uutgt + (q(j)*un(nf)+Ou)/(1+q(j));
            vvtgt = vvtgt + (q(j)*vn(nf)+Ov)/(1+q(j));
            
            rhobb = rhobb + f(i,ns);
            rhos  = rhos  + 6*W(j)*rhon(n)*(cx(j)*Ou+cy(j)*Ov); % TODO: rho0 ?
         else
            rhobb = rhobb + f(j,n);
         end
      end
      uutgt = uutgt/nDbar;
      vvtgt = vvtgt/nDbar;
      rhotgt = rhobb + rhos;
      
      Pxxeq = rhotgt*cs^2 + rhotgt.*uutgt.^2;
      Pyyeq = rhotgt*cs^2 + rhotgt.*vvtgt.^2;
      Pxyeq =               rhotgt.*uutgt.*vvtgt; % = Pyxeq
      Pxx1  = -rhotgt*cs^2/2/beta * 2*dudx;
      Pyy1  = -rhotgt*cs^2/2/beta * 2*dvdy;
      Pxy1  = -rhotgt*cs^2/2/beta * (dudy + dvdx); % = Pyx1
      Pxx(n) = Pxxeq + Pxx1;
      Pyy(n) = Pyyeq + Pyy1;
      Pxy(n) = Pxyeq + Pxy1; % = Pyx
      
      for i=1:9
         if o(n+cx(i)+cy(i)*Nx)
            j = opp(i);
            ns = n+cx(i)+cy(i)*Nx;
            nf = n+cx(j)+cy(j)*Nx;
            
            f(j,n) = W(j) .* ...
               ( ...
               + rhotgt ...
               + rhotgt.*uutgt * cx(j) / cs^2 ...
               + rhotgt.*vvtgt * cy(j) / cs^2 ...
               + 0.5/cs^4*( (Pxx(n) - rhotgt*cs^2)*(cx(j)*cx(j) - cs^2) ...
               +(Pyy(n) - rhotgt*cs^2)*(cy(j)*cy(j) - cs^2) ...
               +2*(Pxy(n))*(cx(j)*cy(j)) ...
               ) ...
               );
            
         end
      end
   end
   
   % Sum up the force
   Fx = 0; Fy = 0;
   for n=border
      for i=1:9
         if o(n+cx(i)+cy(i)*Nx) % if the ith velocity leads to a solid node
            ns = n+cx(i)+cy(i)*Nx; % that solid node's index
            Fx = Fx + cx(i)*(f(opp(i),n)+f(i,ns));
            Fy = Fy + cy(i)*(f(opp(i),n)+f(i,ns));
         end
      end
   end
   
   % Relaxation/Collision
   rho = reshape(sum(f),Nx,Ny);
   u = reshape((cx * reshape(f,9,Nx*Ny)),Nx,Ny) ./ rho;
   v = reshape((cy * reshape(f,9,Nx*Ny)),Nx,Ny) ./ rho;
   % Inlet
   u(1,2:end-1) = U;
   v(1,2:end-1) = 0;
   rho(1,2:end-1) = 1; % rho before advec -> constant
   % Outlet
   Pxx  = reshape(((cx.*cx) * reshape(f,9,Nx*Ny)),Nx,Ny) ./ rho;
   Pyy  = reshape(((cy.*cy) * reshape(f,9,Nx*Ny)),Nx,Ny) ./ rho;
   Pxy  = reshape(((cx.*cy) * reshape(f,9,Nx*Ny)),Nx,Ny) ./ rho; % = Pyx
   fgr = fgr_next; % grad's population carried from the prev. timestep
   for i = 1:9
      fgr_next(i,:,:) = W(i) .* ...
         ( ...
         + rho ...
         + rho.*u * cx(i) / cs^2 ...
         + rho.*v * cy(i) / cs^2 ...
         + 1/2/cs^4 * ( (Pxx - rho*cs^2)*(cx(i)*cx(i) - cs^2) ...
         +(Pyy - rho*cs^2)*(cy(i)*cy(i) - cs^2) ...
         +2*(Pxy)*(cx(i)*cy(i)) ...
         ) ...
         );
   end
   % Cylinder
   u(bodyNode) = Ou;
   v(bodyNode) = Ov;
   % Domain
   for i = 1:9
      tmp = f(i,:,:); % save a non-collided copy of f
      
      if it >= nPreInitStep
         feq(i,:,:) = rho * W(i) ...
            .* (2 - sqrt(1 + 3*u.^2)) ...
            .* (2 - sqrt(1 + 3*v.^2)) ...
            .* ((2*u+sqrt(1+3*u.^2))./(1-u)).^cx(i) ...
            .* ((2*v+sqrt(1+3*v.^2))./(1-v)).^cy(i);
      end
      f(i,:,:) = f(i,:,:) - alpha*beta*( f(i,:,:) - feq(i,:,:) );
      
      f(i,bodyNode)  = tmp(bodyNode); % no collision in solid
      f(i,wallNode)  = tmp(wallNode); % no collision inside wall boundary
      f(i,1,2:end-1) = feq(i,1,2:end-1); % inlet
      if cx(i) < 0
         f(i,end,2:end-1) = fgr(i,end,2:end-1); % outlet
      end
      %f(i,bb)     = feq(i,bb); %cylinder
   end
   
   % Move Cylinder
   Ox = Nx/4+2+round(Ax*2*R*(1-cos(2*pi*it/T)));
   Oy = Ny/2+2+round(Ay*2*R*sin(2*pi*it/T));
   Ou = 2*pi/T*Ax*2*R*sin(2*pi*it/T);
   Ov = 2*pi/T*Ay*2*R*cos(2*pi*it/T);
   o  = ( (x-Ox).^2 + (y-Oy).^2 ) <= R.^2; % 2D obstacle as 1s and 0s
   bodyNode = find(o);
   oi = o;
   for i=1:9
      oi  = oi | circshift(o,[cx(i) cy(i)]);
   end
   border = find(oi-o); % Linear indexes of obstacle nodes
   % Missing populations?
   
   % Update Entropic Relaxation Parameter
   %alpha = alpha;
   
   % Update Waitbar
   waitbar(it/nIter,w,['Simulating. - ', num2str(100*it/nIter), '% done']);
end
close(w)

%% Display
Cd = FxOut/0.5/U^2/ceil(2*R);
maxrho = max(rhoOut(:));
minrho = min(rhoOut(:));
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
      subplot(2,1,1); hold on
      cyl_plot = cylinder3D([OxOut(i) OyOut(i)],[minrho maxrho],R,10);
      rho_plot = surf(x,y,rhoOut(:,:,i),'edgecolor','none');
      colorbar
      title('Density')
      %zlim([minrho maxrho])
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
