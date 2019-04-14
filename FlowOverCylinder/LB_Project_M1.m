% Lattice Bolzmann Simulation of Flow Around a Cylinder
%
% Milestone 1 :
% No Walls, no obstacle, no inlet/outlet (periodic)

clear; clc
%% Lattice Parameters
nIter = 1000;
Nx = 20;
Ny = 10;
x  = 1:Nx;  % Cells in x-direction
y  = 1:Ny;  % Cells in y-direction
[xx,yy] = meshgrid(x,y);
xx = xx'; yy = yy';
W  = [4/9, 1/9,1/9,1/9,1/9, 1/36,1/36,1/36,1/36];
cx = [  0,   1,  0, -1,  0,    1,  -1,  -1,   1];
cy = [  0,   0,  1,  0, -1,    1,   1,  -1,  -1];

%% Storage Parameters
nInterval = 50;
% Amount of timepoints at which values are recorded
nOut = floor(nIter/nInterval);

%% FLow Parameters
cs   = 1/sqrt(3);
U    = 0.1;
Re   = 50;
visc = 2.*U*Ny/Re;
tau  = visc/cs^2;
beta = 1/(2*tau+1);
alpha = 2; % initial value for the entropic over-relaxation parameter

%% Initialization
rho = ones(Nx,Ny);
u   = U*ones(Nx,Ny);
v   = zeros(Nx,Ny);
feq = zeros(9,Nx,Ny);

% Initialize f and feq
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
   end
   
   % Advection/Free-flight
   for i=1:9
      % periodicity
      f(i,:,:) = circshift(f(i,:,:), [0,cx(i),cy(i)]);
   end
   
   % Relaxation/Collision
   rho = 0; u = 0; v = 0;
   for i = 1:9
      rho = rho + reshape(f(i,:,:),Nx,Ny);
      u = u + reshape(cx(i) * f(i,:,:),Nx,Ny);
      v = v + reshape(cy(i) * f(i,:,:),Nx,Ny);
   end
   u = u./rho; v = v./rho;
   for i = 1:9
      feq(i,:,:) = rho * W(i) ...
         .* (2 - sqrt(1 + 3*u.^2)) ...
         .* (2 - sqrt(1 + 3*v.^2)) ...
         .* ((2*u+sqrt(1+3*u.^2))./(1-u)).^cx(i) ...
         .* ((2*v+sqrt(1+3*v.^2))./(1-v)).^cy(i);
      f(i,:,:) = f(i,:,:) - alpha*beta*( f(i,:,:) - feq(i,:,:) ) ;
   end
   
   % Update Entropic Relaxation Parameter
   %alpha = alpha;
   
   % Update Waitbar
   waitbar(it/nIter,w,['Simulating. - ', num2str(100*it/nIter), '% done']);
end
close(w)


%% Display
uMagOut = sqrt(uOut.^2+vOut.^2);

scrsz = get(0,'ScreenSize');
figure(1)
set(1, 'Position',[1 1 scrsz(3)/2 scrsz(4)/2])

for i = 1:size(rhoOut,3)
   set(1,'Name', ['t = ', num2str(tOut(i)), ...
   ', Re = ', num2str(Re),...
   ' - Press SPACE to advance'],...
   'NumberTitle', 'off')
   if i==1
      subplot(2,1,1)
      % Use surf because contour not rendered for const Zdata
      rho_plot = surf(xx,yy,rhoOut(:,:,i),'edgecolor','none');
      colorbar
      title('Density')
      zlim([0.9 1.1])
      view(2)
      
      subplot(2,1,2)
      u2_plot  = surf(xx,yy,uMagOut(:,:,i),'edgecolor', 'none');
      axis([0 Nx -0.001 Ny])
      view(0,90)
      hold on
      n = Ny/5;
      uv_plot = quiver3(xx(1:n:end,1:n:end),...
         yy(1:n:end,1:n:end),...
         ones(size(xx(1:n:end,1:n:end))),...
         uOut(1:n:end,1:n:end,i),...
         vOut(1:n:end,1:n:end,i),...
         zeros(size(xx(1:n:end,1:n:end))));
      colorbar
      title('Velocity')
      hold off
      
      set(rho_plot, 'ZDataSource', 'rhoOut(:,:,i)')
      set(u2_plot,  'ZDataSource', 'uMagOut(:,:,i)')
      set(uv_plot,  'UDataSource', 'uOut(1:n:end,1:n:end,i)')
      set(uv_plot,  'VDataSource', 'vOut(1:n:end,1:n:end,i)')
   end
   refreshdata(rho_plot)
   refreshdata(u2_plot)
   refreshdata(uv_plot)
   
   pause(0.1);
end
