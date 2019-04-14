classdef Parameters
   %Parameters for 2D LBM simulation.
   
   properties (Constant)
      % Lattice D2Q9
      ex = [1  0 -1  0  1 -1 -1  1  0]
      ey = [0  1  0 -1  1  1 -1 -1  0]
      w  = [1/9 1/9 1/9 1/9 1/36 1/36 1/36 1/36 4/9]
      cs = 1/sqrt(3) % Cs sound speed
      % Matrix used in force calculation in convergence for finding 
      % a=0,1,2,...9
      matrix_f = [7 4 8;3 9 1;6 2 5]

      % Boundary method:
      % 1. bounceback 
      % 2. flippova
      % 3. mei
      % 4. bozidi
      % 5. yu
      method  double {mustBeMember(method,[1,2,3,4,5])} = 3 % default: 1
 
      % Relaxation time for equilibrium, default: 0.8
      Tau = 0.8      
      
      % Velocities are normalized to c = dx/dt
      Uinit   = 0.25	% Uinitial in lattice Units   default: 0.1
      Vinit   = 0.01	% Vinitial in lattice Units   default: 0.0
      Rhoinit = 1.0	% Rhoinitial in lattice Units default: 1.0
      %Rho0In  = 1    % Inlet density?
      Rho0Out = 1    % Outlet density?
 
      % Grid sizes
      nx         double {mustBeInteger} = 400 % default: 200
      ny         double {mustBeInteger} = 100  % default: 50      
      
      % Radius of solid body, [dx] default: 11
      R = 11
      
      % Kinematic viscosity, [m^2/s]
      Nu_SI = 1e-3   % default: 1e-3
      
      % Channel height, [m]
      height_SI = 0.01  % default: 0.01
      
      % Spatial interval
      dx = Parameters.height_SI/Parameters.ny      
      
      % Timestep (uses the viscosity formula to determine timestep)
      dt = Parameters.dx^2/3*(Parameters.Tau-0.5)/Parameters.Nu_SI

      % Reynolds number of the fluid in the pipe
      Re = Parameters.Uinit*(Parameters.dx/Parameters.dt)*...
         Parameters.height_SI/Parameters.Nu_SI
      % Reynolds number of the cyclinder in the fluid
      Re_cylinder = Parameters.Uinit*(Parameters.dx/Parameters.dt)*...
         2*Parameters.R*Parameters.dx/Parameters.Nu_SI
      % 
      FD = Parameters.R*Parameters.Rhoinit*Parameters.Uinit^2 /...
         105.6430/Parameters.Re_cylinder      
      
      % Body forces
      Gx = 0      % Body force in lattice units in x, gravity
      Gy = 0      % Body force in lattice units in y
      
      
      % Stopping criteria
      Utol = 1e-3    % default: 1e-3
      Vtol = 1e-3    % default: 1e-3
      Ures = 100     % default: 100
      Vres = 100     % default: 100
      
      % Number of timesteps
      nStep = 5000
      
      % Plot interval
      nPlot = 100
      nPlotSave = 30000
      % Tecplot output
      DoSaveTecplot = false
      % Log
      nLog = 50
      
   end
   
end