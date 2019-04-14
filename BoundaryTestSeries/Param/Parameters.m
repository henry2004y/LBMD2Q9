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
      method  double {mustBeMember(method,[1,2,3,4,5])} = 1 % default: 1
 
      % Relaxation time for equilibrium, default: 0.8
      Tau = 0.8      
      
      Uinit   = 0.1	% Uinitial in lattice Units   default: 0.1
      Vinit   = 0.0	% Vinitial in lattice Units   default: 0.0
      Rhoinit = 1.0	% Rhoinitial in lattice Units default: 1.0
      %Rho0In  = 1    % Inlet density?
      Rho0Out = 1    % Outlet density?
 
      % Grid sizes
      nx         double {mustBeInteger} = 400 % default: 200
      ny         double {mustBeInteger} = 100  % default: 50      
      
      % Radius of solid body, default: 11
      R = 21
      
      % Viscosity?
      Nu_physical = 1e-3  % default: 1e-3
      
      channel_height = 0.01 % default: 0.01
      
      Nu = (Parameters.Tau - 0.5)/3
      % Reynolds number
      Re = Parameters.Uinit*2*Parameters.ny/Parameters.Nu
      % 
      Re_cylinder = Parameters.Uinit*2*Parameters.R*3/(Parameters.Tau-0.5)
      % 
      FD = Parameters.R*Parameters.Rhoinit*Parameters.Uinit^2 /...
         105.6430/Parameters.Re_cylinder
      %
      t_lattice = Parameters.channel_height^2/Parameters.ny^2/3*...
         (Parameters.Tau-0.5)/Parameters.Nu_physical
      
      % Timestep
      dt = Parameters.t_lattice
      
      % Spatial interval
      dx = Parameters.channel_height/Parameters.ny
      
      % Body forces
      Gx = 0          % Body force in lattice units in x, gravity
      Gy = 0          % Body force in lattice units in y
      
      
      % Stopping criteria
      Utol = 1e-3     % default: 1e-3
      Vtol = 1e-3     % default: 1e-3
      Ures = 100 % default: 100
      Vres = 100 % default: 100
      
      % Number of timesteps
      nStep = 20000
      
      % Plot interval
      nPlot = 500
      % Tecplot output
      DoSaveTecplot = false
      
   end
   
end