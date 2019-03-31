function calc_pressure(Rho,U)
% pressure difference calculations
% hyzhou: I don't understand this part!

u_front_fluid = 0;
u_front_solid = 0;
u_back_fluid  = 0;
u_back_solid  = 0;

for iy=40:48
   u_front_fluid = u_front_fluid + U(91,iy);
end
for iy=40:48
   u_front_solid = u_front_solid + U(92,iy);
end
for iy=40:48
   u_back_fluid = u_back_fluid + U(111,iy);
end
for iy=40:48
   u_back_solid = u_back_solid + U(110,iy);
end

u_front_fluid = u_front_fluid/9;
u_front_solid = u_front_solid/9;
u_back_fluid  = u_back_fluid/9;
u_back_solid  = u_back_solid/9;

r_front_fluid = 0;
r_front_solid = 0;
r_back_fluid  = 0;
r_back_solid  = 0;

for iy=49:57
   r_front_fluid = r_front_fluid + Rho(70,iy);
end
for iy=49:57
   r_front_solid = r_front_solid + Rho(71,iy);
end
for iy=49:57
   r_back_fluid = r_back_fluid + Rho(92,iy);
end
for iy=49:57
   r_back_solid = r_back_solid + Rho(91,iy);
end

r_front_fluid = r_front_fluid/9;
r_front_solid = r_front_solid/9;
r_back_fluid  = r_back_fluid/9;
r_back_solid  = r_back_solid/9;

delta_p = 1/3*(r_front_fluid - r_back_fluid);

% fprintf('%f,%f,%f,%f,%f,%f,%f,%f,%f \n',
% u_front_fluid,u_front_solid,u_back_fluid, u_back_solid,r_front_fluid,r_front_solid,r_back_fluid, r_back_solid,delta_p)
end