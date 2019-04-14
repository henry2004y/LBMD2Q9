function [plot_handle] = cylinder3D(ocoord,zminmax,radius,resolution)
ox = ocoord(1);
oy = ocoord(2);
zmin = zminmax(1);
zmax = zminmax(2);
R = radius;
n = resolution;

theta = linspace(0,2*pi,n);
z = [zmin linspace(zmin,zmax,n) zmax];
[th,zz] = meshgrid(theta,z);
rr = R + 0*zz;
rr([1 end],:) = 0;

X = ox + rr.*cos(th);
Y = oy + rr.*sin(th);
Z = zz;

plot_handle = surf(X,Y,Z,'edgecolor','none');


end

