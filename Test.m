Nx = 180;
Ny = 88;

dy = 111 * 1e3;
dx = dy;

xmin = 0;
xmax = Nx * dx ;
ymin = 0;
ymax = Ny * dy ;

x = (xmin - dx : dx : xmax + dx)';
y = (ymin - dy : dy : ymax + dy)';
[X,Y] = meshgrid(x, y);

i = 2 : Nx + 2;
j = 2 : Ny + 2;

t = 0;
tmax = 3600 * 24 * 30;
dt_output = 3600;

g = 9.8;
Omega = 7.2921 * 1e-5;
a = 6.371 * 10^6;
D = 1;
cfl = 0.5;

phi = linspace(0, 88, Ny + 3)' * pi / 180;
phi0 = median(phi);
beta = 2 * Omega * cos(phi0) / a;
f0 = 2 * Omega * sin(phi0);
f = f0 + beta .* (Y - mean(y));