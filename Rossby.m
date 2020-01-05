Nx = 180;
Ny = 88;

dy = 111 * 1e3;
dx = dy;

xmin = 0;
xmax = Nx * dx;
ymin = 0;
ymax = Ny * dy;

x = (xmin - dx:dx:xmax + dx)';
y = (ymin - dy:dy:ymax + dy)';

[X, Y] = meshgrid(x, y);

i = 2:Nx + 2;
j = 2:Ny + 2;

t = 0
tmax = 3600 * 24 * 30;
dt_output = 3600;

g = 9.8;
Omega = 7.2921 * 1e-5;
a = 6.371 * 10^6;
D = 1;
cfl = 0.5;

%---Topo---%


A = ones(500, 500);
B = A(90:90 + 89, 1:182);

%xtopo = linspace(0, Nx * dx, size(B, 2));
%ytopo = linspace(0, Ny * dy, size(B, 1));

%[Xtopo, Ytopo] = meshgrid(xtopo, ytopo);

%B = interp2(Xtopo, Ytopo, B, X, Y, 'spline');
%B = (B > 0) .* B;

%B = filter2('gaussian', B');
%B = filter2('gaussian', B') * 1.5077e-06;

%B(:, 1) = B(:, end - 2);

%---Topo---%

phi = linspace(0, 88, Ny + 3)' * pi / 180;

phi0 = median(phi);
beta = 2 * Omega * cos(phi0) / a;
f0 = 2 * Omega * sin(phi0);
f = f0 + beta .* (Y - mean(y));
tanphi = tan(phi);

U0 = zeros(Ny + 3, Nx + 3, 3);
total_height = repmat(linspace(5700, 5500, Ny + 3)', 1, Nx + 3);
H = total_height - B;
U0(:, :, 1) = H;

[dh_dx, dh_dy] = gradient(U0(:, :, 1));
U0(:, :, 2) =  - g ./ f .* dh_dy / dx;
U0(:, :, 3) = g ./ f .* dh_dx / dy;

zeta0 = 0;
q = (f + zeta0) ./ H;

U = U0;
n = 1;

global B;

function [F, G] = flux(U)
    global g B
    eta = U(:, :, 1) - B;
    K = 0.5 * (U(:, :, 2).^2 + U(:, :, 3).^2);
    F(:, :, 1) = U(:, :, 1) .* U(:, :, 2);
    F(:, :, 2) = K + g * eta;
    F(:, :, 3) = 0;
    G(:, :, 1) = U(:, :, 1) .* U(:, :, 3);
    G(:, :, 2) = 0;
    G(:, :, 3) = F(:, :, 2);
endfunction

function S = source(U)
    global q a tanphi
    S = zeros(size(U));
    S(:, :, 1) = 0;
    S(:, :, 2) = U(:, :, 1) .* q .* U(:, :, 3)...
     + tanphi .* U(:, :, 2) .* U(:, :, 3) / a;
    S(:, :, 3) =  - U(:, :, 1) .* q .* U(:, :, 2)...
     - tanphi .* U(:, :, 2).^2 / a;
end

function U = bc(U)
    global U0
    U(:, end, 1) = U(:, end - 2, 1);
    U(1, :, 1) = U0(1, :, 1);
    U(end, :, 1) = U0(end, :, 1);
    U(:, end, 2) = U(:, end - 2, 2);
    U(1, :, 2) = U(3, :, 2);
    U(end, :, 2) = U(end - 2, :, 2);
    U(:, end, 3) = U(:, end - 2, 3);
    U(:, 1, :) = U(:, end, :);
    U(1, :, 3) = 0;
    U(end, :, 3) = 0;
end

while t < tmax
    maxu = max(abs(U(:, :, 2)) + sqrt(g * U(:, :, 1)));
    maxv = max(abs(U(:, :, 3)) + sqrt(g * U(:, :, 1)));
    dt = cfl * dx / max([maxu, maxv]);
    
    %U = bc(U);
    
    U(:, end, 1) = U(:, end - 2, 1);
    U(1, :, 1) = U0(1, :, 1);
    U(end, :, 1) = U0(end, :, 1);
    U(:, end, 2) = U(:, end - 2, 2);
    U(1, :, 2) = U(3, :, 2);
    U(end, :, 2) = U(end - 2, :, 2);
    U(:, end, 3) = U(:, end - 2, 3);
    U(:, 1, :) = U(:, end, :);
    U(1, :, 3) = 0;
    U(end, :, 3) = 0;
    
    %[F, G] = flux(U);
    
    eta = U(:, :, 1) - B;
    K = 0.5 * (U(:, :, 2).^2 + U(:, :, 3).^2);
    F(:, :, 1) = U(:, :, 1) .* U(:, :, 2);
    F(:, :, 2) = K + g * eta;
    F(:, :, 3) = 0;
    G(:, :, 1) = U(:, :, 1) .* U(:, :, 3);
    G(:, :, 2) = 0;
    G(:, :, 3) = F(:, :, 2);
    
    %[S] = source(U);
    
    S = zeros(size(U));
    S(:, :, 1) = 0;
    S(:, :, 2) = U(:, :, 1) .* q .* U(:, :, 3)...
     + tanphi .* U(:, :, 2) .* U(:, :, 3) / a;
    S(:, :, 3) =  - U(:, :, 1) .* q .* U(:, :, 2)...
     - tanphi .* U(:, :, 2).^2 / a;
    
    U_iph(j, i, :) = 0.5 * (U(j, i + 1, :) + U(j, i, :)) - ...
    (0.5 * dt / dx) * (F(j, i + 1, :) - F(j, i, :));
    U_jph(j, i, :) = 0.5 * (U(j + 1, i, :) + U(j, i, :)) - ...
    (0.5 * dt / dy) * (G(j + 1, i, :) - G(j, i, :));

    [F_iph, ~] = flux(U_iph);
    [~, G_jph] = flux(U_jph);
    U(j, i, :) = U(j, i, :) + dt * S(j, i, :) - ...
    (dt / dx) * (F_iph(j, i, :) - F_iph(j, i - 1, :)) - ...
    (dt / dy) * (G_jph(j, i, :) - G_jph(j - 1, i, :));

    t = t + dt;
    
    %U  =  bc(U);
    %[F, G] = flux(U);
    %S = source(U);
    %U(j, i, :) = 0.5 * (U(j, i + 1, :) + U(j, i - 1, :))...
    % - 0.5 * dt / dx * (F(j, i + 1, :) - F(j, i - 1, :));
    %U(j, i, :) = 0.5 * (U(j + 1, i, :) + U(j - 1, i, :))...
    % - 0.5 * dt / dy * (G(j + 1, i, :) - G(j - 1, i, :)) + dt * S(j, i, :);

end