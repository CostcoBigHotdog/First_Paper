%%
%  Version 8_26
%  Correct the forgetting of dt in 8_25 ihx part
%


%% Kinetics parameters
D = 1.43;       % Diffusion coefficient 
v = 88.46e6;
sum_a = 0.004835; 
sum_f = 0.001782; 
rho_0 = 1.7141;
rho = rho_0;
lifetime = 4.48274e-7;
fisnum = 2.92;
ln = 2.95630e-7;

lamb1 = 0.0127095;
lamb2 = 0.0301051; 
lamb3 = 0.112269; 
lamb4 = 0.327525; 
lamb5 = 1.23099;
lamb6 = 8.1784;

beta1 = 6.75e-5;
beta2 = 6.66e-4;
beta3 = 5.11e-4;
beta4 = 1.43e-3; 
beta5 = 6.43e-4; 
beta6 = 1.7e-4;
beta = beta1 + beta2 + beta3 + beta4 + beta5 + beta6;

% reactivity feedback
% a_d = -834e-5;  % Doppler
a_d = -400e-5;  % Doppler

% a_f = -0.2625e-5; % Fuel expansion
a_f = -0.02625e-5; % Fuel expansion

a_ext = -1e-5;


%% Core geometry & net
L = 47 * 2.54;          % Core height
nx = 100;               % Number of spatial points
dx = L / (nx - 1);     % Spatial step size
dt = 1e-3;              % Time step size
nt = 40000;             % Number of time points
T = dt *  (nt - 1);     % Time domain length
x = linspace(0, L, nx);     % Spatial grid
t = linspace(0, T, nt);     % Time grid

% Fuel and cladding size
r = 94.2 / 2 * 2.54;
A_core = pi * r * r;
r_f = 0.272;
r_cl = 0.328;
k_cl = 0.11;
A_f = pi * r_f * r_f;
A_cl = pi * r_cl * r_cl - pi * r_f * r_f;
n_pin = 11382;
W_na = 5.4e6;
v_na = 600;
A_na = W_na/v_na /n_pin;


% Hotpool 
% nedd h >= a*dt, 
L_hotPool = 792;
nx_hotPool = 100;
dx_hotPool = L_hotPool / (nx_hotPool-1);
k_na = 1.42;

% IHX
L_ihx = 536;
nx_ihx = 100;
dx_ihx = L_ihx / (nx_ihx-1);

A_ihxP = 4961.882;
A_ihxS = 4233.779;
A_ihxH = 10667.358;

%% Initial condition
% Core kinetics
p = zeros(nx, nt);  % relative power
p(:, 1) = 1;
load('p_initial.mat');
% p(:, 1) = p_initial;
p_0 = zeros(nx, 1);
% p_0(:, 1) = 425e6 / n_pin / L;
p_0(:, 1) = 840e6 / n_pin / L;
p_real = zeros(nx, nt);
p_real(:, 1) = p(:, 1) .* p_0(:, 1) .* n_pin;
p_fractional = zeros(1, nt);
p_fractional(1, 1) = 1;


c1 = zeros(nx, nt); 
c1(:, 1) = beta1 * fisnum * sum_f / lamb1 * p(:, 1); 
c2 = zeros(nx, nt); 
c2(:, 1) = beta2 * fisnum * sum_f / lamb2 * p(:, 1); 
c3 = zeros(nx, nt); 
c3(:, 1) = beta3 * fisnum * sum_f / lamb3 * p(:, 1); 
c4 = zeros(nx, nt); 
c4(:, 1) = beta4 * fisnum * sum_f / lamb4 * p(:, 1); 
c5 = zeros(nx, nt); 
c5(:, 1) = beta5 * fisnum * sum_f / lamb5 * p(:, 1); 
c6 = zeros(nx, nt); 
c6(:, 1) = beta6 * fisnum * sum_f / lamb6 * p(:, 1); 

% Fuel temperature
T_f = zeros(nx, nt);
load('T_fInitial.mat');
% T_f(:, 1) = 810;
T_f(:, 1) = T_fInitial;

% Cladding temperature
T_cl = zeros(nx, nt);
load('T_clInitial.mat');
% T_cl(:, 1) = linspace(320, 470, nx);
T_cl(:, 1) = T_clInitial;

% Sodium coolant temperature
T_na = zeros(nx, nt);
load('T_naInitial.mat');
% T_na(:, 1) = linspace(320, 470, nx);
T_na(:, 1) = T_naInitial;

% Hotpool temperature
T_hotPool = zeros(nx_hotPool, nt);
T_hotPool(:, 1) = 468;

% Ihx temperature
T_ihxP = zeros(nx_ihx, nt);
load('T_ihxPInitial.mat');
T_ihxP(:, 1) = T_ihxPInitial;

T_ihxS = zeros(nx_ihx, nt);
load('T_ihxSInitial.mat');
T_ihxS(:, 1) = T_ihxSInitial;


%% PDE parameters and Boundary conditions

% Core kinetics

% Fuel heat transfer
rhoCA_fuel = 16.2 * 0.34 * A_f *(840/425);

% Cladding heat transfer
rhoCA_clad = 6.551 * 0.33 * A_cl *(840/425);
% h_clad = 2 * pi * r_f * 0.2 *(840/425);
h_clad = 2 * pi * r_f * 0.4 *(840/425);
A2 =   (1 - 2 * k_cl * A_cl * (1/rhoCA_clad) * dt / dx^2) * eye(nx) ...
     + (k_cl * A_cl  * (1/rhoCA_clad) * dt / dx^2) * diag(ones(nx-1, 1), 1) ...
     + (k_cl * A_cl  * (1/rhoCA_clad) * dt / dx^2) * diag(ones(nx-1, 1), -1);

% Core sodium coolant heat transfer
rhoCA_na = 0.927 * 1.228 * A_na ;
h_na = 2 * pi * r_cl * 3.9 ;    % A * h
A3 =   (1 - (v_na * dt /dx)) * eye(nx) ...
     + (v_na * dt /dx) * diag(ones(nx-1, 1), -1);
core_inlet = zeros(nx, 1); 
core_inlet(1) = 360 * (v_na * dt /dx);

% Hotpool heat conduction
W_hotPool = 500 ;
A4 =   (1 - (W_hotPool * dt /dx_hotPool)) * eye(nx_hotPool) ...
     + (W_hotPool * dt /dx_hotPool) * diag(ones(nx_hotPool-1, 1), -1);
hotPool_inlet = zeros(nx_hotPool, 1); 

% IHX heat transfer
rhoCA_ihxP = 0.927 * 1.228 * A_ihxP *(840/425);
rhoCA_ihxS = 0.927 * 1.228 * A_ihxS *(840/425);
WC_ihxP = -1126420.13 * 1.228 *(840/425);
% WC_ihxS = 1152879.67 * 1.228 *(840/425);
WC_ihxS = 1272879.67 * 1.228 *(840/425);
h_ihx = 0.2815 * A_ihxH *(840/425) * 1000;

A5 =   (1 + (WC_ihxP / rhoCA_ihxP * dt / dx_ihx)) * eye(nx_hotPool) ...
     - (WC_ihxP / rhoCA_ihxP * dt / dx_ihx) * diag(ones(nx_hotPool-1, 1), 1);
ihxP_inlet = zeros(nx_ihx, 1);

A6 =   (1 - (WC_ihxS / rhoCA_ihxS * dt / dx_ihx)) * eye(nx_hotPool) ...
     + (WC_ihxS / rhoCA_ihxS * dt / dx_ihx) * diag(ones(nx_hotPool-1, 1), -1);
ihxS_inlet = zeros(nx_ihx, 1);


%% Time iteration
% Start timing
tic;
% Time stepping using the backward Euler method
for k = 2:nt
    
    % Core kinetics update
    A = (1 - (( rho - beta ) / lifetime - (1/ln) - sum_a * v) * dt + 2 * D * v * dt / dx^2) * eye(nx) ...
            - (D * v * dt / dx^2) * diag(ones(nx-1, 1), 1) ...                                                                                                
            - (D * v * dt / dx^2) * diag(ones(nx-1, 1), -1);
    % A(1, 1)= A(1, 1) + (-D * v * dt / dx^2) / (1 - (-1 / (3 * 0.7104 * D)) * dx );
    % A(nx, nx)= A(nx, nx) + (-D * v * dt / dx^2) / (1 - ((-1 / (3 * 0.7104 * D)) * dx));
    
    f_c1 = beta1 * fisnum * sum_f .* p(:, k-1) - lamb1 * c1(:, k-1);
    c1(:, k) = c1(:, k-1) + dt * f_c1;
    f_c2 = beta2 * fisnum * sum_f .* p(:, k-1) - lamb2 * c2(:, k-1);
    c2(:, k) = c2(:, k-1) + dt * f_c2;
    f_c3 = beta3 * fisnum * sum_f .* p(:, k-1) - lamb3 * c3(:, k-1);
    c3(:, k) = c3(:, k-1) + dt * f_c3;
    f_c4 = beta4 * fisnum * sum_f .* p(:, k-1) - lamb4 * c4(:, k-1);
    c4(:, k) = c4(:, k-1) + dt * f_c4;
    f_c5 = beta5 * fisnum * sum_f .* p(:, k-1) - lamb5 * c5(:, k-1);
    c5(:, k) = c5(:, k-1) + dt * f_c5;
    f_c6 = beta6 * fisnum * sum_f .* p(:, k-1) - lamb6 * c6(:, k-1);
    c6(:, k) = c6(:, k-1) + dt * f_c6;

    f_p = v * (lamb1 * c1(:, k) + lamb2 * c2(:, k) + lamb3 * c3(:, k) + lamb4 * c4(:, k) + lamb5 * c5(:, k) + lamb6 * c6(:, k));
    b = p(:, k-1) + dt * f_p;
    p(:, k) = A \ b;
    p_real(:, k) = p(:, k) .* p_0 .* n_pin;

    current_power = sum(p(:,k));
    p_fractional(1, k) = current_power/nx;

    % rho_ext = a_ext * (current_power- 0.5 * nx);
    % rho = rho + rho_ext;


    % Fuel temperature update
    T_f(:,k) =  T_f(:,k-1) + dt .* ( p(:,k-1) .* p_0 + h_clad .* (T_cl(:, k-1) - T_f(:, k-1)) ) ./ rhoCA_fuel ;

    % Reactivity change
    if k <= 50000     
         rho_ext = a_ext * (current_power- nx);
         rho = rho + rho_ext;
    elseif (k > 50000) && (k <= 60000)
         % rho = rho_0 + 5e-8 * (k-50000) + a_d * log(T_f(50, k) / T_f(50, 50000)) + a_f * (T_f(50, k) - T_f(50, 50000));
    elseif k > 60000
         % rho = rho_0 + 50e-5 + a_d * log(T_f(50, k) / T_f(50, 50000)) + a_f * (T_f(50, k) - T_f(50, 50000));
    end
    
    % Cladding temperature update
    f_cl = h_clad / rhoCA_clad .* (T_f(:, k-1) - T_cl(:, k-1)) + ...
           h_na / rhoCA_clad .* (T_na(:, k-1) - T_cl(:, k-1));
    T_cl(:, k) = A2 * T_cl(:, k-1) + dt .* f_cl;

    % Sodium coolant temperature update
    f_na = h_na/rhoCA_na .* (T_cl(:, k-1) - T_na(:, k-1)); 
    T_na(:, k) = A3 * T_na(:, k-1) +  dt .* f_na + core_inlet;

    % Hotpool inlet/all position temperature update
    hotPool_inlet(1) = (W_hotPool * dt / dx_hotPool) * T_na(nx, k-1);
    T_hotPool(:, k) = A4 * T_hotPool(:, k-1) + hotPool_inlet;

    % IHX temperature update
    % Primary side
    ihxP_inlet(nx_ihx) = - (WC_ihxP / rhoCA_ihxP * dt / dx_ihx) * T_hotPool(nx_hotPool, k-1);
    T_ihxP(:, k) = A5 *  T_ihxP(:, k-1) + ihxP_inlet + dt * h_ihx / rhoCA_ihxP .* (T_ihxS(:, k-1) - T_ihxP(:, k-1));
    % Secondary side
    ihxS_inlet(1) =  (WC_ihxS / rhoCA_ihxS * dt / dx_ihx) * 326;
    T_ihxS(:, k) = A6 *  T_ihxS(:, k-1) + ihxS_inlet + dt * h_ihx / rhoCA_ihxS .* (T_ihxP(:, k-1) - T_ihxS(:, k-1));

    % core_inlet(1) = T_ihxP(1, k) * (v_na * dt /dx); 
  
end
% Stop timing
elapsedTime = toc;
fprintf('Elapsed time: %.6f seconds\n', elapsedTime);


%% Plot
figure(1)
y = p_fractional(1,:);
x_ms = 1:length(y);     
x_s = x_ms / 1000;
plot(x_s, y);
ylim([0.3 2]);
xlabel('Time (s)');
ylabel('Fractional Power');

figure(2)
hold on

% Plot each line with distinct colors and line styles
plot(x_s, T_na(100, :), 'r-', 'DisplayName', 'Core Coolant outlet Temperature'); % 红色实线
plot(x_s, T_f(50, :), 'b--', 'DisplayName', 'Fuel Temperature (center)'); % 蓝色虚线
% plot(x_s, T_ihxS(100, :), 'c-.', 'DisplayName', 'IHX Secondary Sodium Outlet Temperature'); % 绿色点划线

% Add labels and title
xlabel('Time (s)');
ylabel('Temperature (°C)');
ylim([200 1500]);

% Add text annotations to the plot
text(x_s(end), T_na(50, end), 'Core Coolant outlet Temperature', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'Color', 'r');
text(x_s(end), T_f(50, end), 'Fuel Temperature (center)', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'Color', 'b');
% text(x_s(end), T_ihxS(100, end), 'IHX Secondary Sodium Outlet Temperature', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'Color', 'c');

% Optionally, add a title
% title('Temperature')
hold off
