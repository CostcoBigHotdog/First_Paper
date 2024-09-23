%% Time & space steps
dt = 1e-3;              % Time step size
nt = 40000;             % Number of time points
T = dt *  (nt - 1);     % Time domain length

% TES 
L_hx = 1000;          % Hx length
nx_hx = 50;               % Number of spatial points
dx_hx = L_hx / (nx_hx - 1);     % Spatial step size

% Steam generator (SG)
L_sg = 1588;          % Hx length
nx_sg = 100;               % Number of spatial points
dx_sg = L_sg / (nx_sg - 1);     % Spatial step size

%% TES geometry
% TES
C_ms = 1.52;
A_hx = 4995e4;
C_na = 1.228;

A_tank = 1460e4;    % section of molten salt tank
rho_ms = 1.804;

% Steam generator (SG)
C_ms = 1.52;
C_water = 4.94;
C_steam = 4.51;
A_sg = 4995e4;

%% Initial conditions
% TES
T_hxP = zeros(nx_hx, nt);
T_hxP(:, 1) = linspace(326, 478, nx_hx);
T_hxTube = zeros(nx_hx, nt);
T_hxTube(:, 1) = linspace(310, 429, nx_hx);
T_hxS = zeros(nx_hx, nt);
T_hxS(:, 1) = linspace(295, 380, nx_hx);

L_hotTank = zeros(1, nt);   % height
L_hotTank(1,1) = 0.5e2;
L_coldTank = zeros(1, nt);
L_coldTank(1,1) = 6e2;
M_hotTank = zeros(1, nt);   % mass
M_hotTank(1,1) = L_hotTank(1,1) * A_tank * rho_ms;
M_coldTank = zeros(1, nt);
M_coldTank(1,1) = L_coldTank(1,1) * A_tank * rho_ms;
T_hotTank = zeros(1, nt);   % temp
T_hotTank(1,1) = 380;
T_coldTank = zeros(1, nt);
T_coldTank(1,1) = 295;

% Steam generator (SG)
T_sgP = zeros(nx_sg, nt);
% T_sgP(:, 1) = linspace(295, 380, nx_sg);
load('T_sgP1.mat');
T_sgP(:, 1) = T_sgP1;

T_sgTube = zeros(nx_sg, nt);
% T_sgTube(:, 1) = linspace(265, 340, nx_sg);
load('T_sgTube1.mat');
T_sgTube(:, 1) = T_sgTube1;

H_sgS = zeros(nx_sg, nt);
% H_sgS(:, 1) = linspace(1014, 2874, nx_sg);
load('H_sgS1.mat');
H_sgS(:, 1) = H_sgS1;
T_sgS = zeros(nx_sg, nt);
load('T_sgS1.mat');
T_sgS(:, 1) = T_sgS1;

%% PDE parameters
% Heat exchanger
% WC_hxP = -2 * 1272879.67 * 1.228 *(840/425);    % primary side sodium 
WC_hxP = - 840e6/(478 - 326);
WC_hxS = 840e6/(380 - 295);  % second sides molten salt, mass flow rate from power output
rhoCA_hxP = 0.927 * 11.46e6 / 1000 * C_na;  % comes from reference volume
rhoCA_hxS = 1.84 * 13.41e6 / 1000 * C_ms;   % comes from reference volume 
rhoCA_hxTube = 5.4245e3 * 840 / 10 ;   % comes from yimeng's code, make a scale regarding to power

h_hxP = 0.005 * A_hx * 10e-2;
h_hxS = 0.0042 * A_hx * 10e-2;

A7 =   (1 + (WC_hxP / rhoCA_hxP * dt / dx_hx)) * eye(nx_hx) ...
     - (WC_hxP / rhoCA_hxP * dt / dx_hx) * diag(ones(nx_hx-1, 1), 1);
hxP_inlet = zeros(nx_hx, 1);
hxP_inlet(nx_hx) = - (WC_hxP / rhoCA_hxP * dt / dx_hx) * 478;

A8 =   (1 - (WC_hxS / rhoCA_hxS * dt / dx_hx)) * eye(nx_hx) ...
     + (WC_hxS / rhoCA_hxS * dt / dx_hx) * diag(ones(nx_hx -1, 1), -1);
hxS_inlet = zeros(nx_hx, 1);
hxS_inlet(1) = (WC_hxS / rhoCA_hxS * dt / dx_hx) * T_coldTank(1,1);

% Steam generator (SG)
WC_sgP = - 840e6/(380 - 295);
W_sgS = 840e6/(2874 - 1014);  % second sides molten salt, mass flow rate from power output

rhoCA_sgP = 1.84 * 13.41e6 / 1000 * C_ms;  % comes from reference volume
rhoCA_sgTube = 5.4245e3 * 840 / 10 ;   % comes from yimeng's code, make a scale regarding to power
rhoA_sgS = W_sgS / 2000;  

h_sgP = 1 * A_sg  * 10e-2;
h_sgS = 0.002 * A_sg * 10e-2;

A7 =   (1 + (WC_sgP / rhoCA_sgP * dt / dx_sg)) * eye(nx_sg) ...
     - (WC_sgP / rhoCA_sgP * dt / dx_sg) * diag(ones(nx_sg-1, 1), 1);
sgP_inlet = zeros(nx_sg, 1);
sgP_inlet(nx_sg) = - (WC_sgP / rhoCA_sgP * dt / dx_sg) * 380;

A8 =   (1 - (W_sgS / rhoA_sgS * dt / dx_sg)) * eye(nx_sg) ...
     + (W_sgS / rhoA_sgS * dt / dx_sg) * diag(ones(nx_sg -1, 1), -1);
sgS_inlet = zeros(nx_sg, 1);
sgS_inlet(1) = (W_sgS / rhoA_sgS * dt / dx_sg) * 1014;

%% Time iteration
% Start timing
tic;
% Time stepping using the backward Euler method
for k = 2:nt
    
    % TES hx temperature update
    % ihxP_inlet(nx_ihx) = - (WC_ihxP / rhoCA_ihxP * dt / dx_ihx) * T_hotPool(nx_hotPool, k-1);
    T_hxP(:, k) = A7 * T_hxP(:, k-1) + hxP_inlet + h_hxP * dt / rhoCA_hxP .* (T_hxTube(:, k-1) - T_hxP(:, k-1));
    T_hxTube(:, k) = T_hxTube(:, k -1) + h_hxP * dt / rhoCA_hxTube .* (T_hxP(:, k-1) - T_hxTube(:, k-1))...
                                   + h_hxS * dt / rhoCA_hxTube .* (T_hxS(:, k-1) - T_hxTube(:, k-1));
    % ihxS_inlet(1) =  (WC_ihxS / rhoCA_ihxS * dt / dx_ihx) * 326;
    T_hxS(:, k) = A8 * T_hxS(:, k-1) + hxS_inlet + h_hxS * dt / rhoCA_hxS .* (T_hxTube(:, k-1) - T_hxS(:, k-1));

    % SG temperature update
    % ihxP_inlet(nx_ihx) = - (WC_ihxP / rhoCA_ihxP * dt / dx_ihx) * T_hotPool(nx_hotPool, k-1);
    T_sgP(:, k) = A7 * T_sgP(:, k-1) + sgP_inlet + h_sgP * dt / rhoCA_sgP .* (T_sgTube(:, k-1) - T_sgP(:, k-1));
    T_sgTube(:, k) = T_sgTube(:, k -1) + h_sgP * dt / rhoCA_sgTube .* (T_sgP(:, k-1) - T_sgTube(:, k-1))...
                                       + h_sgS * dt / rhoCA_sgTube .* (T_sgS(:, k-1) - T_sgTube(:, k-1));
    H_sgS(:, k) = A8 * H_sgS(:, k-1) + sgS_inlet + h_sgS * dt / rhoA_sgS .* (T_sgTube(:, k-1) - T_sgS(:, k-1));
    % T_sgS(:, k) = region1(H_sgS(:, k), nx_sg);
    for i = 1:nx_sg
        if H_sgS(i, k) < 1228.95
            T_sgS(i, k) = 235 + (H_sgS(i, k) - 1014) / C_water;
        elseif H_sgS(i, k) >= 1228.95 && H_sgS(i, k) <= 2781.48
            T_sgS(i, k) = 278.52;
        else
            T_sgS(i, k) = 278.52 + (H_sgS(i, k) - 2781.48) / C_steam;
        end
    end


  
end
% Stop timing
elapsedTime = toc;
fprintf('Elapsed time: %.6f seconds\n', elapsedTime);


