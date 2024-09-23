%% Time & space steps
% Steam generator (SG)
L_sg = 1588;          % Hx length
nx_sg = 100;               % Number of spatial points
dx_sg = L_sg / (nx_sg - 1);     % Spatial step size
dt = 1e-3;              % Time step size
nt = 500000;             % Number of time points
T = dt *  (nt - 1);     % Time domain length

%% Geometry
% Steam generator (SG)
C_ms = 1.52;
C_water = 4.94;
C_steam = 4.51;
A_sg = 4995e4;

%% Initial conditions
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
% for i = 1:nx_sg
%     if H_sgS(i, 1) < 1228.95
%         T_sgS(i, 1) = 235 + (H_sgS(i, 1) - 1014) / C_water;
%     elseif H_sgS(i, 1) >= 1228.95 && H_sgS(i, 1) <= 2781.48
%         T_sgS(i, 1) = 278.52;
%     else
%         T_sgS(i, 1) = 278.52 + (H_sgS(i, 1) - 2781.48) / C_steam;
%     end
% end
load('T_sgS1.mat');
T_sgS(:, 1) = T_sgS1;

%% PDE parameters
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


