% ==================================================================
% SMALL OPEN ECONOMY: IMPORTED GREENFLATION SHOCK
% Section 6.6: Impact of Global Battery Price Shock
% ==================================================================
clear all; close all; clc;

if exist('thesis_dtc', 'dir')
    rmdir('thesis_dtc', 's'); % 's' removes all subdirectories/files
end
if exist('+thesis_dtc', 'dir')
    rmdir('+thesis_dtc', 's');
end

%run dynare
dynare thesis_dtc noclearall nolog;

global M_ oo_ options_

% --- CRITICAL SAFETY PATCH ---
% Enforce the same installation friction as Figure 1
% This ensures the "SOE Risk Shift" comparison is fair.
i_adj = strmatch('phi_adj_bat', M_.param_names, 'exact');
M_.params(i_adj) = 50.0; 

% 2. Isolate the Global Price Shock
options_.irf = 40;
options_.order = 1;
options_.noprint = 1;

% Manually zero out the e_ren shock (Domestic)
M_.xo_variance = zeros(M_.exo_nbr, M_.exo_nbr);

% Set e_price shock (Global) to 10% 
% (A 10% spike in Lithium prices is a standard stress test)
i_price = strmatch('e_price', M_.exo_names, 'exact');
M_.xo_variance(i_price, i_price) = (0.10)^2; 

% 3. Run Simulation
[info, oo_, options_] = stoch_simul(M_, options_, oo_, []);

% 4. PLOTTING FIGURE
figure('Name', 'Figure 6: Imported Greenflation');

% Panel 1: The Global Shock
subplot(1,2,1);
plot(oo_.irfs.P_tech_e_price * 100, 'r-', 'LineWidth', 2);
title('Global Battery Price (P_{tech})');
ylabel('% Increase');
xlabel('Quarters');
grid on;

% Panel 2: The Domestic Fallout
subplot(1,2,2);
plot(oo_.irfs.Y_e_price * 100, 'k-', 'LineWidth', 2);
title('Impact on Domestic GDP');
ylabel('% Deviation');
xlabel('Quarters');
grid on;
