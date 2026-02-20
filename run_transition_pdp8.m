% ==================================================================
% SCENARIO ANALYSIS: THE PDP8 TRANSITION (2025-2050)
% Open Economy Simulation: Renewable Ramp + Trade Balance Effects
% ==================================================================

% Add this at the top of your script
if exist('thesis_pdp8', 'dir')
    rmdir('thesis_pdp8', 's'); % 's' removes all subdirectories/files
end
if exist('+thesis_pdp8', 'dir')
    rmdir('+thesis_pdp8', 's');
end

clear all; close all; clc;

% 1. Run the Deterministic Model
% We use 'noclearall' to preserve data for plotting
dynare thesis_pdp8 noclearall nolog;

% 2. Setup Plotting Data
global M_ oo_ options_
T = 1:100; % 100 Quarters (25 Years)

figure('Name', 'Chapter 6: PDP8 Transition Path (Open Economy)');

% --- Panel 1: The Policy Target (Renewable Capacity) ---
subplot(3,2,1);
% Check if E_ren is exogenous (simulated in oo_.exo_simul)
id_E = strmatch('E_ren', M_.exo_names, 'exact');
% exo_simul usually has extra periods for initial/terminal conditions
path_E = oo_.exo_simul(2:101, id_E); 
plot(T, path_E, 'b-', 'LineWidth', 2);
title('Renewable Capacity Target (PDP8)');
ylabel('Capacity'); grid on;

% --- Panel 2: The Reliability Valley ---
subplot(3,2,2);
plot(T, oo_.endo_simul(strmatch('u', M_.endo_names, 'exact'), 1:100), 'r-', 'LineWidth', 2);
title('Grid Reliability (u)');
ylabel('Utilization'); grid on;
yline(0.98, 'k--', 'Reliability Threshold');

% --- Panel 3: Investment Agility ---
subplot(3,2,3);
id_Ib = strmatch('I_b', M_.endo_names, 'exact');
id_Ig = strmatch('I_g', M_.endo_names, 'exact');
plot(T, oo_.endo_simul(id_Ib, 1:100), 'b-', 'LineWidth', 2); hold on;
plot(T, oo_.endo_simul(id_Ig, 1:100), 'k--', 'LineWidth', 1.5);
legend('Private Batteries (Imported)', 'Public Grid (Domestic)');
title('Investment Response');
grid on;

% --- Panel 4: Innovation Response ---
subplot(3,2,4);
plot(T, oo_.endo_simul(strmatch('A_bat', M_.endo_names, 'exact'), 1:100), 'g-', 'LineWidth', 2);
title('Battery Innovation (A_{bat})');
grid on;

% --- Panel 5: THE NEW OPEN ECONOMY RESULT ---
subplot(3,2,5:6); % Spans both columns at the bottom
% We plot the Trade Balance to show the cost of importing tech
id_TB = strmatch('TradeBal', M_.endo_names, 'exact');
plot(T, oo_.endo_simul(id_TB, 1:100), 'm-', 'LineWidth', 2);
title('External Balance (Trade Deficit)');
ylabel('Net Exports'); xlabel('Quarters (2025-2050)');
grid on;
yline(0, 'k-', 'Balanced Trade');
