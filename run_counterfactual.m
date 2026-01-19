% ==================================================================
% COUNTERFACTUAL ANALYSIS: THE VALUE OF INNOVATION
% Compares "Agile Innovation" vs. "Stagnant Technology"
% ==================================================================
clear all; close all; clc;

% 1. CLEANUP (Prevent File Lock Errors)
% Forcefully remove previous Dynare output folders to avoid overwrite issues
if exist('thesis_dtc', 'dir')
    rmdir('thesis_dtc', 's'); 
end
if exist('+thesis_dtc', 'dir')
    rmdir('+thesis_dtc', 's');
end

% 2. INITIALIZE DYNARE
% Run the model once to create the global structures
dynare thesis_dtc noclearall nolog;

% 3. SETUP EXPERIMENT
global M_ oo_ options_

% Find parameter indices
i_eta = strmatch('eta_dtc', M_.param_names, 'exact');
i_adj = strmatch('phi_adj_bat', M_.param_names, 'exact');

% --- CRITICAL: Enforce Consistency ---
% We must use the same friction as Figure 1 (phi_adj_bat = 50.0).
% This ensures the "Innovation Rescue" comparison is fair and matches the Stress Test.
M_.params(i_adj) = 50.0; 

% ------------------------------------------------------------------
% RUN 1: NO INNOVATION (Counterfactual)
% ------------------------------------------------------------------
fprintf('Simulating Counterfactual (No Innovation)...\n');

% Manually set eta_dtc to 0 (Technology is fixed/exogenous)
M_.params(i_eta) = 0; 

% Force a re-simulation (Quietly)
options_.noprint = 1; 
options_.irf = 40;
options_.order = 1;

% Use [] as the last argument to ensure all variables are calculated
[info, oo_, options_] = stoch_simul(M_, options_, oo_, []); 

% Save the results
if isfield(oo_.irfs, 'Y_e_ren')
    Y_no_innov = oo_.irfs.Y_e_ren;
else
    error('IRF for Y_e_ren not found in Run 1.');
end

% ------------------------------------------------------------------
% RUN 2: WITH INNOVATION (Thesis Baseline)
% ------------------------------------------------------------------
fprintf('Simulating Baseline (With Innovation)...\n');

% Manually set eta_dtc back to 0.5 (Innovation is active/endogenous)
M_.params(i_eta) = 0.5;

% Force a re-simulation
[info, oo_, options_] = stoch_simul(M_, options_, oo_, []);

% Save the results
if isfield(oo_.irfs, 'Y_e_ren')
    Y_with_innov = oo_.irfs.Y_e_ren;
else
    error('IRF for Y_e_ren not found in Run 2.');
end

% ==================================================================
% PLOTTING FIGURE 4
% ==================================================================
figure('Name', 'Figure 4: The Innovation Rescue Effect');
hold on;
grid on;

% 1. Define X-axis
T_steps = length(Y_no_innov);
x_axis = 1:T_steps;

% 2. Ensure vectors are Row Vectors (1xN) for 'fill' command
Y_no_innov   = reshape(Y_no_innov, 1, []); 
Y_with_innov = reshape(Y_with_innov, 1, []);

% 3. Plot the Shaded Area (The "Welfare Gain")
% Creates a grey polygon between the two curves
fill([x_axis, fliplr(x_axis)], ...
     [Y_with_innov, fliplr(Y_no_innov)], ...
     [0.9 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.5); 

% 4. Plot the Lines
p1 = plot(x_axis, Y_with_innov, 'b-', 'LineWidth', 3);       % Baseline
p2 = plot(x_axis, Y_no_innov,   'r--', 'LineWidth', 3);      % Counterfactual

% 5. Labels & Formatting
xlabel('Quarters');
ylabel('% Deviation from SS');
title('The Innovation Rescue Effect');
legend([p1 p2], 'Thesis Model (Endogenous Innovation)', 'Counterfactual (Fixed Tech)', ...
       'Location', 'SouthEast');
axis tight;
hold off;
