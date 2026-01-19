% ==================================================================
% DRIVER SCRIPT: SENSITIVITY ANALYSIS (ELASTICITY OF SUBSTITUTION)
% Tests if the "Agility Gap" holds when Grid/Batteries are Substitutes.
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

% 2. Define Scenarios
rho_values = [0.4, 0.55]; % 0.4 = Complements (Vietnam), 0.9 = Substitutes
scenario_names = {'Baseline (Complements \rho=0.4)', 'Weak Substitution (\rho=0.5)'};
line_styles = {'b-', 'r--'}; 

% Find parameter index
i_rho = strmatch('rho', M_.param_names, 'exact');

% --- FIX: PRE-ALLOCATE MATRICES WITH NaNs ---
% We expect 40 periods (rows) and 2 scenarios (columns).
% Using NaNs ensures that if a run fails, the plot doesn't crash.
results_Ig = nan(40, length(rho_values));
results_Ib = nan(40, length(rho_values));

% 3. Run Loop
fprintf('Running Sensitivity Analysis...\n');
options_.irf = 40;
options_.order = 1;
options_.noprint = 1;

for i = 1:length(rho_values)
    % Update rho
    M_.params(i_rho) = rho_values(i);
    
    % Run Simulation
    [info, oo_, options_] = stoch_simul(M_, options_, oo_, []);
    
    if info(1) == 0
        % --- FIX: DIRECT INDEXING ---
        % Store results in the specific column 'i'
        if isfield(oo_.irfs, 'I_g_e_ren')
             results_Ig(:, i) = oo_.irfs.I_g_e_ren';
             results_Ib(:, i) = oo_.irfs.I_b_e_ren';
        else
             fprintf('Warning: IRFs not found for rho=%.2f\n', rho_values(i));
        end
    else
        fprintf('Error: Simulation failed for rho = %.2f (Code: %d)\n', rho_values(i), info(1));
    end
end

% 4. Plot Comparison
figure('Name', 'Sensitivity Analysis: Agility Gap');
t = 1:40;

% Panel A: Private Battery Investment (I_b)
subplot(1,2,1);
plot(t, results_Ib(:,1), line_styles{1}, 'LineWidth', 2); hold on;
% This line will now be invisible (but safe) if the second run failed
plot(t, results_Ib(:,2), line_styles{2}, 'LineWidth', 2);
title('Private Response (I_b)');
ylabel('Deviation from SS');
xlabel('Quarters');
grid on;
legend(scenario_names, 'Location', 'NorthEast');

% Panel B: Public Grid Investment (I_g)
subplot(1,2,2);
plot(t, results_Ig(:,1), line_styles{1}, 'LineWidth', 2); hold on;
plot(t, results_Ig(:,2), line_styles{2}, 'LineWidth', 2);
title('Public Response (I_g)');
xlabel('Quarters');
grid on;

% Add Note
sgtitle('Figure A.2: Sensitivity of Agility Gap to Substitution Elasticity');
