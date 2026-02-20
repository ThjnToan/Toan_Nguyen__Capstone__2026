% ==================================================================
% DRIVER SCRIPT: SENSITIVITY ANALYSIS (ELASTICITY OF SUBSTITUTION)
% Tests if the "Agility Gap" holds when Grid/Batteries are Substitutes.
% ==================================================================
clear all; close all; clc;

% 1. Folder Cleanup
if exist('thesis_dtc', 'dir')
    rmdir('thesis_dtc', 's'); 
end
if exist('+thesis_dtc', 'dir')
    rmdir('+thesis_dtc', 's');
end

% 2. Run Dynare (Load the Model Structure)
dynare thesis_dtc noclearall nolog;
global M_ oo_ options_

% 3. Define Scenarios
rho_values = [0.4, 0.9]; % 0.4 = Complements (Vietnam), 0.9 = High Substitution
scenario_names = {'Baseline (Complements \rho=0.4)', 'High Substitution (\rho=0.9)'};
line_styles = {'b-', 'r--'}; 

% Find parameter index
i_rho = strmatch('rho', M_.param_names, 'exact');

% --- PRE-ALLOCATE MATRICES ---
results_Ig = nan(40, length(rho_values));
results_Ib = nan(40, length(rho_values));

% 4. Run Loop
fprintf('\nRunning Sensitivity Analysis...\n');
options_.irf = 40;
options_.order = 1;
options_.noprint = 1;

for i = 1:length(rho_values)
    % Update rho
    M_.params(i_rho) = rho_values(i);
    
    % Run Simulation
    [info, oo_, options_] = stoch_simul(M_, options_, oo_, []);
    
    if info(1) == 0
        % The shock in your updated mod file is 'e_ren_shock'
        if isfield(oo_.irfs, 'I_g_e_ren_shock')
             results_Ig(:, i) = oo_.irfs.I_g_e_ren_shock';
             results_Ib(:, i) = oo_.irfs.I_b_e_ren_shock';
        elseif isfield(oo_.irfs, 'I_g_e_ren')
             % Fallback for older mod versions
             results_Ig(:, i) = oo_.irfs.I_g_e_ren';
             results_Ib(:, i) = oo_.irfs.I_b_e_ren';
        else
             fprintf('Warning: IRF fields not found. Check shock name in .mod file.\n');
        end
    else
        fprintf('Error: Simulation failed for rho = %.2f (Code: %d)\n', rho_values(i), info(1));
    end
end

% 5. Plot Comparison
figure('Name', 'Sensitivity Analysis: Agility Gap');
t = 1:40;

% Panel A: Private Battery Investment (I_b)
subplot(1,2,1);
plot(t, results_Ib(:,1), line_styles{1}, 'LineWidth', 2); hold on;
plot(t, results_Ib(:,2), line_styles{2}, 'LineWidth', 2);
title('Private Response (I_b)');
ylabel('Deviation from SS');
xlabel('Quarters');
grid on;
legend(scenario_names, 'Location', 'SouthEast'); % Moved Legend to avoid covering lines

% Panel B: Public Grid Investment (I_g)
subplot(1,2,2);
plot(t, results_Ig(:,1), line_styles{1}, 'LineWidth', 2); hold on;
plot(t, results_Ig(:,2), line_styles{2}, 'LineWidth', 2);
title('Public Response (I_g)');
xlabel('Quarters');
grid on;

% Add Note
sgtitle('Figure A.2: Sensitivity of Agility Gap to Substitution Elasticity');
