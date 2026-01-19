% ==================================================================
% FINDING THE OPTIMAL FLEXIBILITY MIX
% Sensitivity of Welfare to Battery Share (Mu)
% ==================================================================
clear all; close all; clc;

% 1. CLEANUP (Prevent File Lock Errors)
if exist('thesis_dtc', 'dir')
    rmdir('thesis_dtc', 's'); 
end
if exist('+thesis_dtc', 'dir')
    rmdir('+thesis_dtc', 's');
end

% 2. INITIALIZE DYNARE
% Run the model once to create the global structures
dynare thesis_dtc noclearall nolog;

global M_ oo_ options_

% 3. SETUP EXPERIMENT
mu_values = 0.05:0.05:0.60; 
output_drop = zeros(length(mu_values), 1);

% Find indices
i_mu   = strmatch('mu', M_.param_names, 'exact');
i_adj  = strmatch('phi_adj_bat', M_.param_names, 'exact');
i_eren = strmatch('e_ren', M_.exo_names, 'exact');

% --- CRITICAL CONSISTENCY CHECK ---
% We enforce the same "Stress Test" settings as Figure 1 to ensure
% Figure 5 is logically consistent with the rest of the thesis.

% 1. Force High Adjustment Costs (The "Friction")
M_.params(i_adj) = 50.0; 

% 2. Force the 50% Shock (The "Dunkelflaute")
M_.xo_variance(i_eren, i_eren) = (0.50)^2; 

fprintf('Running Optimization Loop with Friction = %.1f...\n', M_.params(i_adj));

% 4. RUN THE LOOP
options_.irf = 40; 
options_.order = 1;
options_.noprint = 1; % Keep command window clean

for i = 1:length(mu_values)
    % Update Mu
    M_.params(i_mu) = mu_values(i);
    
    % Run Simulation
    % Use [] to calculate all variables
    [info, oo_, options_] = stoch_simul(M_, options_, oo_, []);
    
    % Capture the "Max Drop" (Recession Depth)
    if isfield(oo_.irfs, 'Y_e_ren')
        % We multiply by 100 to get percentage
        output_drop(i) = min(oo_.irfs.Y_e_ren) * 100; 
    else
        error('IRF for Y_e_ren not found. Check shock name.');
    end
    
    fprintf('  Mu = %.2f | Max Output Drop = %.6f%%\n', mu_values(i), output_drop(i));
end

% 5. PLOT FIGURE 5
figure('Name', 'Figure 5: The Transition Trap');
hold on; 

% Plot the Frontier Curve
plot(mu_values, output_drop, 'b-o', 'LineWidth', 2, 'MarkerFaceColor', 'b');

% Highlight Vietnam (0.16)
xline(0.16, 'r--', 'LineWidth', 2);
text(0.17, min(output_drop)*0.95, 'Vietnam (Current)', 'Color', 'r', 'FontSize', 10);

% Highlight The Optimal Point (Closest to zero drop)
[min_loss, min_idx] = max(output_drop); 
optimal_mu = mu_values(min_idx);
xline(optimal_mu, 'g--', 'LineWidth', 2);
text(optimal_mu + 0.02, min_loss, 'Optimal Mix', 'Color', 'g', 'FontSize', 10, 'FontWeight', 'bold');

% Formatting
xlabel('Battery Share of Flexibility (\mu)', 'FontSize', 12);
ylabel('Output Response to Shock (%)', 'FontSize', 12);
title('The Transition Trap: Optimal Flexibility Mix');
grid on; box on;

% Add Annotation Box
dim = [.15 .2 .3 .3]; 
str = {'The "Reliability Trap"', 'Low \mu makes the grid', 'too brittle.'};
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on', ...
           'BackgroundColor', 'white', 'EdgeColor', 'k');

hold off;
