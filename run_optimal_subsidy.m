% ==================================================================
% FINDING THE OPTIMAL R&D SUBSIDY
% This script runs 'thesis_dtc.mod' 21 times to find the peak.
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

% 2. Setup The Experiment
global M_ oo_ options_

% We test subsidy aggressiveness from 0 (None) to 20 (High)
sub_values = 0:1:20; 
welfare_results = zeros(length(sub_values), 1);
innovation_gain = zeros(length(sub_values), 1);

% Find where 'phi_sub' is located in the parameter list
i_phi = strmatch('phi_sub', M_.param_names, 'exact');

% --- CONSISTENCY PATCH ---
% Enforce Stress Test Calibration
i_adj = strmatch('phi_adj_bat', M_.param_names, 'exact');
M_.params(i_adj) = 50.0; 

i_eren = strmatch('e_ren', M_.exo_names, 'exact');
M_.xo_variance(i_eren, i_eren) = (0.50)^2;

fprintf('Searching for Optimal Subsidy...\n');
options_.irf = 40;
options_.order = 1;
options_.noprint = 1; % Don't print 20 separate reports

for i = 1:length(sub_values)
    % A. Update the Parameter
    M_.params(i_phi) = sub_values(i);
    
    % B. Run the Model (Quietly)
    [info, oo_, options_] = stoch_simul(M_, options_, oo_, []);
    
    % C. Store Results
    if info(1) == 0 % If simulation succeeded
        % Metric 1: Total Consumption (Welfare Proxy)
        welfare_results(i) = sum(oo_.irfs.C_e_ren); 
        
        % Metric 2: Total Technology Gain
        innovation_gain(i) = sum(oo_.irfs.A_bat_e_ren);
    else
        welfare_results(i) = NaN;
    end
    
    fprintf('  Phi = %d | Welfare Change = %.4f\n', sub_values(i), welfare_results(i));
end

% 3. Plot Figure 4.7
figure('Name', 'Optimal R&D Subsidy');

% Left Axis: Innovation Speed (Benefit)
yyaxis left
plot(sub_values, innovation_gain, 'b-o', 'LineWidth', 2);
ylabel('Innovation Gain (A_{bat})', 'Color', 'b');
xlabel('Subsidy Aggressiveness (\phi_{sub})');

% Right Axis: Consumption Welfare (Net Effect)
yyaxis right
plot(sub_values, welfare_results, 'r-s', 'LineWidth', 2);
ylabel('Net Consumption Impact (Welfare)', 'Color', 'r');

title('Figure 4.7: The Optimal R&D Subsidy');
grid on;
legend('Tech Gain', 'Welfare Cost', 'Location', 'SouthWest');
