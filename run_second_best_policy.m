% ==================================================================
% DISTINCTION UPGRADE: SECOND-BEST POLICY ANALYSIS
% "Does a broken price signal justify government subsidies?"
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

subsidy_range = 0:2:30; % Test subsidies from 0% to 30%
market_scenarios = [1.0, 0.3]; % 1.0 = Perfect Market, 0.3 = Regulated (Sticky)
scenario_labels = {'Perfect Market (\chi=1.0)', 'Regulated Market (\chi=0.3)'};
colors = {'b-o', 'r-s'};

welfare_results = zeros(length(subsidy_range), length(market_scenarios));

% Parameters indices
i_sub = strmatch('phi_sub', M_.param_names, 'exact');
i_chi = strmatch('chi_price', M_.param_names, 'exact');
i_welfare = strmatch('Welfare', M_.endo_names, 'exact');

fprintf('Running Policy Robustness Check...\n');
options_.irf = 40;
options_.order = 1;
options_.noprint = 1;

% 2. The Double Loop
for s = 1:length(market_scenarios)
    M_.params(i_chi) = market_scenarios(s); % Set Market Regime
    
    for i = 1:length(subsidy_range)
        M_.params(i_sub) = subsidy_range(i); % Set Subsidy
        
        % Run
        [info, oo_, options_] = stoch_simul(M_, options_, oo_, []);
        
        if info(1) == 0
            % Capture Welfare Loss (IRF Drop)
            % Note: We look for the "Least Bad" outcome (Max Welfare)
            W_drop = min(oo_.irfs.Welfare_e_ren); 
            welfare_results(i, s) = W_drop;
        else
            welfare_results(i, s) = NaN;
        end
    end
end

% 3. The Distinction Graph
figure('Name', 'Second-Best Policy Analysis');
hold on;

% Plot Scenario 1: Perfect Market (Blue)
plot(subsidy_range, welfare_results(:,1), colors{1}, 'LineWidth', 2, 'MarkerFaceColor', 'b');

% Plot Scenario 2: Regulated Market (Red)
plot(subsidy_range, welfare_results(:,2), colors{2}, 'LineWidth', 2, 'MarkerFaceColor', 'r');

% Find Peaks (Optimal Policies)
[max_w1, idx1] = max(welfare_results(:,1));
[max_w2, idx2] = max(welfare_results(:,2));

% Annotate Optimal Points
plot(subsidy_range(idx1), max_w1, 'k*', 'MarkerSize', 15, 'LineWidth', 2);
plot(subsidy_range(idx2), max_w2, 'k*', 'MarkerSize', 15, 'LineWidth', 2);

text(subsidy_range(idx1), max_w1+0.00005, [' Optimal: \phi=' num2str(subsidy_range(idx1))], 'HorizontalAlignment', 'center');
text(subsidy_range(idx2), max_w2+0.00005, [' Optimal: \phi=' num2str(subsidy_range(idx2))], 'HorizontalAlignment', 'center');

xlabel('R&D Subsidy Aggressiveness (\phi_{sub})');
ylabel('Welfare Impact (Utils)');
legend(scenario_labels{1}, scenario_labels{2}, 'Optimal Point', 'Location', 'SouthWest');
grid on;
hold off;
