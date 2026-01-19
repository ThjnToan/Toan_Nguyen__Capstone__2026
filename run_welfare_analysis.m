% ==================================================================
% WELFARE ANALYSIS: CONSUMPTION EQUIVALENT VARIATION (CEV)
% Calculates the welfare cost of the renewable reliability crisis.
% ==================================================================
clear all; close all; clc;

% 1. Folder Cleanup
if exist('thesis_dtc', 'dir')
    rmdir('thesis_dtc', 's'); 
end
if exist('+thesis_dtc', 'dir')
    rmdir('+thesis_dtc', 's');
end

% 2. Run Dynare
dynare thesis_dtc noclearall nolog;
global oo_ M_ options_

% 3. SHOCK SETUP (Dynare 6.5+)
% We analyze the Domestic Renewable Shock (Dunkelflaute) in isolation.
M_.Sigma_e = zeros(M_.exo_nbr, M_.exo_nbr);

% A. Set Renewable Shock to 50%
i_ren = strmatch('e_ren', M_.exo_names, 'exact');
M_.Sigma_e(i_ren, i_ren) = (0.50)^2; 

% B. Ensure Price Shock is OFF
i_price = strmatch('e_price', M_.exo_names, 'exact');
M_.Sigma_e(i_price, i_price) = (0.00)^2; 

% 4. PARAMETER CONSISTENCY
% Ensure we are using the intended friction from the mod file.
% (Optional: explicit override to ensure safety)
i_adj = strmatch('phi_adj_bat', M_.param_names, 'exact');
M_.params(i_adj) = 50.0; 

% 5. Run Simulation
options_.irf = 40;
options_.order = 1;
options_.noprint = 1;
[info, oo_, options_] = stoch_simul(M_, options_, oo_, []);

% 6. Extract Values
W_ss_index = strmatch('Welfare', M_.endo_names, 'exact');
W_ss = oo_.steady_state(W_ss_index);

% Get Welfare Drop
if isfield(oo_.irfs, 'Welfare_e_ren')
    W_impact_drop = min(oo_.irfs.Welfare_e_ren); 
    W_crisis_event = W_ss + W_impact_drop;
else
    error('Welfare IRF not found.');
end

% 7. Calculate CEV (Lambda)
beta_idx = strmatch('beta', M_.param_names, 'exact');
beta = M_.params(beta_idx);
lambda = 100 * (1 - exp( (W_crisis_event - W_ss) * (1-beta) ));

% 8. Print Final Results
fprintf('\n=================================================\n');
fprintf('   WELFARE COST ANALYSIS (Vietnam Energy Crisis)   \n');
fprintf('=================================================\n');
fprintf('Steady State Welfare:               %.4f Utils\n', W_ss);
fprintf('Crisis Impact (Total Drop):         %.6f Utils\n', W_impact_drop);
fprintf('-------------------------------------------------\n');
fprintf('CONSUMPTION EQUIVALENT VARIATION (CEV)\n');
fprintf('Equivalent permanent consumption cut:\n');
fprintf('          %.6f%% \n', abs(lambda)); 
fprintf('=================================================\n');
