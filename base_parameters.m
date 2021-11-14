parameters = containers.Map();

% This file contains all of the base parameter values as determined from parameter fitting and literature
% Fit parameters are all fit to Chylinski et al data in parameter_fitting.m
% This should be the starting point for any simulation run

% Therapy parameters - cannot be modified
parameters('Cas_degradation') = 10^(-2.0);      % fraction/hour, fit, but lowered for stable cells (vs. fast-dividing mES in Chylinski); (assuming bound, unbound, and post-edit are equal)
parameters('Cas_gRNA_binding') = 10^(-4.2577);  % fit; assuming identical for all gRNAs
parameters('alpha_r_sgRNA2') = 10^(3.3090);     % molecules/hour, fit
parameters('alpha_p_Cas9') = 10^(3.0415);       % molecules/hour, fit
parameters('delta_Cas9') = parameters('Cas_degradation');
parameters('delta_g') = 10^(0.0003);            % fraction/hour, fit; assuming identical for all gRNAs
parameters('k_cat') = 10^(-4.7518);             % fit

% Base kill switch parameters
parameters('alpha_r_sgRNA1') = parameters('alpha_r_sgRNA2'); % assume not readily modified 

% Delay mechanism parameters
parameters('alpha_p_TF') = parameters('alpha_p_Cas9'); % Can be modulated - explore values
parameters('delta_TF') = parameters('Cas_degradation'); % Assume reasonably similar and stable proteins
parameters('K_A') = 2.34 * 10^6;                % from Calin Belta paper; Cannot be readily modulated
parameters('K_R') = parameters('K_A');          % Assume same transition; Cannot be readily modulated
parameters('n') = 0.92;                         % from Calin Belta paper; Cannot be readily modulated

parameters('alpha_p_Cre') = parameters('alpha_p_Cas9'); % Can be modulated - explore values
parameters('delta_Cre') = parameters('Cas_degradation'); % Assume reasonably similar and stable proteins
parameters('k_cre') = 10^(-7.1535);             % fit; Cannot be readily modulated

% Second copies of delay mechanisms
parameters('alpha_p_TF2') = parameters('alpha_p_TF'); % Can be modulated - explore values
parameters('delta_TF2') = parameters('delta_TF'); % Assume reasonably similar and stable proteins
parameters('alpha_p_CreH') = parameters('alpha_p_Cre'); % Can be modulated - explore values
parameters('delta_CreH') = parameters('delta_Cre'); % Assume reasonably similar and stable proteins
