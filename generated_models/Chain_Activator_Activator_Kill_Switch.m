function [time_interval, y_out, y] = Chain_Activator_Activator_Kill_Switch(time_span, parameters, initial, step)
% time_span is the hours values [start, stop]
% parameters is a Map of names to numbers (e.g., rate constants, decay rates, Hill coefficients)
% initial is a Map of variable names to initial values
% step is the number of hours between samples in output; defaults to 1
% Returns vector of time, matrix of output levels at those time points, matrix of all species
    if nargin < 4, step = 1; end
    
    % Define names for input/output variable indexes
    AAV = 1;
	genome = 8;

    % Set initial values
    y0=zeros(1,12);
    y0(AAV) = initial('AAV');
	y0(genome) = initial('genome');
    
    % Run ODE
    solution = ode15s(@(t,x) diff_eq(t, x, parameters), time_span, y0);
    
    % Evaluate species levels at given times
    time_interval = time_span(1):step:time_span(end);
    y = deval(solution, time_interval);
    y_out = y([AAV],:);
end

% ODE differential function
function dx=diff_eq(t, x, parameters)
    % Unpack parameters from parameter map
    Cas_degradation = parameters('Cas_degradation');
	Cas_gRNA_binding = parameters('Cas_gRNA_binding');
	K_A = parameters('K_A');
	alpha_p_Cas9 = parameters('alpha_p_Cas9');
	alpha_p_TF = parameters('alpha_p_TF');
	alpha_p_TF2 = parameters('alpha_p_TF2');
	alpha_r_sgRNA1 = parameters('alpha_r_sgRNA1');
	alpha_r_sgRNA2 = parameters('alpha_r_sgRNA2');
	delta_Cas9 = parameters('delta_Cas9');
	delta_TF = parameters('delta_TF');
	delta_TF2 = parameters('delta_TF2');
	delta_g = parameters('delta_g');
	k_cat = parameters('k_cat');
	n = parameters('n');
    
    % Unpack individual species from x
    x = max(1e-12,real(x)); % Truncate values just above zero
    AAV = x(1);
	Cas9 = x(2);
	Cas9_sgRNA1 = x(3);
	Cas9_sgRNA2 = x(4);
	TF = x(5);
	TF2 = x(6);
	edited_genome = x(7);
	genome = x(8);
	postedit_Cas9_sgRNA1 = x(9);
	postedit_Cas9_sgRNA2 = x(10);
	sgRNA1 = x(11);
	sgRNA2 = x(12);
    
    % Compute derivative for each species
    d_AAV = - k_cat*AAV*Cas9_sgRNA1;
	d_Cas9_sgRNA1 =  Cas_gRNA_binding*Cas9*sgRNA1 - Cas_degradation*Cas9_sgRNA1 - k_cat*AAV*Cas9_sgRNA1;
	d_Cas9_sgRNA2 =  Cas_gRNA_binding*Cas9*sgRNA2 - Cas_degradation*Cas9_sgRNA2 - k_cat*Cas9_sgRNA2*genome;
	d_genome = - k_cat*Cas9_sgRNA2*genome;
	d_postedit_Cas9_sgRNA1 =  k_cat*AAV*Cas9_sgRNA1 - Cas_degradation*postedit_Cas9_sgRNA1;
	d_postedit_Cas9_sgRNA2 =  k_cat*Cas9_sgRNA2*genome - Cas_degradation*postedit_Cas9_sgRNA2;
	d_edited_genome =  k_cat*Cas9_sgRNA2*genome;
	d_TF =  alpha_p_TF*AAV - delta_TF*TF;
	d_TF2 =  alpha_p_TF2*(TF^n)/(K_A^n + TF^n)*AAV - delta_TF2*TF2;
	d_Cas9 =  alpha_p_Cas9*AAV - Cas_gRNA_binding*Cas9*sgRNA1 - Cas_gRNA_binding*Cas9*sgRNA2 - delta_Cas9*Cas9;
	d_sgRNA1 =  alpha_r_sgRNA1*(TF2^n)/(K_A^n + TF2^n)*AAV - Cas_gRNA_binding*Cas9*sgRNA1 - delta_g*sgRNA1;
	d_sgRNA2 =  alpha_r_sgRNA2*AAV - Cas_gRNA_binding*Cas9*sgRNA2 - delta_g*sgRNA2;
    
    % Pack derivatives for return, ensuring none are complex or go below zero
    dx = max(-x,real([d_AAV, d_Cas9, d_Cas9_sgRNA1, d_Cas9_sgRNA2, d_TF, d_TF2, d_edited_genome, d_genome, d_postedit_Cas9_sgRNA1, d_postedit_Cas9_sgRNA2, d_sgRNA1, d_sgRNA2])');
end
