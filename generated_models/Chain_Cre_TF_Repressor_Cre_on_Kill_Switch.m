function [time_interval, y_out, y] = Chain_Cre_TF_Repressor_Cre_on_Kill_Switch(time_span, parameters, initial, step)
% time_span is the hours values [start, stop]
% parameters is a Map of names to numbers (e.g., rate constants, decay rates, Hill coefficients)
% initial is a Map of variable names to initial values
% step is the number of hours between samples in output; defaults to 1
% Returns vector of time, matrix of output levels at those time points, matrix of all species
    if nargin < 4, step = 1; end
    
    % Define names for input/output variable indexes
    AAV = 1;
	Cre_regulated_region = 6;
	genome = 10;

    % Set initial values
    y0=zeros(1,14);
    y0(AAV) = initial('AAV');
	y0(Cre_regulated_region) = initial('Cre_regulated_region');
	y0(genome) = initial('genome');
    
    % Run ODE
    solution = ode45(@(t,x) diff_eq(t, x, parameters), time_span, y0);
    
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
	K_R = parameters('K_R');
	alpha_p_Cas9 = parameters('alpha_p_Cas9');
	alpha_p_Cre = parameters('alpha_p_Cre');
	alpha_p_TF = parameters('alpha_p_TF');
	alpha_r_sgRNA1 = parameters('alpha_r_sgRNA1');
	alpha_r_sgRNA2 = parameters('alpha_r_sgRNA2');
	delta_Cas9 = parameters('delta_Cas9');
	delta_Cre = parameters('delta_Cre');
	delta_TF = parameters('delta_TF');
	delta_g = parameters('delta_g');
	k_cat = parameters('k_cat');
	k_cre = parameters('k_cre');
	n = parameters('n');
    
    % Unpack individual species from x
    AAV = x(1);
	Cas9 = x(2);
	Cas9_sgRNA1 = x(3);
	Cas9_sgRNA2 = x(4);
	Cre = x(5);
	Cre_regulated_region = x(6);
	TF = x(7);
	edited_Cre_regulated_region = x(8);
	edited_genome = x(9);
	genome = x(10);
	postedit_Cas9_sgRNA1 = x(11);
	postedit_Cas9_sgRNA2 = x(12);
	sgRNA1 = x(13);
	sgRNA2 = x(14);
    
    % Compute derivative for each species
    d_AAV = - k_cat*AAV*Cas9_sgRNA1;
	d_Cas9_sgRNA1 =  Cas_gRNA_binding*Cas9*sgRNA1 - Cas_degradation*Cas9_sgRNA1 - k_cat*AAV*Cas9_sgRNA1;
	d_Cas9_sgRNA2 =  Cas_gRNA_binding*Cas9*sgRNA2 - Cas_degradation*Cas9_sgRNA2 - k_cat*Cas9_sgRNA2*genome;
	d_genome = - k_cat*Cas9_sgRNA2*genome;
	d_postedit_Cas9_sgRNA1 =  k_cat*AAV*Cas9_sgRNA1 - Cas_degradation*postedit_Cas9_sgRNA1;
	d_postedit_Cas9_sgRNA2 =  k_cat*Cas9_sgRNA2*genome - Cas_degradation*postedit_Cas9_sgRNA2;
	d_edited_genome =  k_cat*Cas9_sgRNA2*genome;
	d_TF =  alpha_p_TF*(edited_Cre_regulated_region/AAV)*AAV - delta_TF*TF;
	d_Cre =  alpha_p_Cre*AAV - delta_Cre*Cre;
	d_Cre_regulated_region = - k_cre*Cre_regulated_region*Cre^4 + (Cre_regulated_region/AAV)*d_AAV;
	d_edited_Cre_regulated_region =  k_cre*Cre_regulated_region*Cre^4 + (edited_Cre_regulated_region/AAV)*d_AAV;
	d_Cas9 =  alpha_p_Cas9*AAV - Cas_gRNA_binding*Cas9*sgRNA1 - Cas_gRNA_binding*Cas9*sgRNA2 - delta_Cas9*Cas9;
	d_sgRNA1 =  alpha_r_sgRNA1*(K_R^n)/(K_R^n + TF^n)*AAV - Cas_gRNA_binding*Cas9*sgRNA1 - delta_g*sgRNA1;
	d_sgRNA2 =  alpha_r_sgRNA2*AAV - Cas_gRNA_binding*Cas9*sgRNA2 - delta_g*sgRNA2;
    
    % Pack derivatives for return
    dx = [d_AAV, d_Cas9, d_Cas9_sgRNA1, d_Cas9_sgRNA2, d_Cre, d_Cre_regulated_region, d_TF, d_edited_Cre_regulated_region, d_edited_genome, d_genome, d_postedit_Cas9_sgRNA1, d_postedit_Cas9_sgRNA2, d_sgRNA1, d_sgRNA2]';
end
