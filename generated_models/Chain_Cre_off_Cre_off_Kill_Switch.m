function [time_interval, y_out, y] = Chain_Cre_off_Cre_off_Kill_Switch(time_span, parameters, initial, step)
% time_span is the hours values [start, stop]
% parameters is a Map of names to numbers (e.g., rate constants, decay rates, Hill coefficients)
% initial is a Map of variable names to initial values
% step is the number of hours between samples in output; defaults to 1
% Returns vector of time, matrix of output levels at those time points, matrix of all species
    if nargin < 4, step = 1; end
    
    % Define names for input/output variable indexes
    AAV = 1;
	CreH_regulated_region = 7;
	Cre_regulated_region = 8;
	genome = 12;

    % Set initial values
    y0=zeros(1,16);
    y0(AAV) = initial('AAV');
	y0(CreH_regulated_region) = initial('CreH_regulated_region');
	y0(Cre_regulated_region) = initial('Cre_regulated_region');
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
	alpha_p_Cas9 = parameters('alpha_p_Cas9');
	alpha_p_Cre = parameters('alpha_p_Cre');
	alpha_p_CreH = parameters('alpha_p_CreH');
	alpha_r_sgRNA1 = parameters('alpha_r_sgRNA1');
	alpha_r_sgRNA2 = parameters('alpha_r_sgRNA2');
	delta_Cas9 = parameters('delta_Cas9');
	delta_Cre = parameters('delta_Cre');
	delta_CreH = parameters('delta_CreH');
	delta_g = parameters('delta_g');
	k_cat = parameters('k_cat');
	k_cre = parameters('k_cre');
    
    % Unpack individual species from x
    x = max(1e-12,real(x)); % Truncate values just above zero
    AAV = x(1);
	Cas9 = x(2);
	Cas9_sgRNA1 = x(3);
	Cas9_sgRNA2 = x(4);
	Cre = x(5);
	CreH = x(6);
	CreH_regulated_region = x(7);
	Cre_regulated_region = x(8);
	edited_CreH_regulated_region = x(9);
	edited_Cre_regulated_region = x(10);
	edited_genome = x(11);
	genome = x(12);
	postedit_Cas9_sgRNA1 = x(13);
	postedit_Cas9_sgRNA2 = x(14);
	sgRNA1 = x(15);
	sgRNA2 = x(16);
    
    % Compute derivative for each species
    d_AAV = - k_cat*AAV*Cas9_sgRNA1;
	d_Cas9_sgRNA1 =  Cas_gRNA_binding*Cas9*sgRNA1 - Cas_degradation*Cas9_sgRNA1 - k_cat*AAV*Cas9_sgRNA1;
	d_Cas9_sgRNA2 =  Cas_gRNA_binding*Cas9*sgRNA2 - Cas_degradation*Cas9_sgRNA2 - k_cat*Cas9_sgRNA2*genome;
	d_genome = - k_cat*Cas9_sgRNA2*genome;
	d_postedit_Cas9_sgRNA1 =  k_cat*AAV*Cas9_sgRNA1 - Cas_degradation*postedit_Cas9_sgRNA1;
	d_postedit_Cas9_sgRNA2 =  k_cat*Cas9_sgRNA2*genome - Cas_degradation*postedit_Cas9_sgRNA2;
	d_edited_genome =  k_cat*Cas9_sgRNA2*genome;
	d_Cre =  alpha_p_Cre*AAV - delta_Cre*Cre;
	d_Cre_regulated_region = - k_cre*Cre_regulated_region*Cre^4 + (Cre_regulated_region/AAV)*d_AAV;
	d_edited_Cre_regulated_region =  k_cre*Cre_regulated_region*Cre^4 + (edited_Cre_regulated_region/AAV)*d_AAV;
	d_CreH =  alpha_p_CreH*(Cre_regulated_region/AAV)*AAV - delta_CreH*CreH;
	d_CreH_regulated_region = - k_cre*CreH_regulated_region*CreH^4 + (CreH_regulated_region/AAV)*d_AAV;
	d_edited_CreH_regulated_region =  k_cre*CreH_regulated_region*CreH^4 + (edited_CreH_regulated_region/AAV)*d_AAV;
	d_Cas9 =  alpha_p_Cas9*AAV - Cas_gRNA_binding*Cas9*sgRNA1 - Cas_gRNA_binding*Cas9*sgRNA2 - delta_Cas9*Cas9;
	d_sgRNA1 =  alpha_r_sgRNA1*(CreH_regulated_region/AAV)*AAV - Cas_gRNA_binding*Cas9*sgRNA1 - delta_g*sgRNA1;
	d_sgRNA2 =  alpha_r_sgRNA2*AAV - Cas_gRNA_binding*Cas9*sgRNA2 - delta_g*sgRNA2;
    
    % Pack derivatives for return, ensuring none are complex or go below zero
    dx = max(-x,real([d_AAV, d_Cas9, d_Cas9_sgRNA1, d_Cas9_sgRNA2, d_Cre, d_CreH, d_CreH_regulated_region, d_Cre_regulated_region, d_edited_CreH_regulated_region, d_edited_Cre_regulated_region, d_edited_genome, d_genome, d_postedit_Cas9_sgRNA1, d_postedit_Cas9_sgRNA2, d_sgRNA1, d_sgRNA2])');
end
