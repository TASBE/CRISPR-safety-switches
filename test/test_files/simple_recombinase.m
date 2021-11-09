function [time_interval, y_out, y] = Recombinase(time_span, parameters, initial, step)
% time_span is the hours values [start, stop]
% parameters is a Map of names to numbers (e.g., rate constants, decay rates, Hill coefficients)
% initial is a Map of variable names to initial values
% step is the number of hours between samples in output; defaults to 1
% Returns vector of time, matrix of output levels at those time points, matrix of all species
    if nargin < 4, step = 1; end
    
    % Define names for input/output variable indexes
    AAV = 1;
	Cas9 = 2;
	Cre_regulated_region = 4;

    % Set initial values
    y0=zeros(1,5);
    y0(AAV) = initial('AAV');
	y0(Cre_regulated_region) = initial('Cre_regulated_region');
    
    % Run ODE
    solution = ode45(@(t,x) diff_eq(t, x, parameters), time_span, y0);
    
    % Evaluate species levels at given times
    time_interval = time_span(1):step:time_span(end);
    y = deval(solution, time_interval);
    y_out = y([Cas9],:);
end

% ODE differential function
function dx=diff_eq(t, x, parameters)
    % Unpack parameters from parameter map
    alpha_p_Cas9 = parameters('alpha_p_Cas9');
	alpha_p_Cre = parameters('alpha_p_Cre');
	delta_Cas9 = parameters('delta_Cas9');
	delta_Cre = parameters('delta_Cre');
	k_cre = parameters('k_cre');
    
    % Unpack individual species from x
    AAV = x(1);
	Cas9 = x(2);
	Cre = x(3);
	Cre_regulated_region = x(4);
	edited_Cre_regulated_region = x(5);
    
    % Compute derivative for each species
    d_AAV = 0;
	d_Cas9 =  alpha_p_Cas9*(edited_Cre_regulated_region/AAV)*AAV - delta_Cas9*Cas9 - delta_Cas9*Cas9;
	d_Cre =  alpha_p_Cre*AAV - delta_Cre*Cre;
	d_Cre_regulated_region = - k_cre*Cre_regulated_region*Cre^4 + (Cre_regulated_region/AAV)*d_AAV;
	d_edited_Cre_regulated_region =  k_cre*Cre_regulated_region*Cre^4 + (edited_Cre_regulated_region/AAV)*d_AAV;
    
    % Pack derivatives for return
    dx = [d_AAV, d_Cas9, d_Cre, d_Cre_regulated_region, d_edited_Cre_regulated_region]';
end
