function [time_interval, y_out, y] = simple_repression(time_span, parameters, initial, step)
% time_span is the hours values [start, stop]
% parameters is a Map of names to numbers (e.g., rate constants, decay rates, Hill coefficients)
% initial is a Map of variable names to initial values
% step is the number of hours between samples in output; defaults to 1
% Returns vector of time, matrix of output levels at those time points, matrix of all species
    if nargin < 4, step = 1; end
    
    % Define names for input/output variable indexes
    AAV = 1;
	Cas9 = 2;

    % Set initial values
    y0=zeros(1,3);
    y0(AAV) = initial('AAV');
    
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
    K_R = parameters('K_R');
	alpha_p_Cas9 = parameters('alpha_p_Cas9');
	alpha_p_TF = parameters('alpha_p_TF');
	delta_Cas9 = parameters('delta_Cas9');
	delta_TF = parameters('delta_TF');
	n = parameters('n');
    
    % Unpack individual species from x
    AAV = x(1);
	Cas9 = x(2);
	TF = x(3);
    
    % Compute derivative for each species
    d_AAV = 0;
	d_Cas9 =  alpha_p_Cas9*(K_R^n)/(K_R^n + TF^n)*AAV - delta_Cas9*Cas9;
	d_TF =  alpha_p_TF*AAV - delta_TF*TF;
    
    % Pack derivatives for return
    dx = [d_AAV, d_Cas9, d_TF]';
end
