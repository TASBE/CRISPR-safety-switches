% Hardcode the data
xFig1 = [2 2 2 3 3 3 4 4 4 5 5 5 6 6 6 10 10 10 13 13 13];
yFig1 = [0.53 0.54, 0.59 5.62 5.83 6.26 36.9 36.2 38.7 67.8 65.6 66.9 79 76 78.7 96 97.2 96.7 95.2 96.8 96];
% Get the average for each unique time point
[G, days] = findgroups(xFig1);
yMeansFig1 = splitapply(@mean, yFig1, G);
% Convert to hours
hours = days .* 24;

% Plot just to check the data (i.e. for typos)
plot(xFig1.*24, yFig1, 'ko')
hold on;
plot(hours, yMeansFig1, 'p')

% Create the model, see the end of the file

% Fit the model using a starting point: Not sure what values to use ASK
% JAKE
parameters0 = [1 1 1 1 1 1 1];
x = lsqcurvefit(@fit_Cre_on_Kill_Switch, parameters0, hours, yMeansFig1);

% Run the model with the fit parameters
res = fit_Cre_on_Kill_Switch(x, [0 312]);

% Plot the model results
plot(xFig1.*24, yFig1, 'ko')

%% Model Function
function yResults = fit_Cre_on_Kill_Switch(parameters, t)

    % Intial values (Same number of elements as x in the DiffEq function)
    y0=ones(1, 10); % Set everything to 1? % ASK JAKE
    
    % Run ODE
    [T, Cv] = ode45(@diff_eq, t, y0); % Cv is a name I got from a function on stack exchange, change it?
    
    % ODE differential function
    function dx=diff_eq(t, x)
        % lsqcurve fit won't let me pass x0 as a map, only as a double
    	alpha_p_Cre = parameters(1);
        alpha_p_Cas9 = parameters(1);
        alpha_p_GFP = parameters(1);
        
    	delta_Cre = parameters(2);
        delta_Cas9 = parameters(2);
        delta_GFP = parameters(2);
        Cas_degradation = parameters(2);
        
        alpha_r_sgRNA1 = parameters(3);
        delta_g = parameters(4);
        
    	k_cre = parameters(5);
        Cas_gRNA_binding = parameters(6);
        k_cat = parameters(7);       
        
        % Unpack individual species from x
        AAV = x(1);
    	Cas9 = x(2);
    	Cas9_sgRNA1 = x(3);
    	Cre = x(4);
    	Cre_regulated_region = x(5);
    	edited_Cre_regulated_region = x(6);
    	postedit_Cas9_sgRNA1 = x(7);
    	sgRNA1 = x(8);
        % Add GFP
        gfpGene = x(9);
        gfp = x(10);
        
        % Compute derivative for each species
        % Nothing is targting AAV genome
        d_AAV = 0; % Does this need a different name? 
        % Cre gets expressed, degrades
        d_Cre =  alpha_p_Cre*AAV - delta_Cre*Cre;
        % Cas gets expressed, binds, degrades
        d_Cas9 =  alpha_p_Cas9*AAV - Cas_gRNA_binding*Cas9*sgRNA1 - ...
            delta_Cas9*Cas9;
        % Cre modifies the genome
        d_Cre_regulated_region = - k_cre*Cre_regulated_region*Cre^4 + ...
            (Cre_regulated_region/AAV)*d_AAV;
    	d_edited_Cre_regulated_region =  k_cre*Cre_regulated_region*Cre^4 + ...
            (edited_Cre_regulated_region/AAV)*d_AAV;
        % sgRNA is expressed
        d_sgRNA1 =  alpha_r_sgRNA1*(edited_Cre_regulated_region/AAV)*AAV - ...
            Cas_gRNA_binding*Cas9*sgRNA1 - delta_g*sgRNA1;
        % sgRNA binds to Cas and edits  
    	d_Cas9_sgRNA1 =  Cas_gRNA_binding*Cas9*sgRNA1 - ...
            Cas_degradation*Cas9_sgRNA1 - k_cat*gfpGene*Cas9_sgRNA1;
        d_postedit_Cas9_sgRNA1 =  k_cat*gfpGene*Cas9_sgRNA1 - ...
            Cas_degradation*postedit_Cas9_sgRNA1;
        % gfp has its own gene
        d_gfpGene = -k_cat*gfpGene*Cas9_sgRNA1;
        % gfp gets expressed
        d_gfp = alpha_p_GFP*gfpGene - delta_GFP*gfp;
    
        % Pack derivatives for return
        dx = [d_AAV, d_Cas9, d_Cas9_sgRNA1, d_Cre, d_Cre_regulated_region, ...
            d_edited_Cre_regulated_region, d_postedit_Cas9_sgRNA1, d_sgRNA1, d_gfpGene, d_gfp]';
    end
    
    
    % Solve for just the percent of GFP (what the paper is reporting)
    yResults = (Cv(:, 10) / 1)'; % Change 1 to be the original amount, ASK JAKE
end



