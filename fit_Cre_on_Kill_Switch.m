%% Model Function
function yResults = fit_Cre_on_Kill_Switch(parameters, t)
    
    t = [0 t];

    % Intial values (Same number of elements as x in the DiffEq function)
    % Postitions 1, 5, & 9 (all the genes) start at 1, everything else is 0
    y0 = zeros(1, 10);
    y0(1) = 1; % Lentivirus (LV)
    y0(5) = 1; % Cre_regulated_region
    y0(9) = 1; % GFP gene
    
    % Print parameters
    %parameters
    
    % Run ODE
    [T, Cv] = ode15s(@diff_eq, t, y0); % Cv is a name I got from a function on stack exchange, change it?
    
    % ODE differential function
    function dx=diff_eq(t, x)
        % lsqcurve fit won't let me pass x0 as a map, only as a double
    	alpha_p_Cre = 10^parameters(1);
        alpha_p_Cas9 = 10^parameters(1);
        alpha_p_GFP = 10^parameters(1);
        
    	delta_Cre = 10^parameters(2);
        delta_Cas9 = 10^parameters(2);
        delta_GFP = 10^parameters(2);
        Cas_degradation = 10^parameters(2);
        
        alpha_r_sgRNA1 = 10^parameters(3);
        delta_g = 10^parameters(4);
        
    	k_cre = 10^parameters(5);
        Cas_gRNA_binding = 10^parameters(6);
        k_cat = 10^parameters(7);
        
        % Unpack individual species from x
        LV = x(1);
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
        d_Cre =  alpha_p_Cre*LV - delta_Cre*Cre;
        % Cas gets expressed, binds, degrades
        d_Cas9 =  alpha_p_Cas9*LV - Cas_gRNA_binding*Cas9*sgRNA1 - ...
            delta_Cas9*Cas9;
        % Cre modifies the genome
        d_Cre_regulated_region = - k_cre*Cre_regulated_region*Cre^4 + ...
            (Cre_regulated_region/LV)*d_AAV;
    	d_edited_Cre_regulated_region =  k_cre*Cre_regulated_region*Cre^4 + ...
            (edited_Cre_regulated_region/LV)*d_AAV;
        % sgRNA is expressed
        d_sgRNA1 =  alpha_r_sgRNA1*(edited_Cre_regulated_region/LV)*LV - ...
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
    
%     parameters
    % Calculate max GFP: Solve for GFP when gfp_gene is equal to 1 and
    % d_gfp is equal to 0
    % GFP = alpha_p_GFP / delta_GFP
    max_GFP = 10^parameters(1) / 10^parameters(2);
    % Calculate percent negative GFP at each time point
    %yResults = (1 - (Cv(2:end, 10) / max_GFP))' * 100; % 
    yResults = max(0,(1 - (Cv(2:end, 10) / (max_GFP*0.1))))' * 100; % alternate version: fraction at <10% max
end

