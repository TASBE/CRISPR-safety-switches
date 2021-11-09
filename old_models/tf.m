function [tint, yint] = tf(tspan, vars)
    % FUNCTION NAME:
    %   tf
    %
    % DESCRIPTION:
    %   ODE model of safety-switch with time delay controled by TF
    %
    % INPUT:
    %   tspan - (double) Vector of the min and max time
    %   vars - (cell) List containing all variables used for all models
    %
    % OUTPUT:
    %   tint - (double) Vector of all time points
    %   yint - (double) Vecotr of all calculated values (Cas9 concentration
    %   in nm)
    %
    % ASSUMPTIONS AND LIMITATIONS:
    %   None
    %
    % REVISION HISTORY:
    %   07/11/2021 - HS
    %       * Initial implementation, cleaned "GeneTherapySystemElimination.m"
    %
    
    T1 = vars(1); % TF trasncription rate
    T2 = vars(2); % TF translation rate
    Tg = vars(3); % gRNA transcription rate
    Dtfr = vars(4); % TF mRNA degradation rate
    Dtfp = vars(5); % TF protein degradation rate
    Dg = vars(6); % gRNA degradation rate
    n = vars(7); % Hill coefficience 
    Ka = vars(8); % parameter for TF binding promoter
    Tcr = vars(9); % Cas9 transcription rate
    Tcp = vars(10); % Cas9 translation rate
    Dcr = vars(11); % Cas9 mRNA degradation rate
    Dcp = vars(12); % Cas9 protein degradation rate
    Kcg = vars(13); % Cas9-gRNA binding rate
    k1 = vars(14); % parameter for Cas9-gRNA cutting
    Km1 = vars(15); % parameter for Cas9-gRNA cutting

    y0 = [0, 0, 0, 0, 0, 0, 0, 0, 46, 50, 50, 50, 0, 0, 5, 10]'; % Time delayed entire system elimination
    sol = ode45(@(t,x) f1(t, x, T1, T2, Tg, Dtfr, Dtfp, Dg, n, Ka, Tcr, Tcp, Dcr, Dcp, Kcg, k1, Km1), tspan, y0);
    tint = 1:tspan(end);
    y = deval(sol, tint);
    yint = y(10, :);
end

function dx=f1(t, x, T1, T2, Tg, Dtfr, Dtfp, Dg, n, Ka, Tcr, Tcp, Dcr, ...
        Dcp, Kcg, k1, Km1)
   % Variable names
    TFr = x(1); % [TFr]
    TFp = x(2); % [TFp]
    g1 = x(3); % [g1]
    g2 = x(4); % [g2]
    Cr = x(5); % [Cr]
    Cp = x(6);  % [Cp]
    Cg1 = x(7); % [Cg1]
    Cg2= x(8); % [Cg2]
    n1 = x(9); % [n1]
    n3 = x(10); % [n3]
    n2 = x(11); % [n2]
    n4 = x(12); % [n4]
    g3 = x(13); % [g3]
    Cg3 = x(14); % [Cg3]
    n5 = x(15); % [n5]
    n6 = x(16); % [n6]

    % Equations
    dTFrdt = n1*T1-TFr*Dtfr; % TF Eqn 1
    dTFpdt = TFr*T2-TFp*Dtfp; % TF Eqn 2
    dg1dt = n2*Tg*(1-TFp^n/((Ka*(n2+n4))^n+TFp^n))-g1*Cp*Kcg-g1*Dg; % gRNA1 Eqn 11
    dg2dt = n4*Tg*(1-TFp^n/((Ka*(n2+n4))^n+TFp^n))-g2*Cp*Kcg-g2*Dg; % gRNA2 Eqn 12
    dCrdt = n3*Tcr-Cr*Dcr; % Cas9 Eqn 4
    dCpdt = Cr*Tcp-g1*Cp*Kcg-g2*Cp*Kcg-g3*Cp*Kcg-Cp*Dcp; % Cas9 Eqn 14 (but expanded?)
    dCg1dt = g1*Cp*Kcg-Cg1*Dcp; % Cas9-gRNA1 Eqn 7
    dCg2dt = g2*Cp*Kcg-x(8)*Dcp; % Cas9-gRNA2 Eqn 15 (But should it be Dcg?)
    dn1dt = -k1*Cg1*n1/(Km1+n1+n3+n2+n4+n5)-k1*Cg2*n1/(Km1+n1); % TF Eqn 17 (but missing minus nx times d at the end?)
    dn3dt = -k1*Cg1*n3/(Km1+n1+n3+n2+n4+n5); % Cas9 Eqn 19 (but missing minus nx times d at the end?)
    dn2dt = -k1*Cg1*n2/(Km1+n1+n3+n2+n4+n5); % gRNA1 Eqn 18 (but missing minus nx times d at the end?)
    dn4dt = -k1*Cg1*n4/(Km1+n1+n3+n2+n4+n5); % gRNA2 Eqn 20 (but missing minus nx times d at the end?)
    dg3dt = n5*Tg-g3*Cp*Kcg-g3*Dg; % gRNA3 Eqn 13
    dCg3dt = g3*Cp*Kcg-Cg3*Dcp; % Cas9-gRNA3 Eqn 16 Uses Dcg in the paper
    dn5dt = -k1*Cg1*n5/(Km1+n1+n3+n2+n4+n5); % gRNA3 plasmid Eqn 21 (but missing minus nx times d at the end?)
    dn6dt = -k1*Cg3*n6/(Km1+n6); % genome Eqn 22 (but missing minus nx times d at the end?)
    
    % Asign results to output
    dx = [dTFrdt;
        dTFpdt;
        dg1dt;
        dg2dt;
        dCrdt;
        dCpdt;
        dCg1dt;
        dCg2dt;
        dn1dt;
        dn3dt;
        dn2dt;
        dn4dt;
        dg3dt;
        dCg3dt;
        dn5dt;
        dn6dt];
end
