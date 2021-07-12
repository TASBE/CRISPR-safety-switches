function [tint, yint] = cre(tspan, vars)
    % FUNCTION NAME:
    %   cre
    %
    % DESCRIPTION:
    %   ODE model of safety-switch with no time delay
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
    %       * Initial implementation, cleaned "creRecombination.m"
    %
    
    % Unpack variables
    Tc = vars(1);
    T1 = vars(2);
    Dr = vars(3);
    k1 = vars(4);
    Dp = vars(5);
    k2 = vars(6);
    Km2 = vars(7);
    n = vars(8);
    Ka = vars(9);
    k3 = vars(10);
    Km3 = vars(11);
    
    options = [];
    lag = [44.2974, 27.4716];
    sol = dde23(@RepressDelay, lag, history1, tspan, options, Tc, T1, ...
        Dr, k1, Dp, k2, Km2, n, Ka, k3, Km3);
    tint = 1:tspan(end);
    y = deval(sol, tint);
    yint = y(11, :)/10;
end

function dydt=RepressDelay(t, y, Z, Tc, T1, Dr, k1, Dp, k2, Km2, n, Ka, ...
        k3, Km3)
    ylag2=Z(:,2); % I think I can simplify this, and just have one number in lag. The other was being used for DNA repair, and both were used for the comb model.
    
    % Variable Names
    unknown = y(1);
    unknown = y(2);
    unknown = y(3);
    unknown = y(4);
    unknown = y(5);
    unknown = y(6);
    unknown = y(7);
    unknown = y(8);
    unknown = y(9);
    unknown = y(10);
    unknown = y(11);
    
    % Equations    
    sth = y(10)*0.666394-y(1)*0.1; % Cre Anyway to remove these magic numbers?
    sth = y(11)*Tc-k1*y(2)*y(3)-k1*y(2)*y(8)-y(2)*Dp; % Cas9
    sth = T1*y(6)*(1-y(1)^n/((Ka*(y(5)*2+y(6)))^n+y(1)^n))-k1*y(2)*y(3)-y(3)*Dr; % gRNA1 Eqn 25? Lots of things are different
    sth = k1*y(2)*y(3)-y(4)*Dp; % Cas9-gRNA1
    sth = -k3*y(1)*y(5)/(Km3+y(5)); % loxp-polyT-loxp-gRNA1
    sth = k3*ylag2(1)*ylag2(5)/(Km3+ylag2(5)); % loxp-gRNA1
    sth = -k2*y(9)*y(7)/(Km2+y(7)); % Genome
    sth = y(12)*T1-k1*y(2)*y(8)-y(8)*Dr; % gRNA2
    sth = k1*y(2)*y(8)-y(9)*Dp; % Cas9-gRNA2
    sth = -k2*y(4)*y(10)/(Km2+y(11)+y(10)); % Cre plasmid
    sth = -k2*y(4)*y(11)/(Km2+y(11)+y(10)); % Cas9 plasmid
    sth = 0; % gRNA2 plasmid (Can I remove this entirely?)
    
    % Asign results to output
    dydt=[
        ];
end

function s = history1(t)
        s=[0, 0, 0, 0, 10, 0, 10, 0, 0, 10, 10, 0]';
end