function [tint, yint] = tfDelayNoSgrna2(tspan, vars)
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

%     y0 = [0, 0, 0, 0, 0, 0, 0, 0, 46, 50, 50, 50, 0, 0, 5, 10]'; % Time delayed entire system elimination
    y0 = [0, 0, 0, 0, 0, 0, 46, 50, 50, 0, 0, 5, 10]'; % Edited to remove the sgRNA 2 numbers
    sol = ode45(@(t,x) f1(t, x, T1, T2, Tg, Dtfr, Dtfp, Dg, n, Ka, Tcr, Tcp, Dcr, Dcp, Kcg, k1, Km1), tspan, y0);
    tint = 1:tspan(end);
    y = deval(sol, tint);
    yint = y(10, :);
end

function dx=f1(t, x, T1, T2, Tg, Dtfr, Dtfp, Dg, n, Ka, Tcr, Tcp, Dcr, ...
        Dcp, Kcg, k1, Km1)
    dx=[x(7)*T1-x(1)*Dtfr; % TF
        x(1)*T2-x(2)*Dtfp; % TF
        x(9)*Tg*(1-x(2)^n/((Ka*(x(9)))^n+x(2)^n))-x(3)*x(5)*Kcg-x(3)*Dg; % gRNA1
%         x(12)*Tg*(1-x(2)^n/((Ka*(x(11)+x(12)))^n+x(2)^n))-x(4)*x(6)*Kcg-x(4)*Dg; % gRNA2
        x(8)*Tcr-x(4)*Dcr; % Cas9
        x(4)*Tcp-x(3)*x(5)*Kcg*x(5)*Kcg-x(10)*x(5)*Kcg-x(5)*Dcp; % Cas9
        x(3)*x(5)*Kcg-x(6)*Dcp; % Cas9-gRNA1
%         x(4)*x(6)*Kcg-x(8)*Dcp; % Cas9-gRNA2
        -k1*x(6)*x(7)/(Km1+x(7)+x(8)+x(9)+x(12))-k1*x(7)/(Km1+x(7)); % TF
        -k1*x(6)*x(8)/(Km1+x(7)+x(8)+x(9)+x(12)); % Cas9
        -k1*x(6)*x(9)/(Km1+x(7)+x(8)+x(8)+x(12)); % gRNA1
%         -k1*x(7)*x(12)/(Km1+x(9)+x(10)+x(11)+x(12)+x(15)) % gRNA2
        x(12)*Tg-x(10)*x(5)*Kcg-x(10)*Dg; % gRNA3
        x(10)*x(5)*Kcg-x(11)*Dcp; % Cas9-gRNA3
        -k1*x(6)*x(12)/(Km1+x(7)+x(8)+x(9)+x(12)); % gRNA3 plasmid
        -k1*x(11)*x(13)/(Km1+x(13))]; % genome
end