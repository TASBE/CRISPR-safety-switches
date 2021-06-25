function [tint, yint] = creAndDna(tspan, vars)
    % Unpack variables
    Tc = vars(1); % Cas9 generation rate. transcription and translation is simplified as 1 step here (Taken from FitCas9assumeDelay)
    T1 = vars(2); % gRNA transcription rate. This can also be written as n*T1.
    Dr = vars(3); % gRNA degradation rate.
    k1 = vars(4); % Cas9 & gRNA binding rate.
    Dp = vars(5); % Cas9 protein degradation rate. Cas9-gRNA complex also use this dagradation rate,
% it's assumed that gRNA is stabilized by binding with Cas9.
    k2 = vars(6); % "Parameter describing the activity of Cre recombination"
    Km2 = vars(7); % "Parameter describing the activity of Cre recombination"
    n = vars(8); % Hill coefficient
    Ka = vars(9); % Promoter binding
    k3 = vars(10); % ?
    Km3 = vars(11); % ?

    options = [];
    lag = [44.2974, 27.4716];
    sol = dde23(@CreDNArepair, lag, history3, tspan, options, Tc, T1, ...
        Dr, k1, Dp, k2, Km2, n, Ka, k3, Km3);
    tint = 1:tspan(end);
    y = deval(sol, tint);
    yint = y(11, :)/10;
end

function dydt=CreDNArepair(t,y,Z, Tc, T1, Dr, k1, Dp, k2, Km2, n, Ka, ...
        k3, Km3)
    ylag1=Z(:,1);
    ylag2=Z(:,2);

    dydt=[y(11)*Tc-k1*y(1)*y(2)-k1*y(1)*y(9)-y(1)*Dp;%Cas9
        y(13)*T1*(1-y(6)^n/((Ka*(y(12)*2+y(13)))^n+y(6)^n))-k1*y(1)*y(2)-y(2)*Dr;%gRNA1
        k1*y(1)*y(2)-y(3)*Dp;%Cas9-gRNA1
        -k2*y(10)*y(4)/(Km2+y(4)+y(11));%Cre plasmid
        0;%nothing
        y(4)*0.666394-y(6)*0.1;%Cre
        -k2*y(3)*y(7)/(Km2+y(7));%Terminator-gRNA2 plasmid
        k3*ylag1(3)*ylag1(7)/(Km2+ylag1(7));%gRNA2 plasmid
        T1*y(8)-k1*y(1)*y(9)-y(9)*Dr;%gRNA2
        k1*y(1)*y(9)-y(10)*Dp;%Cas9-gRNA2
        -k2*y(10)*y(11)/(Km2+y(4)+y(11));%Cas9 plasmid
        -k3*y(6)*y(12)/(Km3+y(12));%loxp-polyT-loxp gRNA1 plasmid
        k3*ylag2(6)*ylag2(12)/(Km3+ylag2(12))];%loxp-gRNA1 plasmid
end

function s = history3(t)
        s=[0, 0, 0, 10, 0, 0, 100, 0, 0, 0, 10, 100, 0]';
end
