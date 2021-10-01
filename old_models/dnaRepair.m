function [tint, yint] = dnaRepair(tspan, vars)
    % Unpack variables
    Tc = vars(1);
    T1 = vars(2);
    Dr = vars(3);
    k1 = vars(4);
    Dp = vars(5);
    k2 = vars(6);
    Km2 = vars(7);
    
    options = [];
    lag = [44.2974, 27.4716];
    sol = dde23(@DNArepair, lag, history2, tspan, options, Tc, T1, ...
        Dr, k1, Dp, k2, Km2);
    tint = 1:tspan(end);
    y = deval(sol, tint);
    yint = y(8, :)/10;
end

function dydt=DNArepair(t,y,Z, Tc, T1, Dr, k1, Dp, k2, Km2);
    ylag1=Z(:,1);

    dydt=[y(8)*Tc-k1*y(1)*y(2)-k1*y(1)*y(6)-y(1)*Dp;%Cas9
        y(9)*T1-k1*y(1)*y(2)-y(2)*Dr;%gRNA1
        k1*y(1)*y(2)-y(3)*Dp;%Cas9-gRNA1
        -k2*y(3)*y(4)/(Km2+y(4));%Terminator-gRNA2 plasmid
        k2*ylag1(3)*ylag1(4)/(Km2+ylag1(4));%gRNA2 plasmid
        y(5)*T1-k1*y(1)*y(6)-y(6)*Dr;%gRNA2
        k1*y(1)*y(6)-y(7)*Dp;%Cas9-gRNA2
        -k2*y(7)*y(8)/(Km2+y(8));%Cas9 plasmid
        0];%gRNA1 plasmid
end

function s = history2(t)
        s=[0, 0, 0, 10, 0, 0, 0, 10, 10]';
end