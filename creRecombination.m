function [tint, yint] = creRecombination(tspan, vars)
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

    %Cre recombination:
%     for i=1:30
%         a=i*24,
%         tspan=[0,a];
%         sol=dde23(@RepressDelay, lag, @history1, tspan);
%         e(2,i+1)=sol.y(11,end)/10;
%     end
    
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
    ylag2=Z(:,2);
        
    dydt=[y(10)*0.666394-y(1)*0.1;%Cre
        y(11)*Tc-k1*y(2)*y(3)-k1*y(2)*y(8)-y(2)*Dp;%Cas9
        T1*y(6)*(1-y(1)^n/((Ka*(y(5)*2+y(6)))^n+y(1)^n))-k1*y(2)*y(3)-y(3)*Dr;%gRNA1
        k1*y(2)*y(3)-y(4)*Dp;%Cas9-gRNA1
        -k3*y(1)*y(5)/(Km3+y(5));%loxp-polyT-loxp-gRNA1
        k3*ylag2(1)*ylag2(5)/(Km3+ylag2(5));%loxp-gRNA1
        -k2*y(9)*y(7)/(Km2+y(7));%Genome
        y(12)*T1-k1*y(2)*y(8)-y(8)*Dr;%gRNA2
        k1*y(2)*y(8)-y(9)*Dp;%Cas9-gRNA2
        -k2*y(4)*y(10)/(Km2+y(11)+y(10));%Cre plasmid
        -k2*y(4)*y(11)/(Km2+y(11)+y(10));%Cas9 plasmid
        0];%gRNA2 plasmid;
end

function s = history1(t)
        s=[0, 0, 0, 0, 10, 0, 10, 0, 0, 10, 10, 0]';
end