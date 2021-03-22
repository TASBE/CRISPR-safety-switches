function sol = perturbableFitCas9AssumeDelay(tspan, lag, Tc, T1)
        options = [];
        sol=dde23(@ddefun, lag, [3.8039, 0, 0, 10, 0], tspan, options, Tc, T1);% Had to hardcodde history, otherwise I get "too many arguements error", not sure why

    function dydt = ddefun(t,y,Z,Tc,T1)
        ylag1 = Z(:,1);
        n=10;
        Dr=0.5358;%1.15138
        k1=0.2379;%0.3119
        Dp=0.5785;%0.2557
        k2=0.18;%0.1308
        Km2=4.7509;%11.7129

        dydt=[n*Tc-k1*y(1)*y(2)-y(1)*Dp;%Cas9 protein
            T1-k1*y(1)*y(2)-y(2)*Dr;%gRNA
            k1*y(1)*y(2)-y(3)*Dp;%Cas9-gRNA
            -k2*y(3)*y(4)/(Km2+y(4));%GFP DNA
            k2*ylag1(3)*ylag1(4)/(Km2+ylag1(4))];%mutated GFP
    end

% function s = history(t)
% s=[3.8039, 0, 0, 10, 0]';%[103.7, 0, 0, 10, 0]
end

