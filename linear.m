function [tint, yint] = linear(tspan, vars)
    Tc = vars(1);
    T1 = vars(2);
    Dr = vars(3);
    k1 = vars(4);
    Dp = vars(5);
    k2 = vars(6);
    Km2 = vars(7);
    
    %Linear:
%     y0=[0, 0, 0, 10, 10];
%     for i=1:30
%         a=i*24,
%         t=[0,a];
%         [x,y]=ode45(@f1,t,y0);
%         e(1,i+1)=y(end,4)/10;
%     end
   
    y0=[0, 0, 0, 10, 10];
    sol = ode45(@(t,x) f1(t, x, Tc, T1, Dr, k1, Dp, k2, Km2), tspan, y0);
    tint = 1:tspan(end);
    y = deval(sol, tint);
    yint = y(4, :)/10;
end

function dx=f1(t, x, Tc, T1, Dr, k1, Dp, k2, Km2)
    dx=[x(4)*Tc-k1*x(1)*x(2)-x(1)*Dp; % Cas9
        T1*x(5)-k1*x(1)*x(2)-x(2)*Dr; % gRNA1
        k1*x(1)*x(2)-x(3)*Dp; % Cas9-gRNA1
        -k2*x(3)*x(4)/(Km2+x(4)); % Cas9 plasmid
        0]; % gRNA1 plasmid
end
