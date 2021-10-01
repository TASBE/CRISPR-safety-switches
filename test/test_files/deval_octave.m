function y=deval(sol,t)
       %% y=deval(sol,t)
       %% sol= solution of ode y'(t)=f(t,y(t))  y(0)=y0
       %%       obtained with sol=ode45(@f,[t0,T],y0)
       %%   t=  vector of times [t(1) ...t(n)]
       %%   y= vector of solution values [y(t(1)) ... y(t(n))]
       %%      interpolated  from sol data

       Y=sol.y;
       x=sol.x;
       n=length(t);
       l=size(Y,1);
       y=zeros(l,n);
       K=5;
       for k=1:l
              x_tmp=x(1:3);
              ind=find((t<x(2)));
              P=polyfit(x_tmp,Y(k,1:3),K);
              new_y=polyval(P,t(ind)) ;
              for p=2:length(x)-2
                     x_tmp=x(p-1:p+2);
                     ind=find((t>=x(p))&(t<x(p+1)));
                     P=polyfit(x_tmp,Y(k,p-1:p+2),K);
                     new_y=[new_y polyval(P,t(ind))];
              end
              x_tmp=x(end-2:end);
              ind=find((t>=x(end-1)));
              P=polyfit(x_tmp,Y(k,end-2:end),K);
              new_y=[new_y polyval(P,t(ind))];
              y(k,:)=new_y;
       end
