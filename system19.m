clear
a=-0.1;
epsilon=0.01;
alpha=1.8;
sigma=1;

dt=0.01;
t=[0 : dt : 100];
TN=size(t);
x=zeros(TN);
int_s= zeros(TN);
int_t= zeros(TN);
int_Lr= zeros(TN);
int_theta= zeros(TN);

Lr= SDE_StableLevyMotion1(TN(2), alpha);
location_t0= find(t==0);
for i=1:TN(2)
    if t(i)<0
        int_s(i)= dt*sum(exp(-2.0*(t(i) - t(1:i))) * dt^(1.0/alpha).*Lr(1:i));
    elseif t(i)==0
        int_s(i)=0;
    else
        int_s(i)= dt*sum(exp(-2.0*(t(i) - t(location_t0:i))) * dt^(1.0/alpha).*Lr(location_t0:i));
    end
end
for i= 1: TN(2)
    int_t(i)=  dt*sum( exp(-2.0*(t(i)-t(1:i)) ) .* int_s(1:i) );
    int_theta(i)= dt*sum(exp(-2.0*(t(i)-t(1:i)) ) * dt^(1.0/alpha).*Lr(1:i));
end

x0=[0, 0.5, 2, 5, 10, -10, -5, -2, -0.5];
for k=1:9
    x(1)=x0(k);
    for i=1:TN(2)-1
        x(i+1) = x(i) - epsilon*a*x(i) - epsilon/(1+x(i)^2)*( sigma*int_theta(i) - sin(x(i))/2.0 + epsilon*( -a*x(i)*cos(x(i))/4.0 +...
            sin(x(i))*cos(x(i))/8.0/(1.0+x(i)^2) ) + epsilon*sigma*cos(x(i))/(1.0+x(i)^2)*int_t(i) );
    end

    hold on
    plot(t, x, 'linewidth', 2)
%   set(gca,'XLIM', [0 1000])
    box on
end
title('Solutions of the reduced system (19) with a = -0.1','FontSize',10,'FontWeight','bold');
