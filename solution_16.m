   dt=0.01;
T=[dt:dt:10000];
N=size(T);
x=zeros(N);
y=zeros(N);

a=-0.1;
epsilon = 0.01;
alpha = 1.8;
sigma = 1;

x0=[0, 0.5, 2, 5, 10, -10, -5, -2, -0.5];
for k=1:9
    x(1)=x0(k);
    y(1)=0;
    Ln=SDE_StableLevyMotion1(N(2),alpha);
    for i=1:N(2)-1
        x(i+1)= x(i)- dt*epsilon*(a*x(i)+y(i)/(1.0+x(i)^2));
        y(i+1)= y(i)- dt*(2*y(i)+sin(x(i)))+sigma*dt^(1.0/alpha)*Ln(i);
    end
    plot3(T, x, zeros(N),'linewidth',1)
    hold on
    plot3(zeros(N),x, y, 'linewidth',1)
end 

% plot3(T, x, zeros(N))
% hold on
% plot3(zeros(N),x, y)
xlabel('t')
ylabel('x') 
zlabel('y')
title('Solutions of the system (16) with a = -0.1','FontSize',10,'FontWeight','bold');

