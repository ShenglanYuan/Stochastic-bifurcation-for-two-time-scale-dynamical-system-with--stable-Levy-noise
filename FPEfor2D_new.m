clc;clear;
close all
tic
%system parameters
global alpha
alpha = 1.8; %alpha-stable index
eps = 0.01; %scale
sigma = 1; a =-0.1; %coefficients
x0 = 0; y0 = 0; %initial state
C_alpha = (alpha/(2^(1-alpha)*sqrt(pi))) * (gamma((1+alpha)/2)/gamma(1-alpha/2));

%discretization
global Lx1 Lx2 Ly1 Ly2
Lx1 = -1.; Lx2 = 1.; %x in (Lx1,Lx2)
Ly1 = -1.; Ly2 = 1.; %y in (Ly1,Ly2)
dx = 0.05; dy = 0.05;
dt = 0.0001;
x = Lx1:dx:Lx2; y = Ly1:dy:Ly2;
nx = length(x); ny = length(y);
P = zeros((nx-2)*(ny-2),1);

%initial 
mu=[x0,y0];% 均值向量
Sigma=[2*dt*sigma/(eps^(1/alpha)) 0;0 2*dt*sigma/(eps^(1/alpha))];% 协方差矩阵
[X,Y]=meshgrid(Lx1:dx:Lx2,Ly1:dy:Ly2);%在XOY面上，产生网格数据
p=mvnpdf([X(:) Y(:)],mu,Sigma);%求取联合概率密度，相当于Z轴 
p=reshape(p,size(X));
figure
surf(X,Y,p),axis tight,title('二维正态分布图')

pp = p(2:end-1,2:end-1);
P = pp(:);
%differential part
count = 1;
tmp = sigma^alpha * C_alpha / (eps*alpha);
for i = 2:nx-1
    for j = 2:ny-1
        s1(count) = myfun(x(i),y(j),a,1) + 0*2/eps - tmp*myfun(x(i),y(j),a,4);
        s2(count) =  myfun(x(i),y(j),a,3)/(2*dy*eps);
        s3(count) =  -myfun(x(i),y(j),a,3)/(2*dy*eps);
        s4(count) =  myfun(x(i),y(j),a,2)/(2*dx);
        s5(count) =  -myfun(x(i),y(j),a,2)/(2*dx);
%         tmp = sigma^alpha * C_alpha / (eps*alpha);
%         B(count) = tmp*(1/(y(j)-Ly1)^alpha + 1/(Ly2-y(j))^alpha); 
        if j == ny-1
            s2(count) = 0;
        end
        if j == 2
            s3(count) = 0;
        end
        count = count + 1;    
    end
end
A = zeros(nx-2,ny-2);
A =diag(s1) + diag(s2(1:end-1),1) + diag(s3(2:end),-1) + diag(s4(1:end-(nx-2)),nx-2) + diag(s5(nx-2+1:end),-(nx-2)); %square matrix

%integral part
C1 = zeros((nx)*(ny)); C2 = C1;BB = C1;
count1 = 1; count2 = 1; count3 = 0;

for i = 2:nx-1
    for j = 2:ny-1
        for k = -(j-1):ny-j
            if k == 0
                continue
            end
            if k >= -(j-2) && k <= ny-j-1
                C1(i+(nx-1)*count2+1+count3,j+(ny-1)*count1+i-1+k) = dy/(abs(k*dy)^(1+alpha));
                C2(i+(nx-1)*count2+1+count3,j+(ny-1)*count1+i-1) = C2(i+(nx-1)*count2+1+count3,j+(ny-1)*count1+i-1) + dy/(abs(k*dy)^(1+alpha));
            end
            
        end 
        count3 = count3 + 1;
    end
    count3 = 0;
    count1 = count1 + 1;
    count2 = count2 + 1;
end
CC = sigma^alpha*C_alpha*(C1 - C2)/eps;
for i = 2:nx-1
    for j = 2:ny-1
        C((i-2)*(nx-2)+1:(i-2)*(nx-2)+nx-2,(j-2)*(ny-2)+1:(j-2)*(ny-2)+ny-2) = ...
            CC((i-1)*(nx)+2:(i-1)*(nx)+nx-1,(j-1)*(ny)+2:(j-1)*(ny)+ny-1);
    end
end


%Time iteration
L=500;
v = 1; w = 1;
for t = 1:L
    P1 = P + w*(A + v*C)*dt*P;
%     P = P1;
    P2 = 3/4*P + 1/4*P1 + 1/4*w*(A + v*C)*dt*P1;
    P = 1/3*P + 2/3*P2 + 2/3*w*(A + v*C)*dt*P2;
    if t == 1 || t== 10 || t == 20 || t == 30|| t == 40 ||t == 50 || t ==100 || t ==300 || t ==500
    figure
    Q = zeros(nx,ny);
    b = reshape(P,nx-2,ny-2)';
    Q(2:nx-1,2:ny-1) = (b);
    surf(X,Y,Q)
    xlabel('x')
    ylabel('y')
    zlabel('p(x,y,t)')
    %title('Probability density function of non-local Fokker-Planck equation (15) with t=0','FontSize',10,'FontWeight','bold');
%     axis([-2 2 -2 ])
    end
end

%Plot
% figure
% Q = zeros(nx,ny);
% b = reshape(P,nx-2,ny-2)';
% Q(2:nx-1,2:ny-1) = (b);
% surf(X,Y,Q)
% xlabel('x')
% ylabel('y')
% myfun(x(4),y(4),a,3)/eps
toc
function f = myfun(x,y,a,index)
global Ly1 Ly2 alpha
    if index == 1
        f = a - 2*x*y/(1+x^2)^2;
%         f = -0;
    else if index == 2
            f = a*x + y/(1+x^2);
%             f = 0;
        else if index == 3
            f = 2*y + sin(x);
%             f = 0;
            else
                f = 1/(y-Ly1)^alpha + 1/(Ly2-y)^alpha;
            end
        end
    end
end

% eps = 0.01;
% syms xx
% f=int(1/(-xx)^(1+alpha),xx,-1,-eps)  + int(1/(xx)^(1+alpha),xx,eps,1)



