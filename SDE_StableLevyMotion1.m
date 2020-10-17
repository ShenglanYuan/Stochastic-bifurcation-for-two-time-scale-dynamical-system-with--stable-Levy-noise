% SDE_StableLevyMotion.m

% Simulate sample solution paths of the following SDE drien by \alpha-stable Levy process
%dX_t=(X_t-X^3_t)dt + dL^{\alpha}_t,     where L^{\alpha}_t is the \alpha-stable Levy process and t \in[0,1]

% In the code, I use the Euler Scheme to get the algorithm
%X_{t_{i+1}}=(X_{t_i}-X^3_{t_i})\Delta t+\Delta L^{\a}_t
%and I let alpha=1.9 which can be chaged to any real number in [0.2]

% Reference BOOK: A Janocki, A Weron, Simulation and chaotic 
% behavior of alpha-stable stochastic processes. Marcel Dekker, Inc. 1994. 
% All the following citations are about this book.

% The levy simulation here is about S_alpha(1, 0, 0), and we can get other
% simulations S_alpha(A, 0, 0) by the properties1-2 on page24.
function X=SDE_StableLevyMotion1(M,a)

%initial setting
alpha = a;                           %index of stability,  
T=10;                                   %the time interval                                 %number of step to take

t=T*(0:M)/M;                             %the time steps
Gamma = rand(1, M)*pi-pi/2;            %generate uniformly distributed random number on (-\pi/2,\pi/2)             
ISENorm = -log(rand(1, M));            
%%% it can be proved easily via the direct derivation of the Distribution function
deltaT=T/M;                            %time size
%main program
DeltaL=zeros(1,M);
X=zeros(1, M+1);
for i=1:M
    DeltaL(i)=(t(i+1)-t(i))^(1/alpha)*sin(alpha*Gamma(i))*(cos((1-alpha)*...
           Gamma(i))/ISENorm(i))^(1/alpha-1)/(cos(Gamma(i)))^(1/alpha);
       %%% 1--The formula is from (3.5.1) on Page 48 of the above book. I check 
       %%%    this expression due to (3.5.1)
       %%% 2-- The first term "(t(i+1)-t(i))^(1/alpha)" is obtained by of Levy motion
       %%%  from definition 2.5.2(3) on Page30. i.e., X(t)-X(s)\sim
       %%%  S_alpha((t-s)^(1/alpha),beta,0).
end
X=DeltaL;