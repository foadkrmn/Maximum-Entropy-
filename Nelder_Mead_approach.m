% Applying Nelder-Mead method to find Lagrange Multiplyers to generate distribution with hieghst
% entrpy with constraint. 
% This code is developed by Foad Karimian.
% To verify the code, following data can be used based on example provided
% in Mohammad-Djafari, A. (1992). A Matlab program to calculate the maximum entropy 
% distributions: if x=[-1:0.01:1], mu=[0 0.3 0 0.15], then independent of
% initial value for lambda, final lambda values should be lambda=[0.9392 0
% -3.3414 0 4.6875].

clear all;
clc;

xmin=-6;
xmax=6;
dx=0.01;
x=[xmin:dx:xmax];
mu=[0 1];

mu=mu(:);                    %import mu and make a vector
x=x(:);                      %make a vector of x
M=length(mu);                %determines sumation over indicies
phi=ones(length(x),M);       %function to generate moments (mean, variance,...)
phi(:,1)=phi(:,1).*x;        %first column is x

for i=2:M
    phi(:,i)=phi(:,i-1).*x;  %generate x^i
end

phmu=zeros(length(x),M);

for i=1:M
    phmu(:,i)=phi(:,i)-mu(i);       %generates x^i - mu(i)
end

l0=zeros(M,1);

Q = @(l) sum(exp(-phmu*l).*dx)

options = optimset('Display','iter','PlotFcns',@optimplotfval);
lambda = fminsearch(Q,l0,options);  %minimizing potential function to find Lagrangian multiplyers

q = sum(exp(-phmu*lambda).*dx);     %calculate potential value

pdf=exp(-phmu*lambda)./q;             %generate distribution  

lambda0 = log(q.*exp(-lambda.'*mu));    %find lambda0, normalizing factor

lambda=-lambda;

lambda=[lambda0;lambda]

figure(1)
plot(x,pdf)

cdf=zeros(length(x),1);
cdf(1)=pdf(1);

for i=2:length(x)
    cdf(i)=pdf(i)+cdf(i-1);
end

figure(2)
plot(x,cdf)
