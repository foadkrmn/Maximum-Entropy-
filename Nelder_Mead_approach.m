% Applying Nelder-Mead method to find Lagrange Multiplyers to generate distribution with hieghst
% entrpy with constraint. 
% This code is developed by Foad Karimian.

clear all;
clc;

xmin=-3;
xmax=3;
dx=0.001;
x=[xmin:dx:xmax];
mu=[0 1 0.2884 1.9006];

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

p=exp(-phmu*lambda)./q;             %generate distribution

lambda0 = log(sum(exp(-(phi*lambda)).*dx));  %find lambda0, normalizing factor

lambda=[lambda0;lambda];
lambda

plot(x,p)
