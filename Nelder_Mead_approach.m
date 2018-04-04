% Applying Nelder-Mead method to find lamda for MaxEnt 
clear all;

xmin=-1.58;
xmax=1.78;
dx=0.01;
x=[xmin:dx:xmax];
mu=[0 1 0.2884 1.9006];

mu=mu(:)                     %import mu and make a vector
x=x(:);                      %make a vector of x
dx=x(2)-x(1);                %find dx
M=length(mu);                %determines sumation over indicies
phi=ones(length(x),M);       %function to generate moments (mean, variance,...)
phi(:,1)=phi(:,1).*x;        %first column is x

for i=2:M
    phi(:,i)=phi(:,i-1).*x;  %generate x^i
end


for k=1:M
    
    Q = @(lambda) sum(dx.*exp(-sum(lambda(k).*(phi(:,k)-mu(k)),2))); %Potential function 

end

lambda0=[-0.5,0.6, 0.2, -1];
options = optimset('Display','iter','PlotFcns',@optimplotfval);
l = fminsearch(Q,lambda0,options)  %minimizing potential function to find Lagrangian multiplyers


for k=1:M
    q = sum(dx.*exp(-sum(l(k).*(phi(:,k)-mu(k)),2))); %find potential function value
end

f=ones(length(x),1);

for k=1:M
    f=(exp(-sum(l(k).*(phi(:,k)-mu(k)),2)))/q;
end
s=sum(f)*dx
plot(x,f)

    
