% Applying Nelder-Mead method to find lamda for MaxEnt 
clear all;

% xmin=-1.58;
% xmax=1.78;
% dx=0.01;
% x=[xmin:dx:xmax];
% mu=[0 1 0.3];
% M=length(mu);

lambda0=[-0.1,0.6,0.3,-0.3]
l = fminsearch(@maxent,lambda0)

% for i=1:M
%     phi=zeros(length(x),M);       
%     phi(:,1)=ones(size(x)).*x;
%     for i=2:M
%          phi(:,M)=x.*phi(:,M-1);   %next columns generates x^i
%     end
%     l0=sum(dx.*exp(-sum(l.'*phi(:,i))));
% end


function f = maxent(lambda)     %fuction to create Q

xmin=-1.58;
xmax=1.78;
dx=0.01;
x=[xmin:dx:xmax];
mu=[0 1 0.2884 1.9006];

mu=mu(:);                     %import mu and make a vector
x=x(:);                       %make a vector of x
dx=x(2)-x(1);                 %find dx
M=length(mu);                 %determines sumation over indicies
phi=zeros(length(x),M);       %function to generate moments (mean, variance,...)
phi(:,1)=ones(size(x)).*x;    %first column is x

for i=2:M
    phi(:,M)=x.*phi(:,M-1);   %next columns generates x^i
end

for k=1:M
    
    f = sum(dx.*exp(-sum(lambda(k).*(phi(:,k)-mu(k)))))
end

end