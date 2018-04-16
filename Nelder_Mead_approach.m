% This code is developed by Foad Karimian
% This part of code calculates MLE estimations of parameters for Lognormal
% and 3-paramter Weibull distributions

clear all;
clc;

%Actual data in form of Log(cycles) is entered in variable "data"
data=[4.615287063	4.6218	4.6295	4.6360	4.6380	4.6578	4.6725	4.6747	4.6748	4.6800	4.6839	4.6874	4.6906	4.6928	4.6945	4.6977	4.7256	4.7313	4.7603	4.7705	4.7798	4.7808	4.7998	4.8177	4.8289	4.8304	4.8320	4.8363	4.8552	4.8600];

Lognorm=mle(data);      %MLE estimation of mean and std
Lognorm=Lognorm.';      %make a vector of dist. parameters

xmin=-5.*Lognorm(2)+Lognorm(1);     %set x interval with respect to std, extending from 4*std before mean to 4*std after mean
xmax=5.*Lognorm(2)+Lognorm(1);
dx=0.0001;
x=[xmin:dx:xmax];
x=x(:);

normp = normpdf(x,Lognorm(1),Lognorm(2));    %generate pdf of lognormal distribution
normc = normcdf(x,Lognorm(1),Lognorm(2));    %generate cdf of lognormal distribution

figure(3);                          %draw pdf lognormal distribution
plot(x,normp)
title('Lognormarl PDF')
xlabel('data')
ylabel('PDF')

figure(4);                          %draw cdf lognormal distribution
plot(x,normc)
title('Lognormarl CDF')
xlabel('data')
ylabel('CDF')

movegui(figure(3),'northwest')      %place lognormal pdf on top-left on screen 
movegui(figure(4),'north')      %place lognormal cdf on top-center on screen 

custpdf = @(x,a,b,c) (x>c).*(b/a).*(((x-c)/a).^(b-1)).*exp(-((x-c)/a).^b);  %PDF for 3-parameter Weibull
opt = statset('MaxIter',1e5,'MaxFunEvals',1e5,'FunValCheck','off');
weibull= mle(data,'pdf',custpdf,'start',[0.2 2.4 4],'Options',opt,...       %MLE estimation of parameters
    'LowerBound',[0 0 -Inf],'UpperBound',[Inf Inf min(data)]);

weibull=weibull.';                  %make a vector of parameters

weipdf=custpdf(x,weibull(1),weibull(2),weibull(3));    %generate PDF of 3-parameter weibull distribtion

weicdf=zeros(length(x),1);
weicdf(1)=weipdf(1).*dx;

for i=2:length(x)
    weicdf(i)=weipdf(i).*dx+ weicdf(i-1);                   %generate CDF of 3-parameter weibull distribtion
end

figure(5)           %draw pdf 3-parameter weibull distribution
plot(x,weipdf)
title('3 Parameter Weibull PDF')
xlabel('data')
ylabel('PDF')

figure(6)           %draw cdf 3-parameter weibull distribution
plot(x,weicdf)
title('3 Parameter Weibull CDF')
xlabel('data')
ylabel('CDF')

movegui(figure(5),'west')      %place figure on left-center on screen
movegui(figure(6),'center')      %place figure on center on screen


% Applying Nelder-Mead method to find Lagrange Multipliers to generate distribution with highst
% entropy with constraint. 
% This code is developed by Foad Karimian.
% To verify the code, following data can be used based on example provided
% in Mohammad-Djafari, A. (1992). A Matlab program to calculate the maximum entropy 
% distributions: if x=[-1:0.01:1], mu=[0 0.3 0 0.15], then independent of
% initial value for lambda, final lambda values should be lambda=[0.9392 0
% -3.3414 0 4.6875].

data = data(:);
S = std(data);
Mean = mean(data);
data_st=(data-Mean)/S;

Mst=mean(data_st);
Sst=std(data_st);
sk=skewness(data_st);
kr=kurtosis(data_st);
mu=[Mst Sst sk kr];

xmin=-5;
xmax=5;
dxst=(xmax-xmin)/(length(x)+1);
x=[xmin+dxst:dxst:xmax-dxst];

mu=mu(:);                    %import mu and make a vector
x=x(:);                      %make a vector of x
M=length(mu);                %determines summation over indicies
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
    
lambda0 = log(sum(exp(-phi*lambda).*dx));   %find lambda0, normalizing factor

cdf=zeros(length(x),1);
cdf(1)=pdf(1).*dx;

for i=2:length(x)
    cdf(i)=pdf(i).*dx+cdf(i-1);     %generate cdf of MaxEnt
end


figure(1);          %draw pdf MaxEnt distribution
plot(x,pdf)
title('MaxEnt PDF')
xlabel('Data')
ylabel('PDF')


figure(2);            %draw pdf MaxEnt distribution
plot(x,cdf)
title('MaxEnt CDF')
xlabel('Data')
ylabel('CDF')

movegui(figure(1),'southwest')      %place figure on lower-left on screen
movegui(figure(2),'south')      %place figure on lower-center on screen

% Following part of code calculates MxEnt pdf and cdf for unstandardized
% data

x=(x.*std(data))+mean(data);        %Changing X interval to match interval for Lognormal and Weibull

phmun=zeros(length(x),M);           %generating x^i-mu(i) for new x's and mu's
mun=[Mean;S;sk;kr];                 %new mu is actual mean and std of data

for i=1:M
    phmun(:,i)=phi(:,i)-mun(i);       %generates x^i - mu(i) for actual data (unstandardiz)
end


qin = sum(exp(-(phmun*lambda)).*dx);     %calculate potential function value

pdfin=exp(-(phmun*lambda))./qin;         %generate PDF  

cdfin=zeros(length(x),1);                %generate CDF 
cdfin(1)=pdfin(1).*dx;

for i=2:length(x)
    cdfin(i)=pdfin(i).*dx+cdfin(i-1);     %generate cdf of MaxEnt
end

cdf_data=zeros(length(data),1);

for i=1:length(data)
    cdf_data(i)=(i-0.3)/(length(data)+0.4);
end

figure(7)                                   %Lognormal, weibull and MaxEnt PDF in one figure
plot(x,normp,':',x,weipdf,'--',x,pdfin)
title('Lognormal vs 3-p Weibull vs MaxEnt PDF')
legend('Lognormal','3 Paramter Weibull','MaxEnt')

figure(8)                                   %Lognormal, weibull and MaxEnt CDF in one figure
plot(x,normc,':',x,weicdf,'--',x,cdfin,data,cdf_data,'s')
title('Lognormal vs 3-p Weibull vs MaxEnt CDF')
legend({'Lognormal','3 Parameter Weibull','MaxEnt'},'Location','northwest')

movegui(figure(7),'northeast')      %place figure on upper-right on screen
movegui(figure(8),'east')      %place figure on right on screen

disp("MLE parameters for Lognormal distribution are:");
disp(["mean:",num2str(Lognorm(1));      %display Lognormal Distribution parameters
     "std:",num2str(Lognorm(2))]);
 
disp("MLE parameters for 3-Parameter Weibull distributions are:");
disp(["alpha:",num2str(weibull(1));    %display 3-Parameter Weibull Distribution parameters
     "beta:",num2str(weibull(2));
     "Log(n0)",num2str(weibull(3))]);
 
disp("Lagrangian multipliers (lambda0, lambda1, ...) are:");       %display Lagrangian Multipliers parameters
disp(num2str(lambda0));
disp(num2str(-lambda));

% This part of the code find r.m.s. of all distributions from data points

d_lognorm=zeros(length(data),1);        % for lognormal distribution

for i=1:length(data)
    d_lognorm(i,1)=(i-0.3)/(length(data)+0.4) - interp1(x,normc,data(i)); 
end

rms_lognorm=(sum(d_lognorm.^2)/length(data))^0.5;

d_wei=zeros(length(data),1);        % for weibull distribution

for i=1:length(data)
    d_wei(i,1)=(i-0.3)/(length(data)+0.4) - interp1(x,weicdf,data(i)); 
end

rms_wei=(sum(d_wei.^2)/length(data))^0.5;

d_maxent=zeros(length(data),1);        % for maxent distribution

for i=1:length(data)
    d_maxent(i,1)=(i-0.3)/(length(data)+0.4) - interp1(x,cdfin,data(i)); 
end

rms_maxent=(sum(d_maxent.^2)/length(data))^0.5;

disp("Tail fit for Lognormal distribution are (for two first data points):");
disp([num2str(d_lognorm(1)),'   ',num2str(d_lognorm(2))]);
disp("R.M.S. for Lognormal distribution is:");
disp(num2str(rms_lognorm));

disp("Tail fit for 3 Parameter Weibull distribution are (for two first data points):");
disp([num2str(d_wei(1)),'   ',num2str(d_wei(2))]);
disp("R.M.S. for 3 Parameter Weibull distribution is:");
disp(num2str(rms_wei));

disp("Tail fit for MaxEnt distribution are (for two first data points):");
disp([num2str(d_maxent(1)),'   ',num2str(d_maxent(2))]);
disp("R.M.S. for MaxEnt distribution is:");
disp(num2str(rms_maxent));
