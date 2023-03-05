clc
clear all
close all

ResultFolder = [cd '\Set 3\Weak'];
rf = dir(ResultFolder);
load([ResultFolder '\' rf(3).name]);
resolution = 201;
facit = [1/120000 1/40000 1/600000]; % weak
%facit = [1/120000 1/20000 1/100000]; % strong

rl0 = round(log(facit(1)));
rl1 = round(log(facit(2)));
plot_range = 3;
l0linspace = linspace(rl0-2*plot_range,rl0+plot_range,resolution);
l1linspace = linspace(rl1-2*plot_range,rl1+plot_range,resolution);
[L0,L1] = meshgrid(l0linspace,l1linspace);

%scatter(ML_lambda(:,1),ML_lambda(:,2))
Z = zeros(resolution);
S = cov(ML_lambda(:,1), ML_lambda(:,2));
h = S*(1/(length(ML_lambda))^(1/3));
hrho = corr(ML_lambda(:,1), ML_lambda(:,2));
for e = 1:length(ML_lambda)
    mu_e = [ML_lambda(e,1); ML_lambda(e,2)];
    Z = Z + exp(-((L0-mu_e(1)).^2/h(1,1) + (L1-mu_e(2)).^2)/h(2,2) - 2./(h(1,2)*h(2,1))*((L0-mu_e(1)).*(L1-mu_e(2))));
    si_e = h*eye(2);
end
surf(L0,L1,Z,'EdgeColor','none')
view(0,90)
xlim([rl0-2*plot_range rl0+plot_range])
ylim([rl1-2*plot_range rl1+plot_range])