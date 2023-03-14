clc
clear all
close all

facit = [1/120000 1/40000 1/600000]; % weak
%facit = [1/120000 1/20000 1/100000]; % strong
resolution = 201;
h = figure('units','centimeters','position',[0 0 16.8 8.5]);
rl0 = round(log(facit(1)));
rl1 = round(log(facit(2)));
plot_range = 3;
l0linspace = linspace(rl0-2*plot_range,rl0+plot_range,resolution);
l1linspace = linspace(rl1-2*plot_range,rl1+plot_range,resolution);
[L0,L1] = meshgrid(l0linspace,l1linspace);

% WRONG WRONG WRONG WRONG WRONG WRONG WRONG WRONG
ResultFolder = [cd '\Set 4\Weak\N0 512\'];
rf           = dir(ResultFolder);
Z1 = 0;
Z2 = 0;
for r = 3:length(rf)
    if contains(rf(r).name,'.mat')
        load([ResultFolder '\' rf(r).name])
        Z1 = Z1 + sum(Z,1);
        Z2 = Z2 + sum(Z,2);
    end
end
Z1 = (Z1-min(Z1));
Z2 = (Z2-min(Z2));

subplot(1,2,1)
plot(l0linspace,Z1,'k')
hold on
plot([log(facit(1)) log(facit(1))],[-1 max(Z1)*2],'k--')
xlim([l0linspace(1) l0linspace(end)])
ylim(max(Z1)*[0.5 1.1])
xlabel('\lambda_0')
subplot(1,2,2)
plot(l1linspace,Z2,'k')
hold on
plot([log(facit(2)) log(facit(2))],[-1 max(Z2)*2],'k--')
xlim([l1linspace(1) l1linspace(end)])
ylim(max(Z2)*[0.5 1.1])
xlabel('\lambda_1')