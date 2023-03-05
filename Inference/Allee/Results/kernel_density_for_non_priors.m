clc
clear all
close all

ResultFolder = [cd '\Set 3\Weak'];
rf = dir(ResultFolder);
load([ResultFolder '\' rf(3).name]);
%%
resolution = 201;
facit = [1/120000 1/40000 1/600000]; % weak
%facit = [1/120000 1/20000 1/100000]; % strong
h = figure('units','centimeters','position',[0 0 16.8 8.5]);
rl0 = round(log(facit(1)));
rl1 = round(log(facit(2)));
plot_range = 3;
l0linspace = linspace(rl0-2*plot_range,rl0+plot_range,resolution);
l1linspace = linspace(rl1-2*plot_range,rl1+plot_range,resolution);
[L0,L1] = meshgrid(l0linspace,l1linspace);
subplot(1,2,1)
%Z = KDE_no_prior(L0,L1,ML_lambda);
Z = Zmod;
surf(L0,L1,Z,'EdgeColor','none')
view(0,90)
axis equal
xlim([rl0-1*plot_range rl0+plot_range])
ylim([rl1-1*plot_range rl1+plot_range])
xlabel('log(\lambda_0)')
ylabel('log(\lambda_1)')
title('Posterior for \lambda_0, \lambda_1')
subplot(1,2,2)
        [m1,i1] = max(Z);
        [m2,i2] = max(m1);

        best_loc_index = [i2 i1(i2)]; % x-koordinat y-koordinat
        ML = [l0linspace(best_loc_index(1)) l1linspace(best_loc_index(2))];
[c,level] = contour(L0,L1,Z - min(min(Z)),(0.1:0.2:0.9)*max(max(Z - min(min(Z)))));
%clabel(c,level)
axis equal
xlim([rl0-0.75*plot_range rl0+0.25*plot_range])
ylim([rl1-0.5*plot_range rl1+0.5*plot_range])
hold on
scatter(log(facit(1)),log(facit(2)),'red','filled')
scatter(ML(1),ML(2),'blue','filled')
plot(ML_lambda(:,1), ML_lambda(:,2),'b-o')
xlabel('log(\lambda_0)')
ylabel('log(\lambda_1)')
title('Confidence contours of posterior')
sgtitle('Weak Allee effect')
%%
%facit = [1/120000 1/40000 1/600000]; % weak
facit = [1/120000 1/20000 1/100000]; % strong
h = figure('units','centimeters','position',[0 0 16.8 8.5]);
rl0 = round(log(facit(1)));
rl1 = round(log(facit(2)));
plot_range = 3;
l0linspace = linspace(rl0-2*plot_range,rl0+plot_range,resolution);
l1linspace = linspace(rl1-2*plot_range,rl1+plot_range,resolution);
[L0,L1] = meshgrid(l0linspace,l1linspace);
subplot(1,2,1)
%Z = KDE_no_prior(L0,L1,ML_lambda);
Z = Zmod;
surf(L0,L1,Z,'EdgeColor','none')
view(0,90)
axis equal
xlim([rl0-1*plot_range rl0+plot_range])
ylim([rl1-1*plot_range rl1+plot_range])
xlabel('log(\lambda_0)')
ylabel('log(\lambda_1)')
title('Posterior for \lambda_0, \lambda_1')
subplot(1,2,2)
        [m1,i1] = max(Z);
        [m2,i2] = max(m1);

        best_loc_index = [i2 i1(i2)]; % x-koordinat y-koordinat
        ML = [l0linspace(best_loc_index(1)) l1linspace(best_loc_index(2))];
[c,level] = contour(L0,L1,Z - min(min(Z)),(0.1:0.2:0.9)*max(max(Z - min(min(Z)))));
%clabel(c,level)
axis equal
xlim([rl0-0.75*plot_range rl0+0.25*plot_range])
ylim([rl1-0.5*plot_range rl1+0.5*plot_range])
hold on
scatter(log(facit(1)),log(facit(2)),'red','filled')
scatter(ML(1),ML(2),'blue','filled')
plot(ML_lambda(:,1), ML_lambda(:,2),'b-o')
xlabel('log(\lambda_0)')
ylabel('log(\lambda_1)')
title('Confidence contours of posterior')
sgtitle('Strong Allee effect')




function z = KDE_no_prior(x,y,ML_lambda)
z = 0;
S = cov(ML_lambda(:,1), ML_lambda(:,2));
h = S*(1/(length(ML_lambda))^(1/3));
for e = 1:length(ML_lambda)
    mu_e = [ML_lambda(e,1); ML_lambda(e,2)];
    zz = (x-mu_e(1)).^2/max(max(h)) + (y-mu_e(2)).^2/max(max(h));
    z = z + exp(-zz);
end
end
