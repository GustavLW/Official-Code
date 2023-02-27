clc
clear all
close all

TopFolder  = fileparts(fileparts(pwd));
SimFolder  = [TopFolder '\Simulation'];
DataFolder = [SimFolder '\Datasets'];
%dir(SimFolder)
addpath(([SimFolder, filesep]))
addpath(([fileparts(pwd) '\Diffusion']))
addpath(([fileparts(pwd) '\Allee']))
addpath(([fileparts(pwd) '\Interaction']))

df = dir(DataFolder);
df = df(3:end);

load([DataFolder '\' df(end).name])

cfile = '24jan23.mat';
load(cfile)
T = size(repeated_trials{1,1},1);
Q = size(repeated_trials,1);
A = length(repeated_trials{1,4});
%% get the parameter winners
format long
for q = 1:Q
    A_param = exp(repeated_trials{q,1}(end,:))
end
%%
ttt = [];
for q = 1:Q
tt = log(repeated_trials{q,3}(end,:))';
ttt = [ttt tt];
end

boxplot(ttt,'Notch','on')
xlim([0.5 3.5])
ylim([-5.5 -4])
hold on
plot([0 4],log(observed_cells{end-1}(end))*[1 1],'k--')
grid on
xlabel('\Delta t (minutes)')
xticklabels([10 20 40])
ylabel('$$\log(\hat{\sigma})$$','interpreter','latex')
%% COMPARE WINNING POTENTIAL COMPARED TO UNDERLYING
B_cell  = cell(Q,1);
t = T;
FROM_HERE = 150;
for q = 1:size(repeated_trials)
    B_array = [];
    for a = 1:A
        b_tmp = repeated_trials{q,4}{a}.Bfitness;
        B_array = [B_array;b_tmp (1:length(b_tmp))' a*ones(length(b_tmp),1)];
    end
    B_cell{q} = B_array;
end

for q = 1:Q
    B_array =  B_cell{q};
for b = 1:length(B_array)-1
    if B_array(length(B_array)-b,1) == B_array(length(B_array) - (b-1),1)
        B_array(length(B_array)-b+1,1) = -10^(10);
    end
end
 B_cell{q} = B_array;
end
%plot(B_array(:,1))
LOSERS = 20;
get_the_best = zeros(Q,LOSERS,2);
best_spread  = zeros(Q,LOSERS);
loser_distance = zeros(Q+1,size(get_the_best,2))+10^(-16);
winner_distance = zeros(Q+1,1)+10^(-16);
for q = 1:Q
    x     = linspace(0,5,1001);
    [B,I] = maxk(B_cell{q}(:,1),LOSERS);
    B_param = exp(repeated_trials{q,1}(t,:));
    facit = observed_cells{end-1}(5:end-(length(observed_cells)-2));
    %subplot(3,1,q)
    %histogram(B*2^(q-1)/length(observed_cells{1}.location))
    best_spread(q,:) = B*2^(q-1)/length(observed_cells{1}.location);
    get_the_best(q,:,:) = [B_cell{q}(I,3), B_cell{q}(I,2)];

    mposdiff = max(U_pot(x,B_param)-U_pot(x,facit));
    devidevi = find(U_pot(x,B_param)-U_pot(x,facit)<-2*mposdiff);
    devidevi = [FROM_HERE devidevi];
    FROM_HERE = devidevi(end);


        for b = 1:size(get_the_best,2)
            A_param = exp(repeated_trials{q,4}{get_the_best(q,b,1)}.U_param(get_the_best(q,b,2),:));
            loser_distance(q+1,b) = norm(U_pot(x(FROM_HERE:end),facit)-U_pot(x(FROM_HERE:end),A_param));
        end
        winner_distance(q+1) = norm(U_pot(x(FROM_HERE:end),facit)-U_pot(x(FROM_HERE:end),B_param));
   big_mean = mean(loser_distance');
   big_var = var(loser_distance');     
end
big_mean
big_var
%%

clc


h = figure('units','centimeters','position',[0 0 20*0.8 25*0.8]);
for q = 1:Q
    filename = ['particle_swarm_result_' cfile(1:end-4) '.png'];
    x       = linspace(0,5,1001);
 
    
    for t = T:T
        subplot(3,2,q)
        hold off
        B_param = exp(repeated_trials{q,1}(t,:));
        hold on
        plot(x,U_pot(x,facit),'k','LineWidth',1.5)
        plot(x,U_pot(x,B_param),'b','LineWidth',1.5)
        for b = 1:size(get_the_best,2)
            A_param = exp(repeated_trials{q,4}{get_the_best(q,b,1)}.U_param(get_the_best(q,b,2),:));
            plot(x,U_pot(x,A_param),'magenta')
        end
        plot(x,U_pot(x,facit),'k','LineWidth',1.5)
        plot(x,U_pot(x,B_param),'b','LineWidth',1.5)        

        plot([0.5 4],[0 0],'k')
        grid on
        axis([0.5 3 min(U_pot(x,facit))*[2 -5]])
        xlabel('r')
        ylabel('U(r;\theta)')
        yticks([])
        title([num2str(size(get_the_best,2)) ' best propsed potentials for r_\theta=' num2str(q/4) '.'])
        legend('Underlying potential','Highest scoring potential')
        drawnow;
        
    end
    

end



subplot(3,2,5:6)
rt  = 0:1/4:1;
m1  = big_mean;
m21  = m1 + sqrt(big_var);
m22  = m1 - sqrt(big_var);
hold off
plot(rt,(winner_distance),'b-o','LineWidth',2)
hold on
plot(rt,m1,'m-o','LineWidth',1.5)
plot(rt,m21,'m--')
plot(rt,m22,'m--')

ylim([(1-0.1*sign(min(m22)))*min(m22) 1.1*max(m21)])
xticks(0:1/4:1)
xticklabels([0 0.25 0.5 0.75 1])
xlabel('r_\theta')
ylabel('Error')
legend('Winner',['Average for top ' num2str(LOSERS)],'1 standard deviation','Location','northwest')
grid on
%saveas(h,filename)

%%

%% PERFORMANCE OF FACIT POTENTIAL FOR DIFFERENT TIME DILATIONS.
clc
Q = 12;
obs_array = cell(Q,2);
for q = 1:Q
obs_cop = observed_cells;
        for i = 1:length(observed_cells)-2
            obs_cop{i}.location = observed_cells{i}.location(:,q:q:end);
        end
        obs_cop{end-1}(1) = observed_cells{end-1}(1)*q;
        obs_cop{end}      = observed_cells{end}(q:q:end);
        obs_array{q,1}    = obs_cop;
end

scsc = zeros(Q,1);
tic
for q = 1:Q
    agent = agents{1};
    agent.U_param = log(facit);
    S = 30;
    L = obs_array{q,1}{end-1}(1)/12;
    [agent,~] = score_agent(RD(q:q:end),agent,1,obs_array{q,1},L,S);
    scsc(q) = agent.fitness(1);
    disp(q)
    toc
end
%%
K_dt = floor(144./(1:12))';
plot(K_dt,scsc./K_dt,'k')
xlabel('Number of observations')
ylabel('Score of underlying potential')
axis([min(K_dt) max(K_dt) 0.5 4])
grid on
%%
save('facit_score','scsc')



%% GRAFIK FÖR ANIMERA LÖSNINGAR
clc
close all
q = 3;
agents = repeated_trials{q,4};
best_ever_location = repeated_trials{q,1};

figure('units','centimeters','position',[0 0 25 25]);
filename = 'particle_swarm_result.gif';
c = linspace(1,10,A);
zoom = 3;
for t = 2:T
    current_params = zeros(A,4*pot_type - 2);
    current_best   = current_params;
    for a = 1:A
        current_params(a,:) = agents{a}.U_param(t,:);
        current_best(a,:)   = agents{a}.U_best(t,:);
    end
    if size(agents{a}.U_param(t,:)) < 4
        hold off
        scatter(current_params(:,1),current_params(:,2),[],c,'*')
        hold on
        scatter(current_best(:,1),current_best(:,2),[],c,'o')
        scatter(best_ever_location(t,1),best_ever_location(t,2),25,'red','o')
        scatter(log(facit(1)),log(facit(2)),25,'red','s')
        axis([-15 -0 0 5])
        grid on
    else
        subplot(2,2,1)
        hold off
        scatter(current_params(:,1),current_params(:,4),[],c,'*')
        hold on
        scatter(current_best(:,1),current_best(:,4),[],c,'o')
        plot(best_ever_location(2:t,1),best_ever_location(2:t,4),'r--o')
        scatter(log(facit(1)),log(facit(4)),100,'red','s')
        axis([log(facit(1))-zoom log(facit(1))+zoom log(facit(4))-zoom log(facit(4))+zoom])
        grid on
        xlabel('k_1')
        ylabel('k_2')
        subplot(2,2,2)
        hold off
        scatter(current_params(:,2),current_params(:,5),[],c,'*')
        hold on
        scatter(current_best(:,2),current_best(:,5),[],c,'o')
        plot(best_ever_location(2:t,2),best_ever_location(2:t,5),'r--o')
        scatter(log(facit(2)),log(facit(5)),100,'red','s')
        axis([log(facit(2))-zoom log(facit(2))+zoom log(facit(5))-zoom log(facit(5))+zoom])
        grid on
        xlabel('l_1')
        ylabel('l_2')
        subplot(2,2,3)
        hold off
        scatter(current_params(:,3),current_params(:,6),[],c,'*')
        hold on
        scatter(current_best(:,3),current_best(:,6),[],c,'o')
        plot(best_ever_location(2:t,3),best_ever_location(2:t,6),'r--o')
        scatter(log(facit(3)),log(facit(6)),100,'red','s')
        axis([log(facit(3))-zoom log(facit(3))+zoom log(facit(6))-zoom log(facit(6))+zoom])
        grid on
        xlabel('a_1')
        ylabel('a_2')
        subplot(2,2,4)
        hold off
        x       = linspace(0,5,1001);
        xs      = linspace(0,5,21);
        facit = observed_cells{end-1}(5:end-(length(observed_cells)-2));
        B_param = exp(repeated_trials{q,1}(t,:));
        plot(xs,U_pot(xs,facit),'ro')
        hold on
        plot(xs+0.125,U_pot(xs+0.125,B_param),'bd')
        plot(x,U_pot(x,facit),'r')
        plot(x,U_pot(x,B_param),'b')
        plot([0.5 4],[0 0],'k')
        grid on
        axis([0.5 4 min(U_pot(x,facit))*[2 -5]])
        xlabel('Distance')
        ylabel('Potential energy')
    end
    drawnow;
    frame = getframe(1);
      im = frame2im(frame);
      [imind,cm] = rgb2ind(im,256);
    if t == 1
         imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append');

    end
end