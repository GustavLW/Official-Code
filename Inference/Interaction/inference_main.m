function inference_main(are_we_local,EXPERIMENT_NR,observed_cells)



time_dilate = 0;
A = 1*feature('numcores');  % number of optimization agents
Q = 3;
T = 4;  % number of generations
S = 1;
pot_type = 2;
RD = calculate_radial_distribution(observed_cells,1);
obs_array = cell(Q,2);

if are_we_local == 1
    if time_dilate ==1
        facit = observed_cells{end-1}(5:end-(length(observed_cells)-2));
        for q = 1:Q
            obs_cop = observed_cells;
            for i = 1:length(observed_cells)-2
                obs_cop{i}.location = observed_cells{i}.location(:,2^(q-1):2^(q-1):end);
            end
            obs_cop{end-1}(1) = observed_cells{end-1}(1)*2^(q-1);
            obs_cop{end}      = observed_cells{end}(2^(q-1):2^(q-1):end);
            obs_array{q,1}    = obs_cop;
            s1 = zeros(A,length(facit));
            for a = 1:A
                s1(a,:) = log(facit) +  rand_sphere_shell(length(facit),0.75);
            end
            obs_array{q,2} = s1;
        end
    else
        try
            facit = observed_cells{end-1}(5:end-(length(observed_cells)-2));
            s1 = zeros(A,Q,length(facit));
            for q = 1:Q
                for a = 1:A
                    s1(a,q,:) = log(facit) +  rand_sphere_shell(length(facit),q/4);
                end
                obs_array{q,1} = observed_cells;
                obs_array{q,2} = squeeze(s1(:,q,:));
            end
        catch
            facit = [10 0.55 4*10^(-3) 4 1.2 6*10^(-5)];
            para_ranges = log(facit') + 1*[-1 1; -1 1; -1 1; -1 1; -1 1; -1 1];
            s0 = hyper_cube_sampling(Q,para_ranges);
            s1 = zeros(A,Q,length(facit));
            for q = 1:Q
                for a = 1:A
                    s1(a,q,:) = s0(q,:) + rand_sphere_shell(length(facit),0.5);
                end
                obs_array{q,1} = observed_cells;
                obs_array{q,2} = squeeze(s1(:,q,:));
            end
        end
    end
elseif are_we_local == 0
    Q           = 1;
    obs_array   = cell(Q,2);
    facit       = [10 0.55 4*10^(-3) 4 1.2 6*10^(-5)];
    para_ranges = log(facit') + 1*[-1 1; -1 1; -1 1; -1 1; -1 1; -1 1];
    s0 = hyper_cube_sampling(64,para_ranges);
    s1 = zeros(A,64,length(facit));
    for q = 1:64
        for a = 1:A
            s1(a,q,:) = s0(q,:) + rand_sphere_shell(length(facit),0.5);
        end
    end
    obs_array{Q,1} = observed_cells;
    obs_array{Q,2} = squeeze(s1(:,EXPERIMENT_NR,:));
    disp('test')
end
repeated_trials = cell(Q,4); % location evolution, fitness, sigma for best location, all agents for good measure

clc
tic
for q = 1:Q
    [best_ever_location,best_ever_fitness,best_ever_sigma,agents] =...
    particle_swarm_optimization(q,A,pot_type,T,facit,obs_array{q,2},RD,obs_array{q,1},obs_array{q,1}{end-1}(1)/12,S);
    repeated_trials{q,1} = best_ever_location;
    repeated_trials{q,2} = best_ever_fitness;
    repeated_trials{q,3} = best_ever_sigma;
    repeated_trials{q,4} = agents;
    disp(['Repeated trial ' num2str(q) '/' num2str(Q) ' completed.'])
    toc
end

if are_we_local == 1
    c = clock;
    c1 = num2str(c(1)); % år
    c2 = num2str(c(2)); % månad
    tmp = length(c2);
    if tmp == 1
        c2 = strcat('0',c2);
    end
    c3 = num2str(c(3)); % dag
    tmp = length(c3);
    if tmp == 1
        c3 = strcat('0',c3);
    end
    c4 = num2str(c(4)); % timme
    tmp = length(c4);
    if tmp == 1
        c4 = strcat('0',c4);
    end
    c5 = num2str(c(5)); % sekund
    tmp = length(c5);
    if tmp == 1
        c5 = strcat('0',c5);
    end
    filename  = strcat(c1,c2,c3,'_',c4,c5,'_result','.mat');
elseif are_we_local == 0
    filename  = strcat('real_',num2str(EXPERIMENT_NR),'_result','.mat');  
end
save(strcat('Results/',filename),'repeated_trials')

function point = rand_sphere_shell(n,r_radii)
rang = [pi*rand(1,n-2) 2*pi*rand];
point    = ones(1,n);
point(2) = sin(rang(1));
point(end) = prod(sin(rang));
if n > 3
    for i = 3:n-1
        point(i) = point(i-1).*sin(rang(i-1));
    end
end
point(1:end-1) = point(1:end-1).*(cos(rang(1:end)));
point = r_radii*point;
end

end






