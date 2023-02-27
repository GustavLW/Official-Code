
clc
clear all
close all
%rng(921111)
dset = 1;
n1m  = [1 1];   % governs initial cell number
n2m  = [3 3];   % governs initial density
n3m  = [1 1];   % governs strength of allee effect (1 = no proliferation)
n4m  = [1 1];   % governs bimodality
n5m  = [2 2];   % governs interaction parameters
for n1 = n1m(1):n1m(2)                              
    for n2 = n2m(1):n2m(2)                          
        for n3 = n3m(1):n3m(2)                      
            for n4 = n4m(1):n4m(2)                  
                for n5 = n5m(1):n5m(2)              
                    dt    = 1;                      % one time unit is one hour; one time step is a second
                    freq  = 600.0;                 % frequency of observations in seconds
                    h     = 24;                     % total number of hours the simulation will be run
                    Kobs  = h*3600/freq;            % number of of observations
                    K     = round(h*3600);          % total number of discrete time steps [seconds]
                    scale = 50;                     % one length unit is one average cell radius. domain is scale X scale
                    N     = 32*2^n1;                % size of initial population
                    Vmin  = 0.00021;                % depth of potential
                    alpha = 3.50000;                % steepness of potential
                    s_bas = exp(-4.5);              % diffusion coefficients
                    sigma = s_bas*ones(1,N);  
                    if n4 == 2
                        sigma(4:4:end) = s_bas*exp(-2);
                    end
                    l0 = [0.000 0.020 0.250]/3600;        % base division rate
                    l1 = [0.000 0.200 2.500]/3600;        % allee parameter
                    om = [0.000 0.003 0.600]/3600;        % base death rate
                    lambda0 = l0(n3);
                    lambda1 = l1(n3);
                    omega   = om(n3);

                    population = cell(N,1);
                    observed_cells = cell(N,1);
                    Omega = scale*[0 0;1 0; 1 1; 0 1]';
                    for i = 1:N
                        population{i}     = create_cell(0,[],[],scale/2+scale*(1+n2)*(2*rand(2,1)-1)/10,0,0);
                        observed_cells{i} = create_cell(0,[],[],zeros(2,Kobs),[],[]);
                    end
                    observed_neighbours = cell(Kobs,1);
                    population = calculate_neighbours(population,1);
                    tagged_particle = population{1};
                    tagged_particle.location = zeros(2,K);
                    tagged_particle.local_density = zeros(1,K);
                    U_param = [Vmin alpha];
                    if n5 == 2
                        % ALTERNATIVE PARAMETRISATION - a fucking mess
                        a1 = 0.004;
                        a2 = 0.00006;
                        k1 = 10;
                        k2 = 4;
                        l1 = 0.55;
                        l2 = 1.2;
                        U_param = [k1 l1 a1 k2 l2 a2];
                    end
                    tic
                    N_old = N;
                    tot_dead = 0;
                    %%% %%% BEGIN SIMULATION %%% %%%
                    for k = 1:K
                        %%% %%% photograph the situation %%% %%%
                        if mod(k,freq) == 0
                            obs = k/(freq);
                            for i = 1:N
                                observed_cells{i}.b_time          = ceil(population{i}.b_time/(freq));
                                observed_cells{i}.d_time          = ceil(population{i}.d_time/(freq));
                                observed_cells{i}.parent          = population{i}.parent;
                                observed_cells{i}.children        = population{i}.children;
                                observed_cells{i}.location(:,obs) = population{i}.location;
                                observed_neighbours{obs}          = [observed_neighbours{obs};population{i}.neighbours];
                            end
                        end

                        % add cells loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        if (lambda0 + lambda1 + omega) > 0
                            newcomers   = cell(0);
                            newcomers_o = cell(0);
                            exiters     = cell(0);
                            newindex    = N + 1;
                            any_bd   = 0;
                            any_b    = 0;
                            for i = 1:min(N,N_old)
                                tmp      = population{i}.birth_evolution;
                                tmp2     = population{i}.birth_barrier;
                                tmpD     = population{i}.death_evolution;
                                tmpD2    = population{i}.death_barrier;
                                % compute local density %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                loc_dens = local_density(population,i,N,N_old,3.068891070360678);
                                population{i}.local_density = loc_dens;
                                % update internal states
                                tmp  = tmp + dt*(1-loc_dens)*(lambda0+lambda1*loc_dens)*(1-tmp);
                                tmpD = tmpD + dt*omega*(1-tmpD);
                                % check if birth will occur %%%%%%%%%%%%%%%%
                                if tmp2 < tmp
                                    any_bd = any_bd + 1;
                                    any_b  = any_b + 1;
                                    crand     = 2*pi*rand;
                                    placetmp  = [cos(crand);sin(crand)]*(1 + 0.2*abs(randn));
                                    newcomers = [newcomers;...
                                        create_cell(k,i,[],population{i}.location + placetmp,0,0)];
                                    tmp3      = population{i}.children;
                                    tmp3      = [tmp3 newindex];
                                    newindex  = newindex + 1;
                                    population{i}.children        = tmp3;
                                    population{i}.birth_evolution = 0;
                                    population{i}.birth_barrier   = rand;

                                    %%% newcomer copy for observation loop
                                    newcomers_o = [newcomers_o;...
                                        create_cell(0,[],[],(population{i}.location + placetmp).*ones(2,Kobs),[],[])];
                                    %%% newcomer copy for observation loop
                                else
                                    population{i}.birth_evolution = tmp;
                                end
                                % check if death will occur %%%%%%%%%%%%%%%%
                                if tmpD2 < tmpD
                                    any_bd                        = any_bd + 1;
                                    tot_dead                      = tot_dead + 1;
                                    population{i}.d_time          = k;
                                    population{i}.location        = 10^30*[1;1];
                                    population{i}.birth_evolution = -10^20;
                                    population{i}.birth_barrier   = 10^20;
                                    population{i}.death_evolution = -10^20;
                                else
                                    population{i}.death_evolution = tmpD;
                                end
                            end
                            population = [population; newcomers];
                            observed_cells = [observed_cells; newcomers_o];
                            if any_bd > 0 % recalculate neighbours if anyone got born and shit
                                sig_prob = 1 - (n4-1)/2;
                                n_sig = sigma(1)*(0.1 + 0.9*binornd(1,sig_prob,1,any_b));
                                population = calculate_neighbours(population,k);
                                N_old = N;
                                sigma = [sigma n_sig];
                            end
                        end

                        % propagate loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        N = length(population); % to account for newcomers
                        if Vmin > 0 % if no interaction takes place, let's skip this part.
                            F = interactions(population,U_param);
                        else
                            F = zeros(2,N);
                        end 
                        for i = 1:N
                            if population{i}.d_time > k
                                population{i}.location(1) = max(min(population{i}.location(1)...
                                    + dt*F(1,i) + sqrt(dt)*sigma(i)*randn(1,1),scale),0);
                                population{i}.location(2) = max(min(population{i}.location(2)...
                                    + dt*F(2,i) + sqrt(dt)*sigma(i)*randn(1,1),scale),0);
                                % addition of the randomness at this stage is ok since the
                                % interactions have already been calculated.
                            end
                        end
                        if mod(k,3600) == 0
                            disp(['Hour: ' num2str(k/3600) '. Dataset: ' num2str(dset) '/' num2str(span(n1m)*span(n2m)*span(n3m)*span(n4m)*span(n5m)) '. Current population is: ' num2str(N-tot_dead) '.'])
                            %%% RECALCULATE NEIGHBOURS HERE
                            population = calculate_neighbours(population,k);
                            N_old      = N;
                            %%% RECALCULATE NEIGHBOURS HERE
                            toc
                        end

                        tagged_particle.location(:,k)    = population{1}.location;
                        tagged_particle.local_density(k) = population{1}.local_density;
                    end
                    %%% %%% %%% SIMULATION DONE   %%% %%% %%%
                    %%% %%% %%% SAVE OBSERVATIONS %%% %%% %%%
                    for i = 1:N
                        observed_cells{i} = rmfield(observed_cells{i},{'local_density','birth_evolution','birth_barrier','death_evolution','death_barrier'});
                    end
                    observed_cells{i+1} = [freq [lambda0 lambda1 omega] U_param sigma];
                    fin_pop = size(observed_neighbours{end},1);

                    observed_cells{i+2} = cell(Kobs,2);
                    for k = 1:Kobs
                        alive   = zeros(1,N);
                        for i = 1:N
                            if and(k >= observed_cells{i}.b_time,k < observed_cells{i}.d_time)
                                alive(i) = 1;
                            else
                                observed_cells{i}.location(:,k) = [NaN;NaN];
                            end
                        end
                        tmp_pop                      = observed_neighbours{k};
                        cur_pop                      = size(tmp_pop,1);
                        tmp_fin                      = zeros(fin_pop);
                        tmp_fin(1:cur_pop,1:cur_pop) = tmp_pop;
                        observed_neighbours{k}       = sparse(tmp_fin);
                        observed_cells{i+2}{k,1} = sparse(tmp_fin);
                        observed_cells{i+2}{k,2} = alive;
                    end

                    

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
                    filename  = strcat(c1,c2,c3,'_',c4,c5,'_dataset','.mat');
                    save(strcat('Datasets\',filename),'observed_cells')
                    % observed_cells is the object we shall save
                    dset = dset + 1;
                end
            end
        end
    end
end

function s = span(v)
    s = v(2) - v(1) + 1;
end
