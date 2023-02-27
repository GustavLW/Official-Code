clc
clear all
close all

load extracted_tracks.mat

N = size(ECAET,2);
K = size(ECAET,3)-2;
%%
observed_cells = cell(N+2,1);
observed_cells{N+1} = 20*60;
for i = 1:N
    observed_cells{i} = struct('b_time',[],'d_time',[],'parent',[],'children',[],'neighbours',[],'location',[]);
end
observed_cells{N+2} = cell(K,2);

for i = 1:N
    loc_tmp = squeeze(ECAET(:,i,:));
    s_tmp   = find(isnan(loc_tmp(1,:)) == 0);
    observed_cells{i}.location = loc_tmp;
    observed_cells{i}.b_time = s_tmp(1);
    observed_cells{i}.d_time = s_tmp(end);
end

for k = 1:K
    pic_tmp = squeeze(ECAET(:,:,k));
    d_tmp   = zeros(N);
    alive   = zeros(1,N);
    for i = 1:N-1
        if ~isnan(pic_tmp(1,i))
            alive(i) = 1;
        end
        for j = i+1:N 
            d_tmp(i,j) = sign(max(0, 3 - norm(pic_tmp(:,i)-pic_tmp(:,j)) ));
            d_tmp(j,i) = d_tmp(i,j);
        end
    end
    if ~isnan(pic_tmp(1,N))
        alive(N) = 1;
    end
    d = sparse(d_tmp);
    observed_cells{N+2}{k,1} = d;
    observed_cells{N+2}{k,2} = alive;
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
filename  = strcat(c1,c2,c3,'_',c4,c5,'_dataset_real','.mat');
save(strcat(filename),'observed_cells')