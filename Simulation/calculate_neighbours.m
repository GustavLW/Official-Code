function population = calculate_neighbours(population,k)
has_happened = 0; %checks if entire population has been wiped out
if isstruct(population{end}) == false
   appendix = population{end};
   population = population(1:end-1);
   has_happened = 1;
end
N      = length(population);
x      = zeros(2,N);
N_diff = N - size(population{1}.neighbours,2);

% this code assures that the neighbour-matrix is of the correct width,
% and unpacks all cell locations
for i = 1:N
    if size(population{i}.location,2) == 1
    x(:,i) = population{i}.location;
    tmp2   = population{i}.neighbours;
    if sum(size(tmp2))==0 && k > 1
        tmp2 = population{i-1}.neighbours;
        tmp2 = zeros(size(tmp2));
    else
        tmp2   = [tmp2, zeros(1,N_diff)];
    end
    population{i}.neighbours = tmp2;
    else
    x(:,i) = population{i}.location(:,k);
    tmp2   = population{i}.neighbours(k,:);
    if sum(size(tmp2))==0 && k > 1
        tmp2 = population{i-1}.neighbours(k,:);
        tmp2 = zeros(size(tmp2));
    else
        tmp2   = [tmp2, zeros(1,N_diff)];
    end
    population{i}.neighbours(k,:) = tmp2;
    end
end
% this code assures that the neighbour-matrix is of the correct width
% and unpacks all cell locations


% this code calculate which cells are our neighbours at time k
for i = 1:N
    if population{i}.d_time >= k
        tmp = zeros(1,N);
        for j = 1:N
            if j~=i && population{j}.d_time >= k
                if norm(x(:,i)-x(:,j)) < 3
                    tmp(j) = 1;
                end
            end
        end
        if size(population{i}.location,2) == 1
            population{i}.neighbours = tmp;
        else
            disp('happened?')
            population{i}.neighbours(k,:) = tmp;
        end
    end
end

if has_happened == 1
    population{end + 1} = appendix;
end


