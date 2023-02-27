function F = interactions(population,U_param,neiA,aliA)
switch nargin
    case 2
        % for data generation simulation
        N_nei = size(population{1}.neighbours,2);
        N     = length(population);
        x     = zeros(2,N);
        F     = zeros(2,N);
        for i = 1:N
            x(:,i) = population{i}.location;
        end

        for i = 1:N
            xi  = x(:,i);
            nei = population{i}.neighbours;
            if sum(size(nei)) > 0
                for j = 1:min(N,N_nei)
                    if nei(j) == 1
                        xj     = x(:,j);
                        F(:,i) = F(:,i) + pair_potential(xi,xj,U_param);
                    end
                end
            end
        end
    case 4
        % for inference simulation
        % population should be a 2xN matrix, neiA a symmetric NxN matrix
        % stage 1 testing cleared pepelaugh
        N     = length(population);
        x     = population;
        F     = 0*x;
        for n = 1:length(aliA)
            i = aliA(n);
            nei_tmp = neiA(i,:);
            nei_tmp = find(nei_tmp>0);
            xi      = x(:,i);
            if and(~isempty(nei_tmp),~isnan(xi(1)))
                for nn = 1:length(nei_tmp)
                    j  = nei_tmp(nn);
                    xj = x(:,j);
                    F(:,i) = F(:,i) + pair_potential(xi,xj,U_param);
                end
            end
        end
end
