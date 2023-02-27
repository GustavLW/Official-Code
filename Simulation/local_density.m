function rho = local_density(population,i,N,N_old,a)
rho = 0;
xi  = population{i}.location;
nei = population{i}.neighbours;

for j = 1:min(N,N_old)
    if nei(j) == 1
        xj  = population{j}.location;
        rho = rho + exp(-a*norm(xi-xj))*a;
    end
end

rho = min(rho,1);