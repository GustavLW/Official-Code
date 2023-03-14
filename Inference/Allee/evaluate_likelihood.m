function likel = evaluate_likelihood(rho,beta,alive,dt,l0,l1)
% takes in one lambda0 and one lambda1 and return that pair's likelihood

N = size(rho,1);
K = size(rho,2);
likel  = 0;
for k = 1:K
    a_tmp = find(alive(:,k)==1);
    for i = a_tmp(1):a_tmp(end)
        if sum(rho(i,k:K)) ~= 0 && sum(rho(i,1:k)) ~=0 && rho(i,k) < 1
            hik = dt*(1 - rho(i,k))*(l0+l1*rho(i,k));
            likel = likel + beta(i,k)*log(hik) - hik - log(factorial(beta(i,k)));
        end
    end
end
