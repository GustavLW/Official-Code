function pop_tmp = sample_from_next_state(pop_tmp,nei_tmp,ali_tmp,U_param,sigma,freq,L)
% has passed stage 1 testing yee always wins
% input:    population, neighbours (matrix), theta, sigma
% freq:     seconds between observations
% K:        number of simulation steps
dt          = freq/L;
N           = size(pop_tmp,2);

for k = 1:L-1
    F = interactions(pop_tmp,U_param,nei_tmp,ali_tmp);
    for n = 1:length(ali_tmp)
        i = ali_tmp(n);
        pop_tmp(:,i) = pop_tmp(:,i) + dt*F(:,i) + sqrt(dt)*sigma(i)*[randn;randn];
    end
end


