function pop_tmp = unpack_image(observed_cells,k)
ali_tmp   = find(observed_cells{end}{k,2} == 1);
N = length(observed_cells) - 2;
pop_tmp = NaN*ones(2,N);

for n = 1:length(ali_tmp)
   i = ali_tmp(n);
   pop_tmp(:,i) = observed_cells{i}.location(:,k); 
end

