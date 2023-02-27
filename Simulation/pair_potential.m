function f = pair_potential(xi,xj,U_param) 
% calculate the influence of cell j on cell i 

r    = norm(xi-xj);
if r < 3 % extremely arbitrary
    f = -V(xi,xj,U_param);  
else
    f = [0;0];
end

