function new_cell = create_cell(time,parent,neighbours,location,local_density,birth_evolution)

new_cell = struct('b_time',time,'d_time',10^128,...                          % times
                  'parent',parent,'children',[],'neighbours',neighbours,...  % family structure
                  'location',location,'local_density',local_density,...      % dynamics, spatial
                  'birth_evolution',birth_evolution,'birth_barrier',rand,... % dynamics, birth
                  'death_evolution',0,'death_barrier',rand);   % dynamics, death