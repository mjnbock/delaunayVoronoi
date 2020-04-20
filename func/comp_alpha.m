function retval = comp_alpha( t )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% comp_alpha.m
% computes alpha anisotropy parameter
% interaction==3: alpha indicates dcrit anisotropy
%
%
% t  simulation time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% make sure that the values here pass sanity check in initialize.m
if t > 8*3600  &&  t < 16*3600
  alpha = 1.4;
else
  alpha = 1.0;
end

retval = alpha;

