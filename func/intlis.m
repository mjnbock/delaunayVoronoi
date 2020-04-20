function retval = intlis( v,u, M,R, mindist )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% intlis.m
%%%%%%%%
%%%%%%%% intersect line w/ sphere in 2d
%%%%%%%% M center sphere (row vector [x,y]), R radius sphere
%%%%%%%% v starting point straight line (row vector)
%%%%%%%% u normalized direction of straight line (row vector)
%%%%%%%% mindist: minimal distance (when to consider two points equal)
%%%%%%%%
%%%%%%%% Remark:
%%%%%%%% solving
%%%%%%%% (x-M)^2 = R^2
%%%%%%%% x = v + g*u    w/ g real running parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vM = v - M;

% discriminant is 4*(b^2-c)
b = u*vM';
c = vM*vM' - R^2;

D = b^2 - c;

if abs(D) < mindist
  % D approx. 0
  g = -b;
  x = v + g*u;
  retval = [ x(1), x(2) ];
elseif D > 0
  g = -b;
  g1 = g + sqrt(D);
  g2 = g - sqrt(D);
  x1 = v + g1*u;
  x2 = v + g2*u;
  retval = [ x1(1),x1(2); x2(1),x2(2) ];
else
  retval = [];
end

return;

