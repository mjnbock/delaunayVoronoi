function retval = intlip( v1,u1, v2,u2, mindist )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% intlip.m
%%%%%%%%
%%%%%%%% intersect two straight lines in 2d
%%%%%%%% v1 starting point straight line 1 (row vector)
%%%%%%%% v2 starting point straight line 2 (row vector)
%%%%%%%% u1 (normalized) direction straight line 1 (row vector)
%%%%%%%% u2 (normalized) direction straight line 2 (row vector)
%%%%%%%% mindist: minimal distance (when to consider two points equal)
%%%%%%%%
%%%%%%%% Remark:
%%%%%%%% straight lines constructed via
%%%%%%%% x1 = v1 + lmd1*u1
%%%%%%%% x2 = v2 + lmd2*u2
%%%%%%%% and solving
%%%%%%%% /  |   |  \ / lmd1 \   / (v1-v2)(1) \
%%%%%%%% | u1  -u2 | |      | = |            | = / a b \ / lmd1 \ = / c \
%%%%%%%% \  |   |  / \ lmd2 /   \ (v1-v2)(2) /   \ d e / \ lmd2 /   \ f /    
%%%%%%%% operator \ or so not used: no certainty about # solutions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

w = v2 - v1;

a =  u1(1);
b = -u2(1);
c =   w(1);
d =  u1(2);
e = -u2(2);
f =   w(2);

% consider determinants
d0 = a*e - d*b;
d1 = a*f - d*c;
d2 = c*e - f*b;
if abs(d0) < mindist
  if abs(d1) < mindist  &&  abs(d2) < mindist
    retval = ones(2,2);
  else
    retval = [];
  end
  return;
end

lmd1 = d2/d0;
lmd2 = d1/d0;

pt1 = v1 + lmd1*u1;
pt2 = v2 + lmd2*u2;

vdist = pt1 - pt2;
dist = sqrt( vdist*vdist' );

if dist < mindist
  retval = pt1;
else
  error('stop intlip.m: solutions not coinciding');
end

