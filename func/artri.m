function retval = artri(p1,p2,p3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculates area of triangle given by points p1, p2, p3
% according to Heron formula
%
% p1,p2,p3    triangle points (any dimension)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sizp1 = size(p1);
sizp2 = size(p2);
sizp3 = size(p3);

if  sizp1(1)~=sizp2(1) || sizp1(2)~=sizp2(2) ||...
    sizp2(1)~=sizp3(1) || sizp2(2)~=sizp3(2)
  error('stop artri.m: triangle points of different dimension');
end

p1p2 = p1 - p2;
p2p3 = p2 - p3;
p3p1 = p3 - p1;

a = sqrt( p1p2*p1p2' );
b = sqrt( p2p3*p2p3' );
c = sqrt( p3p1*p3p1' );

s = (a+b+c)/2;

retval = sqrt( s*(s-a)*(s-b)*(s-c) );

