function retval = vFpol( r,rho,bi,pol, fpol,fpolnum,kappol )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% vFpol.m
% computes polarity force
%
%
% r        radius of cell body
% rho      filament density of cell
% bi       relative amount of bound integrin
% pol      cell polarity vector
%
% fpol     strength of polarity force
% fpolnum  number of points to approximate cell circumference
% kappol   strength of polarity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% computing approximation
dph = 2*pi/fpolnum;
phi = ( 0:dph:(2*pi+dph/2) )';


% radial unit vectors and polarity direction
cosphi = cos(phi);
sinphi = sin(phi);

phip = phi - atan2( pol(2),pol(1) );
ppol = sqrt( pol(1)^2 + pol(2)^2 );


% polarity mode
polscal = ppol * ( (1-kappol) + kappol*cos(phip) );


% sum contributions
scosphi = sum( polscal .* cosphi );
ssinphi = sum( polscal .* sinphi );


% prefactors and return
Fpolx = fpol* bi*rho*r * scosphi * dph; % TODO think about areal
Fpoly = fpol* bi*rho*r * ssinphi * dph; % influence of integrin

retval = [ Fpolx, Fpoly ];

