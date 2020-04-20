function retval = vFloc( r,rho,bi, Rmax,vph, floc,flocnum )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% vFloc.m
% computes locomotion force
%
%
% r        radius of cell body
% rho      filament density of cell
% bi       relative amount of bound integrin
%
% Rmax     Pmax-radius of cell
% vph      [phm,php] phi range of marginal arc
%
%
% floc     strength of locomotion force
% flocnum  number of points to approximate marginal arc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% computing approximation
php = vph(2);
phm = vph(1);
if php < phm
  %warning('vFloc.m: unexpected angle range');
  error('stop vFloc.m: unexpected angle range');
end
deltph = php - phm;
dph = deltph/flocnum;
phi = ( (phm+dph/2):dph:(php-dph/2) )';


% summing up
scosphi = sum( cos(phi) );
ssinphi = sum( sin(phi) );


% prefactors and return
%Flocx = floc * (rho*r/Rmax) * scosphi * Rmax*dph;
%Flocy = floc * (rho*r/Rmax) * ssinphi * Rmax*dph;
Flocx = floc* bi*rho*r * scosphi * dph; % TODO think about areal
Flocy = floc* bi*rho*r * ssinphi * dph; % influence of integrin

retval = [ Flocx, Flocy ];

