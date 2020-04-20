function retval = arcdens( vxb,vxs, Rij,vMij,vth,...
                           rhob,rhos,rb,rs,cdb,cds,...
                           uh,uv, myfintnum,mindist )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arcdens.m
% computes neighboring arc filament densities according to B158,B159
% for spherical cell-cell contacts; acting on smaller cell!!!
%
%
% vxb      position vector center bigger  cell
% vxs      position vector center smaller cell
%
% Rij      radius contact sphere
% vMij     position vector center contact sphere / start vertex postion
% vthm     start and stop angle contact arc / stop vertex position
%
% rhob     filament density bigger  cell
% rhos     filament density smaller cell
% rb       radius bigger  cell body
% rs       radius smaller cell body
% cdb      relative amount cadherin bigger cell
% cds      relative amount cadherin smaller cell
%
% mindist  equality threshold distance
% fintnum  number of points to approximate contact arc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% computing approximation
%arcapx1 = ceil( fintnum * deltth/2/pi ); % works only for arcs
%arcapx2 = ceil( fintnum / 10          ); % scale for straight line unclear
%arcapx  = max([arcapx1,arcapx2]);
arcapx = myfintnum;
if isnan(Rij)
  v1 = vMij;
  v2 = vth;
  vdeltv = v2 - v1;
  dgamx = vdeltv(1)/arcapx;
  dgamy = vdeltv(2)/arcapx;
  if abs( vdeltv(1) ) < mindist % A approx in global coords
    Ax = ones(arcapx,1)*v1(1);
  else
    Ax = ( (v1(1)+dgamx/2):dgamx:(v2(1)-dgamx/3) )'; % -dgx/3 for fixed size(Ax)
  end
  if abs( vdeltv(2) ) < mindist
    Ay = ones(arcapx,1)*v1(2);
  else
    Ay = ( (v1(2)+dgamy/2):dgamy:(v2(2)-dgamy/3) )'; % -dgx/3 for fixed size(Ay)
  end
  costhet =  ones(arcapx,1); % in local coords
  sinthet = zeros(arcapx,1);
else
  thp = vth(2);
  thm = vth(1);
  if thp < thm
    error('stop arcdens.m: unexpected angle range');
  end
  deltth = thp - thm;
  dth = deltth/(arcapx);
  theta = ( (thm+dth/2):dth:(thp-dth/2) )';
  %dgam = Rij*dth;
  uAx = - cos(theta);
  uAy = - sin(theta);
  Ax = vMij(1) - Rij*uAx; % A approx in global coords
  Ay = vMij(2) - Rij*uAy;
  costhet = - (  uAx*uh(1) +  uAy*uh(2) ); % as in cosphib: > 0, in local coords
  sinthet = - (  uAx*uv(1) +  uAy*uv(2) );
end


% find filament directions
Rbx = vxb(1) - Ax;
Rby = vxb(2) - Ay;
Rsx = vxs(1) - Ax;
Rsy = vxs(2) - Ay;
absRb = sqrt( Rbx.*Rbx + Rby.*Rby );
absRs = sqrt( Rsx.*Rsx + Rsy.*Rsy );
uRbx = Rbx ./ absRb;
uRby = Rby ./ absRb;
uRsx = Rsx ./ absRs;
uRsy = Rsy ./ absRs;


% finding trigonometric functions in LOCAL coordinates
cosphib =   ( uRbx*uh(1) + uRby*uh(2) );
sinphib = - ( uRbx*uv(1) + uRby*uv(2) );

cosphis = - ( uRsx*uh(1) + uRsy*uh(2) ); % as in cosphib: > 0
sinphis = - ( uRsx*uv(1) + uRsy*uv(2) );


% find filament and cadherin density
% in local coords we have
%dphibdth = Rij*cos(theta+phib)./Rb;
%dphisdth = Rij*cos(theta+phis)./Rs;
%rhobb = rb*rhob*dphibdth;
%rhoss = rs*rhos*dphisdth;
%rho = sqrt( rhobb .* rhoss )/Rij;
% which can be rewritten to
dphibdth = ( cosphib.*costhet - sinphib.*sinthet )./absRb;
dphisdth = ( cosphis.*costhet + sinphis.*sinthet )./absRs;
rhobb = rb*rhob*dphibdth;
rhoss = rs*rhos*dphisdth;
rhopair = sqrt( rhobb .* rhoss );
cdpair  = sqrt(   cdb  * cds   ) * ones(length(rhopair),1);

% return densities in columns
retval = [ cosphib, sinphib, ...
           cosphis, sinphis, ...
           rhopair, cdpair   ];

