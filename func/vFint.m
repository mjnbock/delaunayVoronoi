function retval = vFint( delta,deltcri,deltmin, rb,rs, dgam,finth,fintv, ...
                         uh,uv, mindist, mydensities );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% vFint.m
% computes interaction force according to B158,B159
% for spherical cell-cell contacts; acting on smaller cell!!!
%
%
% delta    distance of cell body surfaces 
% deltcri  optimal distance of cell surfaces
% deltmin  minimal distance of cell surfaces
%
% mindist  minimal floating point distance
%
% cosphib  cos(phi) filament bigger cell               \
% sinphib  sin(phi) filament bigger cell               |
% cosphis  cos(phi) filament smaller cell              | 
% sinphis  sin(phi) filament smaller cell               > mydensities(n,:)
%                                                      |
% rhopair  filament pairing density on neighboring arc |
% cdpari   cadherin  paring density on neighboring arc /
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% setup densities
cosphib = mydensities(:,1);
sinphib = mydensities(:,2);
cosphis = mydensities(:,3);
sinphis = mydensities(:,4);
rhopair = mydensities(:,5);
cdpair  = mydensities(:,6);

% compute log scaling term, (Pmax-1) scaling done in 
% initialize.m for performance
ddeltcri = deltcri * (rb+rs);
ddeltmin = deltmin * (rb+rs);
if delta < ddeltmin
  warning('vFint.m: cells closer than deltmin, working arround');
  %error('stop vFint.m: cells closer than deltmin, decrease dt');
  llog = log(          mindist/(ddeltcri-ddeltmin) );
else
  llog = log( (delta-ddeltmin)/(ddeltcri-ddeltmin) );
end


% integrate force density
% TODO think about integrin vertical component
FFinth = sum( rhopair.*(cosphib+cosphis)/2 )*dgam*mean(cdpair);
FFintv = sum( rhopair.*(sinphib+sinphis)/2 )*dgam*mean(cdpair);


% return
FFinth = FFinth*finth*llog;
FFintv = FFintv*fintv*llog;

Fintxs = [   FFinth*uh(1), FFintv*uv(1) ];
Fintys = [   FFinth*uh(2), FFintv*uv(2) ];
Fintxb = [ - FFinth*uh(1), FFintv*uv(1) ];
Fintyb = [ - FFinth*uh(2), FFintv*uv(2) ];

retval = [ sum(Fintxs), sum(Fintys); ... 
           sum(Fintxb), sum(Fintyb); ...
            Fintxs(1) ,  Fintys(1) ; ...
            Fintxs(2) ,  Fintys(2) ; ...
            Fintxb(1) ,  Fintyb(1) ; ...
            Fintxb(2) ,  Fintyb(2)       ];

