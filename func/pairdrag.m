function retval = pairdrag( fvis, dgam, mydensities );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pairdrag.m
% computes viscous interaction drag coefficient
% pairdrag(i,j) = fvis * \int_{arc} rhopair*cdpair \ddd gamma
%
% fvis     strength of viscous interaction (per one filament-cadherin pair)
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

% computing drag
retval = fvis * sum(rhopair.*cdpair) * dgam;

