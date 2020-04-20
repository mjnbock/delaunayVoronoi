%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% number of cells
numcells = 111;
dillute = 0.0;         % Fraction of cells to be removed (dilluted) 
                       % after creation.
                       % Mind that with lattices of any type, a bunch of 
                       % cells may be removed to arrive at the desired 
                       % numcells.


% tissue parameters
Pmax   = 3;
Delcri = 0.5;
r0     = 2.5;


% save which tissues
save_full_tissue = 1;  % full tissue, numcells2 > numcells
save_base_tissue = 1;  % basic tissue, numcells
save_dill_tissue = 1;  % dilluted tissue, < numcells


% lattice settings
%lattice = 'square';
lattice = 'square';
spacing = 2*r0*Pmax*Delcri;


% positions
x_pspd = 0.00*spacing; % Poisson prefactor
x_gspd = 0.05*spacing; % Gaussian prefactor

y_pspd = 0.00*spacing;
y_gspd = 0.05*spacing;


% radii and weights
r_pspd  = 0.00*r0;
r_gspd  = 0.15*r0; % 0.08 ... 0.12

%w0        = r0;
%w_pspd   = r_spd;
%w_gpread  = r_spd;


% cell properties, discriminated per type
frac = [ 1.00, 0.00, 0.00, 0.00, 0.00, 0.00 ]; % percentage of cell type

Tp0     = [  120,  120,  120,  120,  120,  120 ];
Tp_pspd = [ 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 ].*Tp0;
Tp_gspd = [ 0.10, 0.10, 0.10, 0.10, 0.10, 0.10 ].*Tp0;

a0      = [ 8,    8,    8,    8,    8,    8    ]; % formerly rho
a_pspd  = [ 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 ].*a0;
a_gspd  = [ 0.05, 0.05, 0.05, 0.05, 0.05, 0.05 ].*a0;

b0      = [ 1.0,  1.0,  1.0,  1.0,  1.0,  1.0  ];
b_pspd  = [ 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 ].*b0;
b_gspd  = [ 0.05, 0.05, 0.05, 0.05, 0.05, 0.05 ].*b0;

c0      = [ 1.0,  3.0,  1.0,  1.0,  1.0,  1.0  ];
c_pspd  = [ 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 ].*c0;
c_gspd  = [ 0.05, 0.05, 0.05, 0.05, 0.05, 0.05 ].*c0;

m0      = [ 1.0,  1.0,  1.0,  1.0,  1.0,  1.0  ];
m_pspd  = [ 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 ].*m0;
m_gspd  = [ 0.05, 0.05, 0.05, 0.05, 0.05, 0.05 ].*m0;

p0      = [ 1.0,  1.0,  1.0,  1.0,  1.0,  1.0  ];
p_pspd  = [ 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 ].*p0;
p_gspd  = [ 0.05, 0.05, 0.05, 0.05, 0.05, 0.05 ].*p0;


% artificial constants
smalldist = 1e-6;     % floats with lower difference are considered equal



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% POSITIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%%%%%%%%%%%%%%%% hexagonal lattice %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if lattice(1) == 'h'
  %------------- determine maximal # of y=0 cells --------------
  imax = -1;
  numcells2 = -1;
  while numcells2<numcells
    imax = imax+2;
    numcells2 = imax;
    numcellrow = imax-1;
    jmax = 0;
    while numcellrow >= (imax+1)/2
      numcells2 = numcells2 + 2*numcellrow;
      numcellrow = numcellrow-1;
      jmax = jmax+1;
    end
  end %while
  
  %------------- compute cell positions ------------------------
  x = zeros(1,imax);
  y = zeros(1,imax);
  xspacing = spacing;
  yspacing = spacing*sqrt(3)/2;

  j=0; % ---------------
  xstart = ( - floor((imax-1)/2) - mod(j,2)/2 )*xspacing;
  xend   = (   floor((imax-1)/2) + mod(j,2)/2 )*xspacing;

  xx = xstart:xspacing:xend+smalldist;
  yy = j*yspacing*ones(1,imax);
  xx = xx + x_pspd*rand(1,imax) + x_gspd*randn(1,imax);
  yy = yy + y_pspd*rand(1,imax) + y_gspd*randn(1,imax);
  x = xx;
  y = yy;

  imax = imax-1; % -----

  for j=1:jmax
    xstart = ( - floor((imax-1)/2) - mod(j,2)/2 )*xspacing;
    xend   = (   floor((imax-1)/2) + mod(j,2)/2 )*xspacing;

    xx = xstart:xspacing:xend+smalldist;
    yy = j*yspacing*ones(1,imax);      % upper row
    xx = xx + x_pspd*rand(1,imax) + x_gspd*randn(1,imax);
    yy = yy + y_pspd*rand(1,imax) + y_gspd*randn(1,imax);
    x = [ x, xx ];
    y = [ y, yy ];
    
    xx = xstart:xspacing:xend+smalldist;
    yy = -j*yspacing*ones(1,imax);     % lower row
    xx = xx + x_pspd*rand(1,imax) + x_gspd*randn(1,imax);
    yy = yy + y_pspd*rand(1,imax) + y_gspd*randn(1,imax);
    x = [ x, xx ];
    y = [ y, yy ];
    
    imax = imax-1;
  end

%%%%%%%%%%%%%%%% square lattice %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif lattice(1) == 's'
  imax = ceil(sqrt(numcells));
  jmax = imax;
  numcells2 = imax*jmax;
  x = [];
  y = [];
  for i=1:imax
    xx = 0:spacing:imax*spacing-smalldist;
    yy = i*spacing*ones(1,imax);
    xx = xx + x_pspd*rand(1,imax) + x_gspd*randn(1,imax);
    yy = yy + y_pspd*rand(1,imax) + y_gspd*randn(1,imax);
    x = [ x, xx ];
    y = [ y, yy ];
  end
  
%%%%%%%%%%%%%%%% other lattice %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
  error('stop createconf.m: unknown lattice type');
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% CELL PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%%%%%%%%%%%%%%%% determine radius and weights %%%%%%%%%%%%%%%%%%
r = r0 + r_pspd*rand(1,numcells2) + r_gspd*randn(1,numcells2);
w = r;


%%%%%%%%%%%%%%%% determine type %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initial check
fsum = sum(frac);
ftyp = find(frac);
sizftyp = size(ftyp);
if sizftyp(1)==0 || sizftyp(2)==0
  error('stop createconf.m: no celltypes defined in frac');
end
if abs( sum(frac) - 1.00 ) > smalldist
  frac = frac/fsum;
end

% compute # cells per type according to frac
[ fracmax, fmind ] = max(frac);
sizfrac = size(frac);
numtyp = sizfrac(2);
ctnum = zeros(1,sizfrac(2));
for i=1:numtyp
  if i ~= fmind
    ctnum(i) = ceil( frac(i)*numcells2 );
  end
end
ctnum(fmind) = numcells2 - sum(ctnum);

% set type of cells
typ = zeros(1,numcells2);
cellpool = 1:numcells2;
for i=1:numtyp
  for j=1:ctnum(i)
    sizcellpool = size(cellpool);
    k = ceil( sizcellpool(2)*rand(1,1) );
    typ(cellpool(k)) = i;
    cellpool = setdiff( cellpool, [cellpool(k)] );
  end
end


%%%%%%%%%%%%%%%% determine other cell properties %%%%%%%%%%%%%%%
Tp = zeros(1,numcells2);
a  = zeros(1,numcells2);
b  = zeros(1,numcells2);
c  = zeros(1,numcells2);
m  = zeros(1,numcells2);
p  = zeros(1,numcells2);

for i=1:numcells2
  Tp(i) = Tp0(typ(i)) + Tp_pspd(typ(i))*rand(1,1) + Tp_gspd(typ(i))*randn(1,1);
  a(i)  =  a0(typ(i)) +  a_pspd(typ(i))*rand(1,1) +  a_gspd(typ(i))*randn(1,1);
  b(i)  =  b0(typ(i)) +  b_pspd(typ(i))*rand(1,1) +  b_gspd(typ(i))*randn(1,1);
  c(i)  =  c0(typ(i)) +  c_pspd(typ(i))*rand(1,1) +  c_gspd(typ(i))*randn(1,1);
  m(i)  =  m0(typ(i)) +  m_pspd(typ(i))*rand(1,1) +  m_gspd(typ(i))*randn(1,1);
  p(i)  =  p0(typ(i)) +  p_pspd(typ(i))*rand(1,1) +  p_gspd(typ(i))*randn(1,1);
end

% TODO chech that Tp a b c m p are positive



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% CELL REMOVAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%%%%%%%%%%%%%%%% restrict to numcells %%%%%%%%%%%%%%%%%%%%%%%%%%
cellpool = 1:numcells2;
for i=1:numcells2-numcells
  sizcellpool = size(cellpool);
  k = ceil( sizcellpool(2)*rand(1,1) );
  cellpool = setdiff( cellpool, [cellpool(k)] ); 
end

x1   =   x(cellpool);
y1   =   y(cellpool);
r1   =   r(cellpool);
w1   =   w(cellpool);
typ1 = typ(cellpool);
Tp1  =  Tp(cellpool);
a1   =   a(cellpool);
b1   =   b(cellpool);
c1   =   c(cellpool);
m1   =   m(cellpool);
p1   =   p(cellpool);


%%%%%%%%%%%%%%%% remove desired fraction %%%%%%%%%%%%%%%%%%%%%%%
cellpool = 1:numcells;
for i=1:floor(dillute*numcells)
  sizcellpool = size(cellpool);
  k = ceil( sizcellpool(2)*rand(1,1) );
  cellpool = setdiff( cellpool, [cellpool(k)] ); 
end
sizcellpool = size(cellpool);
numcells3 = sizcellpool(2);

x2   =   x1(cellpool);
y2   =   y1(cellpool);
r2   =   r1(cellpool);
w2   =   w1(cellpool);
typ2 = typ1(cellpool);
Tp2  =  Tp1(cellpool);
a2   =   a1(cellpool);
b2   =   b1(cellpool);
c2   =   c1(cellpool);
m2   =   m1(cellpool);
p2   =   p1(cellpool);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% PRINT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%% save full tissue (numcells2 cells) %%%%%%%%%%%%
if save_full_tissue == 1
  ft_filnam = sprintf('%ic_%s_full.txt', numcells, lattice )
  [ ftf, msg ] = fopen( ft_filnam, 'w' );
  if ftf == -1
    error('stop createconf.m: could not open %s / %s\n', ft_filnam, msg );
  end

  fprintf( ftf, '%s\n', '% x,y position' );
  fprintf( ftf, '%s\n', '% r,w radius and weight' );
  fprintf( ftf, '%s\n', '% typ type of cell' );
  fprintf( ftf, '%s\n', '% Tp persistence time' );
  fprintf( ftf, '%s\n', '% rho density of radial actin filaments on cell body # / \mu m' );
  fprintf( ftf, '%s\n', '% bi bound integrin within whole cell (to be gauged)' );
  fprintf( ftf, '%s\n', '% cd cadherin on cell-cell contacts   (to be gauged)' );
  fprintf( ftf, '%s\n', '% my myosin level                     (to be gauged)' );
  fprintf( ftf, '%s\n', '% ps protrusivity, scales stochastic perturbation' );
  fprintf( ftf, '%s\n', '%' );
  fprintf( ftf, '%s\n', '% x       y      r     w     typ  Tp     rho   bi    cd     my    ps' );

  for i=1:numcells2
    fprintf( ftf, '% 6.2f, % 6.2f,  %4.2f, %4.2f,  ', x(i), y(i), r(i),w(i) );
    fprintf( ftf, '%1.i, %4.0f,  %5.2f, ',          typ(i),Tp(i), a(i)      );
    fprintf( ftf, '%4.2f, %4.2f,  %4.2f, %4.2f\n',    b(i), c(i), m(i),p(i) );
  end
  
  fclose(ftf);
end
 

%%%%%%%%%%%%%%%% save basic tissue (numcells cells) %%%%%%%%%%%%
if save_base_tissue == 1
  bt_filnam = sprintf('%ic_%s_basic.txt', numcells, lattice )
  [ btf, msg ] = fopen( bt_filnam, 'w' );
  if btf == -1
    error('stop createconf.m: could not open %s / %s\n', bt_filnam, msg );
  end
  
  fprintf( btf, '%s\n', '% x,y position' );
  fprintf( btf, '%s\n', '% r,w radius and weight' );
  fprintf( btf, '%s\n', '% typ type of cell' );
  fprintf( btf, '%s\n', '% Tp persistence time' );
  fprintf( btf, '%s\n', '% rho density of radial actin filaments on cell body # / \mu m' );
  fprintf( btf, '%s\n', '% bi bound integrin within whole cell (to be gauged)' );
  fprintf( btf, '%s\n', '% cd cadherin on cell-cell contacts   (to be gauged)' );
  fprintf( btf, '%s\n', '% my myosin level                     (to be gauged)' );
  fprintf( btf, '%s\n', '% ps protrusivity, scales stochastic perturbation' );
  fprintf( btf, '%s\n', '%' );
  fprintf( btf, '%s\n', '% x       y      r     w     typ  Tp     rho   bi    cd     my    ps' );

  for i=1:numcells
    fprintf( btf, '% 6.2f, % 6.2f,  %4.2f, %4.2f,  ', x1(i),y1(i),r1(i),w1(i) );
    fprintf( btf, '%1.i, %4.0f,  %5.2f, ',         typ1(i),Tp1(i),a1(i)       );
    fprintf( btf, '%4.2f, %4.2f,  %4.2f, %4.2f\n',    b1(i),c1(i),m1(i),p1(i) );
  end
  
  fclose(btf);
end


%%%%%%%%%%%%%%%% save dilluted tissue %%%%%%%%%%%%%%%%%%%%%%%%%%
if save_dill_tissue == 1
  dt_filnam = sprintf('%ic_%s_dilluted.txt', numcells, lattice )
  [ dtf, msg ] = fopen( dt_filnam, 'w' );
  if dtf == -1
    error('stop createconf.m: could not open %s / %s\n', dt_filnam, msg );
  end
  
  fprintf( dtf, '%s\n', '% x,y position' );
  fprintf( dtf, '%s\n', '% r,w radius and weight' );
  fprintf( dtf, '%s\n', '% typ type of cell' );
  fprintf( dtf, '%s\n', '% Tp persistence time' );
  fprintf( dtf, '%s\n', '% rho density of radial actin filaments on cell body # / \mu m' );
  fprintf( dtf, '%s\n', '% bi bound integrin within whole cell (to be gauged)' );
  fprintf( dtf, '%s\n', '% cd cadherin on cell-cell contacts   (to be gauged)' );
  fprintf( dtf, '%s\n', '% my myosin level                     (to be gauged)' );
  fprintf( dtf, '%s\n', '% ps protrusivity, scales stochastic perturbation' );
  fprintf( dtf, '%s\n', '%' );
  fprintf( dtf, '%s\n', '% x       y      r     w     typ  Tp     rho   bi    cd     my    ps' );

  for i=1:numcells3
    fprintf( dtf, '% 6.2f, % 6.2f,  %4.2f, %4.2f,  ',x2(i),y2(i),r2(i),w2(i) );
    fprintf( dtf, '%1.i, %4.0f,  %5.2f, ',        typ2(i),Tp2(i),a2(i)       );
    fprintf( dtf, '%4.2f, %4.2f,  %4.2f, %4.2f\n',   b2(i),c2(i),m2(i),p2(i) );
  end
  
  fclose(dtf);
end

