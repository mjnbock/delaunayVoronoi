%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% modifications for quality-assurance tests in ./qa
if exist('testing','var');
  disp('initializing testing environment:');
  if testing==2
    ttmax = test_tmax
    plotincr   = test_plotincr
    datincr    = test_datincr
    expdatincr = test_expdatincr
  end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% computing/loading initial values
if simcontin == 1
  finst;
  ttmax = ttmax + t;
  randn( 'state', state );
  sffn = strcat( save_path, save_fin_file );
  inidat = load( sffn );
  opmod = 'at';
  if movie == 1
    movtime = t/movtdiv;
  end
else
  t = 0;
  simstep = 0;
  opmod = 'wt';
  if movie == 1
    movtime = 0;
  end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% graphics
if plotting == 1
  if runfromoct == 1
    toggle_octplot;
  end
  fig=figure();
  if runfromoct ~= 1
    %axis off
    set(fig,'DoubleBuffer','on');
    % mov = avifile('D:\Simulations\Multi_cell\multicell.avi');
    % mov.quality = 100;
    %     hold on;
    %     axis equal
    %     axis([-5 15 -5 17]);
    %     axis tight
    %     axis off
  end
end
if static_scale == 1
  axis( statscal );
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% initial data structure setup
dum = size(inidat);
cm = dum(1);

% cell
sizinidat = size(inidat);
if sizinidat(2) ~= 11
  error('main.m: wrong # columns in inidat array, need x y r w typ Tp rho bi cd my ps');
end
x   = inidat(:, 1);     % center x-coord
y   = inidat(:, 2);     % center y-coord
r   = inidat(:, 3);     % current radii
w   = inidat(:, 4);     % current weights
r0  = inidat(:, 3);     % initial radii
w0  = inidat(:, 4);     % initial weights
typ = inidat(:, 5);     % type of cell
Tp  = inidat(:, 6);     % persistence time
rho = inidat(:, 7);     % filament densities
% TODO 2009-08-04: gauge and clarify interpretation of these quantities
bi  = inidat(:, 8);     % relative amount of bound integrin ~ 1
cd  = inidat(:, 9);     % relative amount of cadherin ~ 1
my  = inidat(:,10);     % relative amount of myosin ~ 1
ps  = inidat(:,11);     % relative amount of perturbation ~ 1 (protrusiveness)

% dynamics of cells
Fint  = zeros(cm,2);  % interaction force acting on cell
Finth = zeros(cm,2);  % fint horizontal
Fintv = zeros(cm,2);  % fint vertical
Floc  = zeros(cm,2);  % locomotion force acting on cell
Fpol  = zeros(cm,2);  % polarity   force acting on cell
Fvis  = zeros(cm,2);  % viscous    force acting on cell
v   = zeros(cm,2);    % velocity of cell
pol = zeros(cm,2);    % polarity vector of cell
%pol = randn(cm,2);    % polarity vector of cell

drag = zeros(1,cm);
for i = 1:cm
  % all in all, drag(i) should be ~ 5e5 ... 1.5e6

  % scale by mass computed from cell volume
  %cell_mass_density = 1e-9; % 1 pg / micron^3, but need mg as units
  %cell_mass = cell_mass_density * 4*( Pmax*r(i) )^3*pi/3
  % does not work for small cells (unstable);
  % how is the areal radius related to the volume radius of a cell?


  if friction == 1
    % do not scale drag coefficient, take constant mass 5e-6 mg per cell
    cell_mass = 5e-6;
    drag(i) = 2*drag_coeff1 * cell_mass; % TODO why 2?
    % = 1.5e6 w/ drag_coeff=1.5e10
  elseif friction == 2
    % scale by cell area
    % Note: 2011-03-02
    %       We assume a mean cell body radius r of
    %       approximately 2 um. 
    %       Only for this radius a drag_coeff2=1326.3
    %       will yield the expected migration velocity.
    %       Bigger cells will be slower, smaller cells
    %       will be faster.
    %       ------------------------------------------------
    %       From Drosophila wing development data one estimates 
    %       mean(r) ~ 1.04 um @ dcrit=0.50, Pmax=3.0
    %       mean(r) ~ 0.89 um @ dcrit=0.66, Pmax=3.0
    %       [Current Biology 19:1950 (2010)].
    %       ------------------------------------------------
    %       During Drosophila germ band elongation, we have
    %       mean(r) ~ 4.00 um @ dcrit=0.33, Pmax=3.0
    %       mean(r) ~ 3.33 um @ dcrit=0.50, Pmax=3.0
    %       [Nature Cell Biology 10:1401 (2008)]
    %       ------------------------------------------------
    %       MDCK cells in almost polarized epithelia would be
    %       mean(r) ~ 6.77 um @ dcrit=0.33, Pmax=3.0 (Acell=400um^2)
    %       In migrating, mesenchymal state we have
    %       mean(r) ~ 6.77 um @ dcrit=0.68, Pmax=3.0 (Acell=800um^2)
    %       [PRL 104:168104  (2010)]
    %       ------------------------------------------------
    %       see also notes D51-D53, ~ 2010-03-01
    drag(i) =  drag_coeff2 * (Pmax*r(i))^2*pi;
  else
    error('initialize.m: unknown friction model');
  end
end

% growth, mitosis and cellcycle
Tgrow1=zeros(cm,1);
Tmito1=zeros(cm,1);
for i=1:cm
  % setup time interval of cell growth (between two mitotic divisions)
  Tgrow1(i) = Tgrow0(typ(i)) + PTgrow(typ(i))*rand(1,1);
  % setup time interval of cell mitosis (aka cytokinesis)
  Tmito1(i) = Tmito0(typ(i)) + PTmito(typ(i))*rand(1,1);
end
Tgrow=ones(cm,1)*(-dt);
Tmito=ones(cm,1)*(max(Tmito1)+dt);

Gam0 = zeros(1,cm);        % length marginal boundary of cell
Gamp = zeros(1,cm);        % length interaction pair boundary of cell

if perturb_vel ~= 1
  evel = 0.0;
end


% short explanation of the spirit of variable naming:
%
% c Cell
% p (cell) Pair
% t (cell) Tripel (pairwise overlapping)
% a contact arc candidate
% f free closure arc
%
% i,j,k,l  % cell indices
% m,mm,... % overlapping cell pair indices
% n,nn,... % overlap triple (or vertex) indices
% o,oo,... % pseudo vertex indices
% p,pp,... % Voronoi arc indices
% q,qq,... % free, marginal arc indices


% relating and linking the various involved objects
% Cell -> Pairs, Triples
%cm                        % # cells       next two: XXX dangerous /2
c2a = zeros(cm,nol/2);     % indices of Voronoi arcs arround cell {p}
c2p = zeros(cm,nol/2);     % indices of pairs containing cell {m}
c2am = zeros(1,cm);        % # of Voronoi arcs arround cell
c2pm = zeros(1,cm);        % # of pairs containing cell


% Pair -> Cells, Triples
pm = 0;                    % # pairs
p2s = zeros(nol*cm,7);     % contact sphere  1:Rij  2:Mijx  3:Mijy
                           %   4:cos(thmax)  5:thmax
                           %   6:thij  orientation angle of vx_big - vx_small
                           %   7:delta cell body surface distance
                           % or plane: 1:NaN  2:x0   3:y0  (local origin)
			   %   4:slpx 5:slpy  normalized slope of straight line
                           %   6: thij  orientation angle of vx_i - vx_j
                           %   7:delta cell body surface distance
p2c = zeros(nol*cm,2);     % indices of cells i,j forming pair
p2t = zeros(nol*cm,nol);   % indices of tripels containing pair n,nn,...
p2v = zeros(nol*cm,nol);   % vertices of pair n,nn,...
p2tm = zeros(1,nol*cm);    % # tripels containing pair
p2vm = zeros(1,nol*cm);    % # vertices of pair
% TODO p2d wether pair delaunay neighbors


% Tripel -> Vertices, Cells, Pairs
tm = 0;                    % # triples
t2c = zeros(nol*cm,3);     % indices of cells forming tripel i,j,k
t2p = zeros(nol*cm,3);     % indices of pairs forming tripel m,mm,...


% Vertices, Vertices -> Cells
vm = 0;                    % # vertices
v2v = zeros(nol*cm,3);     % vertex coords 1:x 2:y, 3:P_distance
v2c = zeros(nol*cm,3);     % indices of cells forming vertex i,j,k


% pseudo Wertices arising from Pmax
wm = 0;                    % # pseudo vertices
w2w = zeros(2*nol*cm,3);   % pseudo vertex coordinates, 1:x 2:y, 3:P_distance
p2w = zeros(nol*cm,2);     % pseudo vertices of pair o,oo


% contact arcs
am = 0;                    % number of contact arc candidates XXX dangerous /2
a2a = zeros(nol*cm/2,14);  % 1:m  pair index
                           % 2:R  radius of contact sphere
                           % 3:Mx 4:My  center of sphere
                           % 5:theta1 6:theta2  angle range of contact surface
                           % 7,8 + 9,10:     start and end vertex
                           % 11,12  type  of start and and vertex 0w/1t
                           % 13,14  index of start and and vertex  o/n
                           % or flat contact surface
                           % 2: NaN
                           % 3:xi 4:yi  coords of cell i
                           % 5:0  6:0


% free arcs
fm = 0;                    % number of free closure arcs XXX dangerous /5
f2f = zeros(nol*cm/4,6);   % 1: cell index i  2: Pmax*w(i)
                           % 3:xi 4:yi  cell center
                           % 5:thmin 6:thmax w.r.t. cell center


% Voronoi/Delaunay neighbors
nb   = zeros(cm,nol);      % indices of Voronoi neighbors of cell
nbm  = zeros(1,cm);        % number of Voronoi neighbors of cell
nbd  = zeros(1,nol*cm);    % neighbor pair
nbdm = 0;                  % # neighbor pairs


% free cell array
fcm = 0;                   % # of free cells
fc  = zeros(cm,1);         % indices of free cells


% marginal cell array
mcm = 0;                   % # of marginal cells
mc  = zeros(cm,1);         % indices of marginal cells


% actual Vornonoi v/wertices
%vxl = 0;                  % indices of vertices needed in partition n
%wxl = 0;                  % indices of wertices needed in partition o
vxlm = 0;
wxlm = 0;


% temporary variables for voronoi partition
vcl  = zeros(2*nol,6);     % vertex candidate list per pair
                           % 1:orientation angle thv towards Mij or cell i
                           % 1:P_distance w.r.t. cells,
                           % 2:x 3:y  vertex coordinates
                           % 4:index n,o  5:type=0w/1t/2u of v/w/uertex
vl   = zeros(2*nol,6);     % vertex list per pair
                           % 1:orientation angle thv (towards Mij or left cell),
                           % 2:P_distance w.r.t. first cell,
                           % 3:x 4:y  vertex coordinates
                           % 5:index n,o  6:type=0w/1t/2u of v/wertex
svl  = zeros(2*nol,6);     % sorted vertex list per pair
                           % 1:orientation angle thv (towards Mij or left cell),
                           % 2:P_distance w.r.t. first cell,
                           % 3:x 4:y  vertex coordinates
                           % 5:index n,o  6:type=0w/1t/2u of v/wertex
vxcl = zeros(3*nol*cm);    % actual vertex candidate list
wxcl = zeros(3*nol*cm);    % actual wertex candidate list
vxclm = 0;                 % # actual vertices after erosion
wxclm = 0;                 % # actual werteces after erosion
al  = zeros(nol,14);       % per cell arc list, cf. a2a above, 13:phi1 14:phi2
sal = zeros(nol,14);       % per cell arc list, sorted after phi1, the
                           % angular starting point of the arc w.r.t. cell
                           % and phi2 angular ending point of arc w.r.t. cell


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% computing/saving per-run parameter values
% means
rbar   = mean(r);
Tpbar  = mean(Tp);
rhobar = mean(rho);

% relative distances scaled by (Pmax-1), see vFint.m
dmin  =  dmin * (Pmax-1);
dcrit = dcrit * (Pmax-1);

% alpha sanity check
if interaction==3  &&  alpha > 0.0
  if dcrit/alpha < 2*dmin  || dcrit*alpha < 2*dmin
    error('stop initialize.m: dcrit too low, decrease alpha');
  end
  if dcrit*alpha > 2*dcrit || dcrit/alpha > 2*dcrit
    error('stop initialize.m: dcrit too high, decrease alpha');
  end
end

% for plotting
norm_int_force  = 2*pi*rhobar*sqrt( finth^2 + fintv^2 )*rbar / scal_force;
norm_int_forceh = 1.0 * norm_int_force;
norm_int_forcev = 1.0 * norm_int_force;
norm_loc_force  = 1.0 * norm_int_force; 
norm_pol_force  = 1.0 * norm_int_force; 
norm_vis_force  = 1.0 * norm_int_force; 
norm_full_force = 1.0 * norm_int_force;

norm_velocity = 1.0 / scal_velocity;
norm_polarity = 1.0 / scal_polarity;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% opening files and saving of initial configuration
if save_sim == 1
  openfiles;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% initial voronoi tesselation
if plotting == 1  &&  verbose_plot == 1
  plotcraw;
  input('press enter to continue');
end
inivor;    % in octave function calls are expensive; thus the inlining here
           % bear in mind: ALL variables carry over to inivor / runvor
if verbose_print == 1 % debug: printout
  printindices;
  printvertices;
end


% first time plotting
if plotting == 1
  if verbose_plot == 1
    plotcraw;
    input('press enter to continue');
  end
  plotcells;
  if verbose_plot == 1
    input('press enter to continue');
  end
  if plotspic == 1
    filnam = sprintf( '%s%ssingle---%.8i.eps', movpath, movprfx, simstep );
    print( '-depsc', filnam );
    exit(0);
  end
end
% first movie picture
if movie == 1
  if plotting ~= 1
    plotting = 1;
    if runfromoct==1
      toggle_octplot;
      fig=figure();
    else
      fig=figure();
      set(fig,'DoubleBuffer','on');
    end
    plotcells;
  end
  %input('resize plotting window as required, then press enter to continue');
end


% initial data saving
if savdat == 1
  if save_ini == 1

    if simcontin ~= 1
      % save initial cell data
      for i = 1:cm
        fprintf( sif, '%13.10g %13.10g %13.10g %13.10g %13.10g %13.10g\n',...
                        x(i),   y(i),   r(i),   w(i),   rho(i), Tp(i)         );
      end

      % initial rng state and time
      state = randn('state');
      fprintf( sistf, '%% time, counter and final state of rng\n' );
      fprintf( sistf, 't        = %i;\n',   t     );
      fprintf( sistf, 'simstep  = %i;\n', simstep );
      fprintf( sistf, 'state = [ ' );
      sizstat = size(state);
      for stind = 1:sizstat(1)-1
        fprintf( sistf, '%i, ', state(stind) );
      end
      fprintf( sistf, '%i ]\n', state(sizstat(1)) );

      fclose(sif);
      fclose(sistf);
    end

    if movie == 1
      % .eps tissue snapshot 
      if simcontin ~= 1
        filnam = strcat( movpath, movprfx, 'initial.eps' );
      else
        filnam = sprintf( '%s%s%.8i.eps', movpath, movprfx, simstep );
      end
      print( '-depsc', filnam )
    end

  end
end

