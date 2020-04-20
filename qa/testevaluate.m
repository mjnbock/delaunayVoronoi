%%%%%%%%%%%%%%%% parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('/home/mab/delaunay_voronoi/programs/fulloct2d/');
addpath('/home/mab/delaunay_voronoi/programs/fulloct2d/func');

inidatpath = '/home/mab/delaunay_voronoi/programs/fulloct2d/qa/';
%inidatfiles = {  '1c.txt';
%                 '2c_cosalph.txt';
%                 '2c_digest.txt';
%                 '2c_eq_inv.txt';
%                 '2c_eq.txt';
%                 '2c_inv.txt';
%                 '2c.txt';
%                 '3c_stripeq.txt';
%                 '3c_stripe.txt';
%                 '3c_triang.txt';
%                 '3c_tri_eq.txt';
%                 '4c_aj_ul.txt';
%                 '4c_quad.txt';
%                 '4c_sqare.txt';
%                 '4c.txt';
%                 '4c_wide.txt';
%                 '5c_3comnbr.txt';
%                 '5c_digest2.txt';
%                 '5c_digest3.txt';
%                 '5c_digest.txt';
%                 '5c_line_eq.txt';
%                 '5c_line.txt';
%                 '5c.txt';
%                 '6c_close.txt';
%                 '6c.txt';
%                'e5c.txt';        };

%inidatpath = '/home/mab/';
inidatfiles = { '3c_triang.txt' };


Pmax = 3;

nol = 32;       % maximum # of possible overlaps to assume for single cell

mindist = 1e-6; % minimal distance, below points are considered equal
minang  = 1e-6; % minimal distance, below angles are considered equal

check_digestion = 0; % check wether cells are digesting
check_thmax = 0;     % check wether cells are starlike (theta < thmax)

runfromoct = 1;
verbose_plot = 0;
verbose_print = 0;

persistence = 3;

simcontin = 0;
plotting = 1;

simstep = 1;
t = 1;



%%%%%%%%%%%%%%%% file and save settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save options
savdat = 1;          % save simulation data
movie  = 1;          % write out movie pictures
moveps = 1;

% various increments, do what after how many steps
plotincr   = 001;    % plot/picture after how many time steps
datincr    =  01;    % data save increment
expdatincr = 001;    % expensive data save increment

% base path
save_path = './simdata/';  % directory under which everything is saved

% configuration switches
save_sim  = savdat;  % save (any) simulation observables
save_vsim = 0;       % verbosely save simulation observables TODO
save_ini  = 1;       % save initial configuration
save_fin  = 0;       % save  final  configuration
save_vforce = 0;     % verbosely save forces TODO

% simulation data
save_sim_file  = 'sim.dat';    % global observables
save_x_file    = 'x.dat';      % x(t),y(t) for all cells i
save_v_file    = 'v.dat';      % vx(t),vy(t),v(t) for all cells i
save_pol_file  = 'pol.dat';    % polx(t),poly(t),pol(t) for all cells i
save_ftot_file = 'ftot.dat';   % Fx(t),Fy(t),F(t) for all cells i

% observables expensive to measure
save_exp_file  = 'simx.dat';   % global observables

% initial and final configuration
save_ini_file = 'initial.dat'; % initial data: x,y,r,w,rho,Tp
save_fin_file = 'final.dat';   %  final  data: x,y,r,w,rho,Tp
save_ist_file = 'inist.m';     % initial seed and time
save_fst_file = 'finst.m';     %  final  seed and time

% verbose simulation data
save_r_file     = 'r.dat';     % r(t)      for all cells i
save_w_file     = 'w.dat';     % w(t)      for all cells i
save_rho_file   = 'rho.dat';   % rho(t)    for all cells i
save_Tp_file    = 'Tp.dat';    % Tp(t)     for all cells i
save_floc_file  = 'floc.dat';  % Floc(t)   for all cells i
save_fint_file  = 'fint.dat';  % Fint(t)   for all cells i
save_finth_file = 'finth.dat'; % Fint(t)   for all cells i
save_fintv_file = 'fintv.dat'; % Fint(t)   for all cells i

% from initialize.m
opmod = 'at';



%%%%%%%%%%%%%%%% plot settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
toggle_octplot;
%vfig = figure();
fig = figure();

movpath = strcat(save_path,'pics/');   % path to picture directory
movprfx = '';                          % prefix to picture file names

movtime = 0.0;
movtformat = 't = %.1f h'; % format for time in movie
movtdiv    = 60*60;        % divisor of time in movie

static_scale = 0;      % apply static plot scaling
asprat = 4/3;          % aspect ratio for plotting

plot_cell_centnu = 1;  % plot cell center  numbers        1:y else no
plot_cell_center = 0;  % plot cell centers                1:y else no
plot_cell_body   = 1;  % plot cell bodies                 1:y else no
plot_Pmax_sphere = 0;  % plot Pmax spheres                1:y else no
plot_contact_arc = 1;  % plot contact arcs                1:y else no
plot_arc_pair    = 0;  % plot pair indices                1:y else no
plot_margin_arc  = 1;  % plot marginal arcs               1:y else no
plot_margin_num  = 0;  % plot cell of margin              1:y else no
plot_vertex      = 0;  % plot vertex numbers              1:y else no
plot_wertex      = 0;  % plot wertex numbers              1:y else no
plot_Mij         = 0;  % plot Gamma_ij centers            1:y else no
plot_delaunay    = 1;  % plot Delaunay triang             1:y else no
plot_full_force  = 0;  % plot force on cells              1:y else no
plot_int_force   = 0;  % plot interaction force on cells  1:y else no
plot_int_forceh  = 0;  % plot fint horizontal             1:y else no
plot_int_forcev  = 0;  % plot fint vertical               1:y else no
plot_loc_force   = 0;  % plot locomotion  force on cells  1:y else no
plot_pol_force   = 0;  % plot polarity    force on cells  1:y else no
plot_vis_force   = 0;  % plot viscous     force on cells  1:y else no

plot_axis        = 1;  % plot axis                        1:y else no
plot_polarity    = 0;  % plot polarity vector             1:y else no
plot_velocity    = 0;  % plot velocity                    1:y else no

col_cell = 'gbkcmy';   % color of cell type 1                         green
                       % color of cell type 2                          blue
                       % color of cell type 3                         black
                       % color of cell type 4                          cyan
                       % color of cell type 5                       magenta
                       % color of cell type 6                        yellow

col_Pmax_sphere = 'y'; % color of cell Pmax spheres                  yellow
col_contact_arc = 'r'; % color of contact arcs                        green
col_arc_pair    = 'r'; % color of contact arc pair indices            green
col_margin_arc  = 'k'; % color of marginal arcs                yellow/black
col_margin_num  = 'k'; % color of cell number of marginal arc  yellow/black
col_vertex      = 'b'; % color of numbers appearing at vertex point    cyan
col_wertex      = 'c'; % color of numbers appearing at wertex point    blue
col_Mij         = 'g'; % color of Gamma_ij center points              green
col_delaunay    = 'g'; % color of Delaunay triangulation  yellow/green/blue
col_full_force  = 'c'; % color of force on cells                 cyan/black
col_int_force   = 'c'; % color of interaction force on cells           cyan
col_int_forceh  = 'c'; % color of fint horizontal                      cyan
col_int_forcev  = 'm'; % color of fint vertical                     magenta   
col_loc_force   = 'y'; % color of locomotion  force on cells         yellow

lw_cell_body   = 1;    % linewidth cell body     1/2
lw_contact_arc = 1;    % linewidth contact arc   1/1
lw_margin_arc  = 1;    % linewidth margin arc    1/1

sty_cell_center = 'g*'; % color of cell number appearing at cell centers

scal_force = 10.00; % force scaling for plotting

npts_contact_arc = 256; % number of points to approximate cell-cell contacts
npts_margin_arc  = 256; % number of points to approximate cell closure circles
npts_cell_body   = 256; % number of points to approximate cell body circles
npts_Pmax_sphere = 256; % number of points to approximate cell closure circles



%%%%%%%%%%%%%%%% testing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inidatmax = size(inidatfiles);
inidatmax = inidatmax(1);

for i=1:inidatmax
  inidatfile = inidatfiles{i};
  inidatfile = strcat(inidatpath,inidatfiles{i})
  inidat = load(inidatfile);

  dum = size(inidat);
  cm  = dum(1);

  x   = inidat(:, 1);     % center x-coord
  y   = inidat(:, 2);     % center y-coord
  r   = inidat(:, 3);     % radii

  typ = ones(cm,1);
  w   = r;

  v     = [x,y]/10;
  pol   = [x,y]/10;
  Fint  = [x,y]/10;
  Floc  = [x,y]/10;
  Fpol  = [x,y]/10;
  Fvis  = [x,y]/10;

  
  % variables for voronoi partition
  c2a = zeros(cm,nol);       % indices of Voronoi arcs arround cell {p}
  c2f = zeros(cm,nol);       % indices of free, marginal arcs arround cell {q}
  c2p = zeros(cm,nol);       % indices of pairs containing cell {m}
  c2t = zeros(cm,nol);       % indices of vertices containing cell {n}
  c2am = zeros(1,cm);        % # of Voronoi arcs arround cell
  c2fm = zeros(1,cm);        % # of free, marginal arcs arround cell
  c2pm = zeros(1,cm);        % # of pairs containing cell
  c2tm = zeros(1,cm);        % # of tripels containing cell


  % Pair -> Cells, Triples; contact Spheres, pseudo Wertices
  pm = 0;                    % # pairs
  p2s = zeros(nol*cm,7);     % contact sphere  1:Rij  2:Mijx  3:Mijy
                             %   4:cos(thmax)  5:thmax
                             %   6:thij  orientation angle of vx_big - vx_small
                             %   7:delta cell body surface distance
                             % or plane: 1:NaN  2:x0   3:y0  (local origin)
                             %   4:slpx 5:slpy normalized slope of straight line
                             %   6: thij  orientation angle of vx_i - vx_j
                             %   7:delta cell body surface distance
  p2c = zeros(nol*cm,2);     % indices of cells i,j forming pair
  p2t = zeros(nol*cm,nol);   % indices of tripels containing pair n,nn,...
  p2tm = zeros(1,nol*cm);    % # tripels containing pair
  % TODO p2d wether pair delaunay neighbors


  % Tripel -> Cells, Pairs; Vertices, pseudo Wertices
  tm = 0;                    % # triples
  t2v = zeros(nol*cm,3);     % vertex coords, 1:x 2:y, 3:P_distance
  t2c = zeros(nol*cm,3);     % indices of cells forming tripel i,j,k
  t2p = zeros(nol*cm,3);     % indices of pairs forming tripel m,mm,...


  % pseudo Wertices arising from Pmax
  wm = 0;                    % # pseudo vertices
  w2w = zeros(2*nol*cm,3);   % pseudo vertex coordinates, 1:x 2:y, 3:P_distance
  p2w = zeros(nol*cm,2);     % pseudo vertices of pair o,oo


  % contact arcs
  am = 0;                    % number of contact arc candidates
  a2a = zeros(nol*cm,14);    % 1:m  pair index
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
  fm = 0;                    % number of free closure arcs
  f2f = zeros(nol*cm,6);     % 1: cell index i  2: Pmax*w(i)
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
  nbi  = zeros(1,nol);       % neighbors of cell i
  nbj  = zeros(1,nol);       % neighbors of cell j
  vcl  = zeros(2*nol,6);     % vertex candidate list per pair
                             % 1:orientation angle thv (towards Mij or cell i),
                             % 2:P_distance w.r.t. cells,
                             % 3:x 4:y  vertex coordinates
                             % 5:index n,o  6:type=0w/1t/2u of v/w/uertex
  vl   = zeros(2*nol,6);     % vertex list per pair
                             % 1:orientation angle thv (towards Mij or cell i),
                             % 2:P_distance w.r.t. first cell,
                             % 3:x 4:y  vertex coordinates
                             % 5:index n,o  6:type=0w/1t/2u of v/wertex
  svl  = zeros(2*nol,6);     % sorted vertex list per pair
                             % 1:orientation angle thv (towards Mij or cell i),
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


  printf('\ninivor timing:\n');
  asctime( localtime(time()) )
  inivor;
  asctime( localtime(time()) )

  figure(floor(fig));
  printf('\nplotcells timing:\n');
  plotcells;
  asctime( localtime(time()) )

  openfiles;
  printf('\nevaluate timing:\n');
  preevaluate;
  evaluate;
  asctime( localtime(time()) )
  closefiles;

  input('press enter to continue');
end

