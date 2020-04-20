%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Saving setup
% base path and configuration switches from conf.m

% parameters changing in time
save_par_file  = 'param.dat';

% simulation data
save_sim_file  = 'sim.dat';    % global observables
save_x_file    = 'x.dat';      % x(t),y(t) for all cells i
save_v_file    = 'v.dat';      % vx(t),vy(t),v(t) for all cells i
save_A_file    = 'A.dat';      % A(t),Ain(t),Aout(t) for all cells i
save_Per_file  = 'Peri.dat';   % Gam(t),Gamp(t),Gam0(t) for all cells i
save_PA_file   = 'PA.dat';     % perimeter-area ratio for all cells i
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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Plotting setup
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
plot_delaunay    = 0;  % plot Delaunay triang             1:y else no
plot_full_force  = 1;  % plot force on cells              1:y else no
plot_int_force   = 0;  % plot interaction force on cells  1:y else no
plot_int_forceh  = 0;  % plot fint horizontal             1:y else no
plot_int_forcev  = 0;  % plot fint vertical               1:y else no
plot_loc_force   = 0;  % plot locomotion  force on cells  1:y else no
plot_pol_force   = 0;  % plot polarity    force on cells  1:y else no
plot_vis_force   = 1;  % plot viscous     force on cells  1:y else no

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
col_contact_arc = 'r'; % color of contact arcs                          red
col_arc_pair    = 'r'; % color of contact arc pair indices              red
col_margin_arc  = 'k'; % color of marginal arcs                yellow/black
col_margin_num  = 'k'; % color of cell number of marginal arc  yellow/black
col_vertex      = 'c'; % color of numbers appearing at vertex point    cyan
col_wertex      = 'b'; % color of numbers appearing at wertex point    blue
col_Mij         = 'r'; % color of Gamma_ij center points                red
col_delaunay    = 'g'; % color of Delaunay triangulation   green/blue/black
col_full_force  = 'c'; % color of force on cells                       cyan
col_int_force   = 'r'; % color of interaction force on cells            red
col_int_forceh  = 'r'; % color of fint horizontal                       red
col_int_forcev  = 'r'; % color of fint vertical                         red   
col_loc_force   = 'k'; % color of locomotion  force on cells          black
col_pol_force   = 'm'; % color of polarity    force on cells        magenta 
col_vis_force   = 'r'; % color of viscous     force on cells            red 

col_polarity    = 'm'; % color of polarity vector of cells          magenta
col_velocity    = 'g'; % color of velocity of cells        green/blue/black 

lw_cell_body   = 1;    % linewidth cell body     1/2
lw_contact_arc = 1;    % linewidth contact arc   1/1
lw_margin_arc  = 1;    % linewidth margin arc    1/1

sty_cell_center = 'g*'; % color of cell number appearing at cell centers

scal_force    = 10.00; % force    scaling for plotting
scal_velocity =  1e2;  % velocity scaling for plotting
scal_polarity = 00.50; % polarity scaling for plotting

npts_contact_arc = 064; % number of points to approximate cell-cell contacts
npts_margin_arc  = 064; % number of points to approximate cell closure circles
npts_cell_body   = 032; % number of points to approximate cell body circles
npts_Pmax_sphere = 064; % number of points to approximate cell closure circles



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Execution path setup
addpath('/home/mab/projects/delaunay_voronoi/programs/fulloct2d/');
addpath('/home/mab/projects/delaunay_voronoi/programs/fulloct2d/func');
