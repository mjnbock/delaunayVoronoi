%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% program flow
simcontin=0;
if simcontin==0  &&  ~exist('testing','var');
  %close all;
  clear all;
  simcontin=0; % indeed it works like this
end

% cell dynamics
interaction = 3;     % apply interaction forces, cf. fint
                     %     1: isotropic;
                     %     2: anisotropic in arcs, see  alpha
                     %     3: anisotropic in Dcrit, see alpha
                     %     4: line tension + area repulsion EPJE 33:117 TODO

locomotion  = 1;     % apply locomotion  forces, cf. floc

persistence = 3;     % persistent cells, 0: purely overdamped
                     %                   1: phenomenological exponential filter
                     %                   2: overdamped with polarity force
                     %                   3: overdamped with polarity, viscosity

friction = 2;        % friction model 1: same constant for all cells
                     %                2: proportional to free cell area
                     %                3: proportional to actual cell area TODO

cellcycle = 0;       % dynamical cell growth,
                     %      0: off
                     %      1: purely molecular regulation
                     % TODO 2: force-dependent self-organizatoin

perturb_cel = 0;     % perturbation of cell position
perturb_vel = 1;     % perturbation of velocity increment

perturb_gsc = 1;     % scale perturbation by |Gamma|

% plotting options
plotting =  1;       % plot
plotspic =  0;       % stop after first tesselation to record single picture
static_scale = 0;    % apply static plot scaling

% save options
savdat = 1;          % save simulation data
movie  = 1;          % write out movie pictures

% various increments, do what after how many steps
plotincr   = 015;    % plot/picture after how many time steps
datincr    =  15;    % data save increment
expdatincr = 150;    % expensive data save increment

% octave <-> matlab switch
runfromoct = 2;      % run from octave 0:matlab 1:octave<3.0 2:octave>=3.0

% what kind of configurations shall be allowed
check_thmax = 0;     % check wether cells are starlike (theta < thmax)
check_digestion = 1; % check wether cells are digesting 
                     % (overlap w/o pseudo wertices)
check_overlap = 0;   % check wether cell bodies overlap TODO

% force computation mode
force_thmax = 1;     % limit arcs to thmax during force computation
warn_thmax = 1;      % print warning if thmax is violated

% debug options
show_progress = 1;   % displays step and time in console
verbose_print = 0;   % print debugging info
verbose_plot = 0;    % halt after updating plot window



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Initial cell data
if simcontin==0
  inidat = ...
    load('/home/mab/projects/delaunay_voronoi/programs/fulloct2d/data/7c.txt');
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Parameters
% simulation time
dt   =          002;   % time increment (approximately seconds in real time)
tmax =    1*04*3600;   % simulation time

% perturbation
ecel    = 0.00;        % perturbation of cell position
evel    = 5e-3;        % perturbation factor velocity increment ~ 1e-2
                       %   works only for persistence==1
epol    = 0.01;        % perturbation factor polarity increment ~ 1e-1
kapgam  = 0.2;         % 0: perturbation from free boundary Gam0 only
                       % 1: perturbation from contact borders Gamp only

% interaction force
finth =   060.0;       % ~ 40 pN  interaction force per filament bundle (~50)
fintv =   000.0;       % horizontal and vertical component
dcrit =     0.3;       % critical relative distance of cell bodies and
dmin  =     0.1;       % minimum relative distance (interaction force)

alpha =  1.0;          % interaction anisotropy; set alpha<0 to use comp_alpha.m
                       % interaction=3: in dcrit
                       % neutral setting: alpha=1, no anisotropy

fintnum = 032;         % number of points to approximate contact sphere
                       % during force calculation

% locomotion force
floc  =  000.000;      % ~ 10 pN  locomotion force per filament bundle (~50)
flocnum = 032;         % number of points to approximate margin sphere
                       % during force calculation

% polarity force
fpol  =  010.000;      % ~ 10 pN  polarity force per filament bundle (~50)
fpolnum = 70;          % number of points to approximate cell circumference

kappol = 0.5;          % 0.0: symmetric pulling
                       % 0.5: net-forward pulling
                       % 1.0: (unrealistic) forward pulling backward pushing
gampv = 01.00;         % polarity ~ velocity / gampv
polde = 5e-2;          % constant polarity decrease
Tpol =  600.0;         % polarity persistence time

% viscous interaction force
fvis = 200.0;          % viscous cell pair drag per cadherin-filament junction

% velocity dependent drag force
drag_coeff1 = 1.5e10;   % friction==1: ~ 5e10 in 1/s % TODO how is this
drag_coeff2 = 1326.3;   % friction==2: ~ 1326.3 in s^-1 um^-2 per free cell area

% growth cellcycle parameters
Tgrow0  = [ 1*3600, 23*3600, 23*3600, 23*3600, 23*3600, 23*3600 ];
PTgrow  = [     60,    1200,    1200,    1200,    1200,    1200 ];
        % duration and variation of cellular growth phase
growfac = [   1.00,   1.41,   1.41,   1.41,   1.41,   1.41 ];
        % relative size increase during cellular growth phase

% mitosis cellcycle parameters
Tmito0  = [ 1*3600,  1*3600,  1*3600,  1*3600,  1*3600,  1*3600 ];
PTmito  = [     60,      60,      60,      60,      60,      60 ];
        % duration and variation of cellular mitosis phase
mitofac = [   1.41,   1.41,   1.41,   1.41,   1.41,   1.41 ];
        % size decrease during cellular mitosis phase

% Pmax, and equality thresholds
Pmax = 3.00;           % maximal distance function sqrt(P_max)
mindist = 1e-6;        % minimal distance, below points are considered equal
minang  = 1e-6;        % minimal distance, below angles are considered equal

% random numbers
randn('state',sum(100*clock)); % random number generator state

% allocation constants
nol = 20;              % max # of possible overlaps to assume for single cell
                       % utmost upper bound, probably slightly too high
                       % also used as estimate for various other index maxima
                       % >= 20, two shells for each cell, dividable by 4



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Saving setup
% base path
save_path = './simdata/';  % directory under which everything is saved
addpath( save_path );

% configuration switches
save_sim  = savdat;  % save (any) simulation observables
save_vsim = 0;       % verbosely save simulation observables TODO
save_ini  = 1;       % save initial configuration
save_fin  = 1;       % save  final  configuration
save_vforce = 0;     % verbosely save forces TODO



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Movie setup
movpath = strcat(save_path,'pics/');   % path to picture directory
movprfx = '';                          % prefix to picture file names

moveps = 2;                % 1: additionally to jpg, write eps files
                           % 2: additionally to jpg, write eps files
                           %    but only every expdatincr steps
movtformat = 't = %.1f h'; % format for time in movie
movtdiv    = 60*60;        % divisor of time in movie

statscal = [ -15, 20, -10, 15 ];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% run program
base_conf;
if simcontin == 1
  sffn = strcat( save_path, save_fin_file );
  inidat = load( sffn );
end

main;
