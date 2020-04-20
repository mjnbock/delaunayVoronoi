%%%%%%%%%%%%%%%% test config %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
testing=2;       % calculations in evaluate.m, mostly geometric ones
addpath('../');

% various increments, do what after how many steps
test_plotincr   = 001;    % plot/picture after how many time steps
test_datincr    =  01;    % data save increment
test_expdatincr = 001;    % expensive data save increment

% simulation time
test_tmax=0;



%%%%%%%%%%%%%%%% files %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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



%%%%%%%%%%%%%%%% testing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inidatmax = size(inidatfiles);
inidatmax = inidatmax(1);

for i=1:inidatmax
  inidatfile = strcat(inidatpath,inidatfiles{i})


  testing
  conf;

  input('press enter to continue');
end

