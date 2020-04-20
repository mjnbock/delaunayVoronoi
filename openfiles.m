%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% openfiles.m
%%%%%%%%
%%%%%%%% open files to save simulation data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parameters changing in time
sparn = strcat( save_path, save_par_file );
[sparf, msg] = fopen( sparn, opmod );
if sparf == -1
  error( 'openfiles.m: could not open file %s\n  %s\n', sparf, msg );
end
fprintf( sparf, '%% t step  alpha\n' );

% simulation observables
ssfn = strcat( save_path, save_sim_file );
[ssf, msg] = fopen( ssfn, opmod );
if ssf == -1
  error( 'openfiles.m: could not open file %s\n  %s\n', ssfn, msg );
end
fprintf( ssf, '%% t step  size_tissue #cell_mar #cell_int  nb_mea nb_min nb_max  nbm_mea nbm_min nbm_max  nbi_mea nbi_min nbi_max  delt_mea delt_min delt_max  v_mea v_min v_max pol_mea pol_min pol_max\n' );

% coordinates of cell centers
sxfn = strcat( save_path, save_x_file );
[sxf, msg] = fopen( sxfn, opmod );
if sxf == -1
  error( 'openfiles.m: could not open file %s\n  %s\n', sxfn, msg );
end
fprintf( sxf, '%% t step x(1) y(1) x(2) y(2) ... ... x(n) y(n)\n' );

% velocities of cells
svfn = strcat( save_path, save_v_file );
[svf, msg] = fopen( svfn, opmod );
if svf == -1
  error( 'openfiles.m: could not open file %s\n  %s\n', svfn, msg );
end
fprintf( svf, '%% t step vx(1) vy(1) v(1) ... vx(n) vy(n) v(n)\n' );

% Areas of cells
sAfn = strcat( save_path, save_A_file );
[sAf, msg] = fopen( sAfn, opmod );
if sAf == -1
  error( 'openfiles.m: could not open file %s\n  %s\n', sAfn, msg );
end
fprintf( sAf, '%% t step A(1) Ain(1) Aout(1) ... A(n) Ain(n) Aout(n)\n' );

% perimeters of cells
sPerfn = strcat( save_path, save_Per_file );
[sPerf, msg] = fopen( sPerfn, opmod );
if sAf == -1
  error( 'openfiles.m: could not open file %s\n  %s\n', sPerfn, msg );
end
fprintf( sPerf, '%% t step  Gam(1) Gamp(1) Gam0(1) ... Gam(n) Gamp(n) Gam0(n)\n' );

% perimeters of cells
sPAfn = strcat( save_path, save_PA_file );
[sPAf, msg] = fopen( sPAfn, opmod );
if sAf == -1
  error( 'openfiles.m: could not open file %s\n  %s\n', sPAfn, msg );
end
fprintf( sPAf, '%% t step  PAcells(1) ... PAcells(n)\n' );

% polarity vectors of cells
spfn = strcat( save_path, save_pol_file );
[spf, msg] = fopen( spfn, opmod );
if spf == -1
  error( 'openfiles.m: could not open file %s\n  %s\n', spfn, msg );
end
fprintf( spf, '%% t step polx(1) poly(1) pol(1) ... polx(n) poly(n) pol(n)\n' );

% total forces
sftotfn = strcat( save_path, save_ftot_file );
[sftotf, msg] = fopen( sftotfn, opmod );
if sftotf == -1
  error( 'openfiles.m: could not open file %s\n  %s\n', sftotfn, msg );
end
fprintf( sftotf, '%% t step Fmea Fmin Fmax  Fx(1) Fy(1) F(1) ... Fx(n) Fy(n) F(n)\n' );



% expensive simulation observables
sexpfn = strcat( save_path, save_exp_file );
[sexpf, msg] = fopen( sexpfn, opmod );
if sexpf == -1
  error( 'openfiles.m: could not open file %s\n  %s\n', sexpfn, msg );
end
fprintf( sexpf, '%% t step   PAtiss PAcore  AmAt AoAt  Atiss  Ain Acore  Aout Amarg  Ptiss Pcore\n' );



% initial data files
if save_ini == 1  &&  simcontin ~= 1

  % cell data
  sifn = strcat( save_path, save_ini_file );
  [sif, msg] = fopen( sifn, 'wt' );
  if sif == -1
    error( 'openfiles.m: could not open file\n  %s\n', msg );
  end

  % rng state and time
  sistfn = strcat( save_path, save_ist_file );
  [sistf, msg] = fopen( sistfn, 'wt' );
  if sistf == -1
    error( 'openfiles.m: could not open file\n  %s\n', msg );
  end

end



% final data files opened before use
% if save_fin == 1
% 
%   % cell data
%   sffn = strcat( save_path, save_fin_file );
%   [sff, msg] = fopen( sffn, 'wt' );
%   if sff == -1
%     error( 'openfiles.m: could not open file\n  %s\n', msg );
%   end
% 
%   % rng state and time
%   sfstfn = strcat( save_path, save_fst_file );
%   [sfstf, msg] = fopen( sfstfn, 'wt' );
%   if sfstf == -1
%     error( 'openfiles.m: could not open file\n  %s\n', msg );
%   end
% 
% end



% verbose files
if save_vsim == 1

  % radii
  srfn = strcat( save_path, save_r_file );
  [srf, msg] = fopen( srfn, opmod );
  if srf == -1
    error( 'openfiles.m: could not open file %s\n  %s\n', srfn, msg );
  end
  fprintf( srf, '%% t step r(1) r(2) ... r(n)\n' );

  % weights
  swfn = strcat( save_path, save_w_file );
  [swf, msg] = fopen( swfn, opmod );
  if swf == -1
    error( 'openfiles.m: could not open file %s\n  %s\n', swfn, msg );
  end
  fprintf( swf, '%% t step w(1) w(2) ... w(n)\n' );

  % filament densities
  srhofn = strcat( save_path, save_rho_file );
  [srhof, msg] = fopen( srhofn, opmod );
  if srhof == -1
    error( 'openfiles.m: could not open file %s\n  %s\n', srhofn, msg );
  end
  fprintf( srhof, '%% t step rho(1) rho(2) ... rho(n)\n' );

  % persistence times
  sTpfn = strcat( save_path, save_Tp_file );
  [sTpf, msg] = fopen( sTpfn, opmod );
  if sTpf == -1
    error( 'openfiles.m: could not open file %s\n  %s\n', srfn, msg );
  end
  fprintf( sTpf, '%% t step Tp(1) Tp(2) ... Tp(n)\n' );

  % verbose forces
  if save_vforce == 1
    sflocfn = strcat( save_path, save_floc_file );
    [sflocf, msg] = fopen( sflocfn, opmod );
    if sflocf == -1
      error( 'openfiles.m: could not open file %s\n  %s\n', sflocfn, msg );
    end
    fprintf( sflocf, '%% t step Flocx(1) Flocy(1) ... Flocx(n) Flocy(n)\n' );
    sfintfn = strcat( save_path, save_fint_file );
    [sfintf, msg] = fopen( sfintfn, opmod );
    if sfintf == -1
      error( 'openfiles.m: could not open file %s\n  %s\n', sfintfn, msg );
    end
    fprintf( sfintf, '%% t step Fintx(1) Finty(1) ... Fintx(n) Finty(n)\n' );
    sfinthfn = strcat( save_path, save_finth_file );
    [sfinthf, msg] = fopen( sfinthfn, opmod );
    if sfinthf == -1
      error( 'openfiles.m: could not open file %s\n  %s\n', sfinthfn, msg );
    end
    fprintf( sfinthf, ...
             '%% t step Finthx(1) Finthy(1) ... Finthx(n) Finthy(n)\n' );
    sfintvfn = strcat( save_path, save_fintv_file );
    [sfintvf, msg] = fopen( sfintvfn, opmod );
    if sfintf == -1
      error( 'openfiles.m: could not open file %s\n  %s\n', sfintvfn, msg );
    end
    fprintf( sfintvf, ...
             '%% t step Fintvx(1) Fintvy(1) ... Fintvx(n) Fintvy(n)\n' );
  end

end % if verbose files

