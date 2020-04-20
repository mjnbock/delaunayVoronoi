%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Save and quit
if savdat == 1
  if save_fin == 1

    % open files
    sffn = strcat( save_path, save_fin_file ); % cell data
    [sff, msg] = fopen( sffn, 'wt' );
    if sff == -1
      error( 'main.m: could not open file\n  %s\n', msg );
    end
    sfstfn = strcat( save_path, save_fst_file ); % rng state and time
    [sfstf, msg] = fopen( sfstfn, 'wt' );
    if sfstf == -1
      error( 'main.m: could not open file\n  %s\n', msg );
    end

    % save final cell data
    for i = 1:cm
      fprintf( sff, '%13.10g %13.10g %13.10g %13.10g ',...
                       x(i),   y(i),   r(i),   w(i)      );
      fprintf( sff, '%3i %5i ', typ(i), Tp(i) );
      fprintf( sff, '%13.10g %13.10g %13.10g %13.10g %13.10g\n',...
                      rho(i),  bi(i),  cd(i),  my(i),  ps(i)      );
    end
    if plotting == 1
      filnam = strcat( movpath, movprfx, 'final.eps' );
      print( '-depsc', filnam )
    end

    % final rng state and time
    state = randn('state');
    fprintf( sfstf, '%% time, counter and final state of rng\n' );
    fprintf( sfstf, 't        = %i;\n',   t     );
    fprintf( sfstf, 'simstep  = %i;\n', simstep );
    fprintf( sfstf, 'state = [ ' );
    sizstat = size(state);
    for stind = 1:sizstat(1)-1
      fprintf( sfstf, '%i; ', state(stind) );
    end
    fprintf( sfstf, '%i ];\n', state(sizstat(1)) );

  end
end
closefiles;

