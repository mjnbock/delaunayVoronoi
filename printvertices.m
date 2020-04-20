%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% printvertices.m
%%%%%%%%
%%%%%%%% debug: print out verices
%%%%%%%% ATTENTION: all quantities defined and explained in main.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% vertex canditates from t2v
fprintf('v2v vertex n: [ vx, vy ]\n');
for n = 1:tm
  fprintf( stdout, '%i: [ %.16f %.16f ]\n', n, v2v(n,1),v2v(n,2) );
end
fprintf( stdout, '\n');

% wertices from w2w
fprintf('w2w wertex o: [ wx, wy ]\n');
for o = 1:wm
  fprintf( stdout, '%i: [ %.16f %.16f ]\n', o, w2w(o,1),w2w(o,2) );
end
fprintf( stdout, '\n');

fprintf( stdout, '\n' );

