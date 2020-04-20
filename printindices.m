%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% printindices.m
%%%%%%%%
%%%%%%%% debug: print out indices of inital voronoi tesselation
%%%%%%%% ATTENTION: all quantities defined and explained in main.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% pair
fprintf( 1, 'pair: [ cells ]\n' );
for m = 1:pm
  fprintf( 1, '%i: [ %i %i ]\n', m, p2c(m,1), p2c(m,2) );
end
%fprintf( 1, 'pair: [ pseudo vertices ]\n' );
%for m = 1:pm
%  fprintf( 1, '%i: [ %i %i ]\n', m, p2w(m,1), p2w(m,2) );
%end
fprintf( 1, 'cell: [ pairs of cell ]\n' );
for i = 1:cm
  fprintf( 1, '%i: [ ',i );
  for m = 1:c2pm(i)
    fprintf(1, '%i ', c2p(i,m) );
  end
  fprintf( 1, ']\n' );
end

% tripel
fprintf( 1, '\ntripel: [ pairs ]\n' );
for n = 1:tm
  fprintf( 1, '%i: [ %i %i %i ]\n', n, t2p(n,1), t2p(n,2), t2p(n,3) );
end
fprintf( 1, 'tripel: [ cells ]\n' );
for n = 1:tm
  fprintf( 1, '%i: [ %i %i %i ]\n', n, t2c(n,1), t2c(n,2), t2c(n,3) );
end
%fprintf( 1, '\ntripel: [ pseudo vertices ]\n' );
%for n = 1:tm
%  fprintf( 1, '%i: [ %i %i %i ]\n', n, t2w(n,1), t2w(n,2), t2w(n,3) );
%end
fprintf( 1, 'cell: [ tripels of cell ]\n' );
for i = 1:cm
  fprintf( 1, '%i: [ ',i );
  for n = 1:c2tm(i)
    fprintf( 1, '%i ', c2t(i,n) );
  end
  fprintf( 1, ']\n' );
end
fprintf( 1, 'pair: [ tripels of pair ]\n' );
for m = 1:pm
  fprintf( 1, '%i: [ ', m );
  for n = 1:p2tm(m)
    fprintf( 1, '%i ', p2t(m,n) );
  end
  fprintf( 1, ']\n' );
end

fprintf( 1, '\n\n' );

