%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% printindices.m
%%%%%%%%
%%%%%%%% debug: print out indices of inital voronoi tesselation
%%%%%%%% ATTENTION: all quantities defined and explained in main.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% pair
fprintf( 1, 'overlap pair: [ cells ]\n' );
for m = 1:opm
  fprintf( 1, '%i: [ %i %i ]\n', m, op2c(m,1), op2c(m,2) );
end
fprintf( 1, 'cell: [ overlap pairs of cell ]\n' );
for i = 1:cm
  fprintf( 1, '%i: [ ',i );
  for m = 1:c2opm(i)
    fprintf(1, '%i ', c2op(i,m) );
  end
  fprintf( 1, ']\n' );
end
fprintf( 1, '\n\n' );

