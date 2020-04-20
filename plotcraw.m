%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotcraw.m
%
% plots cell configuration
% ATTENTION: all quantities defined and explained in main.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clf;
hold on;



% cell centers bodies closures
if plot_cell_body == 1
  for i = 1:cm
    pl1 = arc( r(i), [x(i),y(i)], [-pi,pi], npts_cell_body, minang );
    line( pl1(1,:), pl1(2,:), 'Color', col_cell(1) );
  end
end
if plot_cell_center == 1 || plot_cell_centnu == 1
  for i = 1:cm
    cnum = num2str(i);
    text( x(i), y(i), cnum, 'Color', col_cell(1) );
  end
end
if plot_Pmax_sphere == 1
  for i = 1:cm
    pl2 = arc( Pmax*r(i), [x(i),y(i)], [-pi,pi], npts_Pmax_sphere, minang );
    line( pl2(1,:), pl2(2,:), 'Color', col_Pmax_sphere );
  end
end



% cell-cell contacts
if plot_contact_arc == 1
  for m = 1:pm
    i = p2c(m,1);
    j = p2c(m,2);
    Rij = p2s(m,1);
    if isnan( Rij )
      x0 = p2s(m,2);
      y0 = p2s(m,3);
      dx = p2s(m,4);
      dy = p2s(m,5);
      ndir  = [ dx, dy ];
      shift = ndir*Pmax*( r(i) + r(j) );
      px = [ x0-shift(1), x0+shift(1) ];
      py = [ y0-shift(2), y0+shift(2) ];
      line( px,py, 'Color',col_contact_arc );
    else
      vMij = p2s(m,2:3);
      thm = p2s(m,5);
      Dphi = p2s(m,6);
      vang = [ Dphi-thm, Dphi+thm ];
      vang = [ -pi, pi ];
      pl = arc( Rij, vMij, vang, npts_contact_arc, minang );
      line( pl(1,:),pl(2,:), 'Color',col_contact_arc );
    end
  end
end
if plot_arc_pair == 1
  for m = 1:pm
    pnum = num2str(m);
    i = p2c(m,1);
    j = p2c(m,2);
    Rij = p2s(m,1);
    if isnan( Rij )
      x0 = p2s(m,2);
      y0 = p2s(m,3);
      text( x0,y0, pnum, 'Color',col_arc_pair );
    else
      vMij = p2s(m,2:3);
      Dphi = p2s(m,6);
      udir = [ cos(Dphi), sin(Dphi) ];
      pt = vMij + Rij*udir;
      text( pt(1),pt(2), pnum, 'Color',col_arc_pair );
    end
  end
end
if plot_Mij == 1
  for m = 1:pm
    tag = strcat( 'M', num2str(m) );
    vMij = p2s(m,2:3);
    text( vMij(1),vMij(2), tag, 'Color',col_Mij );
  end
end



% vertices, wertices
if plot_vertex == 1
  for n = 1:vm
    vnum = num2str(n);
    text( v2v(n,1), v2v(n,2), vnum, 'Color', col_vertex );
  end
end
if plot_wertex == 1
  for o = 1:wm
    wnum = num2str(o);
    text( w2w(o,1), w2w(o,2), wnum, 'Color', col_wertex );
  end
end



% scaling
if plot_Mij == 1
  xx = [ p2s(1:pm,2); x(:) ];
  yy = [ p2s(1:pm,3); y(:) ];
else
  xx = x;
  yy = y;
end

xmin = min(xx);
xmax = max(xx);
ymin = min(yy);
ymax = max(yy);
rmax = max(r)*Pmax;

axis( plotscale( xmin, xmax, ymin, ymax, rmax, asprat, fig ) );



% finishing
if runfromoct==1
  octplot_command('redraw');
else
  drawnow;
end


