#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotcells.m
%
% plots cell configuration
% ATTENTION: all quantities defined and explained in main.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% scaling and preliminaries
clf;

if plot_axis == 0
  axis('off');
end
if runfromoct == 1
  rmax = max(r)*Pmax;
  if static_scale ~= 1
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
    winsiz = plotscale( xmin, xmax, ymin, ymax, rmax, asprat, fig );
    %plotscale( -3,3,  -4,4, 0, 6/8, fig );
  else
    axis( statscal );
    winsiz = statscal;
  end
else
  if static_scale ~= 1
    axis('equal');
  else
    winsiz = statscal;
  end
end
hold on;



% cell center symbols
if plot_cell_center == 1
  plot( x, y, sty_cell_center );
end



% delaunay triangulation: linking of neighbor cells
if plot_delaunay == 1
  for mm = 1:nbdm
    m = nbd(mm);
    i = p2c(m,1);
    j = p2c(m,2);
    plx = [x(i),x(j)];
    ply = [y(i),y(j)];
    line( plx,ply, 'Color',col_delaunay );
  end
end



% cell bodies, Pmax spheres
if plot_cell_body == 1
  for i = 1:cm
    pl1 = arc( r(i), [x(i),y(i)], [-pi,pi], npts_cell_body, minang );
    line( pl1(1,:), pl1(2,:), 'Color',char( col_cell(typ(i)) ), ...
                              'LineWidth', lw_cell_body            );
  end
end
if plot_Pmax_sphere == 1
  for i = 1:cm
    pl2 = arc( Pmax*r(i), [x(i),y(i)], [-pi,pi], npts_Pmax_sphere, minang );
    line( pl2(1,:), pl2(2,:), 'Color',col_Pmax_sphere );
  end
end



% cell-cell contacts
if plot_contact_arc == 1
  for p = 1:am
    Rij = a2a(p,2);
    if isnan( Rij )
      v1 = a2a(p,7:8);
      v2 = a2a(p,9:10);
      plx = [ v1(1), v2(1) ];
      ply = [ v1(2), v2(2) ];
      line( plx,ply, 'Color',col_contact_arc, 'LineWidth',lw_contact_arc );
    else
      Mij = a2a(p,3:4);
      thmin = a2a(p,5);
      thmax = a2a(p,6);
      pl = arc( Rij, Mij, [thmin,thmax], npts_contact_arc, minang );
      line( pl(1,:), pl(2,:), 'Color',col_contact_arc, ...
                              'LineWidth',lw_contact_arc  );
    end
  end
end
if plot_arc_pair == 1
  for p = 1:am
    m = a2a(p,1);
    mm = num2str(m);
    Rij = a2a(p,2);
    if isnan( Rij )
      v1 = a2a(p,7:8);
      v2 = a2a(p,9:10);
      pt = (v1+v2)/2;
      text( pt(1),pt(2), mm, 'Color',col_contact_arc );
    else
      vMij  = a2a(p,3:4);
      thmin = a2a(p,5);
      thmax = a2a(p,6);
      ang = thmin + (thmax-thmin)/2;
      pt = vMij + Rij*[ cos(ang), sin(ang) ];
      text( pt(1),pt(2), mm, 'Color',col_contact_arc );
    end
  end
end
if plot_Mij == 1
  for p = 1:am
    m = a2a(p,1);
    tag = strcat( 'M', num2str(m) );
    vMij = p2s(m,2:3);
    text( vMij(1),vMij(2), tag, 'Color',col_Mij );
  end
end



% free arcs of marginal cells
if plot_margin_arc == 1
  for q = 1:fm
    i   = f2f(q,1);
    R   = f2f(q,2);
    vM  = f2f(q,3:4);
    phi = f2f(q,5:6);
    pl = arc( R, vM, phi, npts_margin_arc, minang );
    line( pl(1,:), pl(2,:), 'Color',col_margin_arc, 'LineWidth',lw_margin_arc );
  end
end
if plot_margin_num == 1
  for q = 1:fm
    i   = f2f(q,1);
    R   = f2f(q,2);
    vM  = f2f(q,3:4);
    ang = f2f(q,5) + ( f2f(q,6) - f2f(q,5) )/2;
    pt = vM + R*[ cos(ang), sin(ang) ];
    ii = num2str(i);
    text( pt(1),pt(2), ii, 'Color',col_margin_num );
  end
end



% vertices
if plot_wertex == 1
  for o = 1:wm
    vnum = num2str(o);
    text( w2w(o,1), w2w(o,2), vnum, 'Color', col_wertex );
  end
end
if plot_vertex == 1
  for n = 1:tm
    vnum = num2str(n);
    text( t2v(n,1), t2v(n,2), vnum, 'Color', col_vertex );
  end
end



% forces
if plot_loc_force == 1
  for i = 1:cm
    vF = Floc(i,:)/norm_loc_force;
    plx = [ x(i), x(i)+vF(1) ];
    ply = [ y(i), y(i)+vF(2) ];
    line( plx,ply, 'Color',col_loc_force );
  end
end
if plot_pol_force == 1  &&  persistence >= 2
  for i = 1:cm
    vF = Fpol(i,:)/norm_pol_force;
    plx = [ x(i), x(i)+vF(1) ];
    ply = [ y(i), y(i)+vF(2) ];
    line( plx,ply, 'Color',col_pol_force );
  end
end
if plot_vis_force == 1  &&  persistence >= 3
  for i = 1:cm
    vF = Fvis(i,:)/norm_vis_force;
    plx = [ x(i), x(i)+vF(1) ];
    ply = [ y(i), y(i)+vF(2) ];
    line( plx,ply, 'Color',col_vis_force );
  end
end
if plot_int_forceh == 1
  for i = 1:cm
    vF = Finth(i,:)/norm_int_forceh;
    plx = [ x(i), x(i)+vF(1) ];
    ply = [ y(i), y(i)+vF(2) ];
    line( plx,ply, 'Color',col_int_forceh );
  end
end
if plot_int_forcev == 1
  for i = 1:cm
    vF = Fintv(i,:)/norm_int_forcev;
    plx = [ x(i), x(i)+vF(1) ];
    ply = [ y(i), y(i)+vF(2) ];
    line( plx,ply, 'Color',col_int_forcev );
  end
end
if plot_int_force == 1
  for i = 1:cm
    vF = Fint(i,:)/norm_int_force;
    plx = [ x(i), x(i)+vF(1) ];
    ply = [ y(i), y(i)+vF(2) ];
    line( plx,ply, 'Color',col_int_force );
  end
end
if plot_full_force == 1
  for i = 1:cm
    vF = ( Floc(i,:) + Fint(i,:) + Fpol(i,:) )/norm_full_force;
    plx = [ x(i), x(i)+vF(1) ];
    ply = [ y(i), y(i)+vF(2) ];
    line( plx,ply, 'Color',col_full_force );
  end
end



% polarity
if plot_polarity == 1
  for i = 1:cm
    vpol = pol(i,:)/norm_polarity;
    plx = [ x(i), x(i)+vpol(1) ];
    ply = [ y(i), y(i)+vpol(2) ];
    line( plx,ply, 'Color',col_polarity );
  end
end



% velocity
if plot_velocity == 1
  for i = 1:cm
    vvel = v(i,:)/norm_velocity;
    plx = [ x(i), x(i)+vvel(1) ];
    ply = [ y(i), y(i)+vvel(2) ];
    line( plx,ply, 'Color',col_velocity );
  end
end



% cell numbers
if plot_cell_centnu == 1
  for i = 1:cm
    cnum = num2str(i);
    text( x(i), y(i), cnum, 'Color',char(col_cell(typ(i))) );
  end
end



% time when doing movie
if runfromoct~=1
  winsiz = axis;
end
if movie==1
  tnum = sprintf(movtformat,movtime);
  xpos = winsiz(1) + ( winsiz(2) - winsiz(1) )*0.80;
  ypos = winsiz(3) + ( winsiz(4) - winsiz(3) )*0.05;
  text( xpos, ypos, tnum, 'Color','k' );
end



% finishing
if runfromoct==1
  axis(winsiz);
  octplot_command('redraw');
else
  axis('equal');
  drawnow;
end

