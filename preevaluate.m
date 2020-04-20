%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% preevaluate.m
%%%%%%%%
%%%%%%%% calculates geometry of cells/tissue while running simulation
%%%%%%%% ATTENTION: all quantities defined and explained in main.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% compute length of free cell boundaries
Gam0 = zeros(1,cm);
for q = 1:fm
  i    = f2f(q,1);
  Rmax = f2f(q,2);
  vph  = f2f(q,5:6);
  Gam0(i) = Gam0(i) + Rmax*( vph(2) - vph(1) );
end



% calculate cell area, perimeter, contact lengths
Ain  = zeros(1,cm);
Aout = zeros(1,cm);
Gamp = zeros(1,cm);
for p = 1:am
  m = a2a(p,1);
  i = p2c(m,1);
  j = p2c(m,2);
  if w(i) > w(j)
    ib = i;
    is = j;
  else
    ib = j;
    is = i;
  end
  vxb = [ x(ib), y(ib) ];
  vxs = [ x(is), y(is) ];
  Rij = a2a(p,2);
  vv1 = a2a(p,7:8);  % start vertex coords
  vv2 = a2a(p,9:10); % stop  vertex coords

  %if i==3 || j==3 % debug
  %  fprintf( stderr, 'contact arc #%i, cells i=%i,j=%i\n', q, i,j );
  %  if ib==3
  %    fprintf( stderr, '  triangle=[%.2f,%.2f],[%.2f,%.2f],[%.2f,%.2f]\n',...
  %                        vxb(1),vxb(2), vv1(1),vv1(2), vv2(1),vv2(2) );
  %    fprintf( stderr, '  Ain old: %f\n', Ain(3)-artri(vxb,vv1,vv2) );
  %  else
  %    fprintf( stderr, '  triangle=[%.2f,%.2f],[%.2f,%.2f],[%.2f,%.2f]\n',...
  %                        vxs(1),vxs(2), vv1(1),vv1(2), vv2(1),vv2(2) );
  %    fprintf( stderr, '  Ain old: %f\n', Ain(3)-artri(vxs,vv1,vv2) );
  %  end
  %  fprintf( stderr, '  Ain+tri: %f\n', Ain(3) );
  %end % debug

  % handle planar <-> spherical contacts
  if isnan(Rij)
    % cell area related to contact surface  given by two triangles
    Ain(ib) = Ain(ib) + artri(vxb,vv1,vv2);
    Ain(is) = Ain(is) + artri(vxs,vv1,vv2);

    % vertex coords set before
    vdist = vv2 - vv1;

    LGam = sqrt( vdist*vdist' );
    Gamp(i) = Gamp(i) + LGam;
    Gamp(j) = Gamp(j) + LGam;
    % area already fully accounted
    %if i==3 || j==3 % debug
    %  fprintf( stderr, '  --> planar contact\n' );
    %  fprintf( stderr, '  Gamp old: %f\n', Gamp(3)-LGam );
    %  fprintf( stderr, '  Gamp new: %f\n', Gamp(3) );
    %end % debug
  else
    vMij   = a2a(p,3:4);  % vMij
    vtheta = a2a(p,5:6);  % [thm,thp] start + stop angle
    dtheta = vtheta(2) - vtheta(1);

    % see D59, holds for both starlike and non-starlike arcs
    Ain(ib) = Ain(ib) + artri(vxb,vMij,vv1) + artri(vxb,vMij,vv2);
    Ain(is) = Ain(is) - artri(vxs,vMij,vv1) - artri(vxs,vMij,vv1);

    Acake = Rij^2*dtheta/2.0;

    Ain(ib) = Ain(ib) - Acake;
    Ain(is) = Ain(is) + Acake;

    LGam = Rij*dtheta;
    Gamp(i) = Gamp(i) + LGam;
    Gamp(j) = Gamp(j) + LGam;

    %if i==3 || j==3 % debug
    %  fprintf( stderr, '  --> circular contact\n' );
    %  if ib==3
    %    fprintf( stderr, '  Ain+tri:  %f\n', Ain(3)+Asphcap );
    %    fprintf( stderr, '  Ain new:  %f\n', Ain(3) );
    %  else
    %    fprintf( stderr, '  Ain+tri:  %f\n', Ain(3)-Asphcap );
    %    fprintf( stderr, '  Ain new:  %f\n', Ain(3) );
    %  end
    %  fprintf( stderr, '  Gamp old: %f\n', Gamp(3)-LGam );
    %  fprintf( stderr, '  Gamp new: %f\n', Gamp(3) );
    %end % debug
  end
end


for q = 1:fm
  i = f2f(q,1);
  R = f2f(q,2);
  %vM = f2f(q,3:4);
  dphi = f2f(q,6) - f2f(q,5);
  Aout(i) = Aout(i) + R^2*dphi/2;
  %if i==3 % debug
  %  fprintf( stderr, 'free arc #%i, cell #%i\n', q, i );
  %  fprintf( stderr, '  Aout old: %f\n', Aout(i)-R^2*dphi/2 );
  %  fprintf( stderr, '  Aout new: %f\n', Aout(i) );
  %end % debug
end

