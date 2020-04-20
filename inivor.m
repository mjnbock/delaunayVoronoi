%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% inivor3.m
%%%%%%%%
%%%%%%%% compute circular voronoi tesselation
%%%%%%%% ATTENTION: all quantities defined and explained in main.m
%%%%%%%% BEWARE: INDEXING and VERTEX MAGIC is NONTRIVIAL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TODO remove superfluous mapping arrays

%%%%%%%% set all counters and neighbor indices to zero
pm = 0;                  % # pairs
tm = 0;                  % # triples
vm = 0;                  % # vertices
wm = 0;                  % # margin vertices
am = 0;                  % # of Voronoi contact arcs
fm = 0;                  % # of free closure arcs
um = 0;                  % # separation vertices

p2w = zeros(nol*cm,2);   % margin vertices of pair o,oo
                         % zeroing needed to allow for digestion
p2tm = zeros(1,nol*cm);  % # tripels containing pair
p2vm = zeros(1,nol*cm);  % # vertices of pair

c2pm = zeros(1,cm);      % # of pairs containing cell
c2am = zeros(1,cm);      % # of arcs of cell

nbm  = zeros(1,cm);      % # of Voronoi neighbors
nbdm = 0;                % # of Voronoi neighbor pairs

swl = zeros(nol,2);      % sorted wertex list for free margins

fcm = 0;                 % # of free cells
mcm = 0;                 % # of marginal cells

vxclm = 0;               % # of candidate vertices
wxclm = 0;               % # of candidate wertices

vxlm = 0;                % # of vertices actually needed in partition
wxlm = 0;                % # of wertices actually needed in partition
%vxl = zeros(1,2*am);    % zeroed
%wxl = zeros(1,2*am);    % below



%%%%%%%% determine pairs of overlapping cells
for i = 1:cm
  jj = i+1:cm;
  distx = x(i) - x(jj);
  disty = y(i) - y(jj);
  dist = sqrt( distx.*distx + disty.*disty );
  cclose = i + find(  dist < Pmax*( w(i) + w(jj) )  )';
  sizcl = size(cclose);
  if sizcl(1) > 0
    for j = cclose
      pm = pm + 1;
      m  = pm;
      % i,j sorted in ascending order
      p2c(m,:) = [ i, j ];
      c2pm(i) = c2pm(i) + 1;
      c2pm(j) = c2pm(j) + 1;
      c2p( i,   c2pm(i) ) = m;
      c2p( j,   c2pm(j) ) = m;
    end
  end
end



%%%%%%%% calculate contact spheres
for m = 1:pm % loop over pairs

  i = p2c(m,1);
  j = p2c(m,2);
  if abs( w(i)-w(j) ) < mindist
    sgn =    1;
    cthm = w(j)/w(i);
  elseif w(i) > w(j)
    sgn  =   1;
    cthm = w(j)/w(i);
  else
    sgn  = - 1;
    cthm = w(i)/w(j);
  end
  dx = ( x(i) - x(j) )*sgn;
  dy = ( y(i) - y(j) )*sgn;
  x0 = ( x(i) + x(j) )/2;
  y0 = ( y(i) + y(j) )/2;
  dij = sqrt( dx^2 + dy^2 );
  delta = dij - r(i) - r(j);

  thij = atan2( dy, dx );
  if abs( w(i)-w(j) ) < mindist
    %
    % The direction of straight line starting from pair-local 
    % origin (x0,y0) is given by (-dy,dx) in global coords.
    %
    vdist = [ dx, dy ];  % record orientation of pair
    udist = vdist/sqrt( vdist*vdist' );
    slp = [ - udist(2), udist(1) ];
    Rij = NaN;
    p2s(m,:) = [ Rij, x0,y0, slp(1),slp(2), thij, delta ];
  else
    % the usual contact sphere
    thm = acos(cthm);        % numerically unstable at wi ~ wj
    wwp =     w(i)^2 + w(j)^2;       % abbreviations similar to
    wwm = ( - w(i)^2 + w(j)^2 )*sgn; % p and m in 2008-03-31 version
    Rij = dij*w(i)*w(j)*sgn/( w(i)^2-w(j)^2 );
    vMij = ShiftRot2d( wwp*dij/wwm/2,0, x0,y0, thij );
    Mijx = vMij(1);
    Mijy = vMij(2);
    p2s(m,:) = [ Rij, Mijx,Mijy, cthm,thm, thij, delta ];
  end

end % loop over pairs



% TODO debug
%printindices;
%plotcraw;
%input('press enter to continue');



%%%%%%%% determine tripels of overlapping cells and vertices 
for m = 1:pm % loop over pairs m
  for mm = m+1:pm % loop over pairs mm
    % due to sortedness comparing these indices is sufficient
    if p2c(m,2) == p2c(mm,1) % if pairs share common cell
      i = p2c( m,  1 ); % i,j,k sorted in ascending order
      j = p2c( mm, 1 ); % pair  m = i,j
      k = p2c( mm, 2 ); % pair mm = j,k
      m3 = intersect( c2p(i,1:c2pm(i)), c2p(k,1:c2pm(k)) ); % pair m3 = i,k
      sizm3 = size(m3); % sizm3 = [1 1] or [0 0], could be checked
      if sizm3 == [1,1] % check if actual tripel
        ii = p2c(m3,1);
        kk = p2c(m3,2);
        if ii ~= i || kk ~= k
          warning('warn inivor.m: strange tripel rejected');
          continue;
        end
      else
        continue;
      end

      % associating overlap tripel
      tm = tm + 1;
      n  = tm;
      t2p(n,:) = [ m, m3, mm ]; % indeed sorted due to loop construction
      t2c(n,:) = [ i, j, k ];
      p2tm( m) = p2tm( m) + 1;
      p2tm(mm) = p2tm(mm) + 1;
      p2tm(m3) = p2tm(m3) + 1;
      p2t(  m, p2tm( m) ) = n;
      p2t( mm, p2tm(mm) ) = n;
      p2t( m3, p2tm(m3) ) = n;

      % associating vertex
      ci = [ x(i), y(i), r(i), w(i) ];
      cj = [ x(j), y(j), r(j), w(j) ];
      ck = [ x(k), y(k), r(k), w(k) ];
      vrtx = vertex( p2s(m,:), p2s(m3,:), p2s(mm,:), ci,cj,ck, mindist );
      sizvrtx = size(vrtx); % 0, 1 or 2 points
      for nn = 1:sizvrtx(1)
        vdisti = [x(i),y(i)]-vrtx(nn,:);
        vdistj = [x(j),y(j)]-vrtx(nn,:);
        vdistk = [x(k),y(k)]-vrtx(nn,:);
        Pi = sqrt( vdisti*vdisti' )/w(i);
        Pj = sqrt( vdistj*vdistj' )/w(j);
        Pk = sqrt( vdistk*vdistk' )/w(k);
        if Pi < Pmax  &&  Pj < Pmax  && Pk < Pmax
          vm = vm + 1;
          v2v(vm,:) = [ vrtx(nn,1),vrtx(nn,2), Pi ];
          v2c(vm,:) = [ i, j, k ];
          p2vm( m) = p2vm( m) + 1;
          p2vm(mm) = p2vm(mm) + 1;
          p2vm(m3) = p2vm(m3) + 1;
          p2v(  m, p2vm( m) ) = vm;
          p2v( mm, p2vm(mm) ) = vm;
          p2v( m3, p2vm(m3) ) = vm;
        end
      end

    end % if pairs share common cell
  end % loop over pairs mm
end % loop over pairs m



%%%%%%%% determine margin vertices of pair
for m = 1:pm % loop over pairs
  i = p2c(m,1);
  j = p2c(m,2);
  si = [ Pmax*w(i), x(i),y(i) ];
  sj = [ Pmax*w(j), x(j),y(j) ];
  pts = wertex( si, sj, mindist );
  sizpts = size(pts);
  if check_digestion == 1 && sizpts(2) < 2 && ...
      ~( Tmito(i)>-dt/2 && Tmito(i)<Tmito1(i)-dt/2 && ...
         abs(Tmito(i)-Tmito(j))>dt/3 || abs(Tmito1(i)-Tmito1(j))>dt/3 )
    error('stop inivor.m: cells overlapping w/o wertex, digestion detected');
  end
  if sizpts(1)==2
    for o=1:2
      pto = pts(o,:);
      wm = wm + 1;
      w2w(wm,:) = [ pto(1),pto(2), Pmax ];
      p2w( m,o) = wm;
    end
  end % else no wertices, even when touching
end



% TODO debug
%figure(floor(vfig)); % octave hack
%plotcraw;
%printindices;
%input('press enter to continue');



%%%%%%%% determine Voronoi arcs and neighbors
for m = 1:pm % loop over pairs

  % initial variables, zero counters
  i = p2c(m,1);
  j = p2c(m,2);
  %printf('pair %i, cells %i,%i\n', m, i, j );
  Rij = p2s(m,1);
  thmax = p2s(m,5);
  thij  = p2s(m,6);
  vclm  = 0; % # vertex candidates of pair
  vlm   = 0; % # w/vertices of pair
  neighbor = 0; % is pair neighbor pair?


  % build neighbor candidate list of pair i,j from their tripels
  nn = p2t(m,1:p2tm(m));
  kk = vec(t2c(nn,:))';
  kk = unique(kk);
  if runfromoct == 1
    cnbc = complement( [i,j], kk );    % cnbc common neighbor candidates
  else
    cnbc = setdiff( kk, [i,j] );       % order indeed reversed
  end


  % build vertex candidate list (trading redundancy in code for performance)
  % 
  % vertex types:
  % -------------
  % 0:   margin     wertex, arising from intersection of Pmax-spheres
  % 1:  ordinary    vertex, situated at triple junctions of cells
  % 2: connection u-vertex, needed for contact circle connection at +/-pi
  % 3: connection w-vertex, needed for   free margin  connection at +/-pi
  %
  if isnan(Rij) % flat contact surface
    vx0 = p2s(m,2:3);
    for no=1:p2vm(m) % vertices
      vclm = vclm + 1;
      n = p2v(m,no);
      vrtx  = v2v(n,1:2);
      Pvrtx = v2v(n, 3 );
      vcl(vclm,2:6) = [ Pvrtx, vrtx(1),vrtx(2), n,1 ];
    end
    for no=1:2 % margin wertices
      o = p2w(m,no);
      if o ~= 0 % checking for digestion earlier
        vclm = vclm + 1;
        wrtx  = w2w(o,1:2);
        Pwrtx = w2w(o, 3 ); % = Pmax
        vcl(vclm,2:6) = [ Pwrtx, wrtx(1),wrtx(2), o,0 ];
      end
    end
    %
    % The connection u-vertex is needed to decide
    % where to draw the arcs comprising the contact border.
    % By separating at the critical angle theta = +/- pi 
    % (think atan2), all angle ranges are contiguous. cf. B152
    % Later on this u-vertex will be deleted, unless it forms
    % the only vertex pair of the arc, like in digesting systems.
    %
    if abs( x(i) - x(j) ) > mindist && x(i) > x(j)
      ux0 = p2s(m,4:5);
      urtx = intlip( vx0,ux0, [x(i),y(i)],[1,0], mindist);
      sizurtx = size(urtx);
      if sizurtx(1) == 1
        vdist = [x(i),y(i)] - urtx;
        Purtx = sqrt( vdist*vdist' )/w(i);
        if Purtx < Pmax+mindist/2 % include in touching case
          vclm = vclm + 1;
          vcl(vclm,:) = [ -pi, Purtx, urtx(1),urtx(2), no,2 ];
          vclm = vclm + 1;
          vcl(vclm,:) = [  pi, Purtx, urtx(1),urtx(2), no,2 ];
        end
      end
    end
  else % spherical contact surface
    vMij = p2s(m,2:3);
    for no=1:p2vm(m)
      vclm = vclm + 1;
      n = p2v(m,no);
      vrtx  = v2v(n,1:2);
      Pvrtx = v2v(n, 3 );
      vcl(vclm,2:6) = [ Pvrtx, vrtx(1),vrtx(2), n,1 ];
    end
    for no=1:2
      o = p2w(m,no);
      if o ~= 0 % checking for digestion earlier
        vclm = vclm + 1;
        wrtx  = w2w(o,1:2);
        Pwrtx = w2w(o, 3 ); % = Pmax
        vcl(vclm,2:6) = [ Pwrtx, wrtx(1),wrtx(2), o,0 ];
      end
    end
    % cf. urtx comment at planar contact above
    urtx = vMij - Rij*[1,0];
    vdist = [x(i),y(i)] - urtx;
    Purtx = sqrt( vdist*vdist' )/w(i);
    if Purtx < Pmax+mindist/2 % include in touching case
      vclm = vclm + 1;
      vcl(vclm,:) = [ -pi, Purtx, urtx(1),urtx(2), no,2 ];
      vclm = vclm + 1;
      vcl(vclm,:) = [  pi, Purtx, urtx(1),urtx(2), no,2 ];
    end
  end % if flat contact surface


  % vertex erosion: find vertices which are not closer to other cells
  for no = 1:vclm
    Pi   = vcl(no,2);
    vrtx = vcl(no,3:4);
    nnoo = vcl(no,5);
    vtyp = vcl(no,6);
    if vtyp == 1 % avoid comparison of equal P-distances
      ijk = v2c(nnoo,:);
      if runfromoct == 1
        ccnbc = complement( ijk, cnbc ); % cnbc common neighbor candidates
      else
        ccnbc = setdiff( cnbc, ijk ); % order indeed reversed
      end
    else
      ccnbc = cnbc;
    end
    bbreak = 0;
    for k = ccnbc
      vdistk = [x(k),y(k)]-vrtx;
      Pk = sqrt(vdistk*vdistk')/w(k);
      if Pk < Pi % vertex closer to other cell
        bbreak = 1;
        break;
      end
    end
    if bbreak == 1 % disregard vertex
      continue;
    else % associate vertex
      if isnan(Rij) 
        if vtyp == 0  ||  vtyp == 1
          vdist = vrtx - [x(i),y(i)];
          vcl(no,1) = atan2( vdist(2), vdist(1) );
        end
      else
        if vtyp == 0  ||  vtyp == 1
          vdist = vrtx - vMij;
          vcl(no,1) = atan2( vdist(2), vdist(1) );
        end
      end
      vlm = vlm + 1;
      vl(vlm,:) = vcl(no,:);
    end
  end
  

  % remove multiply appearing vertices for degenerate case
  if mod(vlm,2) ~= 0
    vclm = 0;
    for no = 1:vlm
      vno = vl(no,3:4);
      vtyp = vl(no,6);
      if vtyp == 2 % urtx always needed
        vclm = vclm + 1;
        vcl(vclm,:) = vl(no,:);
      else
        solo = 1;
        for on = no+1:vlm
          von = vl(on,3:4);
          vdv = vno - von;
          dv = sqrt( vdv*vdv' );
          if dv < mindist
            solo = 0;
          end
        end
        if solo == 1;
          vclm = vclm + 1;
          vcl(vclm,:) = vl(no,:);
        end
      end
    end
    vlm = vclm;
    vl(1:vlm,:) = vcl(1:vclm,:);
  end
  [vvl,no2] = sort( vl(1:vlm,1) ); % sort by thv
  svl(1:vlm,:) = vl(no2,:);


  % determine and associate arcs
  if mod(vlm,2) ~= 0
    error('stop inivor.m: odd number of vertices at i=%i, j=%i',i,j);
    %warning('stop inivor.m: odd number of vertices at i=%i, j=%i',i,j);
    %input('press enter to continue');
  end

  vM = p2s(m,2:3);
  vtyp = svl(1,6);
  if vlm >= 4 && vtyp == 2 % first reconnect splitted arc
    th1 = svl( vlm-1,  1  ); % since first and last one are connection uertices
    th2 = svl(   2  ,  1  ) + 2*pi;
    v1   = svl( vlm-1, 3:4 );
    v2   = svl(   2,   3:4 );
    v1no = svl( vlm-1,  5  );
    v2no = svl(   2,    5  );
    v1t  = svl( vlm-1,  6  );
    v2t  = svl(   2,    6  );
    if check_thmax == 1
      dth1 = th1 - thij;
      dth2 = th2 - thij;
      if abs(dth1) > thmax || abs(dth2) > thmax
        error('stop inivor.m: starlikeness violated');
      end
    end
    neighbor = 1;
    am = am + 1;
    a2a(am,:) = [ m, Rij, vM(1),vM(2), th1,th2, v1(1),v1(2), v2(1),v2(2), ...
                                                  v1t,v2t,    v1no,v2no       ];
    c2am(i) = c2am(i) + 1;
    c2am(j) = c2am(j) + 1;
    c2a( i,   c2am(i) ) = am;
    c2a( j,   c2am(j) ) = am;
    nbm(i) = nbm(i) + 1;
    nbm(j) = nbm(j) + 1;
    nb( i,   nbm(i) ) = j;
    nb( j,   nbm(j) ) = i;
    if verbose_print == 1  ||  verbose_plot == 1
      if v1t == 0
        wxclm = wxclm + 1;
        wxcl(wxclm) = v1no;
      else
        vxclm = vxclm + 1;
        vxcl(vxclm) = v1no;
      end
      if v2t == 0
        wxclm = wxclm + 1;
        wxcl(wxclm) = v2no;
      else
        vxclm = vxclm + 1;
        vxcl(vxclm) = v2no;
      end
    end
    noran = 3:2:vlm-2; % so uertex is deleted
  else
    noran = 1:2:vlm;
  end

  for no=noran % then proceed to associate usual arcs
    v1   = svl(no  ,3:4);
    v2   = svl(no+1,3:4);
    v1no = svl(no,   5 );
    v2no = svl(no+1, 5 );
    v1t  = svl(no  , 6 );
    v2t  = svl(no+1, 6 );
    if isnan(Rij)
      th1 = thij;
      th2 = thij;
    else
      th1 = svl(no  ,1);
      th2 = svl(no+1,1);
    end % if isnan(Rij)
    if check_thmax == 1
      dth1 = th1 - thij;
      dth2 = th2 - thij;
      if abs(dth1) > thmax || abs(dth2) > thmax
        error('stop inivor.m: starlikeness violated');
      end
    end
    neighbor = 1;
    am = am + 1;
    % Rij = Nan if planar contact
    a2a(am,:) = [ m, Rij, vM(1),vM(2), th1,th2, v1(1),v1(2), v2(1),v2(2), ...
                                                  v1t,v2t,    v1no,v2no       ];
    c2am(i) = c2am(i) + 1;
    c2am(j) = c2am(j) + 1;
    c2a( i,   c2am(i) ) = am;
    c2a( j,   c2am(j) ) = am;
    nbm(i) = nbm(i) + 1;
    nbm(j) = nbm(j) + 1;
    nb( i,   nbm(i) ) = j;
    nb( j,   nbm(j) ) = i;
    if verbose_print == 1  ||  verbose_plot == 1
      if v1t == 0
        wxclm = wxclm + 1;
        wxcl(wxclm) = v1no;
      else
        vxclm = vxclm + 1;
        vxcl(vxclm) = v1no;
      end
      if v2t == 0
        wxclm = wxclm + 1;
        wxcl(wxclm) = v2no;
      else
        vxclm = vxclm + 1;
        vxcl(vxclm) = v2no;
      end
    end
  end
  if neighbor == 1
    nbdm = nbdm + 1;
    nbd(nbdm) = m;
  end

end % loop over pairs



%%%%%%%% find free boundary arcs
for i = 1:cm
  atmargin = 0;
  posi = [x(i),y(i)]; % cell center

  
  % erosion for margin connection wertex at phi=+/-pi
  posw  = [x(i),y(i)] - [Pmax*r(i),0]; % ie. phi=+/-pi
  marginw  = 1;
  for k = 1:nbm(i)
    j = nb(i,k);
    vdist = posw - [x(j),y(j)];
    Pj = sqrt(vdist*vdist')/w(j);
    if Pj < Pmax-mindist/2 % marginw=1 in touching case
      marginw = 0;
    end
  end

  
  % construction of wertex list
  wl = zeros(nol,2);
  wlm = 0;

  for p=1:c2am(i)
    alc = a2a( c2a(i,p), : );
    if alc(11)==0
      v1 = alc(7:8);
      vdist1 = v1 - posi;
      phi = atan2( vdist1(2), vdist1(1) );
      wlm = wlm + 1;
      wl(wlm,:) = [phi,0];
    end
    if alc(12)==0
      v2 = alc(9:10);
      vdist2 = v2 - posi;
      phi = atan2( vdist2(2), vdist2(1) );
      wlm = wlm + 1;
      wl(wlm,:) = [phi,0];
    end
  end
  [wwl,oo] = sort( wl(1:wlm,1) );
  swl(1:wlm,:) = wl(oo,:);

  if marginw==1 % cf. erosion above
    wlm = wlm + 1;
    swl(2:wlm,:) = wl(oo,:); 
    swl(1,:) = [-pi,3]; % angle orientation phi and margin wertex type = 3
    wlm = wlm + 1;
    swl(wlm,:) = [ pi,3];
  end

 
  % handle special case of full circle free margin
  if wlm==2 && swl(1,2)==3
    if swl(2,2)~=3
      error('stop inivor.m: unpaired margin wertices at %i',i);
    end
    posww = [x(i),y(i)] + [Pmax*r(i),0];
    marginww = 1;
    for k = 1:nbm(i)
      j = nb(i,k);
      vdist = posww - [x(j),y(j)];
      Pj = sqrt(vdist*vdist')/w(j);
      if Pj < Pmax-mindist/2 % marginww=1 in touching case
        marginww = 0;
      end
    end
    if marginww==1
      phi1 = -pi;
      phi2 =  pi;
      fm = fm + 1;
      f2f(fm,:) = [ i, Pmax*w(i), x(i),y(i), phi1,phi2 ];
      atmargin = 1;      
      mcm = mcm + 1;
      mc(mcm) = i;
      fcm = fcm + 1;
      fc(fm) = i;
    end
    continue;
  end

 
  % safety checks
  if mod(wlm,2)~=0
    error('stop inivor.m: uneven number of type 0 wertices at i=%i',i);
  end
  

  % determine free margin arc with type 3 wertices
  if wlm==0
    ooran = [];
  elseif wlm<2
    error('stop inivor.m: unexpected number of wertices wlm=%i',wlm);
  else
    if swl(1,2)==3 % connect arc at phi=+/-pi
      if wlm<4
        error('stop inivor.m: full free circle expected at %i',i);
      end
      if swl(wlm,2)~=3
        error('stop inivor.m: unpaired margin wertices at %i',i);
      end
      phi1 = swl(wlm-1,1);
      phi2 = swl(  2  ,1) + 2*pi;
      fm = fm + 1;
      f2f(fm,:) = [ i, Pmax*w(i), x(i),y(i), phi1,phi2 ];
      atmargin = 1;
      ooran = 3:2:wlm-2;      
    else
      ooran = 1:2:wlm;
    end
  end

  
  % determine and associate regular free margin arcs
  for o1 = ooran
    o2 = o1 + 1;
    phi1 = swl(o1,1);
    phi2 = swl(o2,1);
    fm = fm + 1;
    f2f(fm,:) = [ i, Pmax*w(i), x(i),y(i), phi1,phi2 ];
    atmargin = 1;
  end

  if atmargin == 1
    mcm = mcm + 1;
    mc(mcm) = i;
  end
end



%%%%%%%% build internal vertex and external wertex list
vxl = zeros(1,2*am);
wxl = zeros(1,2*am);
%vxlm = 0; % zeroed
%wxlm = 0; % above
for p=1:am
  v1t = a2a(p,11);
  v2t = a2a(p,12);
  if v1t==1
    vxlm = vxlm + 1;
    vxl(vxlm) = a2a(p,13);
  end
  if v2t==1
    vxlm = vxlm + 1;
    vxl(vxlm) = a2a(p,14);
  end
  if v1t==2
    wxlm = wxlm + 1;
    wxl(wxlm) = a2a(p,13);
  end
  if v2t==2
    wxlm = wxlm + 1;
    wxl(wxlm) = a2a(p,14);
  end
end

dum = unique(vxl(1,1:vxlm));
vxl = dum;
vxlm = length(dum);

dum = unique(wxl(1,1:wxlm));
wxl = dum;
wxlm = length(dum);



% TODO: debug
%plotcraw;
%plotcells;
%a2a
%input('press enter to continue');

