%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% evaluate.m
%%%%%%%%
%%%%%%%% data and observable evaluation while running simulation
%%%%%%%% ATTENTION: all quantities defined and explained in main.m
%%%%%%%%
% evaluates data and plots / makes movies during simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plotting and movie
if plotting==1 && ( mod(simstep,plotincr)==0 || mod(simstep,expdatincr)==0 )
  movtime = t/movtdiv;
  plotcells;

  if verbose_plot==1
    input('press enter to continue');
  end

  if  movie==1  &&  mod(simstep,plotincr)==0
    if runfromoct >= 1
      filnam = sprintf( '%s%s%.8i.png', movpath, movprfx, simstep );
      print( '-dpng', filnam )
    else
      filnam = sprintf( '%s%s%.8i.jpg', movpath, movprfx, simstep );
      print( '-djpeg', filnam )
    end
    if moveps==1
      filnam = sprintf( '%s%s%.8i.eps', movpath, movprfx, simstep );
      print( '-depsc', filnam )
    end
  end

  if moveps==2  &&  mod(simstep,expdatincr)==0
    filnam = sprintf( '%s%s%.8i.eps', movpath, movprfx, simstep );
    print( '-depsc', filnam )
  end
end



% data collection and saving
if save_sim == 1  &&  mod(simstep,datincr) == 0

  % # marginal cells
  nmargin = mcm/cm;
  nintern = 1 - nmargin;

  % # of neighbors
  nbrmea = mean( nbm ); % global # neighbors
  nbrmin =  min( nbm );
  nbrmax =  max( nbm );
  index = mc(1:mcm);
  nbrmmea = mean( nbm(index) ); % marginal # neighbors
  nbrmmin =  min( nbm(index) );
  nbrmmax =  max( nbm(index) );
  if runfromoct == 1 % TODO earlier, when building marginal cell list
    index = complement( mc(1:mcm), 1:cm );
  else
    index = setdiff( 1:cm, mc(1:mcm) );
  end
  if size(index) > 0
    nbrimea = mean( nbm(index) ); % internal # neighbors
    nbrimin =  min( nbm(index) );
    nbrimax =  max( nbm(index) );
  else
    nbrimea = nan;
    nbrimin = nan;
    nbrimax = nan;
  end

  % normalized pair distances
  mm = a2a(1:am,1);
  ii = p2c(mm,1);
  jj = p2c(mm,2);
  delta = p2s(mm,7) ./ ( r(ii) + r(jj) );
  dltmea = mean(delta);
  [ dltmin, mmin ] = min(delta);
  [ dltmax, mmax ] = max(delta);

  % PA-ratio
  Acells = Ain  + Aout;
  Pcells = Gamp + Gam0;
  PAcells = Pcells ./ sqrt(Acells*pi) / 2;

  % size of tissue
  % TODO actually there should be a cluster analysis; then the definition
  %      of the size of tissue is not so clear; which cluster to take?
  siztiss = 0;
  ist = 0;
  jst = 0;
  if cm == 1 % single cell
    siztiss = 2*Pmax*w(i);
    ist = 1;
    jst = 1;
  else
    for iind = 1:mcm-1
      i  = mc(iind);
      jj = mc(iind+1:mcm);
      xdist = - x(jj) + x(i);
      ydist = - y(jj) + y(i);
      dist = sqrt( xdist.^2 + ydist.^2 );
      [ distmax, j ] = max(dist);
      siztissc = distmax + Pmax*( w(i) + w(j) );
      if siztissc > siztiss
        siztiss = siztissc;
        ist = i;
        jst = j;
      end
    end
  end
  if ist==0 || jst==0
    error('stop evaluate.m: could not determine size of tissue');
  end

  % velocities
  avel = sqrt( v(:,1).^2 + v(:,2).^2 );
  velmea = mean(avel);
  velmin =  min(avel);
  velmax =  max(avel);

  % polarity
  apol = sqrt( pol(:,1).^2 + pol(:,2).^2 );
  polmea = mean(apol);
  polmin =  min(apol);
  polmax =  max(apol);


  % forces
  vF = [ Fint(:,1) + Floc(:,1) + Fpol(:,1) + Fvis(:,1), ...
         Fint(:,2) + Floc(:,2) + Fpol(:,2) + Fvis(:,2) ];
  aF = sqrt( vF(:,1).^2 + vF(:,2).^2 );
  aFmea = mean(aF);
  aFmin =  min(aF);
  aFmax =  max(aF);

  % save parameters
  % XXX update description in openfiles.m base_conf.m after change
  fprintf( sparf, '%10.f %10i ', t, simstep );
  fprintf( sparf, '%10.8f\n',  myalpha );

  % save observables
  % XXX update description in openfiles.m base_conf.m after change
  fprintf( ssf, '%10.f %10i ', t, simstep );
  fprintf( ssf, ' % 10.8f % 10.8f % 10.8f ', siztiss, nmargin, nintern );
  fprintf( ssf, ' % 10.8f % 10.8f % 10.8f ', nbrmea,  nbrmin,  nbrmax );
  fprintf( ssf, ' % 10.8f % 10.8f % 10.8f ', nbrmmea, nbrmmin, nbrmmax );
  fprintf( ssf, ' % 10.8f % 10.8f % 10.8f ', nbrimea, nbrimin, nbrimax );
  fprintf( ssf, ' % 10.8f % 10.8f % 10.8f ', dltmea,  dltmin,  dltmax );
  fprintf( ssf, ' % 10.8f % 10.8f % 10.8f ', velmea,  velmin,  velmax );
  fprintf( ssf, ' % 10.8f % 10.8f % 10.8f',  polmea,  polmin,  polmax );
  fprintf( ssf, '\n' );

  % save positions
  % XXX update description in openfiles.m base_conf.m after change
  fprintf( sxf, '%10.f %10i ', t, simstep );
  for i = 1:cm
    fprintf( sxf, ' % 10.8g % 10.8g ', x(i), y(i) );
  end
  fprintf( sxf, '\n' );

  % save velocities
  % XXX update description in openfiles.m base_conf.m after change
  fprintf( svf, '%10.f %10i ', t, simstep );
  for i = 1:cm
    fprintf( svf, ' % 10.8g % 10.8g % 10.8g ', v(i,1), v(i,2), avel(i) );
  end
  fprintf( svf, '\n' );
  
  % save areas
  % XXX update description in openfiles.m base_conf.m after change
  fprintf( sAf, '%10.f %10i ', t, simstep );
  for i = 1:cm
    fprintf( sAf, ' % 10.8g % 10.8g % 10.8g ', Ain(i)+Aout(i), Ain(i),Aout(i) );
  end
  fprintf( sAf, '\n' );

  % save perimeters
  % XXX update description in openfiles.m base_conf.m after change
  fprintf( sPerf, '%10.f %10i ', t, simstep );
  for i = 1:cm
    fprintf( sPerf, ' % 10.8g % 10.8g % 10.8g ', Gam0(i)+Gamp(i), Gamp(i),Gam0(i) );
  end
  fprintf( sPerf, '\n' );

  % save PA-ratios
  % XXX update description in openfiles.m base_conf.m after change
  fprintf( sPAf, '%10.f %10i ', t, simstep );
  for i = 1:cm
    fprintf( sPAf, ' % 10.8g', PAcells(i) );
  end
  fprintf( sPAf, '\n' );

   % save polarity vectors
  % XXX update description in openfiles.m base_conf.m after change
  fprintf( spf, '%10.f %10i ', t, simstep );
  for i = 1:cm
    fprintf( spf, ' % 10.8g % 10.8g % 10.8g ', pol(i,1), pol(i,2), apol(i) );
  end
  fprintf( spf, '\n' );
  
  % save forces
  % XXX update description in openfiles.m base_conf.m after change
  fprintf( sftotf, '%10.f %10i ', t, simstep );
  fprintf( sftotf, ' % 10.8g % 10.8g % 10.8g ', aFmea, aFmin, aFmax );
  for i = 1:cm
    fprintf( sftotf, ' % 10.8g % 10.8g % 10.8g ', vF(i,1),vF(i,2), aF(i) );
  end
  fprintf( sftotf, '\n' );
end


% expensive data collection and saving
if save_sim == 1  &&  mod(simstep,expdatincr) == 0
  % TODO actually there should be a cluster analysis; then the definition
  %      of the size of tissue is not so clear; which cluster to take?

  % out/in area ratio

  % total area and perimeter of tissue
  Aintiss  = sum(Ain);
  Aouttiss = sum(Aout);
  Atiss = Aintiss + Aouttiss;
  Ptiss = sum(Gam0); % length of free margin of tissue, from preevaluate.m

  % area of tissue margin, exterior to delaunay triangulation
  Acore = 0.0;
  for nn=1:vxlm
    n = vxl(nn);
    i = v2c(n,1);
    j = v2c(n,2);
    k = v2c(n,3);
    Acore = Acore + artri( [x(i),y(i)], [x(j),y(j)], [x(k),y(k)] );
  end
  Amarg = Atiss - Acore;

  % perimeter of delaunay triangulation, and interior vertex list
  Pcore = 0.0;
  for p=1:pm
    i = p2c(m,1);
    j = p2c(m,2);
    v1t = a2a(p,11);
    v2t = a2a(p,12);
    if v1t==0 || v2t==0
      vdist = [x(i),y(i)] - [x(j),y(j)];
      dist = sqrt( vdist*vdist' );
      Pcore = Pcore + dist;
    end
  end

  % area ratios
  AoAt = Aouttiss / Atiss;
  AmAt = Amarg / Atiss;

  % area/perimeter ratios
  PAtiss = Ptiss / sqrt( Atiss * pi ) / 2;
  if Acore > mindist^2
    PAcore = Pcore / sqrt( Acore * pi ) / 2;
  else
    PAcore = NaN; % because when Acore==0 then Pcore==0
  end

  % save areas, circumference, surface tensions 
  % XXX update description in openfiles.m base_conf.m after change
  fprintf( sexpf, '%10.f %10i ', t, simstep );
  fprintf( sexpf, ' % 10.8g % 10.8g ',  PAtiss, PAcore );
  fprintf( sexpf, ' % 10.8g % 10.8g ',  AmAt,   AoAt   );
  fprintf( sexpf, ' % 10.8g ',          Atiss );
  fprintf( sexpf, ' % 10.8g % 10.8g ',  Aintiss,    Acore  );
  fprintf( sexpf, ' % 10.8g % 10.8g ',  Aouttiss,   Amarg  );
  fprintf( sexpf, ' % 10.8g % 10.8g\n', Ptiss,      Pcore  );
end

