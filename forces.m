%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% forces.m
%%%%%%%%
%%%%%%%% compute neighbor pairs while running simulation
%%%%%%%% ATTENTION: all quantities defined and explained in main.m
%%%%%%%%
% computes forces on cells, length of boundaries
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TODO 2009-08-04: - influence of myosin
%                  - influence of integrin bi (per segment, area)
%                  - influence of cadherin cd (per arc, length)
%                  - influence of myosin   my (per segment, area)

% zeroing per time-step arrays
Fint  = zeros(cm,2);
if plot_int_forceh == 1  ||  plot_int_forcev == 1
  Finth = zeros(cm,2);
  Fintv = zeros(cm,2);
end
Floc = zeros(cm,2);
Fpol = zeros(cm,2);

if plot_vis_force == 1
  Fvis = zeros(cm,2);
end
pdrag = zeros(am,1);



% interaction force and |Gamma_ij|
if interaction==1  ||  interaction==2 || interaction==3
  for p = 1:am
    % TODO refine this check to include pathological configs with multiple 
    % contact arcs behind one another as seen from cell center
    % needed, if cells may be digesting
    % TODO make arcs/forces somehow digestion proof
    forcearc = 1;
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


    % ignore cell pairs during mitosis
    if cellcycle >= 1
      if Tmito(i)>-dt/2 && Tmito(i)<Tmito1(i)-dt/2 && ...
          abs(Tmito(i)-Tmito(j))>dt/3 || abs(Tmito1(i)-Tmito1(j))>dt/3
        continue;
      end
    end


    % handle planar <-> spherical contacts 
    if isnan(Rij) % planar
      vv1 = a2a(p,7:8);  % start vertex coords
      vv2 = a2a(p,9:10); % stop  vertex coords
      vdist = vv2 - vv1;
      dgamx = vdist(1)/fintnum;
      dgamy = vdist(2)/fintnum;
      dgam = sqrt( dgamx^2 + dgamy^2 );
    else % spherical
      vv1 = a2a(p,3:4);  % vMij
      vv2 = a2a(p,5:6);  % [thm,thp] start + stop angle
      if force_thmax == 1 % limit angles to be within thmax bound
        thij  = p2s(m,6);
        thmax = p2s(m,5);
        if vv2(2) > pi  &&  thij < 0
          thij = thij + 2*pi;
        end
        thmm = thij - thmax;
        thmp = thij + thmax;

        % thmax fixup
        if vv2(1) < thmm
          vv2(1) = thmm;
          if warn_thmax == 1
            warning('forces.m: thmax violated between %i,%i',i,j);
          end
          if check_thmax == 1 % should be trapped in inivor
            error('forces.m: thmax violated between %i,%i',i,j);
          end
        end
        if vv2(2) > thmp
          vv2(2) = thmp;
          if warn_thmax == 1
            warning('forces.m: thmax violated between %i,%i',i,j);
          end
          if check_thmax == 1 % should be trapped in inivor
            error('forces.m: thmax violated between %i,%i',i,j);
          end
        end
        
        % dgam calculation
        dgam = Rij * (vv2(2)-vv2(1))/fintnum;

        % so the arc is completely out of thmax bounds
        if vv2(2) < vv2(1) % arc w/o force generation
          forcearc = 0;
        end

      end
    end
    
    
    % limit force calculation to arcs actually generating force ...
    if forcearc == 1
      vdist = vxb - vxs;
      dist  = sqrt( vdist*vdist' );
      uh = vdist/dist;
      uv = [-vdist(2),vdist(1)]/dist;
      densities = arcdens( vxb,vxs, Rij, vv1,vv2,...
                           rho(ib),rho(is),r(ib),r(is),cd(ib),cd(is),...
                           uh,uv, fintnum,mindist );

      delta = p2s(m,7);
      if interaction==1 % simple, isotropic
        vF = vFint( delta,dcrit,dmin, r(ib),r(is), dgam,finth,fintv,...
                    uh,uv, mindist, densities );
        Fint(is,:)  = Fint(is,:) + vF(1,:);
        Fint(ib,:)  = Fint(ib,:) + vF(2,:);
      elseif interaction==2 % anisotropic for cells with type > 1
        vF = vFint( delta,dcrit,dmin, r(ib),r(is), dgam,finth,fintv,...
                    uh,uv, mindist, densities );
        if typ(is)>1
          Fint(is,:)  = Fint(is,:) + vF(1,:)*myalpha;
        else
          Fint(is,:)  = Fint(is,:) + vF(1,:);
        end
        if typ(ib)>1
          Fint(ib,:)  = Fint(ib,:) + vF(2,:)/myalpha;
        else
          Fint(ib,:)  = Fint(ib,:) + vF(2,:);
        end
      elseif interaction==3
        dcritib = dcrit;
        dcritis = dcrit;
        if typ(ib) > 1
          dcritib = dcrit / sqrt( myalpha^2 - (myalpha^2-1/myalpha^2)*uh(1)^2 );
        end
        if typ(is) > 1
          dcritis = dcrit / sqrt( myalpha^2 - (myalpha^2-1/myalpha^2)*uh(1)^2 );
        end
        mydcrit = ( dcritib + dcritis )/2;
        vF = vFint( delta,mydcrit,dmin, r(ib),r(is), dgam,finth,fintv,...
                    uh,uv, mindist, densities );
        Fint(is,:)  = Fint(is,:) + vF(1,:);
        Fint(ib,:)  = Fint(ib,:) + vF(2,:);
      end
      if plot_int_forceh == 1 || plot_int_forcev == 1 || save_vforce == 1
        Finth(is,:) = Finth(is,:) + vF(3,:);
        Fintv(is,:) = Fintv(is,:) + vF(4,:);
        Finth(ib,:) = Finth(ib,:) + vF(5,:);
        Fintv(ib,:) = Fintv(ib,:) + vF(6,:);
      end

      % ... also for viscosity, because pairing density rho
      % makes sense within thmax only
      if persistence == 3
        pdrag(m) = pairdrag( fvis, dgam, densities );
        if plot_vis_force == 1
          Fvis(i,:) = Fvis(i,:) - pdrag(m)*( v(i,:) - v(j,:) );
          Fvis(j,:) = Fvis(j,:) - pdrag(m)*( v(j,:) - v(i,:) );
        end
      end
    end
    
  end
end % if interaction==1 || interaction==2 || interaction==3



% modified interaction force based on arc tension and area repulsion
%Ain = zeros(1,cm); % computed in preevaluate.m
if interaction==4
  % cell areas Ain,Aout computed in preevaluate.m
  error('forces.m: interaction==4 not implemented');
end % if interaction == 4



% locomotion force
if locomotion==1
  for q = 1:fm
    i    = f2f(q,1);
    Rmax = f2f(q,2);
    vph  = f2f(q,5:6);
    vF = vFloc( r(i),rho(i),bi(i), Rmax,vph, floc,flocnum );
    Floc(i,:) = Floc(i,:) + vF;
  end
end



% polarity force
if persistence>=2
  for i = 1:cm
    Fpol(i,:) = vFpol( r(i),rho(i),bi(i),pol(i,:), fpol,fpolnum,kappol );
  end
end

