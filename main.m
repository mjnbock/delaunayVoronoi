ttmax = tmax;
% stuuupid hack, matlab somehow looses tmax variable


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% computing/loading initial values
initialize;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Simulate
while t <= ttmax
  if show_progress == 1
    shprg = sprintf('simstep=%i, t=%is', simstep, t );
    disp(shprg);
  end

  % growth and cell division
  if cellcycle==1
    grow;
  end

  % voronoi partition
  runvor;  % can be optimized
  if verbose_print == 1 % debug: printout
    printindices;
    printvertices;
  end
  if plotting == 1
    if verbose_plot == 1
      plotcraw;
      input('press enter to continue');
      plotcells;
      input('press enter to continue');
    end
  end

  % evaluate things needed for force calculation
  preevaluate;

  % compute forces on cells, length of boundaries
  if alpha < 0.0
    if exist('mycomp_alpha.m','file')==2
      myalpha = mycomp_alpha(t);
    else
      myalpha = comp_alpha(t);
    end
  else
    myalpha = alpha;
  end
  forces;
  %warning('checkpoint forces');

  
  % evaluate data, plot cells and make movie
  evaluate;
  %warning('checkpoint evaluate');

  
  % update velocities and positions of cells
  if persistence==0 % purely overdamped
    for i = 1:cm
      vF = Fint(i,:) + Floc(i,:);
      v(i,:) = vF/drag(i);
      x(i) = x(i) + v(i,1)*dt;
      y(i) = y(i) + v(i,2)*dt;
    end
  end

  if persistence==1 % phenomenological persistence time in velocity
    for i = 1:cm
      vvi = v(i,:);
      x(i) = x(i) + vvi(1)*dt;
      y(i) = y(i) + vvi(2)*dt;
      vF = Fint(i,:) + Floc(i,:);
      dvi = ( vF/drag(i) - vvi )*dt/Tp(i);
      if perturb_vel == 1
        if perturb_gsc == 1
          scalgam = kapgam*Gam0(i) + (1-kapgam)*Gamp(i);
          %scalgam = sqrt( (Gam0(i)+Gamp(i))/4/pi/Pmax/rbar );
          %scalgam = sqrt( (Gam0(i)+Gamp(i))/4/pi/Pmax/r(i) );
          scalgam = sqrt( scalgam/4/pi/Pmax/r(i) );
        else
          scalgam = 1/sqrt(2);
        end
        %dvi = dvi + evel*randn(1,2)*scalgam*sqrt( dt/Tp(i) );
        dvi = dvi + evel*randn(1,2)*ps(i)*scalgam*sqrt( dt/Tpbar );
      end
      v(i,:) = v(i,:) + dvi;
    end
  end

  if persistence==2 % phenomenological cell polarity vector
    for i = 1:cm
      vF = Fint(i,:) + Floc(i,:) + Fpol(i,:);
      v(i,:) = vF/drag(i);
      x(i) = x(i) + v(i,1)*dt;
      y(i) = y(i) + v(i,2)*dt;

      dpi = ( v(i,:)/gampv - pol(i,:) )*dt/Tpol;
      scalgam = kapgam*Gam0(i) + (1-kapgam)*Gamp(i);
      %scalgam = sqrt( (Gam0(i)+Gamp(i))/4/pi/Pmax/rbar );
      %scalgam = sqrt( (Gam0(i)+Gamp(i))/4/pi/Pmax/r(i) );
      scalgam = sqrt( scalgam/4/pi/Pmax/r(i) );
      dpi = dpi + epol*randn(1,2)*ps(i)*scalgam*sqrt(dt);

      polabs = sqrt( pol(i,1)^2 + pol(i,2)^2 );
      if polabs > polde*dt/Tpol
        dpi = dpi - polde*pol(i,:)/polabs * dt/Tpol;
      end

      pol(i,:) = pol(i,:) + dpi;
    end
  end

  if persistence==3 % phenomenological cell polarity + viscous pair interaction
    % compute velocities
    vFact = zeros(cm,2);
    mDrag = zeros(cm,cm);
    for i = 1:cm
      vFact(i,:) = Fint(i,:) + Floc(i,:) + Fpol(i,:);
      mDrag(i,i) = drag(i);
    end
    for p = 1:am
      m = a2a(p,1);
      i = p2c(m,1);
      j = p2c(m,2);
      mDrag(i,j) =   pdrag(m);
      mDrag(j,i) = - pdrag(m);
    end
    v(:,1) = mDrag\vFact(:,1);
    v(:,2) = mDrag\vFact(:,2);

    % increment coords and polarity vector
    for i = 1:cm
      x(i) = x(i) + v(i,1)*dt;
      y(i) = y(i) + v(i,2)*dt;

      dpi = ( v(i,:)/gampv - pol(i,:) )*dt/Tpol;
      scalgam = kapgam*Gam0(i) + (1-kapgam)*Gamp(i);
      scalgam = sqrt( scalgam/4/pi/Pmax/r(i) );
      dpi = dpi + epol*randn(1,2)*ps(i)*scalgam*sqrt(dt);

      polabs = sqrt( pol(i,1)^2 + pol(i,2)^2 );
      if polabs > mindist
        dpi = dpi - polde*pol(i,:)/polabs * dt/Tpol;
      end

      pol(i,:) = pol(i,:) + dpi;
    end
  end

  if perturb_cel==1 % per cell position perturbation
    for i = 1:cm
      x(i) = x(i) + ecel*ps(i)*sqrt(dt/Tpbar/2);
    end
  end

  %input('press enter to resume');
  simstep = simstep + 1;
  t = t + dt;
  %warning('checkpoint time increment');
end % loop over t

%if movie==1 && runfromoct ~= 1
%  mov = close(mov);
%end

finish;

