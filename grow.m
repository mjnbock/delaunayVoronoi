%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% grow.m
%%%%%%%%
%%%%%%%% cell growth and division
% all variables defined and explained in initialize.m
% called from main.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%% initialize mitosis, growth, and then grow %%%%%%%%%%%%%
cmincr = 0;
for i=1:cm
  %%%%%%%% initialize mitosis %%%%%%%%
  if Tgrow(i)>Tgrow1(i)-dt/2 % start mitosis
    % reset time indicators
    Tgrow(i) = -dt;
    Tmito(i) = 0.0; % will be Voronoi pair from next timestep on

    % setup mitosis duration
    Tmito1(i) = Tmito0(typ(i)) + PTmito(typ(i))*randn(1,1);

    % adjust size
    r0(i) = r0(i)/mitofac(typ(i));
    w0(i) = r0(i);

    % determine division plane
    vdist = randn(1,2);
    dist  = vdist*vdist';
    vdist = 10*mindist*vdist/dist;

    % adjust positions, only slightly apart
    x(i) = x(i) - vdist(1);
    y(i) = y(i) - vdist(2);

    % duplicate cell
    x = [ x; x(i) + vdist(1) ];
    y = [ x; y(i) + vdist(2) ];
    r = [ r; r(i) ];
    w = [ w; w(i) ];
    r0 = [ r0; r0(i) ];
    w0 = [ w0; w0(i) ];
    typ = [ typ; typ(i) ];
    Tp  = [ Tp;  Tp(i)  ];
    rho = [ rho; rho(i) ];
    bi  = [ bi;  bi(i)  ];
    cd  = [ cd;  cd(i)  ];
    my  = [ my;  my(i)  ];
    ps  = [ ps;  ps(i)  ];
    Tgrow = [ Tgrow; Tgrow(i) ];
    Tmito = [ Tgrow; Tmito(i) ];
    Tgrow1 = [ Tgrow1; Tgrow1(i) ];
    Tmito1 = [ Tgrow1; Tmito1(i) ];

    v = [ v; v(i,:) ];

    cmincr = cmincr + 1;
  end

  %%%%%%%% initialize growth %%%%%%%%
  if Tmito(i)>Tmito1(i)-dt/2 % start growth
    printf( 'i=%i, growfac(typ(i))=%f\n',i,growfac(typ(i)) );

    Tgrow(i) = 0.0;
    Tmito(i) = -dt;

    % adjust size
    r0(i) = r0(i) * growfac(typ(i));
    w0(i) = r0(i);
    
    % setup growth duration 
    Tgrow1(i) = Tgrow0(typ(i)) + PTmito(typ(i))*randn(1,1);
  end

  %%%%%%%% perform growth %%%%%%%%
  if Tgrow(i)>-dt/2 && Tgrow(i)<Tgrow1(i)-dt/2 % perform growth
    grfi = growfac(typ(i));
    grrr = r0(i)/grfi;
    % linear interpolation from intitial grrr to current r0(i)
    r(i) = grrr*( 1 + (grfi-1)*Tgrow(i)/Tgrow1(i) );
    w(i) = r(i);
    Tgrow(i) = Tgrow(i) + dt;
  end

  % mitosis pairwise, see below
end
cm = cm + cmincr;



%%%%%%%% perform mitosis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for p=1:am
  % During the first timestep of mitosis, this loop is not entered,
  % because the dividing cells do not form a pair just yet.
  m = a2a(p,1);
  i = p2c(m,1);
  j = p2c(m,2);

  %%%%%%%% perform mitosis %%%%%%%%
  if Tmito(i)>-dt/2 && Tmito(i)<Tmito1(i)-dt/2
    % consistency checks
    if abs(Tmito(i)-Tmito(j))>dt/3 || abs(Tmito1(i)-Tmito1(j))>dt/3
      error('grow.m: unconsistent mitosis times, cells %i,%i',i,j);
    end
    if abs(r0(i)-r0(j))>mindist || abs(w0(i)-w0(j))>mindist
      error('grow.m: unconsistent target sizes, cells %i,%i',i,j);
    end
    if abs(r(i)-r(j))>mindist || abs(w(i)-w(j))>mindist
      error('grow.m: unconsistent sizes, cells %i,%i',i,j);
    end

    % prescribe linear angle separation function, see DXXX
    sepphi  = pi*( 1-Tmito(i)/Tmito1(i) );
    separea = pi*r0(i)^2;
    seprad  = sqrt( separea/(2*pi-sepphi-sin(sepphi)) );
    sepdist = 2*seprad*cos(sepphi/2);

    sepcenter = [ x(i)+x(j), y(i)+y(j) ]/2;
    vdist = [ x(j)-x(i), y(j)-y(i) ];
    dist = vdist*vdist';
    udist = vdist/dist;

    % adjust position, size
    if dist > sepdist
      warning('grow.m: small separating distance, cells %i,%i',i,j);
    else
      vdist = udist*d/2;
      x(i) = sepcenter(1) - vdist(1);
      y(i) = sepcenter(2) - vdist(2);
      x(j) = sepcenter(1) + vdist(1);
      y(j) = sepcenter(2) + vdist(2);

      r(i) = seprad;
      r(j) = seprad;

      w(i) = r(i);
      w(j) = r(j);
    end

    % incement time constants
    Tmito(i) = Tmito(i) + dt;
    Tmito(j) = Tmito(j) + dt;
  end
end


% TODO slightly randomize cell sizes

