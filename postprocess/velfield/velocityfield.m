% probably does not work


%%%%%%%% configuration %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xdat='x100.dat'
vdat='v100.dat'

Ncell=23;

outfile='vf.dat'

xsiz=10; %um
ysiz=10; %um

xgrid=0.5; %um
ygrid=0.5; %um



%%%%%%%% presets %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xx=load(xdat);
vv=load(vdat);

Nx=ceil(xsiz/xgrid); % needs to be even
Ny=ceil(ysiz/ygrid); % needs to be even

vfx=zeros(2*Nx,2*Ny);
vfy=zeros(2*Nx,2*Ny);
vfxnum=zeros(2*Nx,2*Ny);
vfynum=zeros(2*Nx,2*Ny);



%%%%%%%% compute %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:Ncell
  myx  = xx(:,2+2*i-1); % fist two columns counter+time
  myy  = xx(:,2+2*i  ); % fist two columns counter+time
  myvx = vv(:,2+3*i-2);   % fist two columns counter+time
  myvy = vv(:,2+3*i-1); % fist two columns counter+time

  for j=1:Ncell
    if j==i
      continue;
    end
    status = sprintf( 'i=%i, j=%i', i,j );
    disp(status);

    % find position and velocity differences
    deltax  = xx(:,2+2*j-1) - myx;
    deltay  = xx(:,2+2*j  ) - myy;
    deltavx = vv(:,2+3*j-2) - myvx;
    deltavy = vv(:,2+3*j-1) - myvy;

    % which cells are within observed region?
    mynbr = find( deltax>-xsiz && deltax<xsiz && deltay>-ysiz && deltay<ysiz );

    % find bin indices according to position
    xbin = floor(deltax(mynbr))+Nx;
    ybin = floor(deltay(mynbr))+Ny;

    % add velocities
    for k=1:length(mynbr)
       vfx(xbin(k),ybin(k)) = vfx(xbin(k),ybin(k)) + deltavx(mynbr(k));
       vfy(xbin(k),ybin(k)) = vfy(xbin(k),ybin(k)) + deltavy(mynbr(k));
       vfxnum(xbin(k),ybin(k)) = vfxnum(xbin(k),ybin(k)) + 1;
       vfynum(xbin(k),ybin(k)) = vfynum(xbin(k),ybin(k)) + 1;
    end
  end
end

% divide for averaging
for i=1:Nx
  for j=1:Ny
    if vfxnum(i,j)>0
      vfx = vfx ./ vfxnum;
    end
    if vfynum(i,j)>0
      vfy = vfy ./ vfynum;
    end
  end
end



%%%%%%%% print %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ofd = fopen( outfile, 'w' );

fprintf( ofd, '%% mean local velocity field arround single cell\n' );
fprintf( ofd, '%% x y    vx vy v\n' );
for i=1:Nx
  for j=1:Ny
    fprintf( ofd, '%.3f %.3f   ', i*xgrid-xsiz, j*ygrid-ysiz );
    fprintf( ofd, '%.6f %.6f %.6f\n', vfx(i,j), vfy(i,j), ...
                                    sqrt(vfx(i,j)^2+vfy(i,j)^2) );
  end
  fprintf( ofd, '\n' );
end

