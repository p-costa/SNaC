%
% setting up some parameters
%
precision = 'double';     % precision of the real-valued data
r0 = [0.,0.,0.];    % domain origin
non_uniform_grid = true;
%
% read geometry file
%
geofile  = "geometry.out";
data = dlmread(geofile);
ng   = data(1,:);
l    = data(2,:);
dl   = l./ng;
%
% read and generate grid
%
xp = linspace(r0(1)+dl(1)/2.,r0(1)+l(1),ng(1)); % centered  x grid
yp = linspace(r0(2)+dl(2)/2.,r0(2)+l(2),ng(2)); % centered  y grid
zp = linspace(r0(3)+dl(3)/2.,r0(3)+l(3),ng(3)); % centered  z grid
xu = xp + dl(1)/2.;                           % staggered x grid
yv = yp + dl(2)/2.;                           % staggered y grid
zw = zp + dl(3)/2.;                           % staggered z grid
if(non_uniform_grid)
    f   = fopen('grid_x.bin');
    grid_x = fread(f,[ng(3),4],precision);
    fclose(f);
    xp = r0(1) + grid_x(:,2)'; % centered  x grid
    xu = r0(1) + grid_x(:,1)'; % staggered x grid
    f   = fopen('grid_y.bin');
    grid_y = fread(f,[ng(3),4],precision);
    fclose(f);
    yp = r0(2) + grid_y(:,2)'; % centered  y grid
    yv = r0(2) + grid_y(:,1)'; % staggered y grid
    f   = fopen('grid_z.bin');
    grid_z = fread(f,[ng(3),4],precision);
    fclose(f);
    zp = r0(3) + grid_z(:,2)'; % centered  z grid
    zw = r0(3) + grid_z(:,1)'; % staggered z grid
end
%
% read checkpoint binary file
%
filenamei = input("Name of the binary restart file written by CaNS [fld.bin]: ");
if isempty(filenamei)
    filenamei = "fld.bin";
end
data       = zeros([ng(1),ng(2),ng(3),4]); % u(:,:,:),v(:,:,:),w(:,:,:),p(:,:,:)
fldinfo    = zeros([2,1]);
f = fopen(filenamei);
for p = 1:4
    data(:,:,:,p) = reshape(fread(f,ng(1)*ng(2)*ng(3),precision),[ng(1),ng(2),ng(3)]);
end
fldinfo(:) = fread(f,2,precision);
fclose(f);
%
% store data in arrays
%
u = data(:,:,:,1);
v = data(:,:,:,2);
w = data(:,:,:,3);
p = data(:,:,:,4);
time  =       fldinfo(1);
istep = round(fldinfo(2));
