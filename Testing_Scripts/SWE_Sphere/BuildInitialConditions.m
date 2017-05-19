%Setting initial conditions
global imax jmax phi_mean

% Change imax and jmax to change the resolution


setConstants();


p = 4;
q = 1;
ny = jmax;
nx = imax;
npoin = nx*ny;


xmin = 0;
xmax = 2*pi;
ymin = -pi/2;
ymax = pi/2;
xl = xmax-xmin;
yl = ymax-ymin;
dx = xl/nx;
dy = yl/ny;
rade = 6.37E6;
omega = 20;
comega = 7.292E-05;
 jj = 1:jmax;
 ii = 1:imax;
 olatt = ymin+(jj-0.5)*dy;
 olonn = xmin+(ii-0.5)*dx;
 
 ip = 0;
ui = zeros(npoin,1);
vi = zeros(npoin,1);
phii = zeros(npoin,1);
for j=1:jmax
    olat=ymin + (j-0.5)*dy;
    for i=1:imax
        olon=xmin + (i-0.5)*dx;
        ip=ip+1;
        s = sin(olat)*sin(olat)*sin(olat);
        c = cos(olat)*cos(olat);
        ui(ip)=omega*sin(olon)*(s-3*sin(olat)*c);
        vi(ip)=omega*sin(olat)*sin(olat)*cos(olon);
        phii(ip)=phi_mean + 2*comega*rade*omega*s*cos(olat)*sin(olon);
    end
end
y0 = zeros(3*npoin,1);
y0(1:npoin,1) = ui;
y0(npoin+1:2*npoin,1) = vi;
y0(2*npoin+1:3*npoin,1) = phii;

%y0 is the initial conditions required to run the SWE model using ode45.