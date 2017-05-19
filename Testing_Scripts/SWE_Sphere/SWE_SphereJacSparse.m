function jac = SWE_SphereJacSparse(~,xin)
global imax jmax  mxpoi %mxele mxbou nd tol g rk rkd ntimemax istate nrsnapshots mx
setConstants();
nx = imax;
ny = jmax;
node=reshape(1:nx*ny,nx,ny);
npoin =nx*ny;
u0 = xin(1:npoin);
v0 = xin(npoin+1:2*npoin);
phi0 = xin(2*npoin+1:end);
p = 4;
q = 1;
ph = p;
qh = q;
nxh = nx/2;
coord = zeros(mxpoi,2);
xmin = 0;
xmax = 2*pi;
ymin = -pi/2;
ymax = pi/2;
xl = xmax-xmin;
yl = ymax-ymin;
dx = xl/nx;
dy = yl/ny;
rade = 6.37e+06;
alf = 1/3;
%phi_mean = 5.768e4;
%omega = 20.0;
comega = 7.292e-05;
%velmax = -1e5;
alpha_fcor = 0.0;
%alpha = 0.0;
j = 1:ny;
i = 1:nx;
olat = ymin + real(j-0.5)*dy;
olon = xmin + real(i-0.5)*dx;
[XX,YY] = meshgrid(olat,olon);
coord(:,1) = reshape(YY,nx*ny,1);
coord(:,2) = reshape(XX,nx*ny,1);
f = 2*comega*(-cos(coord(:,1)).*cos(coord(:,2))*sin(alpha_fcor) + sin(coord(:,2))*cos(alpha_fcor));
iv = zeros(41*nx*ny,1);
jv=iv;
sv=iv;
%alfh = 0;
%alfu = 0;
alfv = 0;
count = 1;
for i=1:nx
    
    %Longitudinal Nodes
    
    i1 = i-1;
    i2 = i+1;
    i3 = i-p;
    i4 = i+p;
    i3h = i-ph;
    i4h = i+ph;
    
    %Longitudinal periodicity
    
    if(i1 < 1)
        i1=i1+nx;
    end
    if(i2>nx)
        i2 = i2-nx;
    end
    
    %Longitudinal periodicity -P's and +P's
    
    if(i3<1)
        i3=i3+nx;
    end
    if(i4>nx)
        i4 = i4-nx;
    end
    
    %Longitudinal periodicity -P/2's and +P/2's
    
    if (i3h<1)
        i3h = i3h + nx;
    end
    if (i4h>nx)
        i4h = i4h - nx;
    end
    
    for j = 1:ny
        
        %Latitudinal Nodes
        
        j1 = j-1;
        j2 = j+1;
        j3 = j-q;
        j4 = j+q;
        j3h = j-qh;
        j4h = j+qh;
        j1sign = 1;
        j2sign = 1;
        j3sign = 1;
        j4sign = 1;
        j3hsign = 1;
        j4hsign = 1;
        
        %South pole periodicity
        
        ij1 = i;
        if(j1 < 1)
            j1 = 1;
            j1sign = -1;
            ij1 = ij1+nxh;
            if(ij1>nx)
                ij1 = ij1-nx;
            end
        end
        
        %North pole periodicity
        
        ij2 = i;
        if(j2 > ny)
            j2 =ny;
            j2sign = -1;
            ij2 = ij2+nxh;
            if(ij2>nx)
                ij2 = ij2-nx;
            end
        end
        
        %South pole periodicity -Q's
        
        ij3 = i;
        ippj3 = i4;
        impj3 = i3;
        if(j3 < 1)
            j3 = 1-j+q;
            j3sign = -1;
            ij3 = ij3+nxh;
            ippj3 = ippj3 + nxh;
            impj3 = impj3 + nxh;
            if(ij3>nx)
                ij3 = ij3-nx;
            end
            if(ippj3>nx)
                ippj3 = ippj3-nx;
            end
            if(impj3>nx)
                impj3 = impj3 -nx;
            end
        end
        
        %North Pole Periodicity +Q's
        
        ij4 = i;
        ippj4 = i4;
        impj4 = i3;
        if (j4>ny)
            j4 = 2*ny + 1 - j - q;
            j4sign = -1;
            ij4 = ij4 + nxh;
            ippj4 = ippj4 + nxh;
            impj4 = impj4 + nxh;
            if (ij4>nx)
                ij4 = ij4 - nx;
            end
            if (ippj4>nx)
                ippj4 = ippj4 - nx;
            end
            if (impj4>nx)
                impj4 = impj4 - nx;
            end
        end
        
        %South Pole Periodicity -Q/2's
        
        ij3h = i;
        ippj3h = i4h;
        impj3h = i3h;
        if (j3h<1)
            j3h = 1 - j + qh;
            j3hsign = -1;
            ij3h = ij3h + nxh;
            ippj3h = ippj3h + nxh;
            impj3h = impj3h + nxh;
            if (ij3h>nx)
                ij3h = ij3h - nx;
            end
            if (ippj3h>nx)
                ippj3h = ippj3h - nx;
            end
            if (impj3h>nx)
                impj3h = impj3h - nx;
            end
        end
        
        %North Pole Periodicity +Q/2's
        
        ij4h = i;
        ippj4h = i4h;
        impj4h = i3h;
        if (j4h>ny)
            j4h = 2*ny + 1 - j - qh;
            j4hsign = -1;
            ij4h = ij4h + nxh;
            ippj4h = ippj4h + nxh;
            impj4h = impj4h + nxh;
            if (ij4h>nx)
                ij4h = ij4h - nx;
            end
            if (ippj4h>nx)
                ippj4h = ippj4h - nx;
            end
            if (impj4h>nx)
                impj4h = impj4h - nx;
            end
        end
        
        %Centered Diff Grid Points
        
        ip = node(i,j);
        ip1 = node(i1,j);
        ip2 = node(i2,j);
        jp1 = node(ij1,j1);
        jp2 = node(ij2,j2);
        
        %Turkel-Zwas Grid Points
        
        ip3 = node(i3,j);
        ip4 = node(i4,j);
        jp3 = node(ij3,j3);
        jp4 = node(ij4,j4);
        %ip3jp3 = node(impj3,j3);
        %ip4jp3 = node(ippj3,j3);
        %ip3jp4 = node(impj4,j4);
        %ip4jp4 = node(ippj4,j4);
        
        %Staggered Grid Points
        
        ip3h = node(i3h,j);
        ip4h = node(i4h,j);
        jp3h = node(ij3h,j3h);
        jp4h = node(ij4h,j4h);
        ip3hjp3h = node(impj3h,j3h);
        ip4hjp3h = node(ippj3h,j3h);
        ip3hjp4h = node(impj4h,j4h);
        ip4hjp4h = node(ippj4h,j4h);
        
        %Longitudes and Latitudes
        
        olon = coord(ip,1);
        olat = coord(ip,2);
  %      olonpp = olon + p*dx;
  %      olonmp = olon - p*dx;
        olonpq = olon;
        
        if (j4sign==-1)
            olonpq = olonpq + pi;
        end
        olatpq = olat + q*dy;
        olonmq = olon;
        if (j3sign==-1)
            olonmq = olonmq + pi;
        end
        olatmq = olat - q*dy;
        
        olonpqh=olon;
        if (j4hsign==-1)
            olonpqh = olonpqh + pi;
        end
        olatpqh = olat + qh*dy;
        olonmqh = olon;
        if (j3hsign==-1)
            olonmqh = olonmqh + pi;
        end
        olatmqh = olat - qh*dy;
        
        fip = f(ip);
        fip4 = f(ip4);
        fip3 = f(ip3);
        fjp4 = f(jp4);
        fjp3 = f(jp3);
        
        %%%  Jacobian part for the Height  %%%
        
        K1 = 0.5*u0(ip)/(rade*cos(olat)*dx);
        
        K2 = - 0.5/(rade*cos(olat))...
            *((1.0-alf)*( (u0(ip4h)-u0(ip3h))/(ph*dx) +(j4hsign*v0(jp4h)*cos(olatpqh) - ...
            j3hsign*v0(jp3h)*cos(olatmqh))/(qh*dy) ) +alf/2*( (j4hsign*u0(ip4hjp4h)...
            -j4hsign*u0(ip3hjp4h))/(ph*dx) +(j3hsign*u0(ip4hjp3h)-j3hsign...
            *u0(ip3hjp3h))/(ph*dx) + (j4hsign*v0(ip4hjp4h)*cos(olatpqh) -...
            j3hsign*v0(ip4hjp3h)*cos(olatmqh))/(qh*dy) +(j4hsign*v0(ip3hjp4h)*cos(olatpqh)-...
            j3hsign*v0(ip3hjp3h)*cos(olatmqh))/(qh*dy) ) );
        
        K3 = 0.5*v0(ip)/(rade*dy);
        
        K4 = 0.5*phi0(ip)/(rade*cos(olat));
        
        iv(count) = 2*npoin+ip;
        jv(count) = ip;
        sv(count) = (phi0(ip2) - phi0(ip1))*(-0.5)/(dx*rade*cos(olat));
        
       % jac(2*npoin+ip,ip)  = (phi0(ip2) - phi0(ip1))*(-0.5)/(dx*rade*cos(olat));
       count=count+1;
        
       iv(count) = 2*npoin+ip;
        jv(count) = 2*npoin+ip2;
        sv(count) = -K1; 
       
       %jac(2*npoin+ip,2*npoin+ip2) = -K1;
        count=count+1;
        
        iv(count) = 2*npoin+ip;
        jv(count) = 2*npoin+ip1;
        sv(count) = K1; 
        
       % jac(2*npoin+ip,2*npoin+ip1) = K1;
        count=count+1;
        
        iv(count) = 2*npoin+ip;
        jv(count) = 2*npoin+ip;
        sv(count) = K2; 
        
        %jac(2*npoin+ip, 2*npoin+ip) = K2;
        count=count+1;
        
        iv(count) = 2*npoin+ip;
        jv(count) = npoin+ip;
        sv(count) = -0.5*(phi0(jp2)-phi0(jp1))/(dy*rade);
        
        %jac(2*npoin+ip, npoin+ip) = -0.5*(phi0(jp2)-phi0(jp1))/(dy*rade);
        count = count+1;
        iv(count) = 2*npoin+ip;
        jv(count) = 2*npoin+jp2;
        sv(count) = -K3;
        
        %jac(2*npoin+ip,2*npoin+jp2) = -K3;
        count=count+1;
        iv(count) = 2*npoin+ip;
        jv(count) = 2*npoin+jp1;
        sv(count) = K3;
        
        %jac(2*npoin+ip,2*npoin+jp1) = K3;
        count=count+1;
        iv(count) = 2*npoin+ip;
        jv(count) = ip4h;
        sv(count) = -K4*(1-alf)/(ph*dx);
        
        %jac(2*npoin+ip,ip4h) = -K4*(1-alf)/(ph*dx);
        count = count+1;
        iv(count) = 2*npoin+ip;
        jv(count) = ip3h;
        sv(count) = K4*(1-alf)/(ph*dx);
        
        %jac(2*npoin+ip,ip3h) = K4*(1-alf)/(ph*dx);
        count = count+1;
        iv(count) = 2*npoin+ip;
        jv(count) = npoin+jp4h;
        sv(count) = -K4*cos(olatpqh)*j4hsign*(1-alf)/(qh*dy);
        
        %jac(2*npoin+ip,npoin+jp4h) = -K4*cos(olatpqh)*j4hsign*(1-alf)/(qh*dy);
        count = count+1;
        iv(count) = 2*npoin+ip;
        jv(count) = npoin+jp3h;
        sv(count) = K4*(1-alf)*j3hsign*cos(olatmqh)/(qh*dy);
        
        %jac(2*npoin+ip,npoin+jp3h) = K4*(1-alf)*j3hsign*cos(olatmqh)/(qh*dy);
        count = count+1;
        iv(count) = 2*npoin+ip;
        jv(count) = ip4hjp4h;
        sv(count) = -K4*alf*j4hsign/(2*ph*dx);
        
        %jac(2*npoin+ip,ip4hjp4h) = -K4*alf*j4hsign/(2*ph*dx);
        count = count+1;
        iv(count) = 2*npoin+ip;
        jv(count) = ip3hjp4h;
        sv(count) = K4*alf*j4hsign/(2*ph*dx);
        
        %jac(2*npoin+ip,ip3hjp4h) = K4*alf*j4hsign/(2*ph*dx);
        count = count+1;
        iv(count) = 2*npoin+ip;
        jv(count) = ip4hjp3h;
        sv(count) = -K4*alf*j3hsign/(2*ph*dx);
        
        %jac(2*npoin+ip,ip4hjp3h) = -K4*alf*j3hsign/(2*ph*dx);
        count = count+1;
        iv(count) = 2*npoin+ip;
        jv(count) = ip3hjp3h;
        sv(count) = K4*alf*j3hsign/(2*ph*dx);
        
        %jac(2*npoin+ip,ip3hjp3h) = K4*alf*j3hsign/(2*ph*dx);
        count=count+1;
        iv(count) = 2*npoin+ip;
        jv(count) = npoin+ip4hjp4h;
        sv(count) =  -K4*alf*j4hsign*cos(olatpqh)/(2*qh*dy);
        %jac(2*npoin+ip,npoin+ip4hjp4h) = -K4*alf*j4hsign*cos(olatpqh)/(2*qh*dy);
        count = count+1;
        iv(count) = 2*npoin+ip;
        jv(count) = npoin+ip4hjp3h;
        sv(count) = K4*alf*j3hsign*cos(olatmqh)/(2*qh*dy);
        
        %jac(2*npoin+ip,npoin+ip4hjp3h) = K4*alf*j3hsign*cos(olatmqh)/(2*qh*dy);
        count = count+1;
        iv(count) = 2*npoin+ip;
        jv(count) = npoin+ip3hjp4h;
        sv(count) = -K4*alf*j4hsign*cos(olatpqh)/(2*qh*dy);
        %jac(2*npoin+ip,npoin+ip3hjp4h) = -K4*alf*j4hsign*cos(olatpqh)/(2*qh*dy);
        count = count+1;
        iv(count) = 2*npoin+ip;
        jv(count) = npoin+ip3hjp3h;
        sv(count) = K4*alf*j3hsign*cos(olatmqh)/(2*qh*dy);
        
        %jac(2*npoin+ip,npoin+ip3hjp3h) = K4*alf*j3hsign*cos(olatmqh)/(2*qh*dy);
        count = count+1;
        
        %%%  Jacobian part for U  %%%
        iv(count) =ip;
        jv(count) = ip;
        sv(count) =  -0.5*(u0(ip2)-u0(ip1))/(rade*cos(olat)*dx) + ...
            0.5*2*(1.0-alf)*tan(olat)*v0(ip)/rade;
        count = count+1;
        
        %jac(ip,ip) = -0.5*(u0(ip2)-u0(ip1))/(rade*cos(olat)*dx) + ...
        %    0.5*2*(1.0-alf)*tan(olat)*v0(ip)/rade;
        
        iv(count) = ip;
        jv(count) = ip2;
        sv(count) = -0.5*u0(ip)/(rade*cos(olat)*dx);
        count = count+1;
        
        %jac(ip,ip2) = -0.5*u0(ip)/(rade*cos(olat)*dx);
        iv(count) = ip;
        jv(count) = ip1;
        sv(count) = 0.5*u0(ip)/(rade*cos(olat)*dx);
        count = count+1;
        
%       %jac(ip,ip1) = 0.5*u0(ip)/(rade*cos(olat)*dx);

        iv(count) = ip;
        jv(count) = npoin+ip;
        sv(count) = -0.5*(j2sign*u0(jp2) - j1sign*u0(jp1))/(rade*dy)...
            +(1.0-alf)*2*0.5*(fip + u0(ip)*tan(olat)/rade);
        count = count+1;
        
        %jac(ip,npoin+ip) = -0.5*(j2sign*u0(jp2) - j1sign*u0(jp1))/(rade*dy)...
         %   +(1.0-alf)*2*0.5*(fip + u0(ip)*tan(olat)/rade);
        
        iv(count) = ip; 
        jv(count) = jp2;
        sv(count) = -0.5*v0(ip)*j2sign/(rade*dy);
        count = count+1;
        %jac(ip,jp2) = -0.5*v0(ip)*j2sign/(rade*dy);
        
        iv(count) = ip;
        jv(count) = jp1;
        sv(count) =  0.5*v0(ip)*j1sign/(rade*dy);
        count = count+1;
        
        %jac(ip,jp1) = 0.5*v0(ip)*j1sign/(rade*dy);
        
        iv(count) = ip;
        jv(count) = 2*npoin+ip4h;
        sv(count) =  -0.5/(rade*cos(olat)*ph*dx);
        count = count+1;
        
        %jac(ip,2*npoin+ip4h) = -0.5/(rade*cos(olat)*ph*dx);
        
        iv(count) = ip;
        jv(count) = 2*npoin+ip3h;
        sv(count) = 0.5/(rade*cos(olat)*ph*dx);
        count = count+1;
        
        %jac(ip,2*npoin+ip3h) = 0.5/(rade*cos(olat)*ph*dx);
        
        iv(count) = ip;
        jv(count) = ip4;
        sv(count) = 0.5*2*alf*tan(olat)*v0(ip4)/(2*rade);
        count = count+1;
        
        %jac(ip,ip4) = 0.5*2*alf*tan(olat)*v0(ip4)/(2*rade);
        
        iv(count) = ip;
        jv(count) = npoin+ip4;
        sv(count) = 2*0.5*alf*(fip4 + u0(ip4)*tan(olat)/rade)/(2);
        count = count+1;
        
        %jac(ip,npoin+ip4) = 2*0.5*alf*(fip4 + u0(ip4)*tan(olat)/rade)/(2);
        
        iv(count) = ip;
        jv(count) = ip3;
        sv(count) = alf*2*0.5*tan(olat)*v0(ip3)/(2*rade);
        count = count+1;
        
        %jac(ip,ip3) = alf*2*0.5*tan(olat)*v0(ip3)/(2*rade);
        
        iv(count) = ip;
        jv(count) = npoin+ip3;
        sv(count) = alf*2*0.5*(fip3+u0(ip3)*tan(olat)/rade)/2;
        count = count+1;
        
        %jac(ip,npoin+ip3) = alf*2*0.5*(fip3+u0(ip3)*tan(olat)/rade)/2;
        
        %%%  Jacobian part for V  %%%
        
        iv(count) = npoin+ip;
        jv(count) = ip;
        sv(count) = -0.5*(v0(ip2) - v0(ip1))/(rade*cos(olat)*dx)...
            - 2*0.5*(1-alfv)*(fip) - 4*0.5*(1-alfv)*(tan(olat))*u0(ip)/rade;
        count = count+1;
        
        %jac(npoin+ip,ip) = -0.5*(v0(ip2) - v0(ip1))/(rade*cos(olat)*dx)...
        %   -2*(1-alfv)*0.5*(fip + 2*u0(ip)*tan(olat)/rade);
       
        iv(count) = npoin+ip;
        jv(count) = npoin+ip2;
        sv(count) = -0.5*u0(ip)/(rade*cos(olat)*dx);
        count = count+1;
        
        %jac(npoin+ip,npoin+ip2) = -0.5*u0(ip)/(rade*cos(olat)*dx);
        
        iv(count) = npoin+ip;
        jv(count) = npoin+ip1;
        sv(count) = 0.5*u0(ip)/(rade*cos(olat)*dx);
        count = count+1;
        
        %jac(npoin+ip,npoin+ip1) = 0.5*u0(ip)/(rade*cos(olat)*dx);
        
        iv(count) = npoin+ip;
        jv(count) = npoin+ip;
        sv(count) = -0.5*(j2sign*v0(jp2) - j1sign*v0(jp1))/(rade*dy);
        count = count+1;
        
        %jac(npoin+ip,npoin+ip) = -0.5*(j2sign*v0(jp2) - j1sign*v0(jp1))/(rade*dy);
        
        iv(count) = npoin+ip;
        jv(count) = npoin+jp2;
        sv(count) = -0.5*v0(ip)*j2sign/(dy*rade);
        count = count+1;
        
        %jac(npoin+ip,npoin+jp2) = -0.5*v0(ip)*j2sign/(dy*rade);
        
        iv(count) = npoin+ip;
        jv(count) = npoin+jp1;
        sv(count) = 0.5*v0(ip)*j1sign/(dy*rade);
        count = count+1;
        
        %jac(npoin+ip,npoin+jp1) = 0.5*v0(ip)*j1sign/(dy*rade);
        
        iv(count) = npoin+ip;
        jv(count) = 2*npoin+jp4h;
        sv(count) = -0.5/(rade*qh*dy);
        count = count+1;
        
        %jac(npoin+ip,2*npoin+jp4h) = -0.5/(rade*qh*dy);
        
        iv(count) = npoin+ip;
        jv(count) = 2*npoin+jp3h;
        sv(count) = 0.5/(rade*qh*dy);
        count = count+1;
        
        %jac(npoin+ip,2*npoin+jp3h) = 0.5/(rade*qh*dy);
        
        iv(count) = npoin+ip;
        jv(count) = jp4;
        sv(count) = -2*0.5*(alfv/2*fjp4*j4sign) - 4*0.25*alfv*(j4sign^2*u0(jp4)*tan(olatpq)/rade);
        count = count+1;
        
        %jac(npoin+ip,jp4) = -alfv*2*0.5*(j4sign*tan(olatpq))*u0(jp4)/(2*rade)...
        %    -alfv*2*0.5*(fjp4 + j4sign*u0(jp4)*tan(olatmq))*j4sign/rade;
        
        iv(count) = npoin+ip;
        jv(count) = jp3;
        sv(count) = -2*0.25*alfv*fjp3*j3sign - 4*0.25*alfv *(j3sign^2*u0(jp3)/rade*tan(olatmq));
        count = count+1;
        
        %jac(npoin+ip,jp3) = -alfv*0.5*(fjp3 + j3sign*u0(jp3)*tan(olatmq)/rade)*j3sign ...
           % -alfv*0.5*(j3sign*tan(olatmq))*u0(jp3)*j3sign/rade;
        
    end
end
jac = sparse(iv,jv,sv,3*npoin,3*npoin);
return