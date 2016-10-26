function yout = SWE_SphereTZ(t,xin)

global imax jmax mxpoi 
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
comega = 7.292e-05;
alpha_fcor = 0.0;
j = 1:ny;
i = 1:nx;
olat = ymin + real(j-0.5)*dy;
olon = xmin + real(i-0.5)*dx;
[XX,YY] = meshgrid(olat,olon);
coord(:,1) = reshape(YY,nx*ny,1);
coord(:,2) = reshape(XX,nx*ny,1);
f = 2*comega*(-cos(coord(:,1)).*cos(coord(:,2))*sin(alpha_fcor) + sin(coord(:,2))*cos(alpha_fcor));
alfv = 0;

up = zeros(npoin,1);
vp = zeros(npoin,1);
phip = zeros(npoin,1);

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
        
       % olon = coord(ip,1);
        olat = coord(ip,2);
        olatpq=olat + q*dy;
        
        
        
        
        
        olatmq = olat - q*dy;
        
        olatpqh = olat + qh*dy;
        
        olatmqh = olat - qh*dy;
        
        %Coriolis Force
        
%         fip = 2*comega*( -cos(olon)*cos(olat)*sin(alpha) + ...
%             sin(olat)*cos(alpha) );
%         fip4 = 2*comega*( -cos(olonpp)*cos(olat)*sin(alpha) + ...
%             sin(olat)*cos(alpha) );
%         fip3 = 2*comega*( -cos(olonmp)*cos(olat)*sin(alpha) + ...
%             sin(olat)*cos(alpha) );
%         fjp4 = 2*comega*( -cos(olonpq)*cos(olatpq)*sin(alpha) + ...
%             sin(olatpq)*cos(alpha) );
%         fjp3 = 2*comega*( -cos(olonmq)*cos(olatmq)*sin(alpha) + ...
%             sin(olatmq)*cos(alpha) );
        
        
        fip = f(ip);
        fip4 = f(ip4);
        fip3 = f(ip3);
        fjp4 = f(jp4);
        fjp3 = f(jp3);
        
        phip(ip) = - 0.5*u0(ip)/(rade*cos(olat))*(phi0(ip2)-phi0(ip1))/(dx)...
            - 0.5*v0(ip)/(rade)*(phi0(jp2)-phi0(jp1))/(dy) - 0.5*phi0(ip)/(rade*cos(olat))...
            *((1.0-alf)*( (u0(ip4h)-u0(ip3h))/(ph*dx) +(j4hsign*v0(jp4h)*cos(olatpqh) - ...
            j3hsign*v0(jp3h)*cos(olatmqh))/(qh*dy) ) +alf/2*( (j4hsign*u0(ip4hjp4h)...
            -j4hsign*u0(ip3hjp4h))/(ph*dx) +(j3hsign*u0(ip4hjp3h)-j3hsign...
            *u0(ip3hjp3h))/(ph*dx) + (j4hsign*v0(ip4hjp4h)*cos(olatpqh) -...
            j3hsign*v0(ip4hjp3h)*cos(olatmqh))/(qh*dy) +(j4hsign*v0(ip3hjp4h)*cos(olatpqh)-...
            j3hsign*v0(ip3hjp3h)*cos(olatmqh))/(qh*dy) ) );
        
        
        
        
        up(ip) = - 0.5*u0(ip)/(rade*cos(olat))*(u0(ip2)-u0(ip1))/(dx)...
            - 0.5*v0(ip)/rade*(j2sign*u0(jp2)-j1sign*u0(jp1))/(dy)- ...
            0.5/(rade*cos(olat))*(phi0(ip4h)-phi0(ip3h))/(ph*dx) ...
            + 2*0.5*((1.0-alf)*(fip + u0(ip)/rade*tan(olat))*v0(ip)...
            + alf/2*(fip4 + u0(ip4)/rade*tan(olat))*v0(ip4)...
            + alf/2*(fip3 + u0(ip3)/rade*tan(olat))*v0(ip3) );
        
        
          
        vp(ip) =  - 0.5*u0(ip)/(rade*cos(olat))*(v0(ip2)-v0(ip1))/(dx)...
            - 0.5*v0(ip)/rade*(j2sign*v0(jp2)-j1sign*v0(jp1))/(dy)...
            - 0.5/rade*( phi0(jp4h)-phi0(jp3h) )/(qh*dy)...
            - 2*0.5*( (1.0-alfv)*(fip+ u0(ip)/rade*tan(olat))*u0(ip)+ alfv/2*(fjp4...
            + j4sign*u0(jp4)/rade*tan(olatpq))*j4sign*u0(jp4)+alfv/2*(fjp3+j3sign*u0(jp3)...
            /rade*tan(olatmq))*j3sign*u0(jp3) );
    end
end
yout = zeros(3*npoin,1);
yout(1:npoin) = up;
yout(npoin+1:2*npoin) = vp;
yout(2*npoin+1:3*npoin) = phip;

%%% Rescale the time %%%%
%yout = yout/3;

return







