function setConstants()
% global mlat nlon nbound carNbound sbound carSbound eqN eqS dlambda dphi
% global cartXN cartYN dxN dyN mphiN qN f rad cartXS cartYS dxS dyS mphiS qS
% rad = 6381000;
global imax jmax mx mxpoi mxele mxbou nd tol g rk rkd ntimemax istate nrsnapshots
global phi_mean
imax = 72;
jmax = 36;
mx =imax*jmax;
mxpoi = mx;
mxele = mx;
mxbou = mx/5;
nd = 4;
tol = 1e-6;
g = 10.0;
rk = 0.1;
rkd = 0.0;
ntimemax = 150;
istate = mxpoi*3;
nrsnapshots = 37;
phi_mean = 5.768E4;


