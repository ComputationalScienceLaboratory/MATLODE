% z0 = tn , z1 = tn, z2 = tn+1 , z3 = tn+1
function f_tilde = hermiteinterp(t0, t1, y0, y1, y0p, y1p) %#ok<*FNDEF>
% syms x;
% q = x^3 + 3*x^2 + 2*x + 1;
% t0 = 0;
% t1 = 1;
% y0 = subs(q,t0);
% y1 = subs(q,t1);
% yp1 = subs(diff(q,x),t1);
% yp0 = subs(diff(q,x),t0);


z0 = t0; f0 = y0;
z1 = t0; f1 = y0;
z2 = t1; f2 = y1;
z3 = t1; % f3 = y1;

fz0z1   = y0p;
fz1z2   = (f1 - f2)/(z1 - z2); 
fz2z3   = y1p;
fz0z1z2 = (fz0z1 - fz1z2)/(z0 - z2); 
fz1z2z3 = (fz1z2 - fz2z3)/(z1 - z3); 
fz0_z3  = (fz0z1z2 - fz1z2z3)/(z0 - z3);


f_tilde = @(t) (f0 + fz0z1 * (t - t0) + fz0z1z2 * (t - t0)^2 + ...
                fz0_z3 * ((t - t0)^2) * (t - t1));
end