function Av = Mat_Free_Jac(t,y,v, rhs, fcn0, normy)

myeps = 10^-7 * max(10^-5 , normy);

Av = rhs(t,y + myeps*v ) - fcn0;
Av = Av / myeps;

return