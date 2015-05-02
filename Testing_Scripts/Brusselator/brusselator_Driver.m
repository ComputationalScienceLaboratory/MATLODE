%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%           	    COMPUTATIONAL SCIENCE LABORATORY                   %%%
%%%                        Van Der Pol: Driver                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Website: http://csl.cs.vt.edu/                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;

brusselator_Initialize;

options = odeset('Jacobian',jacFun);
[ T, Y ] = ode15s(rhsFun,tspan,x0,options);

plot(Y(:,1),Y(:,2));