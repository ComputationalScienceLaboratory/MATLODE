%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%           	    COMPUTATIONAL SCIENCE LABORATORY                   %%%
%%%                           MOFSET: Driver                           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Website: http://csl.cs.vt.edu/                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
clc;

springs_Initialize;

options = odeset('Stats', 'on','Jacobian',jacFun);
[ T, Y ] = ode45(rhsFun,tspan,x0,options);

plot(T,Y(:,1),T,Y(:,2),T,Y(:,3));