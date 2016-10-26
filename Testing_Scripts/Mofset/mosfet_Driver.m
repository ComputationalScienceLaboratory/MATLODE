%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%           	    COMPUTATIONAL SCIENCE LABORATORY                   %%%
%%%                           MOFSET: Driver                           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Website: http://csl.cs.vt.edu/                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
clc;

mosfet_Initialize;

options = odeset('Stats', 'on', 'Jacobian', jacFun);
[ T, Y ] = ode23s(rhsFun,tspan,Uinit,options);

plot(T,Y);