load initialconditions
y0 = initialconditions;

tspan = [0 3600];

%%%%%Integrate the ODE%%%%%%%%

[tout, yout] = ode45(@SWE_SphereTZ, tspan, y0);


%%%%%% Visualize the final solution %%%%%%

swe_simpleplot(yout(end,:),0);


%%%%%%Integrate it using Forward Euler%%%%%%
return




