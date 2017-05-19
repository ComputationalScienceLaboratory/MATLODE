load initialconditions
y0 = initialconditions;
tspan = [2000 3000];

options = MatlOde_OPTIONS('Jacobian',@SWE_SphereJacSparse,'Y_TLM',eye(7776),'displaySteps',true);

[tout, yout] = MatlOde_RK_FWD_DRIVER_Integrator(@SWE_SphereTZ, tspan, y0, options);

%swe_simpleplot(yout(end,:),0);