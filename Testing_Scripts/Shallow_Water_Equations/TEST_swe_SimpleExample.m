Tspan = [0 3.2];
Ode_Function = @swe_Function_mex;
y0 = swe_Initialize_y0();

Options = MATLODE_OPTIONS('AbsTol', ones(3468,1).*1e-10, ...
    'RelTol', ones(3468,1).*1e-10, 'storeCheckpoint', true);

[ T, Y, Stats ] = MATLODE_ERK_FWD_Integrator( Ode_Function, Tspan, y0, Options );

swe_simpleplot(T,Y);

