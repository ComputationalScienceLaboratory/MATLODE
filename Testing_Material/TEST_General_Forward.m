%% Testing Material: Testing_General_Forward
%
%
%%
function [ steps, errorSolution ] = TEST_General_Forward( MATLODE_Integrator, Ode_Function, Tspan, y0, Option, MaxTolerance, MinTolerance, T_Ref, Y_Ref )

    Npoints = 7;
    steps = zeros(Npoints,1);
    TOLS = logspace(MinTolerance,MaxTolerance,Npoints);

    errorSolution = zeros(Npoints,1);

    for ipt=1:Npoints
        RelTol=TOLS(ipt);
        AbsTol = RelTol;
        Option = MATLODE_OPTIONS(Option, 'AbsTol', AbsTol, 'RelTol', RelTol);

        [ T, Y, Stats ] = ...
        MATLODE_Integrator(Ode_Function, Tspan, y0, Option );

        errorSolution(ipt,1) = rRMS(Y(length(T),:),Y_Ref(length(T_Ref),:));    
        steps(ipt) = Stats.ISTATUS.Nstp;

    end

end

