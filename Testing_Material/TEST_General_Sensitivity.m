%% Testing Material: TEST_General_Sensitivity
%
%%
%  Authored by Tony D'Augustine, Adrian Sandu, and Hong Zhang.
%  Computational Science Laboratory, Virginia Tech.
%  ©2015 Virginia Tech Intellectual Properties, Inc.
%
function [ steps, errorSolution, errorSensitivity, errorQuadrature, errorMu ] = TEST_General_Sensitivity( MATLODE_Integrator, Ode_Function, Tspan, y0, Option, MaxTolerance, MinTolerance, T_Ref, Y_Ref, Sens_Ref, Quad_Ref, Mu_Ref )

    Npoints = 7;
    steps = zeros(Npoints,1);
    TOLS = logspace(MinTolerance,MaxTolerance,Npoints);

    errorSolution = zeros(Npoints,1);
    errorSensitivity = zeros(Npoints,1);
    errorQuadrature = zeros(Npoints,1);
    errorMu = zeros(Npoints,1);
    
    solverName = regexp(func2str(MATLODE_Integrator),'_','split');
    solverFamily = solverName(2);
    solverClass = solverName(3);

    for ipt=1:Npoints
        RelTol=TOLS(ipt);
        AbsTol = RelTol;
        
        Option = MATLODE_OPTIONS(Option, 'AbsTol', AbsTol, 'RelTol', RelTol, ...
            'AbsTol_ADJ',AbsTol,'RelTol_ADJ',RelTol,'AbsTol_TLM',AbsTol,'RelTol_TLM',RelTol);

        % Suppress tolerance warnings for publication
        if ( strcmp(solverClass,'ADJ') ) 
            Option = MATLODE_OPTIONS(Option,'AbsTol_TLM',[],'RelTol_TLM',[]);
            if ( strcmp(solverFamily,'ERK') )
                Option = MATLODE_OPTIONS(Option,'AbsTol_ADJ',[],'RelTol_ADJ',[]);
            end
        else
            Option = MATLODE_OPTIONS(Option,'AbsTol_ADJ',[],'RelTol_ADJ',[]);
        end
       
        if ( strcmp(solverClass,'ADJ') )
            [ T, Y, Sens, Quad, Mu, Stats ] = MATLODE_Integrator(Ode_Function, Tspan, y0, Option );            
            steps(ipt) = Stats.ISTATUS_ADJ.Nstp;
            
            if ( ~isempty(Quad) )
                errorQuadrature(ipt,1) = rRMS(Quad,Quad_Ref);
            end
            if ( ~isempty(Mu) )
                errorMu(ipt,1) = rRMS(Mu,Mu_Ref);
            end
            
        else
            [ T, Y, Sens, Stats ] = MATLODE_Integrator(Ode_Function, Tspan, y0, Option );            
            steps(ipt) = Stats.ISTATUS.Nstp;
        end
        
        errorSolution(ipt,1) = rRMS(Y(length(T),:),Y_Ref(length(T_Ref),:));
        errorSensitivity(ipt,1) = rRMS(Sens,Sens_Ref);        

    end

end

