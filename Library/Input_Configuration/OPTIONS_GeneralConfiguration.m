%%
%
%% Reference
% [1] Tony D'Augustine, Adrian Sandu. MATLODE: A MATLAB ODE Solver and
%     Sensitivity Analysis Toolbox. Submitted to ACM TOMS.
%
function [ OPTIONS, Coefficient ] = OPTIONS_GeneralConfiguration( OPTIONS, family, implementation, y0, tspan )
    
    roundOff = eps/2;

    % RelTol
    if ( OPTIONS.RelTol == 0 )
        OPTIONS.RelTol = ones(size(y0,1),1).*10^-4;     
    elseif ( OPTIONS.RelTol > 0 )
        % DO NOTHING: User supplied or fine tuned value.
    else
        str = [ 'Error: User selected RelTol: ', num2str(OPTIONS.RelTol), '. RelTol must be >= 0.' ];
        error(str);
    end
    
    % RelTol_TLM
    if ( OPTIONS.RelTol_TLM == 0 )
        OPTIONS.RelTol_TLM = ones(size(y0,1),size(y0,1)).*10^-4;
    elseif  ( OPTIONS.RelTol_TLM > 0 )
        % DO NOTHING: User supplied or fine tuned value.
    elseif ( isempty(OPTIONS.RelTol_TLM) == true )
        % DO NOTHING: RelTol_TLM is not required for specified integrator.
    else
        str = [ 'Error: User selected RelTol_TLM: ', num2str(OPTIONS.RelTol_TLM), '. RelTol_TLM must be >= 0.' ];
        error(str);
    end
    
    % RelTol_ADJ
    if ( OPTIONS.RelTol_ADJ == 0 )
        OPTIONS.RelTol_ADJ = OPTIONS.RelTol;
    elseif ( OPTIONS.RelTol_ADJ > 0 )
        % DO NOTHING: User supplied or fine tuned value.
    elseif ( isempty(OPTIONS.RelTol_ADJ) == true )
        % DO NOTHING: RelTol_ADJ is not required for specified integrator.
    else
        str = [ 'Error: User selected RelTol_ADJ: ', num2str(OPTIONS.RelTol_ADJ), '. RelTol_ADJ must be >= 0.' ];
        error(str);
    end
    
    % AbsTol
    if ( OPTIONS.AbsTol == 0 )
        OPTIONS.AbsTol = max(y0)*OPTIONS.RelTol;
    elseif ( OPTIONS.AbsTol > 0 )
        % DO NOTHING: User supplied or fine tuned value.
    else
        str = [ 'Error: User selected AbsTol: ', num2str(OPTIONS.AbsTol), '. AbsTol must be >= 0.' ];
        error(str);
    end
    
    % AbsTol_TLM
    if ( OPTIONS.AbsTol_TLM == 0 )
        OPTIONS.AbsTol_TLM = ones(size(y0,1),size(y0,1)).*10^-4;
    elseif ( OPTIONS.AbsTol_TLM > 0 )
        % DO NOTHING: User supplied or fine tuned value.
    elseif ( isempty(OPTIONS.AbsTol_TLM) == true )
        % DO NOTHING: AbsTol_TLM is not required by specified integrator.
    else
        str = [ 'Error: User selected AbsTol_TLM: ', num2str(OPTIONS.AbsTol_TLM), '. AbsTol_TLM must be >= 0.' ];
        error(str);
    end
    
    % AbsTol_ADJ
    if ( OPTIONS.AbsTol_ADJ == 0 )
        OPTIONS.AbsTol_ADJ = OPTIONS.RelTol;
    elseif ( OPTIONS.AbsTol_ADJ > 0 )
        % DO NOTHING: User supplied or fine tuned value.
    elseif ( isempty(OPTIONS.AbsTol_ADJ) == true ) 
        % DO NOTHING: AbsTol_ADJ is not required by specified integrator.
    else
        str = [ 'Error: User selected AbsTol_ADJ: ', num2str(OPTIONS.AbsTol_ADJ), '. AbsTol_ADJ must be >= 0.' ];
        error(str);
    end

    % Hmin
    if ( OPTIONS.Hmin == 0 )
        OPTIONS.Hmin = 0;
    elseif ( OPTIONS.Hmin > 0 )
        OPTIONS.Hmin = min( abs(OPTIONS.Hmin), abs( tspan(2)-tspan(1) ) );
    else
        str = [ 'Error: User selected Hmin: ', num2str(OPTIONS.Hmin), '. Hmin must be >= 0.' ];
        error(str);
    end
        
    % Hmax
    if ( OPTIONS.Hmax == 0 )
        OPTIONS.Hmax = abs( tspan(2)-tspan(1) );
    elseif ( OPTIONS.Hmax > 0.0 )
        OPTIONS.Hmax = min( abs(OPTIONS.Hmax), abs( tspan(2)-tspan(1) ) );
    else
        str = [ 'Error: User selected Hmax: ', num2str(OPTIONS.Hmax), '. Hmax must be >= 0.' ];
        error(str);
    end

    % FacMin
    if ( OPTIONS.FacMin == 0 )
        OPTIONS.FacMin = 0.2;
    elseif ( OPTIONS.FacMin > 0.0 )
        % DO NOTHING: User supplied or fine tuned value.
    else
        str = [ 'ERROR: User selected FacMin: ', num2str(OPTIONS.FacMin), '. FacMin must be >= 0.' ];
        error(str);
    end
    
    % FacMax
    if ( OPTIONS.FacMax == 0 )
        switch ( family )
            case 'ERK'
                OPTIONS.FacMax = 10.0;
            case 'EXP'
                OPTIONS.FacMax = 6.0;
            case 'EXPK'
                OPTIONS.FacMax = 6.0;
            case 'RK'
                OPTIONS.FacMax = 8.0;
            case 'ROK'
                OPTIONS.FacMax = 6.0;
            case 'ROS'
                OPTIONS.FacMax = 6.0;
            case 'SDIRK'
                OPTIONS.FacMax = 10.0;
        end
    elseif ( OPTIONS.FacMax > 0.0 )
        % DO NOTHING: User supplied or fine tuned value.
    else
        str = ['Error: User selected FacMax: ', num2str(OPTIONS.FacMax), '. FacMin must be >= 0.' ];
        error(str);
    end
    
    % FacRej
    if ( OPTIONS.FacRej == 0 )
        OPTIONS.FacRej = 0.1;
    elseif ( OPTIONS.FacRej > 0 )
        % DO NOTHING: User supplied or fine tuned value.
    else
        str = [ 'Error: User selected FacRej: ', num2str(OPTIONS.FacRej), '. FacMax must be >= 0.' ];
        error(str);
    end
           
    % FacSafe
    if ( OPTIONS.FacSafe == 0 )
        OPTIONS.FacSafe = 0.9;
    elseif ( OPTIONS.FacSafe > 0.0 )
        % DO NOTHING: User supplied or fine tuned value.
    else
        str = [ 'Error: User selected FacSafe: ', num2str(OPTIONS.FacSafe), '. FacSafe must be >= 0.' ];
        error(str);
    end    
    
    % FDincrement 
    if ( OPTIONS.FDIncrement == 0 )
        OPTIONS.FDIncrement = sqrt(max(OPTIONS.RelTol(1),roundOff));
    elseif ( OPTIONS.FDIncrement > 0.0 )
        % DO NOTHING: User supplied or fine tuned value.
    else
        str = [ 'Error: User selected FDincrement: ', num2str(OPTIONS.FDIncrement),'. FDincrement must be > 0.'];
        if ( strcmp(implementation,'ERK') && strcmp(family,'TLM') )
            error(str);
        end
    end
    
    % ThetaMin
    if ( OPTIONS.ThetaMin == 0 )
        OPTIONS.ThetaMin = 1.0d-3;
    elseif ( isempty(OPTIONS.ThetaMin) == true )
        % DO NOTHING: ThetaMin is not required by specified integrator.
    elseif ( OPTIONS.ThetaMin <= 0.0 || OPTIONS.ThetaMin >= 1.0 )
        str = [ 'Error: User selected ThetaMin: ', num2str(OPTIONS.ThetaMin), '. ThetaMin must be <=0 or >= 1.' ];
        error(str);
    else
        % DO NOTHING: User supplied or fine tuned value.
    end    

    % NewtonTol
    if ( OPTIONS.NewtonTol == 0 )
        OPTIONS.NewtonTol = 3.0d-2;
    elseif ( OPTIONS.NewtonTol <= roundOff )
        str = [ 'Error: User selected NewtonTol: ', num2str(OPTIONS.NewtonTol), '. NewtonTol must be >= eps/2.' ];
        error(str);
    else
        % DO NOTHING: User supplied value or fine tuned value.
    end    
    
    % Qmin
    if ( OPTIONS.Qmin == 0 )
        OPTIONS.Qmin = 1.0;
    elseif ( OPTIONS.Qmin > 0 )
        % DO NOTHING: User supplied value or fine tuned value.
    elseif ( isempty(OPTIONS.Qmin) == true )
        % DO NOTHING: Qmin is not required by specified integrator.
    else
        str = [ 'Error: User selected Qmin: ', num2str(OPTIONS.Qmin), '. Qmin must be >= 0.' ];
        error(str);
    end    
    
    % Qmax
    if ( OPTIONS.Qmax == 0 )
        OPTIONS.Qmax = 1.2;
    elseif ( OPTIONS.Qmax > 0 )
        % DO NOTHING: User supplied value or fine tuned value.
    elseif ( isempty(OPTIONS.Qmax) == true )
        % DO NOTHING: Qmax is not required by specified integrator.
    else
        str = [ 'Error: User selected Qmax: ', num2str(RCNTRL_U.Qmax), '. Qmax must be >= 0.' ];
        error(str);
    end    
    
    % Hstart
    if ( OPTIONS.Hstart == 0.0 )
        switch ( family )
            case 'ERK'
                OPTIONS.Hstart = max( OPTIONS.Hmin, roundOff );                
            case 'RK'
                OPTIONS.Hstart = 0.0;
            case 'ROS'
                OPTIONS.Hstart = max( OPTIONS.Hmin, 1d-5 );
            case 'SDIRK'
                OPTIONS.Hstart = max( OPTIONS.Hmin, roundOff );
        end
    elseif ( OPTIONS.Hstart > 0.0 )
        OPTIONS.Hstart = min( abs(OPTIONS.Hstart), abs( tspan(2)-tspan(1) ) );
    else
        str = [ 'Error: User selected Hstart: ', num2str(OPTIONS.Hstat), '. Hstart must be >= 0.' ];
        error(str);
    end    
    
    % ITOL
    if ( OPTIONS.ITOL == 0 )
        switch ( family )
            case 'ROS'
                OPTIONS.ITOL = true; % VectorTol=true
            otherwise
                OPTIONS.ITOL = 1;
        end
    elseif ( OPTIONS.ITOL == 1 )
        switch ( family )
            case 'ROS'
                OPTIONS.ITOL = false; % VectorTol=false
            otherwise
                OPTIONS.ITOL = 0;
        end
    else
        str = [ 'Error: User selected ITOL: ', num2str(ITOL), '. ITOL must be 0 or 1.' ];
        error(str);
    end    
    
    % ERK method number
    RK2 = 1;
    RK3 = 2;
    RK4 = 3;
    RK5 = 4;
    RK6 = 5;
    RK8 = 6;
    
    % RK method number
    R2A = 1;
    R1A = 2;
    L3C = 3;
    GAU = 4;
    L3A = 5;
    
    % ROS method number
    RS2 = 1;
    RS3 = 2;
    RS4 = 3;
    RD3 = 4;
    RD4 = 5;
    
    % SDIRK method number
    S2A = 1;
    S2B = 2;
    S3A = 3;
    S4A = 4;
    S4B = 5;    
    
    % Method
    switch ( family )
        case 'ERK'
            switch ( OPTIONS.Method )
                case 1
                    [ erkMethod, erkELO, erkS, erkName ] = Coefficients_Erk23( RK2 );
                case 2
                    [ erkMethod, erkELO, erkS, erkName ] = Coefficients_Erk3_Heun( RK3 );
                case 3
                    [ erkMethod, erkELO, erkS, erkName ] = Coefficients_Erk43( RK4 );
                case 4
                    [ erkMethod, erkELO, erkS, erkName ] = Coefficients_Dopri5( RK5 );                  
                case 5
                    [ erkMethod, erkELO, erkS, erkName ] = Coefficients_Verme( RK6 );
                case 6
                    [ erkMethod, erkELO, erkS, erkName ] = Coefficients_Dopri853( RK8 );
                case 7
                    [ erkMethod, erkELO, erkS, erkName ] = Coefficients_ForwardEuler( 1 ); % temp
                otherwise
                    [ erkMethod, erkELO, erkS, erkName ] = Coefficients_Dopri5( RK5 );
            end
            Coefficient.Method = erkMethod;
            Coefficient.ELO = erkELO;
            Coefficient.NStage = erkS;
            Coefficient.Name = erkName;
        case 'RK'
            switch ( OPTIONS.Method )
                case 1
                    [ rkMethod, rkELO, rkS, rkName ] = Coefficients_Radau2A( R2A, OPTIONS );
                case 2
                    [ rkMethod, rkELO, rkS, rkName ] = Coefficients_Lobatto3C( L3C, OPTIONS );
                case 3
                    [ rkMethod, rkELO, rkS, rkName ] = Coefficients_Gauss( GAU );
                case 4
                    [ rkMethod, rkELO, rkS, rkName ] = Coefficients_Radau1A( R1A );
                case 5
                    [ rkMethod, rkELO, rkS, rkName ] = Coefficients_Lobatto3A( L3A );
                case 6
                    [ rkMethod, rkELO, rkS, rkName ] = Coefficients_BackwardEuler( 1 ); % temp
                otherwise
                   [ rkMethod, rkELO, rkS, rkName ] = Coefficients_Lobatto3C( L3C, OPTIONS );
            end
            Coefficient.Method = rkMethod;
            Coefficient.ELO = rkELO;
            Coefficient.NStage = rkS;
            Coefficient.Name = rkName;
        case 'ROK'
            switch( OPTIONS.Method )
                otherwise
                    [ rokMethod, rokELO, rokS, rokName ] = Coefficients_ROK4b( 0 );
            end
            Coefficient.Method = rokMethod;
            Coefficient.ELO = rokELO;
            Coefficient.NStage = rokS;
            Coefficient.Name = rokName;
        case 'ROS'
            switch( OPTIONS.Method )
                case 1
                    [ rosMethod, rosELO, rosS, rosName ] = Coefficients_Ros2( RS2 );
                case 2
                    [ rosMethod, rosELO, rosS, rosName ] = Coefficients_Ros3( RS3 );
                case 3
                    [ rosMethod, rosELO, rosS, rosName ] = Coefficients_Ros4( RS4 );
                case 4
                    [ rosMethod, rosELO, rosS, rosName ] = Coefficients_Rodas3( RD3 );
                case 5
                    [ rosMethod, rosELO, rosS, rosName ] = Coefficients_Rodas4( RD4 );
                otherwise % default methods
                    switch ( implementation )
                        case 'FWD'
                            [ rosMethod, rosELO, rosS, rosName ] = Coefficients_Ros4( RS4 );
                        case 'TLM'
                            [ rosMethod, rosELO, rosS, rosName ] = Coefficients_Rodas4( RD4 );
                        case 'ADJ'
                            [ rosMethod, rosELO, rosS, rosName ] = Coefficients_Rodas3( RD3 );
                    end
            end
            Coefficient.Method = rosMethod;
            Coefficient.ELO = rosELO;
            Coefficient.NStage = rosS;
            Coefficient.Name = rosName;
        case 'SDIRK'
            switch ( OPTIONS.Method )
                case 1
                    [ sdirkMethod, sdirkELO, sdirkS, sdirkName ] = Coefficients_Sdirk2A( S2A );
                case 2
                    [ sdirkMethod, sdirkELO, sdirkS, sdirkName ] = Coefficients_Sdirk2B( S2B );
                case 3
                    [ sdirkMethod, sdirkELO, sdirkS, sdirkName ] = Coefficients_Sdirk3A( S3A );
                case 4
                    [ sdirkMethod, sdirkELO, sdirkS, sdirkName ] = Coefficients_Sdirk4A( S4A );
                case 5
                    [ sdirkMethod, sdirkELO, sdirkS, sdirkName ] = Coefficients_Sdirk4B( S4B );
                otherwise
                    [ sdirkMethod, sdirkELO, sdirkS, sdirkName ] = Coefficients_Sdirk4A( S4A );
            end    
            Coefficient.Method = sdirkMethod;
            Coefficient.ELO = sdirkELO;
            Coefficient.NStage = sdirkS;
            Coefficient.Name = sdirkName;
        case 'EXP'
            switch( OPTIONS.Method )
                case 1
                    OPTIONS.OneStepIntegrator = @exp4SingleStep;
                    expName = 'EXP4';
                    expMethod = 1;
                    expELO = 3;
                case 2
                    OPTIONS.OneStepIntegrator = @erow4SingleStep;
                    expName = 'EROW4';
                    expMethod = 2;
                    expELO = 3;
                otherwise
                    OPTIONS.OneStepIntegrator = @exp4SingleStep;
                    expName = 'EXP4';
                    expMethod = 1;
                    expELO = 3;
            end
            Coefficient.Method = expMethod;
            Coefficient.ELO = expELO;
            Coefficient.Name = expName;            
        case 'EXPK'
            switch( OPTIONS.Method )
                otherwise
                    [ expkMethod, expkELO, expkS, expkName ] = Coefficients_EXPK4( 0 );
            end
            Coefficient.Method = expkMethod;
            Coefficient.ELO = expkELO;
            Coefficient.NStage = expkS;
            Coefficient.Name = expkName;
    end
    
    % Max_no_steps
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % (default) ERK     RK      ROS         SDIRK
    %   FWD     10000   200000  200000      200000
    %   TLM     10000   200000  20000       200000
    %   ADJ     10000   20000   200000      10000
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if ( OPTIONS.Max_no_steps == 0 ) % Default Settings
        switch ( family )
            case 'ERK'
                OPTIONS.Max_no_steps = 10000;
            case 'RK'
                switch ( implementation )
                    case 'FWD'
                        OPTIONS.Max_no_steps = 200000;
                    case 'TLM'
                        OPTIONS.Max_no_steps = 200000;
                    case 'ADJ'
                        OPTIONS.Max_no_steps = 20000;
                end
            case 'ROK'
                switch ( implementation )
                    case 'FWD'
                        OPTIONS.Max_no_steps = 200000;
                    case 'TLM'
                       % DO NOTHING
                    case 'ADJ'
                       % DO NOTHING
                end
            case 'ROS'
                switch ( implementation )
                    case 'FWD'
                        OPTIONS.Max_no_steps = 200000; 
                    case 'TLM'
                        OPTIONS.Max_no_steps = 200000;
                    case 'ADJ'
                        OPTIONS.Max_no_steps = 200000;
                end
            case 'SDIRK'
                switch ( implementation )
                    case 'FWD'
                        OPTIONS.Max_no_steps = 200000;
                    case 'TLM'
                        OPTIONS.Max_no_steps = 200000;
                    case 'ADJ'
                        OPTIONS.Max_no_steps = 10000;
                end
            otherwise
                OPTIONS.Max_no_steps = 200000;
        end
    elseif ( OPTIONS.Max_no_steps > 0 )
        % DO NOTHING: User supplied or fine tuned value.
    else
        str = [ 'Error: User selected Max_no_steps: ', num2str(OPTIONS.Max_no_steps), '. Max_no_steps must be > 0.' ];
        error(str);
    end    
    
    % NewtonMaxIt
    if ( OPTIONS.NewtonMaxIt == 0 )
        OPTIONS.NewtonMaxIt = 8;
    elseif ( OPTIONS.NewtonMaxIt > 0 )
        % DO NOTHING: User supplied or fine tuned value.
    elseif ( isempty(OPTIONS.NewtonMaxIt) == true )
        % DO NOTHING: NewtonMaxIt is not required by specified integrator.
    else
        str = [ 'Error: User selected NewtonMaxIt: ', num2str(OPTIONS.NewtonMaxIt), '. NewtonMaxIt must be > 0.' ];
        error(str);
    end  
    
    % StartNewton
    if ( OPTIONS.StartNewton == false )
        OPTIONS.StartNewton = true;
    elseif ( OPTIONS.StartNewton == true )
        OPTIONS.StartNewton = false;
    elseif ( isempty(OPTIONS.StartNewton) == true )
        % DO NOTHING: StartNewton is not required by specified integrator.
    else
        str = [ 'Error: User selected StartNewton: ', num2str(OPTIONS.StartNewton), ', StartNewton must be 0 or 1.' ];
        error(str);
    end    
    
    % DirectTLM
    if ( OPTIONS.DirectTLM == false )
        OPTIONS.DirectTLM = false;
    elseif ( OPTIONS.DirectTLM == true )
        OPTIONS.DirectiTLM = true;
    elseif ( ~isempty( OPTIONS.DirectTLM ) )
        str = [ 'Error: User selected DirectTLM: ', num2str(OPTIONS.DirectTLM), '. DirectTLM must be 0 or 1.' ];
        error(str);
    end   
    
    % SaveLU
    if ( OPTIONS.SaveLU == false )
        OPTIONS.SaveLU = false;
    elseif ( OPTIONS.SaveLU == true )
        OPTIONS.SaveLU = true;
    elseif ( isempty( OPTIONS.SaveLU ) )
        % DO NOTHING: SaveLU is not required by specified integrator.
    else
        str = [ 'Error: User selected SaveLU: ', num2str(OPTIONS.SaveLU), '. SaveLU must be 0 or 1.' ];
        error(str);
    end    
    
    % TLMNewtonEst
    if ( OPTIONS.TLMNewtonEst == false )
        OPTIONS.TLMNewtonEst = false;
    elseif ( OPTIONS.TLMNewtonEst == true )
        OPTIONS.TLMNewtonEst = true;
    elseif ( isempty( OPTIONS.TLMNewtonEst ) )
        % DO NOTHING: TLMNewtonEst is not required by specified integrator.
    else
        str = [ 'Error: User selected TLMNewtonEst: ', num2str(OPTIONS.TLMNewtonEst), '. TLMNewtonEst must be 0 or 1.' ];
        error(str);
    end    
    
    % FDAprox
    if ( OPTIONS.FDAprox == false )
        OPTIONS.FDAprox = false;
    elseif ( OPTIONS.FDAprox == true )
        OPTIONS.FDAprox = true;
    elseif ( isempty( OPTIONS.FDAprox ) )
        % DO NOTHING: FDAprox is not required by specified integrator.
    else
        str = [ 'Error: User selected FDAprox: ', num2str(OPTIONS.FDAprox), '. FDAprox must be 0 or 1.' ];
        error(str);
    end    
    
    
    % AdjointSolve
    if( strcmp(implementation,'ADJ') )
        switch( family )
            case 'ERK'
            case 'RK'              
            case 'ROS'
            case 'SDIRK'
            otherwise
                error('Error: Internal error...');
        end
    end
        
    % DirectADJ
    if ( OPTIONS.DirectADJ == false )
        OPTIONS.DirectADJ = false;
    elseif ( OPTIONS.DirectADJ == true )
        OPTIONS.DirectADJ = true;
    elseif ( ~isempty( OPTIONS.DirectADJ ) )
       str = [ 'Error: User selected DirectADJ: ', num2str(OPTIONS.DirectADJ), '. DirectADJ must be 0 or 1.' ];
       error(str);
    end
    
    % ChunkSize
    if ( OPTIONS.ChunkSize == 0 )
        OPTIONS.ChunkSize = 500;
    elseif ( OPTIONS.ChunkSize > 0 )
        % DO NOTHING: Use user supplied value.        
    else
        str = [ 'Error: User selected ChunkSize: ', num2str(OPTIONS.ChunkSize), '. ChunkSize must be >= 0.' ];
        error(str);
    end    
    
    % storeCheckpoint
    if ( OPTIONS.storeCheckpoint == false )
        OPTIONS.storeCheckpoint = false;
    elseif ( OPTIONS.storeCheckpoint == true )
        OPTIONS.storeCheckpoint = true;
    else
        str = [ 'Error: User selected storeCheckpoint: ', num2str(OPTIONS.storeCheckpoint), '. storeCheckpoint must be 0 or 1.' ];
        error(str);
    end
    
    % TLMTruncErr
    if ( OPTIONS.TLMTruncErr == false )
        OPTIONS.TLMTruncErr = false;
    elseif ( OPTIONS.TLMTruncErr == true )
        OPTIONS.TLMTruncErr = true;
    elseif ( isempty( OPTIONS.TLMTruncErr ) )
        % DO NOTHING: TLMTruncErr is not required by specified integrator.
    else
        str = [ 'Error: User selected TLMTruncErr: ', num2str(OPTIONS.TLMTruncErr), '. TLMTruncErr must be 0 or 1.' ];
        error(str);
    end


end

%%
% <html>
%   <div>
%       <img style="float: right" src="../../../CSL_LogoWithName_1.png" height="50px"></img>
%   </div>
% </html>