%% ROS_ADJ2_DiscreteIntegrator
%
% <html>
%   <div>
%       <img style="float: right" src="../../../MATLODE_LOGO.png" height="150px"></img>
%   </div>
% </html>
%
% <html> Up: <a href="../../../Library/html/Library.html">Library</a> </html>
%
%% Syntax
%    [ Tout, Yout, ISTATUS, RSTATUS, Ierr, stack_ptr, quadrature ] = ROS_ADJ2_DiscreteIntegrator( OdeFunction, Tspan, Y, OPTIONS, Coefficient, stack_ptr )
%
%% Input Parameters
% |OdeFunction|: ODE function function handle
%
% |Tspan|: Time interval
%
% |OPTIONS|: Option struct
%
% |Coefficients|: Constant coefficients associated with method
%
% |stack_ptr|: pointer for global snapshot stack
%
% |adjQuadFlag|: Adjoint quadrature stack flag
%
%% Output Parameters
% |Tout|: Time vector
%
% |Yout|: State vector
%
% |Lambda|: Sensitivity matrix
%
% |ISTATUS|: Integer statistics 
%
% |RSTATUS|: Real statistics
%
% |Ierr|: Error flag
%
%% Description
% Rosenbrock Runge-Kutta discrete 2 adjoint core method.
%
%% Contact Information
%%
% Dr. Adrian Sandu                 | Phone: (540) 231-2193 | Email: sandu@cs.vt.edu
%%
% Tony D'Augustine                 | Phone: (540) 231-6186 | Email: adaug13@vt.edu 
%%
% Computational Science Laboratory | Phone: (540) 231-6186
%% Reference
% [1] Tony D'Augustine, Adrian Sandu. MATLODE: A MATLAB ODE Solver and
%     Sensitivity Analysis Toolbox. Submitted to ACM TOMS.
%
% [2] Hong Zhang, Adrian Sandu. FATODE: a library for forward, adjoint and 
%     tangent linear integration of ODEs, SIAM Journal on Scientific 
%     Computing, 36(5), C504-C523, 2014.
%
function [ Tout, Yout, Lambda, ISTATUS, RSTATUS, Ierr ] = ROS_ADJ2_DiscreteIntegrator( NVAR, OPTIONS, Coefficient, stack_ptr, adjQuadFlag, lambda_Tf )

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   Global Variables
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    global ros_M ros_Gamma ros_Alpha
    global ros_A ros_C

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   Local Variables
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    % note: based off DLAMCH('E') NOT MATLAB's EPS
    Roundoff = 1.11022302462515654E-016;
    
    DeltaMin = 1e-5;
    
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Required
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
    if ( isempty(OPTIONS.Lambda) )
        error('User defined option parameter Lambda is required.');
    end   
    if ( isempty(OPTIONS.Hesstr_vec) )
        error('User defined option parameter Hesstr_vec is required.');
    end
    
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   Initializations
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Lambda = lambda_Tf;

    Direction = 1;
    
    TYindex = 1;
    
    ros_W = zeros(Coefficient.NStage,1);
    
    ISTATUS = ISTATUS_Struct('default');
    RSTATUS = RSTATUS_Struct('default');    
    
    U = zeros(NVAR*Coefficient.NStage,OPTIONS.NADJ);
    V = zeros(NVAR*Coefficient.NStage,OPTIONS.NADJ);
    
    RY0 = zeros(NVAR,OPTIONS.NADJ);
    RY1 = zeros(NVAR,OPTIONS.NADJ);

    % Compute ros_W
    if ( adjQuadFlag ) 
        GammaW = 0;
        for istage=Coefficient.NStage:-1:1
            ros_W(istage) = ros_Gamma(1)*ros_M(istage);
            for j=istage+1:Coefficient.NStage
                HC = ros_C((j-1)*(j-2)/2+istage);
                ros_W(istage) = ros_W(istage)+HC*ros_Gamma(1)*ros_W(j);
            end
            GammaW = GammaW + ros_W(istage)*ros_Gamma(istage);
        end
    end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   Time loop begins
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    while ( stack_ptr > 0 )
        % Recover checkpoints for stage values an vectors
        [ T, H, Ystage, K, stack_ptr ] = rosPop(NVAR,Coefficient.NStage,stack_ptr,OPTIONS.displaySteps);
        Yout(:,TYindex) = Ystage;
        Tout(TYindex,1) = T;        
        TYindex = TYindex + 1;
        
        ISTATUS.Nstp = ISTATUS.Nstp + 1;
        
        % Compute LU decomposition
        if ( ~OPTIONS.SaveLU )
            fjac = OPTIONS.Jacobian(T,Ystage(:));
            ISTATUS.Njac = ISTATUS.Njac + 1;
            
            Tau = 1.0/(Direction*H*ros_Gamma(1));
            
            % Perform decomposition. Note: we do not check singularlity
            e = -fjac + eye(NVAR,NVAR)*Tau;
            ISTATUS.Ndec = ISTATUS.Ndec + 1;
            
        end
        
        if ( adjQuadFlag )
            RY0 = OPTIONS.DRDY(T,Ystage(:));
        end
        
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   Stages
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        for istage=Coefficient.NStage:-1:1
            % Current istage first entry
            istart = NVAR*(istage-1) + 1;
            
            % Compute U
            for m=1:OPTIONS.NADJ
                U(istart:istart+(NVAR-1),m) = Lambda(:,m);
                U(istart:istart+(NVAR-1),m) = U(istart:istart+(NVAR-1),m)*ros_M(istage);
                if ( adjQuadFlag )
                    U(istart:istart+(NVAR-1),m) = U(istart:istart+(NVAR-1),m) + H*ros_W(istage)*RY0(:,m);
                end
            end
            
            for j=istage+1:Coefficient.NStage
                jstart = NVAR*(j-1) + 1;
                HA = ros_A((j-1)*(j-2)/2+istage);
                HC = ros_C((j-1)*(j-2)/2+istage)/(Direction*H);
                for m=1:OPTIONS.NADJ
                    U(istart:istart+(NVAR-1),m) = U(istart:istart+(NVAR-1),m) + HA*V(jstart:jstart+(NVAR-1),m);
                    U(istart:istart+(NVAR-1),m) = U(istart:istart+(NVAR-1),m) + HC*U(jstart:jstart+(NVAR-1),m);
                end
            end
            
            for m=1:OPTIONS.NADJ
                U(istart:istart+(NVAR-1),m) = transpose(e)\U(istart:istart+(NVAR-1),m); %<--------check this...
                ISTATUS.Nsol = ISTATUS.Nsol + 1;
            end
            
            % Compute V
            Tau = T + ros_Alpha(istage)*Direction*H;
            fjac = OPTIONS.Jacobian(Tau,Ystage(istart:istart+(NVAR-1)));
            ISTATUS.Njac = ISTATUS.Njac + 1;
            
            if ( adjQuadFlag )
                RY1 = OPTIONS.DRDY(Tau,Ystage(istart:istart+(NVAR-1)));
            end
            
            for m=1:OPTIONS.NADJ
                V(istart:istart+(NVAR-1),m) = transpose(fjac)*U(istart:istart+(NVAR-1),m);
                if ( adjQuadFlag )
                    V(istart:istart+(NVAR-1),m) = V(istart:istart+(NVAR-1),m) + H*ros_W(istage)*RY1(:,m);
                end
            end
            
        end % stage
        
        % Compute Lambda
        for istage=1:Coefficient.NStage
            istart = NVAR*(istage-1)+1;
            for m=1:OPTIONS.NADJ
                Lambda(:,m) = Lambda(:,m) + V(istart:istart+(NVAR-1),m);
                TMP = OPTIONS.Hesstr_vec(T,Ystage(istart:istart+(NVAR-1)),U(istart:istart+(NVAR-1),m),K(istart:istart+(NVAR-1)));
                Lambda(:,m) = Lambda(:,m) + TMP;
                if ( adjQuadFlag )
                    TMP = OPTIONS.Hesstr_vec_r(T,Ystage,H*ros_W(istage),K(istart:istart+(NVAR-1)));
                    Lambda(:,m) = Lambda(:,m) + TMP;
                end
            end
        end
        
        % NonAutonomous term
        if ( ~OPTIONS.Autonomous )
            Delta = sqrt(Roundoff)*max(DeltaMin,abs(T));
            djdt = OPTIONS.Jacobian(T+Delta,Ystage(:));
            ISTATUS.Njac = ISTATUS.Njac + 1;
            
            % Jacobian time derivatie
            djdt = djdt - fjac;
            djdt = djdt/Delta;
            
            for m=1:OPTIONS.NADJ
                Tmp = zeros(NVAR,1);
                for istage=1:Coefficient.NStage
                    istart = NVAR*(istage-1)+1;
                    Tmp = Tmp + ros_Gamma(istage)*U(istart:istart+(NVAR-1),m);
                end
                Tmp2 = transpose(djdt)*Tmp;
                Lambda(:,m) = Lambda(:,m) + H*Tmp2;
            end
        end 
        
        % Quadautonomous <----- false (temp fix...need to add to ICNTRL)
        if ( adjQuadFlag && true )
            RY1 = OPTIONS.DRDY(T+Delta,Ystage(1));
            RY1 = RY1 - RY0;
            RY1 = RY1/Delta;
            Lambda = Lambda + H*GammaW*RY1;
        end
                
    end % time loop

    % Successful exit
    Ierr = 1;

return;

function [ T, H, Y, K, stack_ptr ] = rosPop( NVAR, rkS, stack_ptr, displaySteps )
    global chk_H chk_T chk_Y chk_K

    if ( stack_ptr <= 0 )
        disp('Pop failed: empty buffer');
    end
        
    H = chk_H(stack_ptr);
    T = chk_T(stack_ptr);
    Y = chk_Y(1:NVAR*rkS,stack_ptr);
    K = chk_K(1:NVAR*rkS,stack_ptr);
    
    % for debugging
    if ( displaySteps == true )
        str = ['Backtracking. Time = ', num2str(T), '; Stepsize = ', num2str(H)];
        disp(str);
    end      
    
    stack_ptr = stack_ptr - 1;

return;
