%% ERK_TLM_Integrator
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
%    [ Tout, Yout, ISTATUS, RSTATUS, Ierr ] = ERK_TLM_Integrator( OdeFunction, Tspan, Y0, OPTIONS, Coefficient)
%
%% Input Parameters
% |OdeFunction|: ODE function function handle
%
% |Tspan|: Time interval
%
% |Y0|: Initial state vector
%
% |OPTIONS|: Option struct
%
% |Coefficients|: Constant coefficients associated with method
%
%% Output Parameters
% |Tout|: Time vector
%
% |Yout|: State vector
%
% |ISTATUS|: Integer statistics 
%
% |RSTATUS|: Real statistics
%
% |Ierr|: Error flag
%
%% Description
% Explicit Runge-Kutta tangent linear core method.
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
function [ Tout, Yout, Y_TLM, ISTATUS, RSTATUS, Ierr ] = ERK_TLM_Integrator( OdeFunction,...
        Tspan, Y, OPTIONS, Coefficient )

    % Force initial value matrix to be N X 1.
    if ( size(Y,2) == 1 )
        % DO NOTHING
    else
        Y = transpose(Y);
    end  

    % Get Problem Size
    NVAR = max(size(Y));
    NTLM = max(size(OPTIONS.Y_TLM));
    
    Tinitial = Tspan(1);
    Tfinal   = Tspan(2);

    Roundoff = 1.11022302462515654E-016;
    
    ISTATUS = ISTATUS_Struct('default');
    RSTATUS = RSTATUS_Struct('default');
    
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Global variables
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    global rkA rkB rkC rkE
    
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Local variables
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    K = zeros( NVAR, Coefficient.NStage );
    K_TLM = zeros( NVAR, Coefficient.NStage, NTLM ); % TLM weights?
    
    Yout = zeros(NVAR,OPTIONS.Max_no_steps);
    Tout = zeros(OPTIONS.Max_no_steps,1);
    TYindex = 1;
    
    Yout(:,TYindex) = Y;
    Tout(TYindex,1) = Tinitial;
             
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Required
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
    if ( isempty(OPTIONS.Y_TLM) )
        error('User defined option parameter Y_TLM is required.');
    end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~        
% ---> Initializations
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    T = Tinitial;
    Tdirection = sign( Tfinal-Tinitial );
    H = max( abs(OPTIONS.Hmin), abs(OPTIONS.Hstart) );
    if ( abs(H) <= 10.0*Roundoff ) 
        H = 1.0d-6;
    end
    H = min( abs(H), OPTIONS.Hmax );
    H = H*sign(Tdirection);
    Reject = false;
    FirstStep = true;
    if( ~OPTIONS.FDAprox )
        FJAC = zeros(NVAR,NVAR); % ALLOCATE( FJAC( NVAR, NVAR ), STAT=i )
    end

    % Determine scaling factor for integration
    SCAL = 1.0 ./ ( OPTIONS.AbsTol + OPTIONS.RelTol .* abs(Y) );
    
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ---> Time loop begins
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    while( (Tfinal-T)*Tdirection - Roundoff > 0.0 ) % Tloop
        
        if ( ISTATUS.Nstp > OPTIONS.Max_no_steps )
            error('Number of steps exceeds maximum bound');
        end
        if ( ( T + 0.1*H == T ) || ( abs(H) <= Roundoff ) )
            error('Step size too small: T + 10*H = T or H < Roundoff');
        end
        
        for istage = 1:Coefficient.NStage % stages           
            K(:,istage) = 0;
            TMP = zeros( NVAR, 1 );

            if ( istage > 1 )
                for j = 1:istage-1
                    TMP = TMP + H*rkA(istage,j)*K(:,j);
                end
            end

            TMP = TMP + Y;
            
            K(:,istage) = OdeFunction( T+rkC(istage)*H,TMP);
            ISTATUS.Nfun = ISTATUS.Nfun + 1;

            if ( ~OPTIONS.FDAprox )
                FJAC = OPTIONS.Jacobian( T+rkC(istage)*H, TMP );
                ISTATUS.Njac = ISTATUS.Njac + 1;
            end

            S_TLM = OPTIONS.Y_TLM;
            for itlm = 1:NTLM % TlmL
                if ( istage > 1 )
                    for j = 1:istage-1
                        S_TLM(:,itlm) = S_TLM(:,itlm) + H*rkA(istage,j)*K_TLM(:,j,itlm);
                    end
                end
            end
            
            if ( ~OPTIONS.FDAprox )
                K_TLM(:,istage,:) = FJAC*S_TLM;
            else
%                [ K_TLM(:,istage,itlm), ISTATUS ] = ...
%                    Forward_FD_JacV( T+rkC(istage)*H, TMP, S_TLM, ...
%                    K(:,istage), OdeFunction, ISTATUS);
%                
%                 [ K_TLM(:,istage,itlm), ISTATUS ] = ...
%                     Complex_FD_JacV( T+rkC(istage)*H, TMP, S_TLM, ...
%                     OdeFunction, ISTATUS);
%                 
%                 [ TMP, K_TLM(:,istage,itlm), ISTATUS ] = ...
%                     Central2_FD_JacV( T+rkC(istage)*H, TMP, S_TLM, ...
%                     OdeFunction, ISTATUS);
%                 
%                 [ TMP, K_TLM(:,istage,itlm), ISTATUS ] = ...
%                     Central4_FD_JacV( T+rkC(istage)*H, TMP, S_TLM, ...
%                     OdeFunction, ISTATUS);
%                 
%                 [ TMP, K_TLM(:,istage,itlm), ISTATUS ] = ...
%                     Left4_FD_JacV( T+rkC(istage)*H, TMP, S_TLM, ...
%                     OdeFunction, ISTATUS);
                error('Finite Difference Approximation is not implemented yet');
            end
            
        end % stages

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ---> Error estimation
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ISTATUS.Nstp = ISTATUS.Nstp + 1;
        TMP = zeros(NVAR,1);
        for i = 1:Coefficient.NStage
            if ( rkE(i) ~= 0.0 )
                TMP = TMP + H*rkE(i)*K(:,i); %<--------prone to cancellation errors
            end
        end
        
        % Perform error norm
        Err = sum((TMP.*SCAL).^2,1);
        Err = max( sqrt( Err/double(NVAR) ), 1.0d-10 );        

        if ( OPTIONS.TLMTruncErr )
            Yerr_TLM = zeros(NVAR,NTLM);
            for itlm = 1:NTLM
                for j = 1:Coefficient.NStage
                    if ( rkE(j) ~= 0.0 )
                        Yerr_TLM(:,itlm) = Yerr_TLM(:,itlm) + H*rkE(j)*K_TLM(:,j,itlm);
                    end
                end
            end
            
            % Perform sensitivity error norm
            SCAL_TLM = 1.0 ./ OPTIONS.AbsTol_TLM+OPTIONS.RelTol_TLM .* abs(Yerr_TLM);
            Err = max([Err, sqrt(sum((OPTIONS.Y_TLM.*SCAL_TLM).^2,1)./NVAR), 1.0d-10] );
        end

        % Computation of new step size Hnew
        Fac = OPTIONS.FacSafe*(Err)^(-1.0/Coefficient.ELO);
        Fac = max( OPTIONS.FacMin, min( OPTIONS.FacMax, Fac ) );
        Hnew = H*Fac;
        
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ---> Accept/Reject step
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % accept
        if ( Err < 1.0 ) % step is accepted
            FirstStep = false;
            ISTATUS.Nacc = ISTATUS.Nacc + 1;

            % Update time
            T = T + H;
            
            % Update solution
            Y_TLM = OPTIONS.Y_TLM;
            for i = 1:Coefficient.NStage
                if ( rkB(i) ~= 0.0 )
                    Y = Y + H*rkB(i)*K(:,i);
                    for itlm = 1:NTLM
                        Y_TLM(:,itlm) = Y_TLM(:,itlm) + H*rkB(i)*K_TLM(:,i,itlm);
                    end                    
                end
            end
            OPTIONS.Y_TLM = Y_TLM;
            
            % Store checkpoint values
            if ( OPTIONS.storeCheckpoint == true )
                Tout(TYindex,1) = T;            
                Yout(:,TYindex) = Y;
                TYindex = TYindex + 1;
            end

            % update scaling coefficients
            SCAL = 1.0 ./ ( OPTIONS.AbsTol + OPTIONS.RelTol .* abs(Y) );

            % Next time step
            Hnew = Tdirection*min( abs(Hnew), OPTIONS.Hmax );
    
            % Last T and H
            RSTATUS.Ntexit = T;
            RSTATUS.Nhexit = H;
            RSTATUS.Nhnew = Hnew;
            
            % No step increase after a rejection
            if ( Reject )
                Hnew = Tdirection*min( abs(Hnew), abs(H) );
            end
            Reject = false;
            
            % COMMENT HERE:
            if( ( T + Hnew/OPTIONS.Qmin - Tfinal )*Tdirection > 0.0 )
                H = Tfinal - T;
            else
                % Hratio = Hnew/H; % Hration is not used.
                H = Hnew;
            end

            % for debug
            if ( OPTIONS.displaySteps == true )
                str = [ 'Accpeted step. Time = ', num2str(T), ';   Stepsize = ', num2str(H) ];
                disp(str);
            end
            
        else % step is rejected
            if ( FirstStep || Reject )
                H = OPTIONS.FacRej*H;
            else
                H = Hnew;
            end
            Reject = true;
            
            % Update Number of rejections
            if ( ISTATUS.Nacc >= 1.0 )
                ISTATUS.Nrej = ISTATUS.Nrej + 1;
            end
            
            %for debug
            if ( OPTIONS.displaySteps == true )
                str = [ 'Rejected step. Time = ', num2str(T), ';   Stepsize = ', num2str(H) ];
                disp(str);
            end
            
        end % accept
    end % time loop
    

    % Successful return
    Ierr = 1;
    
    % Deallocate Memory
    if ( OPTIONS.storeCheckpoint == true )
        Tout(TYindex:OPTIONS.Max_no_steps,:) = [];
        Yout(:,TYindex:OPTIONS.Max_no_steps) = [];
    else
        Tout = T;
        Yout = Y;
    end

    % Temporary Fix
    Yout = transpose(Yout);
            
return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% END: ERK_Integrator_TLM: CORE METHOD %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--------------------------------------------------------------------------









































