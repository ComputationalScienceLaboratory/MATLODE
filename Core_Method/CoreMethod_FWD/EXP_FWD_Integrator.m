%% EXP_FWD_Integrator
%
% <html> Up: <a href="../../../Library/html/Library.html">Library</a> </html>
%
%% Syntax
%
%
%% Input Parameters
%
%
%% Output Parameters
%
%
%% Description
%
%
%% Reference
% [1] Tony D'Augustine, Adrian Sandu. MATLODE
%
% [2] Hong Zhang, Adrian Sandu. FATODE: a library for forward, adjoint and 
%     tangent linear integration of ODEs, SIAM Journal on Scientific 
%     Computing, 36(5), C504-C523, 2014.
%
function [ Tout, Yout, ISTATUS, RSTATUS, Ierr, stack_ptr, quadrature ] = EXP_FWD_Integrator( OdeFunction,...
        Tspan, Y, OPTIONS, Coefficient, adjStackFlag, adjQuadFlag )

    % Force initial value matrix to be N X 1.
    if ( size(Y,2) == 1 )
        % DO NOTHING
    else
        Y = transpose(Y);
    end

    % Get Problem Size
    NVAR = max(size(Y));

    % Initialize time variables
    Tinitial = Tspan(1);
    Tfinal = Tspan(2);

    % note: based off DLAMCH('E') NOT MATLAB's EPS
    Roundoff = 1.11022302462515654E-016;

    DeltaMin = 1e-5;

    if ( OPTIONS.storeCheckpoint == true )
        [ Tout, Yout ] = matlOde_allocateMemory( NVAR, OPTIONS.ChunkSize );
%         Yout = zeros(NVAR,OPTIONS.Max_no_steps);
%         Tout = zeros(OPTIONS.Max_no_steps,1);
    else
        Yout = zeros(NVAR,1);
        Tout = 0;
    end
    TYindex = 1;

    Yout(:,TYindex) = Y;
    Tout(TYindex,1) = Tinitial;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   Global Variables
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    global alpha gamma c gam e b
    % M = 4;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   Initial Settings
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ISTATUS = ISTATUS_Struct('default');
    RSTATUS = RSTATUS_Struct('default');

    stack_ptr = 0;

%     K = zeros(NVAR,Coefficient.NStage);
%     lambda = zeros(OPTIONS.NBasisVectors,Coefficient.NStage);


    quadrature = OPTIONS.Quadrature();

    T = Tinitial;
    H = min( max( abs(OPTIONS.Hmin), abs(OPTIONS.Hstart) ), abs(OPTIONS.Hmax) );
    if ( abs(H) <= 10.0*Roundoff )
        H = DeltaMin;
    end

    if ( Tfinal >= Tinitial )
        Direction = 1;
    else
        Direction = -1;
    end
    H = Direction*H;

    RejectLastH = false;
    RejectMoreH = false;

    % The following variable is used to check the output of
    % embedded error controller/error controller behaviour
    %    Yerr_nrm = zeros(3,1);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   Time loop
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    while ( ( Direction > 0 ) && ( (T-Tfinal)+Roundoff <= 0.0 ) ...
            || ( Direction < 0 ) && ( (Tfinal-T)+Roundoff <= 0.0 ) )

        if ( ISTATUS.Nstp >= OPTIONS.Max_no_steps )
            error('Number of steps exceeds maximum buffer bound \n T= %f;     H= %f', T, H );
        end
        if ( ( ( T+0.1*H ) == T) || ( H <= Roundoff ) )
            error('Step size too small: T + 10*H = T or H < Roundoff \n T= %f;     H= %f', T, H );
        end

        % Limit H if necessary to avoid going beyond Tfinal
        H = min( H, abs( Tfinal-T ) );

        % Compute the function at current time
        Fcn0 = OdeFunction( T, Y );
        ISTATUS.Nfun = ISTATUS.Nfun + 1;

        % Compute the function derivative with respect to T
        if ( ~OPTIONS.Autonomous )
            [ dFdT, ISTATUS ] = fatOde_FunctionTimeDerivative( T, Roundoff, Y, Fcn0, OdeFunction, ISTATUS );
        end

        % Compute the Jacobian at current time
        if ( ~OPTIONS.MatrixFree )
            fjac = OPTIONS.Jacobian(T,Y);
            ISTATUS.Njac = ISTATUS.Njac + 1;
        else
            if( ~isempty(OPTIONS.Jacobian) )
                fjac = @(vee)OPTIONS.Jacobian(T,Y,vee);
            else
                normy = norm(Y);
                fjac = @(vee)Mat_Free_Jac(T,Y,vee,OdeFunction,Fcn0,normy);
            end
        end

        % Repeat step calculation until current step accepted
        accepted = false;
        while ( ~accepted ) % accepted

%             [Varn, Harn, w] = Arnoldi_MF(fjac, Fcn0, dFdT, NVAR, OPTIONS.NBasisVectors, OPTIONS.MatrixFree); % M should be an option!!!!!!!1

%             [ H, ISING, e, ISTATUS ] = fatOde_ROS_PrepareMatrix( OPTIONS.NBasisVectors, H, Direction, gam, Harn, ISTATUS );

%            if ( ISING ~= 0 ) % More than 5 consecutive failed decompositions
%                Ierr = fatOde_ROS_ErrorMessage( -8, T, H );
%                return;
%            end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   Stages
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
           if( ~OPTIONS.Autonomous )

                y = [Y; T];
                f = [Fcn0; 1];
                rhsFun = @(why)([OdeFunction(why(end), why(1:end-1)); 1]);
                epsilon = 1e-8;
                dFdT = OdeFunction(T + epsilon, Y) - Fcn0;

                if( ~OPTIONS.MatrixFree )
                    J = [fjac, dFdT; zeros(1,NVAR+1)];
                else
                    J = @(why)([fjac(why(1:end-1)) + dFdT*why(end); 0]);
                end

            else

                rhsFun = @(why)odeFun(T, why);
                y = Y;
                f = Fcn0;
                if( ~OPTIONS.MatrixFree )
                    J = jacFun(T, Y);
                else
                    J = fjac;
                end
           end

           if ~isempty(OPTIONS.IsJACSymm)
               symmjac = OPTIONS.IsJACSymm;
           else
               symmjac = false;
           end              
           
           [ynew, yerr, ISTATUS] = OPTIONS.OneStepIntegrator(y,H,rhsFun, ...
                    J, f, OPTIONS.MatrixFree, OPTIONS.NBasisVectors, ISTATUS, OPTIONS.AbsTol, OPTIONS.RelTol, OPTIONS.Adaptive_Krylov, symmjac);
        
           % disp(norm(ynew));

            % The following can be used to test the time-step
            % control and embedded method error computation
            % using ODE45 to compute the error.
            % [t45, y45] = ode45(OdeFunction,[T, T+H], Y, ...
            %                    odeset('RelTol',1e-12,'AbsTol',1e-12));

            if( ~OPTIONS.Autonomous )
                Ynew = ynew(1:end-1);
                Yerr = yerr(1:end-1);
            %     Yerr_nrm(1,end + 1) = norm(Yerr, 2);
            %     Yerr = ynew(1:end -1) - y45(end,:)';
            %     Yerr_nrm(2,end) = norm(Yerr, 2);
            %     Yerr_nrm(3,end) = norm(Yerr - yerr(1:end - 1), 2);
            %     Yerr = yerr(1:end-1);
            else
                Ynew = ynew;
                Yerr = yerr;
            %     Yerr_nrm(1,end + 1) = norm(Yerr, 2);
            %     Yerr = ynew - y45(end, :)';
            %     Yerr_nrm(2,end) = norm(Yerr, 2);
            %     Yerr_nrm(3,end) = norm(Yerr - yerr, 2);
            %     Yerr = yerr;
            end

           % display(sprintf('K(20791,:)= %f %f %f %f %f %f %f', K(20791,:)));
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   Error Estimation
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            % Compute error norm
            SCAL = OPTIONS.AbsTol + OPTIONS.RelTol.*max(abs(Y),abs(Ynew));
            Err = max(sqrt(sum((Yerr./SCAL).^2)/NVAR),1d-10);

            % New step size is bounded by FacMin <= Hnew/H <= FacMax
            Fac = min( OPTIONS.FacMax, max( OPTIONS.FacMin, ...
                OPTIONS.FacSafe/Err^(1.0/Coefficient.ELO) ) );
            Hnew = H*Fac;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   Accept/Reject
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            % Check the error magnitude and adjust the step size
            ISTATUS.Nstp = ISTATUS.Nstp + 1;
            if ( ( Err <= 1.0 ) || ( H <= OPTIONS.Hmin ) )
                ISTATUS.Nacc = ISTATUS.Nacc + 1;

                % Update time
                T = T + Direction*H;

                % Update solution
                Y = Ynew;

                if ( OPTIONS.storeCheckpoint == true )
                    Tout(TYindex,1) = T;
                    Yout(:,TYindex) = Y;
                    TYindex = TYindex + 1;
                end

                % Update Memory Allocation
                if ( (ISTATUS.Nacc > OPTIONS.ChunkSize*ISTATUS.Nchk) && (OPTIONS.storeCheckpoint == true) )
                    [ Tout, Yout ] = matlOde_appendMemory(NVAR,Tout,Yout,OPTIONS.ChunkSize);
                    ISTATUS.Nchk = ISTATUS.Nchk + 1;
                end

                Hnew = max( OPTIONS.Hmin, min( Hnew, OPTIONS.Hmax ) );
                if ( RejectLastH ) % No step size increase after a rejected step
                    Hnew = min( Hnew, H );
                end

                % Last T and H
                RSTATUS.Nhexit = H;
                RSTATUS.Nhnew = Hnew;
                RSTATUS.Ntexit = T;

                RejectLastH = false;
                RejectMoreH = false;
                H = Hnew;
                accepted = true;

                % for debugging
                if ( OPTIONS.displaySteps == true )
                    str = ['Accepted step. Time = ', num2str(T), '; Stepsize = ', num2str(H)];
                    disp(str);
                end

            else % Reject step
                if ( RejectMoreH )
                    Hnew = H*OPTIONS.FacRej;
                end

                RejectMoreH = RejectLastH;
                RejectLastH = true;
                H = Hnew;

                if ( ISTATUS.Nacc >= 1 )
                    ISTATUS.Nrej = ISTATUS.Nrej + 1;
                end

                % for debugging
                if ( OPTIONS.displaySteps == true )
                    str = ['Rejected step. Time = ', num2str(T), '; Stepsize = ', num2str(H)];
                    disp(str);
                end

            end
        end
    end

    % Successful return
    Ierr = 1;

    % Decallocate Memory
    if ( OPTIONS.storeCheckpoint == true )
        Tout(ISTATUS.Nacc+1:end,:) = [];
        Yout(:,ISTATUS.Nacc+1:end) = [];
    else
        Tout = T;
        Yout = Y;
    end
    Yout = Yout';

    % Plot norm of errors
    %    figure();
    %semilogy(Yerr_nrm(1, :),'b--o')
    % hold on
    % semilogy(Yerr_nrm(2, :),'g--x')
    % semilogy(Yerr_nrm(3, :),'r*')
    % hold off
    % title(strcat(num2str(OPTIONS.Method),':',num2str(int8(-log10(OPTIONS.AbsTol)))))
return;

%% Major Modification History
% <html>
% <table border=1>
%   <tr>
%       <td><b>Date</b></td>
%       <td>Developer</td>
%       <td>Email</td>
%       <td>Action</td>
%   </tr>
%   <tr>
%       <td>1/1/2014</td>
%       <td>Tony D'Augustine</td>
%       <td>adaug13@vt,edu</td>
%       <td>Release MATLODE_v2.0.00</td>
%   </tr>
% </table>
% </html>
%
