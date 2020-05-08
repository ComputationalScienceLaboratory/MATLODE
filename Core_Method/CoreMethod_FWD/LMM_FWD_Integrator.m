%% LMM_FWD_Integrator
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
% [1] Tony D'Augustine and Adrian Sandu. MATLODE.
%
% [2] P. Tranquilli and A. Sandu, 'Rosenbrock-Krylov Methods for Large
%     Systems of Differential Equations,' SIAM Journal on Scientific 
%     Computing 36(3), A1313–A1338 
%
function [ Tout, Yout, ISTATUS, RSTATUS, Ierr, stack_ptr, quadrature ] = LMM_FWD_Integrator( OdeFunction,...
        Tspan, Y0, OPTIONS, LMM_struct, ~, ~ )
    
    %XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX%
    %OPTIONS.MaxOrder = 5;
    %XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX%
    
    
    % Force initial value matrix to be 1 X N.
    if ( size(Y0,2) == 1 )
        % DO NOTHING
    else
        Y0 = Y0';
    end  

    % Get Problem Size
    NVAR = max(size(Y0));

    % Initialize time variables
    Tinitial = Tspan(1);
    Tfinal = Tspan(2);
        
    % Roundoff error specific to machine
    Roundoff = eps/2;
    
    DeltaMin = 1e-5;
    
    if ( OPTIONS.storeCheckpoint == true )
        [ Tout, Yout ] = matlOde_allocateMemory( NVAR, OPTIONS.ChunkSize );  
    else
        Yout = zeros(NVAR,1);
        Tout = 0;
    end
    TYindex = 1;
    
    Yout(:,TYindex) = Y0;
    Tout(TYindex,1) = Tinitial;    
    
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   Global Variables
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   Initial Settings
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ISTATUS = ISTATUS_Struct('default');
    RSTATUS = RSTATUS_Struct('default');
    
    if ( OPTIONS.storeCheckpoint )
        RSTATUS.stepsizes = zeros(OPTIONS.ChunkSize, 1);
    end
    
    stack_ptr = 0;
    
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
    
    Y = Y0;
    
    % Compute the function at initial time
    F = OdeFunction( T, Y );
    ISTATUS.Nfun = ISTATUS.Nfun + 1;
    
    % Compute the function derivative with respect to T
    if ( ~OPTIONS.Autonomous )
        [ dFdT, ISTATUS ] = fatOde_FunctionTimeDerivative( T, Roundoff, Y, F, OdeFunction, ISTATUS );
    else
        dFdT = zeros(NVAR,1);
    end
    
    Order = 1;
    
    % Preallocate integrator internal state
    Coefficients = LMM_struct.coefficients();
    LMM_state    = LMM_struct.stateInit(OPTIONS.MaxOrder, NVAR, Y, F, H);
    
    SkipJac = false;
    RejectLastH = false;
    RejectMoreH = false;
    NConsecutive = 0; % Count consecutive steps with identical stepsize and order.
        
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   Time loop
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    while ( ( Direction > 0 ) && ( (T-Tfinal)+Roundoff <= 0.0 ) ...
            || ( Direction < 0 ) && ( (Tfinal-T)+Roundoff <= 0.0 ) )
        
        if ( ISTATUS.Nstp >= OPTIONS.Max_no_steps )
            error('Number of steps exceeds maximum buffer bound \n T= %f;     H= %f', T, H );
        end
        if ( ( ( T+0.1*H ) == T) || ( abs(H) <= Roundoff ) )
            error('Step size too small: T + 10*H = T or H < Roundoff \n T= %f;     H= %f', T, H );
        end
        
        % Limit H if necessary to avoid going beyond Tfinal
        if ( abs(Tfinal-T) < H )
            H = abs(Tfinal-T);
            LMM_state = LMM_struct.stateUpdateH(LMM_state, H);
            NConsecutive = 0;
        end
        
        % Compute the Jacobian at current time
        if (~SkipJac)
            [fjac, ISTATUS] = EvaluateJacobian(T, Y, F, OdeFunction, OPTIONS, ISTATUS);
        end
        
        % Setup adjoint vector products
%         if ( OPTIONS.BiOrthogonalLanczos )
%             if ( ~OPTIONS.MatrixFree )
%                 fjact = @(v)fjac'*v;
%             else
%                 if ( ~isempty(OPTIONS.JacobianAdjointVec) )
%                     fjact = @(v)OPTIONS.JacobianAdjointVec(T,Y,v);
%                 elseif ( nargin(OPTIONS.Jacobian) == 2 )
%                     fjact = @(v)Jac'*v;
%                 else
%                     error('Using biorthogonal Lanczos requires some form of Jacobian Adjoint product.');
%                 end
%             end
%         end
        
        % Repeat step calculation until current step accepted
        acceptStep = false;
        while ( ~acceptStep ) % accepted

%   New solutions

            % Ynew: Y(n+1) evaluated at Order, Ypm1: Y(n+1) eval at Order-1, Ypp1: Y(n+1) eval at Order+1
            [Ynew, YE, ELO, H, StepChanged, RejectFlag, LMM_state, ISTATUS] = LMM_struct.onestep(Order, H, Y, F, fjac, dFdT, T, LMM_state, Coefficients, OdeFunction, OPTIONS, ISTATUS);
            if StepChanged
                NConsecutive = 0;
            end
            
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   Error Estimation
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            % Compute error norm
            SCAL = OPTIONS.AbsTol + OPTIONS.RelTol.*max(abs(Y),abs(Ynew));
            Em1 = max(rms(YE(:,1)./SCAL), 1e-10);
            E   = max(rms(YE(:,2)./SCAL), 1e-10);
            Ep1 = max(rms(YE(:,3)./SCAL), 1e-10);
            
            % unscaled timestep factors
            usFACm1 = (1.0/Em1)^(1.0/ELO(1));
            usFAC   = (1.0/E)^(1.0/ELO(2));
            usFACp1 = (1.0/Ep1)^(1.0/ELO(3));
            
%             disp(['RejectFlag = ', num2str(RejectFlag), '.']);
%             disp(['ELO = [', num2str(ELO), '].']);
%             disp(['Order ', num2str(Order-1), ': E = ', num2str(Em1), ', FAC = ', num2str(usFACm1), '.']);
%             disp(['Order ', num2str(Order), ': E = ', num2str(E), ', FAC = ', num2str(usFAC), '.']);
%             disp(['Order ', num2str(Order+1), ': E = ', num2str(Ep1), ', FAC = ', num2str(usFACp1), '.']);

%             opts = odeset('AbsTol', 100*eps, 'RelTol', 100*eps, 'Jacobian', OPTIONS.Jacobian);
%             [~,Y45] = ode45(OdeFunction, [T T+H], Y, opts);
%             disp(['ODE45: E = ', num2str(max(rms((Ynew - Y45(end,:)')./SCAL), 1e-10))])
    
            
            %keyboard
            
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   Accept/Reject 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            % Check the error magnitude and adjust the step size
            ISTATUS.Nstp = ISTATUS.Nstp + 1;
            if (~RejectFlag) && ( ( E <= 1.0 ) || ( H <= OPTIONS.Hmin ) ) 
                ISTATUS.Nacc = ISTATUS.Nacc + 1;
                
                % Update time
                T = T + Direction*H;          
                
                %keyboard
                
                % Update solution
                Y = Ynew;
                
                NConsecutive = min(NConsecutive+1, OPTIONS.MaxOrder+2);
                
                % Change stepsize/order only after enough consecutive steps
                if NConsecutive >= Order+2
                    % New step size is bounded by FacMin <= Hnew/H <= FacMax
                    FACm1 = min( OPTIONS.FacMax, max( OPTIONS.FacMin, OPTIONS.FacSafeLow*usFACm1 ) );
                    FAC   = min( OPTIONS.FacMax, max( OPTIONS.FacMin, OPTIONS.FacSafe*usFAC ) );
                    FACp1 = min( OPTIONS.FacMax, max( OPTIONS.FacMin, OPTIONS.FacSafeHigh*usFACp1 ) );

                    HNm1 = H*FACm1;
                    HN   = H*FAC;
                    HNp1 = H*FACp1;

                    % Select new order and stepsize
                    % Choose the order that corresponds with the maximal H
                    Hnew = HN;
                    Knew = Order;
                    if ( HNm1 > Hnew && Order ~= 1 )
                        Hnew = HNm1;
                        Knew = Order-1;
                    end
                    if ( HNp1 > Hnew && Order < OPTIONS.MaxOrder && Order+1 < ISTATUS.Nacc )
                        Hnew = HNp1;
                        Knew = Order+1;
                    end
                    if ( RejectLastH ) % No step size or order increase after a rejected step
                        Hnew = min( Hnew, H );
                        Knew = min( Knew, Order );
                    end
                    Hnew = max(OPTIONS.Hmin, min(Hnew, OPTIONS.Hmax)); % apply user-defined limits on Hnew
                else
                    Hnew = H;
                    Knew = Order;
                end
                
                % Update Memory Allocation
                if ( (ISTATUS.Nacc > OPTIONS.ChunkSize*ISTATUS.Nchk) && (OPTIONS.storeCheckpoint == true) )
                    [ Tout, Yout ] = matlOde_appendMemory(NVAR,Tout,Yout,OPTIONS.ChunkSize);
                    ISTATUS.Nchk = ISTATUS.Nchk + 1;
                end
                
                if ( OPTIONS.storeCheckpoint == true )
                    Tout(TYindex,1) = T;
                    Yout(:,TYindex) = Ynew;
                    TYindex = TYindex + 1;
                end                
                
                % Last T and H
                RSTATUS.Nhexit = H;
                RSTATUS.Nhnew = Hnew;
                RSTATUS.Ntexit = T;
                
                H = Hnew;
                Order = Knew;
                
                RejectLastH = false;
                RejectMoreH = false;
                acceptStep = true;
                
                % Compute the function at new time
                F = OdeFunction( T, Y );
                ISTATUS.Nfun = ISTATUS.Nfun + 1;
                
                % Compute the function derivative with respect to T
                if ( ~OPTIONS.Autonomous )
                    [ dFdT, ISTATUS ] = fatOde_FunctionTimeDerivative( T, Roundoff, Y, F, OdeFunction, ISTATUS );
                end
                
                % advance LMM internal state
                LMM_state = LMM_struct.stateAdvance(LMM_state, H, Y, F, Order);
                
                % for debugging
                if ( OPTIONS.displaySteps == true )
                    optstr = [' Order = ', num2str(Order), '; E = ', num2str(E), '.'];
                    str = ['Accepted step.', optstr, ' Time = ', num2str(T), '; Stepsize = ', num2str(H)];
                    disp(str);
                end
                
            else % Reject step
                
                if ( RejectMoreH )
                    Hnew = H*OPTIONS.FacRej;
                    if (Order > 1)
                        Order = Order - 1;
                    end
                else
                    FACm1 = min( OPTIONS.FacMax, max( OPTIONS.FacMin, OPTIONS.FacSafeLow*usFACm1 ) );
                    FAC   = min( OPTIONS.FacMax, max( OPTIONS.FacMin, OPTIONS.FacSafe*usFAC ) );
                    
                    HNm1 = H*FACm1;
                    HN   = H*FAC;
                    
                    Hnew = min(H, max(HN, HNm1));
                    if (Order > 1 && HNm1 >= HN)
                       Order = Order - 1;
                    end
                end
                Hnew = max( OPTIONS.Hmin, min( Hnew, OPTIONS.Hmax ) );
                
                NConsecutive = 0;
                RejectMoreH = RejectLastH;
                RejectLastH = true;
                H = Hnew;
                LMM_state = LMM_struct.stateUpdateH(LMM_state, Hnew);
                
                if ( ISTATUS.Nacc >= 1 ) 
                    ISTATUS.Nrej = ISTATUS.Nrej + 1;
                end
                
                % for debugging
                if ( OPTIONS.displaySteps == true )
                    optstr = [' Order = ', num2str(Order), '; E = ', num2str(E), '.'];
                    str = ['Rejected step.', optstr, ' Time = ', num2str(T), '; Stepsize = ', num2str(H)];
                    disp(str);
                end
                
            end
        end
    end
    
    % Successful return
    Ierr = 1;
    
    % Decallocate Memory
    if ( OPTIONS.storeCheckpoint == true )
%         Tout(TYindex:OPTIONS.Max_no_steps) = [];
%         Yout(:,TYindex:OPTIONS.Max_no_steps) = [];      
        Tout(ISTATUS.Nacc+1:end,:) = [];
        Yout(:,ISTATUS.Nacc+1:end) = [];
    else
        Tout = T;
        Yout = Y;
    end
    Yout = transpose(Yout);

return;
end

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
