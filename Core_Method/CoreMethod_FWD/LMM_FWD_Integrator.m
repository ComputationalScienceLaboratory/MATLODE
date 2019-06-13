%% ROK_FWD_Integrator
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
        Tspan, Y0, OPTIONS, LMM_struct, adjStackFlag, adjQuadFlag )
    
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
%         Yout = zeros(NVAR,OPTIONS.Max_no_steps);
%         Tout = zeros(OPTIONS.Max_no_steps,1);
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
        RSTATUS.S1residual = zeros(OPTIONS.ChunkSize, 1);
        RSTATUS.stepsizes = zeros(OPTIONS.ChunkSize, 1);
        if ( OPTIONS.IOArnoldi )
            RSTATUS.blockSizes = zeros(OPTIONS.ChunkSize, 1);
        end
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
    
    dFdT = zeros(NVAR,1);
    
    % Preallocate integrator internal state
    Coefficients = LMM_struct.coefficients();
    k = OPTIONS.MaxOrder + 1;
    Y = [zeros(NVAR, k-1), Y0];
    F = zeros(NVAR, k);
    H = [zeros(1,k-1), H];
    IDX = 1:k; % Y(IDX(k)) is always the current Y, and accepted Ynew are stored in Y(IDX(1)).
    Order = 1;
    
    RejectLastH = false;
    RejectMoreH = false;
        
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   Time loop
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    while ( ( Direction > 0 ) && ( (T-Tfinal)+Roundoff <= 0.0 ) ...
            || ( Direction < 0 ) && ( (Tfinal-T)+Roundoff <= 0.0 ) )
        
        if ( ISTATUS.Nstp >= OPTIONS.Max_no_steps )
            error('Number of steps exceeds maximum buffer bound \n T= %f;     H= %f', T, H(IDX(k)) );
        end
        if ( ( ( T+0.1*H(IDX(k)) ) == T) || ( abs(H(IDX(k))) <= Roundoff ) )
            error('Step size too small: T + 10*H = T or H < Roundoff \n T= %f;     H= %f', T, H(IDX(k)) );
        end
        
        % Limit H if necessary to avoid going beyond Tfinal
        H(IDX(k)) = min( H(IDX(k)), abs( Tfinal-T ) );
        
        % Compute the function at current time
        F(:,IDX(k)) = OdeFunction( T, Y(:,IDX(k)) );
        ISTATUS.Nfun = ISTATUS.Nfun + 1;
        
        % Compute the function derivative with respect to T
        if ( ~OPTIONS.Autonomous )
            [ dFdT, ISTATUS ] = fatOde_FunctionTimeDerivative( T, Roundoff, Y(:,IDX(k)), F(:,IDX(k)), OdeFunction, ISTATUS );
        end
        
        % Compute the Jacobian at current time
        [fjac, ISTATUS] = EvaluateJacobian(T, Y(:,IDX(k)), F(:,IDX(k)), OPTIONS, ISTATUS);
        
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
            [Ynew, YE, ELO, H, RejectFlag, ISTATUS] = LMM_struct.onestep(Order, H, Y, F, k, IDX, fjac, dFdT, T, Coefficients, OdeFunction, OPTIONS, ISTATUS);
            
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   Error Estimation
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            % Compute error norm
            SCAL = OPTIONS.AbsTol + OPTIONS.RelTol.*max(abs(Y(:,IDX(k))),abs(Ynew));
%             Em1 = max(sqrt(sum((YE(:,1)./SCAL).^2)/NVAR), 1e-10);
%             E   = max(sqrt(sum((YE(:,2)./SCAL).^2)/NVAR), 1e-10);
%             Ep1 = max(sqrt(sum((YE(:,3)./SCAL).^2)/NVAR), 1e-10);
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
%             
%             opts = odeset('AbsTol', 100*eps, 'RelTol', 100*eps, 'Jacobian', OPTIONS.Jacobian);
%             [~,Y45] = ode45(OdeFunction, [T T+H(IDX(k))], Y(:,IDX(k)), opts);
%     
%             disp(['ODE45: E = ', num2str(max(rms((Ynew - Y45(end,:)')./SCAL), 1e-10))])
    
            
            %keyboard
            
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   Accept/Reject 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            % Check the error magnitude and adjust the step size
            ISTATUS.Nstp = ISTATUS.Nstp + 1;
            if (~RejectFlag) && ( ( E <= 1.0 ) || ( H(IDX(k)) <= OPTIONS.Hmin ) ) 
                ISTATUS.Nacc = ISTATUS.Nacc + 1;
                
                % Update time
                T = T + Direction*H(IDX(k));          
                
                %keyboard
                
                % Update solution
                Y(:,IDX(1)) = Ynew;
                
                % New step size is bounded by FacMin <= Hnew/H <= FacMax
                FACm1 = min( OPTIONS.FacMax, max( OPTIONS.FacMin, OPTIONS.FacSafeLow*usFACm1 ) );
                FAC   = min( OPTIONS.FacMax, max( OPTIONS.FacMin, OPTIONS.FacSafe*usFAC ) );
                FACp1 = min( OPTIONS.FacMax, max( OPTIONS.FacMin, OPTIONS.FacSafeHigh*usFACp1 ) );
                
                HNm1 = H(IDX(k))*FACm1;
                HN   = H(IDX(k))*FAC;
                HNp1 = H(IDX(k))*FACp1;
            
                % Select new order and stepsize
                % Choose the order that corresponds with the maximal H
                Hnew = HN;
                Knew = Order;
                if ( HNm1 > Hnew && Order ~= 1 )
                    Hnew = HNm1;
                    Knew = Order-1;
                end
                if ( HNp1 > Hnew && Order ~= OPTIONS.MaxOrder && Order+1 < nnz(H) )
                    Hnew = HNp1;
                    Knew = Order+1;
                end
                if ( RejectLastH ) % No step size or order increase after a rejected step
                    Hnew = min( Hnew, H(IDX(k)) );
                    Knew = min( Knew, Order );
                end
                Hnew = max(OPTIONS.Hmin, min(Hnew, OPTIONS.Hmax)); % apply user-defined limits on Hnew
                H(IDX(1)) = Hnew;
                Order = Knew;
                
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
                RSTATUS.Nhexit = H(IDX(k));
                RSTATUS.Nhnew = Hnew;
                RSTATUS.Ntexit = T;
                
                RejectLastH = false;
                RejectMoreH = false;
                acceptStep = true;
                
                % cycle IDX for circular buffer
                tmp = IDX(1);
                for j = 1:k-1
                    IDX(j) = IDX(j+1);
                end
                IDX(k) = tmp;
                
                % for debugging
                if ( OPTIONS.displaySteps == true )
                    optstr = [' Order = ', num2str(Order), '; E = ', num2str(E), '.'];
                    str = ['Accepted step.', optstr, ' Time = ', num2str(T), '; Stepsize = ', num2str(H(IDX(k)))];
                    disp(str);
                end
                
            else % Reject step
                
                if ( RejectMoreH )
                    Hnew = H(IDX(k))*OPTIONS.FacRej;
                    if (Order > 1)
                        Order = Order - 1;
                    end
                else
                    FACm1 = min( OPTIONS.FacMax, max( OPTIONS.FacMin, 0.9*OPTIONS.FacSafeLow*usFACm1 ) );
                    FAC   = min( OPTIONS.FacMax, max( OPTIONS.FacMin, 0.9*OPTIONS.FacSafe*usFAC ) );
                    
                    HNm1 = H(IDX(k))*FACm1;
                    HN   = H(IDX(k))*FAC;
                    
                    Hnew = min(H(IDX(k)), max(HN, HNm1));
                    if (Order > 1 && HNm1 >= HN)
                       Order = Order - 1;
                    end
                end
                Hnew = max( OPTIONS.Hmin, min( Hnew, OPTIONS.Hmax ) );
                
                RejectMoreH = RejectLastH;
                RejectLastH = true;
                H(IDX(k)) = Hnew;
                
                if ( ISTATUS.Nacc >= 1 ) 
                    ISTATUS.Nrej = ISTATUS.Nrej + 1;
                end
                
                % for debugging
                if ( OPTIONS.displaySteps == true )
                    optstr = [' Order = ', num2str(Order), '; E = ', num2str(E), '.'];
                    str = ['Rejected step.', optstr, ' Time = ', num2str(T), '; Stepsize = ', num2str(H(IDX(k)))];
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
        Yout = Y(:,IDX(k));
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
