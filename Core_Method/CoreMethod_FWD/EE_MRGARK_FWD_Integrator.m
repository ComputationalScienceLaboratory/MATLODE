%% ERK_FWD_Integrator
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
%    [ Tout, Yout, ISTATUS, RSTATUS, Ierr, stack_ptr, quadrature ] = ERK_FWD_Integrator( OdeFunction, Tspan, Y, OPTIONS, Coefficient, adjStackFlag, adjQuadFlag, stack_ptr )
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
% |adjStackFlag|: Adjoint snapshot stack flag
%
% |adjQuadFlag|: Adjoint quadrature stack flag
%
% |stack_ptr|: pointer for global snapshot stack
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
% |stack_ptr|: Pointer for global snapshot stack
%
% |quadrature|: Adjoint quadrature
%
%% Description
% Explicit Runge-Kutta forward core method.
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
function [ Tout, Yout, ISTATUS, RSTATUS, Ierr] = MRGARK_ERK_FWD_Integrator( OdeFunctionOne,...
    OdeFunctionTwo, Tspan, Y, OPTIONS, Coefficient, adjStackFlag, adj_QuadFlag, stack_ptr)

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Global Variables
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% global rkA rkB rkC rkE rKAFS rkASF
% Constant computation
[A, BT, BThat, methodOrder, embeddedOrder, Afs, Asf] = coupling(OPTIONS.Method);
AT = A.';
order = min(methodOrder, embeddedOrder);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Local Variables
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Get Problem Size
NVAR = OPTIONS.NVAR;

Tinitial = Tspan(1);
Tfinal = Tspan(end);

Roundoff = max(eps(Tspan))/2;

% Allocate Memory
[ Tout, Yout ] = matlOde_allocateMemory( NVAR, OPTIONS.ChunkSize );

TYindex = 1;

%Force initial value matrix to be N X 1.
if ( size(Y,2) == 1 )
    % DO NOTHING
else
    Y = transpose(Y);
end  

Yout(:, TYindex) = Y;
Tout(TYindex) = Tinitial;

% if ( adjQuadFlag == true ) % initial value of the quadrature
%     quadrature = OPTIONS.Quadrature(Tinitial,Y);
% else
%     quadrature = NaN;
% end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Initial settings
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ISTATUS = ISTATUS_Struct('default');
RSTATUS = RSTATUS_Struct('default');

ISTATUS.Nchk = 1;

T = Tinitial;
Tdirection = sign( Tfinal-Tinitial );
H = max( abs(OPTIONS.Hmin), abs(OPTIONS.Hstart) );
if ( abs(H) <= 10.0*Roundoff)
    H = 1.0e-6;
end

H = min( abs(H), OPTIONS.Hmax );
sign_Tdirection = sign(Tdirection );
H = H*sign_Tdirection;

%Micro Step
M = OPTIONS.InitialM;
h = H/M;

%radius for change in M
Mradius = OPTIONS.RadiusM;

NumStages = numel(BT);
NumVars = numel(Y);

%Building matrices for both fast and slow
KFast = zeros(NumVars, NumStages);
KSlow = zeros(NumVars, NumStages);

Reject = false;
FirstStep = true;

%TEMPORARY VALUES I WILL NEED TO FIX
% scaling Factors when computing next stepsize
facMax = 4;
facMin = 1/4;
fac = 0.9;
counter = 1;
WeightedFastSum = 0;
WeightedFastCount = 0;
WeightedSlowSum = 0;
WeightedSlowCount = 0;
alpha = 0.25;

% Determine scaling factor for integration
% SCAL = 1.0 ./ ( OPTIONS.AbsTol + OPTIONS.RelTol .* abs(Y) );

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Time loop begins
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
while ( T < Tfinal)
    if ( ISTATUS.Nstp > OPTIONS.Max_no_steps )
        error('Number of steps exceeds maximum bound.');
    end
    
    if H == 0
       fprintf('H IS 0');
    end
    
%     case 1: determines if we are close enough to tfinal due to
%             roundoff error
%     case 2: cannot distinguish between tnext and tcurrent
    if abs( Tfinal-T) <= 10*Roundoff*abs(Tfinal)
        T = Tfinal;
        break;
    elseif ( (T+0.1*H == T) || (abs(H) <= Roundoff) )
        fprintf('%d\n', H);
        error('Step size too small; T + 10*H = T or H < eps/2.');
    end
    
    KSlowFast = zeros(NumVars, NumStages);
    
    %Error control vectors
    YLambda = Yout(:, TYindex);
    YHatLambda = Yout(:, TYindex);
    
    %Counter for next slow stage to solve
    ColumnCounter = 1;
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % Stages
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    for lambda = 1 : M
       
        %Extract the fast-slow and slow-fast for current lambda and M
        Afslambda = Afs(lambda, M);
        Asflambda = Asf(lambda, M);
            
        for stage = 1 : NumStages
            
            for col = ColumnCounter : NumStages
                
                % Current fast stage depends on an unsolved slow stage
                if (Afslambda(stage, col) ~= 0)
                    
                    %Time used for error control
                    if lambda == 1
                        TimeSlowInit = tic;
                    end
                    
                    %The time value in which this slow stage will be solved
                    TimeValSlow = T + H * sum(AT(:, col));
                    
                    %Evaluate slow stage
                    KSlow(:,col) = OdeFunctionTwo(TimeValSlow, Yout(:, TYindex)...
                        + H * KSlow(:, 1:(col - 1)) * AT(1 : col -1, col)...
                        + h * KSlowFast(:, col));
                    
                    if lambda == 1
                        TimeSlowEnd = toc(TimeSlowInit);
                    end
                    
                    ColumnCounter = ColumnCounter + 1;
                end
                
            end
            
            %Time for fast stage
            if lambda == 1
                TimeFastInit = tic;
            end
            
            %Time value in which this fast stage will be solved
            TimeValFast = T + h * sum(AT(stage, :)) + ((lambda - 1)/M) * H;
            
            %Solve the current fast stage
            KFast(:, stage) = OdeFunctionOne(TimeValFast, YLambda ...
                + h * KFast(:, 1 : (stage - 1)) * AT(1 : stage - 1, stage)...
                + H * KSlow(:, 1 : (ColumnCounter - 1))...
                * Afslambda(stage, 1 : (ColumnCounter - 1)).');
            
            if (lambda == 1)
                TimeFastEnd = toc(TimeFastInit);
            end
            
            % New Fast stage solved, update the summation matrix to be used
            % for slow stage
            KSlowFast(:,ColumnCounter : end) = KSlowFast(:, ColumnCounter : end)...
                + KFast(:,stage)  * Asflambda(ColumnCounter : end, stage).';
        end
        
        %Update error control vectors
        YLambda = YLambda + h * KFast * BT;
        YHatLambda = YHatLambda + h * KFast * BThat;
    end
    
    % Solve all remaining slow stages
    for col = ColumnCounter : NumStages
        TimeValSlow = T + H * sum(AT(:, col));
        KSlow(:,col) = OdeFunctionTwo(TimeValSlow, Yout(:, TYindex)...
            + H * KSlow(:, 1:(col - 1)) * AT(1 : col -1, col)...
            + h * KSlowFast(:, col));
    end
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % Error estimation
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    YTemp = YLambda + H * KSlow * BT;
    YHat = YHatLambda + H * KSlow * BThat;
    YHatSlow = YLambda + H * KSlow * BThat;
    YHatFast = YHatLambda + H * KSlow * BT;
    sc = OPTIONS.AbsTol + max(abs(YTemp), abs(Yout(:, TYindex))) .* OPTIONS.RelTol;
    
    errSlow = rms((YTemp - YHatSlow) ./ sc);
    errFast = rms((YTemp - YHatFast) ./ sc);
    
    % Arash and Steven's Super Cool POS
    Err = rms((YTemp - YHat) ./ sc);
    if M <= Mradius
        Mnews(1 : (Mradius - M + 1)) = 1 ./ (2 : (Mradius - M + 2));
        Mnews((Mradius - M + 2) : 2*Mradius+1) = (1 : Mradius + M);
    else
        Mnews = max(1, M - Mradius):(M + Mradius);
    end
    % Finding the M value that minimizes the error
    if (~isempty(OPTIONS.TimeValSlow))
        TimeSlowEnd = OPTIONS.TimeValSlow;
    end
    
    if ( ~isempty(OPTIONS.TimeValFast) )
        TimeFastEnd = OPTIONS.TimeValFast;
    end

    WeightedFastSum = TimeFastEnd + (1 - alpha)*WeightedFastSum;
    WeightedSlowSum = TimeSlowEnd + (1 - alpha)*WeightedSlowSum;
    WeightedFastCount = 1 + (1 - alpha)*WeightedFastCount;
    WeightedSlowCount = 1 + (1 - alpha) * WeightedSlowCount;
    TimeSlowAverage = WeightedSlowSum/WeightedSlowCount;
    TimeFastAverage = WeightedFastSum/WeightedFastCount;

%     fprintf('%.5e\n', TimeSlowAverage/TimeFastAverage);
    [~, index] = min ((TimeSlowAverage + Mnews * TimeFastAverage)...
        .* (errSlow + errFast * (M ./ Mnews ).^order)...
        .^(1/(order+1)));
    Mnews = Mnews(index);
    
    %Switch fast and slow functions if number of Micro-steps is below 1
    if Mnews <= 0.5
        counter = counter + 1;
        fprintf('FLIP-FLOP %d\n', counter);
        temp = OdeFunctionOne;
        OdeFunctionOne = OdeFunctionTwo;
        OdeFunctionTwo = temp;
        Mnews = round(1 / Mnews);
    end
    
    Mnews = round(Mnews);
    
    Hfact = (errSlow + errFast * (M / Mnews)) ^ (-1 / (order + 1));
    Hnew =  H * min(facMax, max(facMin, fac * Hfact));
    M = Mnews;
    disp(M);
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % Error estimation
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ISTATUS.Nstp = ISTATUS.Nstp + 1;
%     TMP = zeros(NVAR,1);
%     for i = 1:Coefficient.NStage
%         if (rkE(i) ~= 0.0 )
%             TMP = TMP + H*rkE(i)*K(:,i);
%         end
%     end
%     
%     % Perform error norm
%     Err = sum((TMP.*SCAL).^2,1);
%     Err = max( sqrt( Err/double(NVAR) ), 1.0d-10 );
%     
%     % Computation of new step size Hnew
%     Fac = OPTIONS.FacSafe*(Err)^(-1.0/Coefficient.ELO);
%     Fac = max( OPTIONS.FacMin, min(OPTIONS.FacMax, Fac) );
%     Hnew = H*Fac;
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % Accept/Reject step
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if ( Err <= 1.0 )
        %fprintf('STEP ACCEPTED\n');
        FirstStep = false;
        ISTATUS.Nacc = ISTATUS.Nacc + 1;
        
        % Checkpoints for Adjoint Calculations
        if ( adjStackFlag == true )
            stack_ptr = ERK_Push( NVAR, Coefficient.NStage, T, H, Y, Z, OPTIONS.Max_no_steps, stack_ptr );
        end
        
%         Update Memory Allocation
        if ( (ISTATUS.Nacc > OPTIONS.ChunkSize*ISTATUS.Nchk) && (OPTIONS.storeCheckpoint == true) )
            [ Tout, Yout ] = matlOde_appendMemory(NVAR,Tout,Yout,OPTIONS.ChunkSize);
            ISTATUS.Nchk = ISTATUS.Nchk + 1;
        end
        
        % Update time
        T = T + H;

        
        % Update solution
%         for i = 1:Coefficient.NStage
%             if ( rkB(i) ~= 0.0 )
%                 Y = Y + H*rkB(i)*K(:,i);
%             end
%         end
        
%        Store checkpoint values
%         if ( OPTIONS.storeCheckpoint == true )
% 
%         end
        TYindex = TYindex + 1;
        Tout(TYindex) = T;
        Yout(:, TYindex) = YTemp;
        
        % Update scaling coefficients
        %SCAL = 1.0 ./ ( OPTIONS.AbsTol + OPTIONS.RelTol .* abs(Y) );
        
        % Next time step
        Hnew = Tdirection*min( abs(Hnew), OPTIONS.Hmax);
        
%         Last T and H
        RSTATUS.Ntexit = T;
        RSTATUS.Nhexit = H;
        RSTATUS.Nhnew = Hnew;
        
%         No step increase after a rejection
        if ( Reject )
            Hnew = Tdirection*min( abs(Hnew), abs(H) );
        end
        Reject = false;
        
        
        % Update step size
        if ( (T+Hnew-Tfinal)*Tdirection > 0.0)
            H = Tfinal - T;
        else
            H = max(OPTIONS.Hmin, min(OPTIONS.Hmax, Hnew));
        end
        
        % for debugging
%         if ( OPTIONS.displaySteps == true )
%             str = ['Accepted step. Time = ', num2str(T), '; Stepsize = ', num2str(H)];
%             disp(str);
%         end

        
    else % Step is rejected
        % Update step size
        if ( FirstStep || Reject )
            H = OPTIONS.FacRej*H;
        else
            H = max(OPTIONS.Hmin, min(OPTIONS.Hmax, Hnew));
        end
        
        Reject = true;
%         fprintf('STEP REJECTED\n');

%        Update number of rejections
        if ( ISTATUS.Nacc >= 1.0 )
            ISTATUS.Nrej = ISTATUS.Nrej + 1;
        end
%         
%         %for debugging
%         if ( OPTIONS.displaySteps == true )
%             str = ['Rejected step. Time = ', num2str(T), '; Stepsize = ', num2str(H) ];
%             disp(str);
%         end
        
    end % END accept
    h = H/M;
end % Time loop

% Successful return
Ierr = 1;

% %Deallocate Memory
% if ( OPTIONS.storeCheckpoint == true )
[ Tout, Yout ] = matlOde_deallocateMemory(Tout,Yout,ISTATUS.Nacc);
Tout(ISTATUS.Nacc+1:end,:) = [];
Yout(:,ISTATUS.Nacc+1:end) = [];
% else
%     Tout = T;
%     Yout = Y;
% end

Tout = transpose(Tout);
Yout = transpose(Yout);

return;



% (START) SUBROUTINE: Explcit Runge Kutta Push
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   Saves the next trajectory snapshot for discrete adjoints
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function stack_ptr = ERK_Push( NVAR, rkS, T, H, Y, Z, Max_no_steps, stack_ptr )
global chk_H chk_T chk_Y chk_Z

stack_ptr = stack_ptr + 1;
if ( stack_ptr > Max_no_steps )
    error( 'Push failed: buffer overflow' );
end
chk_H( stack_ptr ) = H;
chk_T( stack_ptr ) = T;
chk_Y( 1:NVAR, stack_ptr ) = Y(1:NVAR);
chk_Z( 1:NVAR, 1:rkS, stack_ptr ) = Z(1:NVAR, 1:rkS);
%chk_Z( 1:rkS, 1:NVAR, stack_ptr ) = Z;

return; % (END) SUBROUTINE: Explicit Runge Kutta Push

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
