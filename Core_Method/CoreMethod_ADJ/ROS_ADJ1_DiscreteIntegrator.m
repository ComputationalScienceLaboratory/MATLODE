function [ Tout, Yout, Lambda, Mu, ISTATUS, RSTATUS, Ierr ] = ROS_ADJ1_DiscreteIntegrator( NVAR, OPTIONS, Coefficient, stack_ptr, adjQuadFlag )

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   Global Variables
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    global ros_Gamma ros_M ros_C ros_Alpha ros_A

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   Local Variables
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % note: based off DLAMCH('E') NOT MATLAB's EPS
    Roundoff = 1.11022302462515654E-016;
    
    DeltaMin = 1e-5;
    
    Mu = OPTIONS.Mu();    
    Lambda = OPTIONS.Lambda;
    
    ISTATUS = ISTATUS_Struct('default');
    RSTATUS = RSTATUS_Struct('default');    
    
    U = zeros(NVAR*Coefficient.NStage,OPTIONS.NADJ);
    V = zeros(NVAR*Coefficient.NStage,OPTIONS.NADJ);    

    RP0 = zeros(OPTIONS.NP,OPTIONS.NADJ);
    RP1 = zeros(OPTIONS.NP,OPTIONS.NADJ);
    
    RY0 = zeros(NVAR,OPTIONS.NADJ);
    RY1 = zeros(NVAR,OPTIONS.NADJ);
    
    FPJAC0 = zeros(NVAR,OPTIONS.NP);
    FPJAC1 = zeros(NVAR,OPTIONS.NP);

    ros_W = zeros(Coefficient.NStage,1);
    
    TYindex = 1;
    
    Direction = 1; % <----- temporary fix
%     if ( Tend >= Tstart )
%         Direction = 1;
%     else
%         Direction = -1;
%     end

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
%   Time loop
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    while ( stack_ptr > 0 )
        % Recover checkpoints for stage values an vectors
        [ T, H, Ystage, K, stack_ptr ] = rosPop(NVAR,Coefficient.NStage,stack_ptr);
        Yout(:,TYindex) = Ystage;
        Tout(TYindex,1) = T;        
        TYindex = TYindex + 1;
        
        ISTATUS.Nstp = ISTATUS.Nstp + 1;
        
        % Compute LU decomposition
        if ( ~OPTIONS.SaveLU )
            fjac = OPTIONS.Jacobian(T,Ystage(:));
            ISTATUS.Njac = ISTATUS.Njac + 1;
            
            Tau = 1/(Direction*H*ros_Gamma(1));
            ISTATUS.Njac = ISTATUS.Njac + 1; %<-------------------????
            
            [ ISING, e ] = lss_decomp( NVAR, Tau, fjac );
            ISTATUS.Ndec = ISTATUS.Ndec + 1;
            
        end
        
        if ( adjQuadFlag )
            RY0 = OPTIONS.DRDY(T,Ystage(:));
        end

       Tmp = 0;
       Tmp3 = 0;
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
                TMP = OPTIONS.Hesstr_vec(T,Ystage,U(istart:istart+(NVAR-1),m),K(istart:istart+(NVAR-1)));
                Lambda(:,m) = Lambda(:,m) + TMP;
                if ( adjQuadFlag )
                    TMP = OPTIONS.Hesstr_vec_r(T,Ystage,H*ros_W(istage),K(istart:istart+(NVAR-1)));
                    Lambda(:,m) = Lambda(:,m) + TMP;
                end
            end
        end        
        
        % Compute Mu
        Beta = 1;
        for istage=1:Coefficient.NStage
            istart = NVAR*(istage-1) + 1;
            Tau = T + ros_Alpha(istage)*Direction*H;
            FPJAC1 = OPTIONS.Jacp(Tau,Ystage(istart:istart+(NVAR-1)));
            if ( adjQuadFlag )
                RP1 = OPTIONS.DRDP(Tau,Ystage(istart:istart+(NVAR-1)));
            end
            for m=1:OPTIONS.NADJ
                Tmp3 = zeros(OPTIONS.NP,1);
                Tmp3 = Tmp3 + transpose(FPJAC1)*U(istart:istart+(NVAR-1),m);
                Mu(:,m) = Mu(:,m) + Tmp3;
                if ( adjQuadFlag )
                    Tmp3 = Tmp3 + H*ros_W(istage)*RP1(:,m);
                    Mu(:,m) = Mu(:,m) + Tmp3;
                end
                Tmp3 = OPTIONS.Hesstr_vec_f_py(T,Ystage,U(istart:istart+(NVAR-1),m),K(istart:istart+NVAR-1));
                Mu(:,m) = Mu(:,m) + Tmp3;
                if ( adjQuadFlag )
                    Tmp3 = OPTIONS.Hesstr_vec_r_py(T,Ystage,H*ros_W(istage),K(istart:istart+NVAR-1));
                    Mu(:,m) = Mu(:,m) + Tmp3;
                end
            end
        end
        
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   NonAutonomous term
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if ( ~OPTIONS.Autonomous )
            Delta = sqrt(Roundoff)*max(DeltaMin,abs(T));
            djdt = OPTIONS.Jacobian(T+Delta,Ystage);
            ISTATUS.Njac = ISTATUS.Njac + 1;
            
            % Jacobian time derivative
            djdt = djdt - fjac;
            djdt = djdt/Delta;
            
            FPJAC0 = OPTIONS.Jacp(T,Ystage);
            FPJAC1 = OPTIONS.Jacp(T+Delta,Ystage);
            FPJAC1 = FPJAC1 - FPJAC0;
            FPJAC1 = FPJAC1/Delta;
            
            for m=1:OPTIONS.NADJ
                Tmp = zeros(NVAR,1);
                for istage=1:Coefficient.NStage
                    istart = NVAR*(istage-1)+1;
                    Tmp = Tmp + ros_Gamma(istage)*U(istart:istart+(NVAR-1),m);
                end
                Alpha = H;
                Mu(:,m) = Mu(:,m) + H*transpose(FPJAC1)*Tmp;
                Tmp2 = transpose(djdt)*Tmp;
                Lambda(:,m) = Lambda(:,m) + H*Tmp2;                                
            end
        end
        
        % QuadAutonomous <----false (temp fix... need to add to ICNTRL )
        if ( adjQuadFlag && true )
            % Approximate d(r_y)/dt and add it to Lambda
            RY1 = OPTIONS.DRDY(T+Delta,Ystage(1)); % <----------------confirm code calls for only first index
            RY1 = RY1 - RY0;
            RY1 = RY1/Delta;
            Lambda = Lambda + H*GammaW*RY1;
            
            % Approximate d(r_p)/dr and add it to Mu
            RP0 = OPTIONS.DRDP(T,Ystage(1));
            RP1 = OPTIONS.DRDP(T+Delta,Ystage(1));
            RP1 = RP1 - RP0;
            RP1 = RP1/Delta;
            Mu = Mu + H*GammaW*RP1;
        end
        
    end % time 
    
    % Successful exit
    Ierr = 1;

return;

