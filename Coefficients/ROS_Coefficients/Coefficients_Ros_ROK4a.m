function [ rosMethod, rosELO, rosS, rosName ] = Coefficients_Ros_ROK4a( ROK4A )

    global ros_Alpha ros_Gamma
    global ros_A ros_C ros_M ros_E
    global ros_NewF

    % Rosenbrock method number
    rosMethod = ROK4A;

    % Number of stages
    rosS = 4;
    
    % rosELO  = estimator of local order - the minimum between the
    %           main and the embedded scheme orders plus 1
    rosELO = 4.0d0;

    % Method name
    rosName = 'ROK4a';

    
    
    gam = 0.572816062482134855408001384976768340931514124329888624090;
    
    gamma = diag(repmat(gam,1,4));
    gamma(2,1) = -1.91153192976055097824558940682133229601582735199337138313;
    gamma(3,1) = 0.3288182406115352215636148889409996289068084425348308862;
    gamma(3,2) = 0;
    gamma(4,1) = 0.0330364423979581129062589491015367687146658980867715225;
    gamma(4,2) = -0.2437515237610823531197130316973934919183493340255052118;
    gamma(4,3) = -0.1706260299199402983371506431639627141896309125577564011;
    
    c = diag(1.0 ./ diag(gamma)) - inv(gamma);
    
    ros_C(1) = c(2,1);
    ros_C(2) = c(3,1);
    ros_C(3) = c(3,2);
    ros_C(4) = c(4,1);
    ros_C(5) = c(4,2);
    ros_C(6) = c(4,3);
    
    ros_Gamma(1) = sum(gamma(1,:));
    ros_Gamma(2) = sum(gamma(2,:));
    ros_Gamma(3) = sum(gamma(3,:));
    ros_Gamma(4) = sum(gamma(4,:));
    
    alpha = zeros(4,4);
    alpha(2,1) = 1;
    alpha(3,1) = 0.1084530016931939175868117432550153454393009345950250308;
    alpha(3,2) = 0.3915469983068060824131882567449846545606990654049749692;
    alpha(4,1) = 0.4345304775600447762471370176270279643577768155816826992;
    alpha(4,2) = 0.1448434925200149254157123392090093214525922718605608997;
    alpha(4,3) = -0.0793739700800597016628493568360372858103690874422435989;
    
    a = alpha * inv(gamma);
    
    ros_A(1) = a(2,1);
    ros_A(2) = a(3,1);
    ros_A(3) = a(3,2);
    ros_A(4) = a(4,1);
    ros_A(5) = a(4,2);
    ros_A(6) = a(4,3);
    
    ros_Alpha(1) = 0;
    ros_Alpha(2) = sum(alpha(2,:));
    ros_Alpha(3) = sum(alpha(3,:));
    ros_Alpha(4) = sum(alpha(4,:));
        
    b = zeros(4,1);
    b(1) = 1/6;
    b(2) = 1/6;
    b(3) = 0;
    b(4) = 2/3;
    
    ros_M = b' * inv(gamma);
    
    bhat = zeros(4,1);
    bhat(1) = 0.5026932257368423534541250307675423348147218619695336514 - b(1);
    bhat(2) = 0.2786755196900585622624861213669585560493517317676223283 - b(2);
    bhat(3) = 0.2186312545730990842833888478654991091359264062628440203 - b(3);
    bhat(4) = 0 - b(4);
    
    ros_E = bhat' * inv(gamma);
    
    ros_NewF = zeros(4,1);
    ros_NewF(1) = true;
    ros_NewF(2) = true;
    ros_NewF(3) = true;
    ros_NewF(4) = true;
    
return;

