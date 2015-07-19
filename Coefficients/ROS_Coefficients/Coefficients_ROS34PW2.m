function [ rosMethod rosELO rosS rosName ] = Coefficients_ROS34PW2( ROW1 )
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Filename: Coefficients_ROS34PW2.m
%
% Original Author:
%
% File Create Date:
%
% Input Arguments:
%   Name        Type
%   ROW1        integer
%
% Output Arguments:
%   Name        Type
%   rosMethod   integer
%   rosELO      double
%   rosS        integer
%
% Modification History:
%   Date        Developer         Email             Action  
%   7/17/2012   Tony D'Augustine  adaug13@vt.edu    Rosenbrock Coefficients Template Created 
%   7/17/2015   Arash Sarshar	  sarshar@vt.edu    Added ROW method based on RANG paper
%  
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Refer to : 
% Rang, Joachim, and L. Angermann.
%"New Rosenbrock W-methods of order 3 for partial differential algebraic equations of index 1.
%" BIT Numerical Mathematics 45.4 (2005): 761-787.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Coefficients_Ros34PW2:
%   Initializes coefficients for ROS34PW2 Rosenbrock-w method.
% 
%       Stages  Order   Stability Property
%       4       3       L-stable, Stiffly accurate
%
% Coefficients_Ros34PW2: INPUT ARGUMENTS
%   ROW (integer):
%
% Coefficients_Ros34PW2: OUTPUT ARGUMENTS
%   rosMethod (integer): Rosenbrock method number.
%   rosELO (double):
%   rosS (integer):
%
% Coefficients_Ros34PW2: GLOBAL VARIABLES
%   ros_A (double):
%   ros_C (double):
%   ros_M (double):
%   ros_E (double):
%   ros_Alpha (double):
%   ros_Gamma (double):
%   ros_NewF (double):
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    global ros_Alpha ros_Gamma
    global ros_A ros_C ros_M ros_E
    global ros_NewF

    rosMethod = ROW1;
	% Method name
    rosName = 'ROS34PW2';
    % Number of stages
    rosS = 4;
    
    % Gamma Matrix
    Gamma_ii=4.3586652150845900e-01;
    Gamma=Gamma_ii*diag(ones(1,rosS));
    Gamma(2,1)=-8.7173304301691801e-01;
    Gamma(3,1)=-9.0338057013044082e-01;
    Gamma(3,2)= 5.4180672388095326e-02;
    Gamma(4,1)= 2.4212380706095346e-01;
    Gamma(4,2)=-1.2232505839045147e+00;
    Gamma(4,3)= 5.4526025533510214e-01;
    
    %%%% Alpha Matrix
    Alpha=zeros(rosS);
    Alpha(2,1)= 8.7173304301691801e-01;
    Alpha(3,1)= 8.4457060015369423e-01;
    Alpha(3,2)=-1.1299064236484185e-01;
    Alpha(4,1)= 0.0000000000000000e+00;
    Alpha(4,2)= 0.0000000000000000e+00;
    Alpha(4,3)= 1.0000000000000000e+00;
    
    %%%% Coefficients for new step solution bi
    b(1) = 2.4212380706095346d-01;
    b(2) =-1.2232505839045147d+00;
    b(3) = 1.5452602553351020d+00;
    b(4) = 4.3586652150845900d-01;
    
    %%%% Coefficients for error estimator \hat b_i
    b_hat(1) = 3.7810903145819369d-01;
    b_hat(2) =-9.6042292212423178d-02;
    b_hat(3) = 5.0000000000000000d-01;
    b_hat(4) = 2.1793326075422950d-01;
    b_hat=b_hat';
    b=b';
    
    %%%% Convert to Implementation Coeefiecnts 
    %%%% Refer to : Hairer, Wanner, BookII, CH IV.7
    A=Alpha/Gamma;
    C=diag(diag(inv(Gamma)))-inv(Gamma);
    M=transpose(inv(Gamma))*b;
    E=transpose(inv(Gamma))*b_hat;
    
%~~~> The coefficient matrices A and C are strictly lower triangular.
%   The lower triangular (subdiagonal) elements are stored in row-wise order:
%   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
%   The general mapping formula is:
%       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )    
%       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )  

    ros_A(1) = A(2,1);
    ros_A(2) = A(3,1);
    ros_A(3) = A(3,2);
    ros_A(4) = A(4,1);
    ros_A(5) = A(4,2);
    ros_A(6) = A(4,3);

    ros_C(1) = C(2,1);
    ros_C(2) = C(3,1);
    ros_C(3) = C(3,2);
    ros_C(4) = C(4,1);
    ros_C(5) = C(4,2);
    ros_C(6) = C(4,3);
    
      %~~~> Y_stage_i ~ Y( T + H*Alpha_i ) 
      %%% Note: non-autonomous ODE was NOT tested ! 
    ros_Alpha(1) = 0.0d0;
    ros_Alpha(2) = A(1);
    ros_Alpha(3) = A(2)+A(3);
    ros_Alpha(4) = A(4)+A(5)+A(6);
    
    %~~~> Gamma_i = \sum_j  gamma_{i,j}
    ros_Gamma(1) = sum(Gamma(1,:));
    ros_Gamma(2) = sum(Gamma(2,:));
    ros_Gamma(3) = sum(Gamma(3,:));
    ros_Gamma(4) = sum(Gamma(4,:));
    
    %~~~> M_i = Coefficients for new step solution bi
    ros_M=M;
    
    %~~~> E_i  = Coefficients for error estimator \hat b_i
    ros_E=E;
    
    %~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
    %   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
    ros_NewF(1)  = true;
    ros_NewF(2)  = true;
    ros_NewF(3)  = true;
    ros_NewF(4)  = true; 
    %~~~> ros_ELO  = estimator of local order - the minimum between the
    %    main and the embedded scheme orders plus 1
    rosELO  = 3.0d0;

return;