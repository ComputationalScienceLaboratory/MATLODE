function [ rosMethod rosELO rosS rosName ] = Coefficients_ROS34PW3( ROW2 )
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Filename: Coefficients_ROS34PW3.m
%
% Original Author:
%
% File Create Date:
%
% Input Arguments:
%   Name        Type
%   ROW2        integer
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
%   Initializes coefficients for ROS34PW3 Rosenbrock-w method.
% 
%       Stages  Order   Stability Property
%       4       4       A-stable 
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

    rosMethod = ROW2;
	% Method name
    rosName = 'ROS34PW3';
    % Number of stages
    rosS = 4;
    
    % Gamma Matrix
    Gamma_ii=  1.0685790213016289;
    Gamma= Gamma_ii*diag(ones(1,rosS));
    Gamma(2,1)=-2.5155456020628817e+00;
    Gamma(3,1)=-8.7991339217106512e-01;
    Gamma(3,2)=-9.6014187766190695e-01;
    Gamma(4,1)=-4.1731389379448741e-01;
    Gamma(4,2)= 4.1091047035857703e-01;
    Gamma(4,3)=-1.3558873204765276e+00;
    
    %%%% Alpha Matrix
    Alpha=zeros(rosS);
    Alpha(2,1)= 2.5155456020628817e+00;
    Alpha(3,1)= 5.0777280103144085e-01;
    Alpha(3,2)= 7.5000000000000000e-01;
    Alpha(4,1)= 1.3959081404277204e-01;
    Alpha(4,2)=-3.3111001065419338e-01;
    Alpha(4,3)= 8.2040559712714178e-01;
    
    %%%% Coefficients for new step solution bi
    b(1) = 2.2047681286931747e-01;
    b(2) = 2.7828278331185935e-03;
    b(3) = 7.1844787635140066e-03;
    b(4) = 7.6955588053404989e-01;
    
    %%%% Coefficients for error estimator \hat b_i
    b_hat(1) = 3.1300297285209688e-01;
    b_hat(2) =-2.8946895245112692e-01;
    b_hat(3) = 9.7646597959903003e-01;
    b_hat(4) = 0.0000000000000000e+00;
    b_hat=b_hat';
    b=b';
    
    %%%% Convert to Implementation Coeefiecnts 
    %%%% Refer to : Hairer, Wanner, BookII, CH IV.7
    A=Alpha/Gamma;
    C=diag(diag(inv(Gamma)))-inv(Gamma);
    M=transpose(inv(Gamma))*b;
    E=transpose(inv(Gamma))*(b-b_hat);
    
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
    rosELO  = 4.0d0;

return;
