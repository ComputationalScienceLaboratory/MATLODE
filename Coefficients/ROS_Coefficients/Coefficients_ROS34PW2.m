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
%   RS4         integer
%
% Output Arguments:
%   Name        Type
%   rosMethod   integer
%   rosELO      double
%   rosS        integer
%
% Modification History:
%   Date        Developer         Email             Action  
%   7/17/2012   Tony D'Augustine  adaug13@vt.edu    Translated Fortran90 Ros4() method into MATLAB.
%   6/16/2015   Arash Sarshar	  sarshar@vt.edu    Added ROW method based on RANG paper
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Coefficients_Ros4:
%   Initializes coefficients for ROS34PW2 Rosenbrock-w method.
% 
%       Stages  Order   Stability Property
%       4       3       L-stable
%
% Coefficients_Ros4: INPUT ARGUMENTS
%   ROW (integer):
%
% Coefficients_Ros4: OUTPUT ARGUMENTS
%   rosMethod (integer): Rosenbrock method number.
%   rosELO (double):
%   rosS (integer):
%
% Coefficients_Ros4: GLOBAL VARIABLES
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

    % rosMethod = RS4
    rosMethod = ROW1;
    
	% Method name
    rosName = 'ROS34PW2';

    % Number of stages
    rosS = 4;

     
%~~~> The coefficient matrices A and C are strictly lower triangular.
%   The lower triangular (subdiagonal) elements are stored in row-wise order:
%   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
%   The general mapping formula is:
%       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )    
%       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )  
    ros_A(1) =  8.7173304301691801d-01;
    ros_A(2) =  8.4457060015369423d-01;
    ros_A(3) = -1.1299064236484185d-01;
    ros_A(4) = 0.0d0;
    ros_A(5) = 0.0d0;
    ros_A(6) = 1.0d0;

    ros_C(1) =-8.7173304301691801d-01;
    ros_C(2) =-9.0338057013044082d-01;
    ros_C(3) = 5.4180672388095326d-02;
    ros_C(4) = 2.4212380706095346d-01;
    ros_C(5) =-1.2232505839045147d+00;
    ros_C(6) = 5.4526025533510214d-01;
    
      %~~~> Y_stage_i ~ Y( T + H*Alpha_i ) \sum(j=1,j=i-1) A(i,j)
    ros_Alpha(1) = 0.0d0;
    ros_Alpha(2) = ros_A(1);
    ros_Alpha(3) = ros_A(2)+ros_A(3);
    ros_Alpha(4) = ros_A(4)+ros_A(5)+ros_A(6);
    
    %~~~> Gamma_i = \sum_j  gamma_{i,j} is it ????    
    ros_Gamma(1) = 4.3586652150845900d-01;
    ros_Gamma(2) = ros_C(1)+ros_Gamma(1);
    ros_Gamma(3) = ros_C(2)+ros_C(3)+ros_Gamma(1);
    ros_Gamma(4) = ros_C(4)+ros_C(5)+ros_C(6)+ros_Gamma(1);
    
    %~~~> M_i = Coefficients for new step solution bi
    ros_M(1) = 2.4212380706095346d-01;
    ros_M(2) =-1.2232505839045147d+00;
    ros_M(3) = 1.5452602553351020d+00;
    ros_M(4) = 4.3586652150845900d-01;
    
    %~~~> E_i  = Coefficients for error estimator \hat b_i
    ros_E(1) = 3.7810903145819369d-01;
    ros_E(2) =-9.6042292212423178d-02;
    ros_E(3) = 5.0000000000000000d-01;
    ros_E(4) = 2.1793326075422950d-01;
    
    
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ros_A(1) =  2.218787467653286;
    ros_A(2) =  0.0;
    ros_A(3) =  0.0;
    ros_A(4) = 1.208587690772214;
    ros_A(5) = 7.511610241919324e-2;
    ros_A(6) = 5.0e-1;
    
    ros_C(1) =-2.218787467653286;
    ros_C(2) =-9.461966143940745e-2;
    ros_C(3) =-7.913526735718213e-3;
    ros_C(4) =-1.870323744195384;
    ros_C(5) =-9.624340112825115e-2;
    ros_C(6) = 2.726301276675511e-1;
    
    ros_Gamma(1) = 4.358665215084590e-1;
    
        %~~~> M_i = Coefficients for new step solution bi
    ros_M(1) = 3.285609536316354e-1;
    ros_M(2) =-5.785609536316354e-1;
    ros_M(3) = 2.5e-1;
    ros_M(4) = 1.0;
    
    %~~~> E_i  = Coefficients for error estimator \hat b_i
    ros_E(1) = -2.5e-1;
    ros_E(2) = 0.0;
    ros_E(3) = 2.5e-1;
    ros_E(4) = 1.0;
    
    
    
    
    %~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
    %   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
    ros_NewF(1)  = true;
    ros_NewF(2)  = true;
    ros_NewF(3)  = true;
    ros_NewF(4)  = true; %%test tes

    %~~~> ros_ELO  = estimator of local order - the minimum between the
    %    main and the embedded scheme orders plus 1
    rosELO  = 4.0d0;

return;
