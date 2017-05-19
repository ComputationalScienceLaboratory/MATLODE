function [ rosMethod rosELO rosS rosName ] = Coefficients_Ros2( RS2 )
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Filename: Coefficients_Ros2.m
%
% Original Author:
%
% File Create Date:
%
% Input Arguments:
%   Name        Type
%   RD3         integer
%
% Output Arguments:
%   Name        Type
%   rosMethod   integer
%   rosELO      double
%   rosS        integer
%
% Modification History:
%   Date        Developer         Email             Action  
%   7/17/2012   Tony D'Augustine  adaug13@vt.edu    Translated Fortran90 Ros2() method into MATLAB.
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Coefficients_Ros2:
%   Initializes coefficients for Ros2 Rosenbrock method.
%
% Coefficients_Ros2:
%   RS2 (integer): Rosenbrock method number.
%
% Coefficients_Ros2: INPUT ARGUMENTS
%   rosMethod (integer):
%   rosELO (double):
%   rosS (integer):
%
% Coefficients_Ros2: GLOBAL VARIABLES
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

    g = 1.0d0 + 1.0d0/sqrt(2.0d0);
   
    rosMethod = RS2;

    % Method name
    rosName = 'Ros2';

%~~~> Number of stages
    rosS = 2;
   
%~~~> Y_stage_i ~ Y( T + H*Alpha_i )
    ros_Alpha(1) = 0.0d0;
    ros_Alpha(2) = 1.0d0 ;
    
%~~~> Gamma_i = \sum_j  gamma_{i,j}
    ros_Gamma(1) = g;
    ros_Gamma(2) =-g;
    
%~~~> The coefficient matrices A and C are strictly lower triangular.
%   The lower triangular (subdiagonal) elements are stored in row-wise order:
%   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
%   The general mapping formula is:
%       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )    
%       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )    
    ros_A(1) = (1.d0)/g;
    ros_C(1) = (-2.d0)/g;
    
%~~~> M_i = Coefficients for new step solution
    ros_M(1)= (3.d0)/(2.d0*g);
    ros_M(2)= (1.d0)/(2.d0*g);
    
% E_i = Coefficients for error estimator    
    ros_E(1) = 1.d0/(2.d0*g);
    ros_E(2) = 1.d0/(2.d0*g);
    
%~~~> ros_ELO = estimator of local order - the minimum between the
%    main and the embedded scheme orders plus one
    rosELO = 2.0d0;
    
%~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
%   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
    ros_NewF(1) = true;
    ros_NewF(2) = true;
    
return;
