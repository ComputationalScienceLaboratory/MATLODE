function [ rosMethod rosELO rosS rosName ] = Coefficients_Rodas3( RD3 )
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Filename: Coefficients_Rodas3.m
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
%   7/17/2012   Tony D'Augustine  adaug13@vt.edu    Translated Fortran90 Rodas3() method into MATLAB.
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Coefficients_Rodas3:
%   Initializes coefficients for Rodas3 Rosenbrock method.
%
%       Stages  Order   Stability Property
%       4       3       Stiffly-stable
%
% Coefficients_Rodas3: INPUT ARGUMENTS
%   RD3 (integer): Rosenbrock method number.
%
% Coefficients_Rodas3: OUTPUT ARGUMENTS
%   rosMethod (integer): Rosenbrock method number.
%   rosELO (double):
%   rosS (integer):
%
% Coefficients_Rodas3: GLOBAL VARIABLES
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

    
    rosMethod = RD3;
    
%~~~> Name of the method
	rosName = 'Rodas3';
    
%~~~> Number of stages
    rosS = 4;
    
%~~~> Y_stage_i ~ Y( T + H*Alpha_i )
    ros_Alpha(1) = 0.0d+00;
    ros_Alpha(2) = 0.0d+00;
    ros_Alpha(3) = 1.0d+00;
    ros_Alpha(4) = 1.0d+00;
%~~~> Gamma_i = \sum_j  gamma_{i,j}    
    ros_Gamma(1) = 0.5d+00;
    ros_Gamma(2) = 1.5d+00;
    ros_Gamma(3) = 0.0d+00;
    ros_Gamma(4) = 0.0d+00;

%~~~> The coefficient matrices A and C are strictly lower triangular.
%   The lower triangular (subdiagonal) elements are stored in row-wise order:
%   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
%   The general mapping formula is:
%       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )    
%       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )  
    ros_A(1) = 0.0d+00;
    ros_A(2) = 2.0d+00;
    ros_A(3) = 0.0d+00;
    ros_A(4) = 2.0d+00;
    ros_A(5) = 0.0d+00;
    ros_A(6) = 1.0d+00;

    ros_C(1) = 4.0d+00;
    ros_C(2) = 1.0d+00;
    ros_C(3) =-1.0d+00;
    ros_C(4) = 1.0d+00;
    ros_C(5) =-1.0d+00;
    ros_C(6) =-(8.0d+00/3.0d+00);
    
%~~~> M_i = Coefficients for new step solution
    ros_M(1) = 2.0d+00;
    ros_M(2) = 0.0d+00;
    ros_M(3) = 1.0d+00;
    ros_M(4) = 1.0d+00;
%~~~> E_i  = Coefficients for error estimator    
    ros_E(1) = 0.0d+00;
    ros_E(2) = 0.0d+00;
    ros_E(3) = 0.0d+00;
    ros_E(4) = 1.0d+00;

%~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
%   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
    ros_NewF(1)  = true;
    ros_NewF(2)  = false;
    ros_NewF(3)  = true;
    ros_NewF(4)  = true;

%~~~> ros_ELO  = estimator of local order - the minimum between the
%    main and the embedded scheme orders plus 1
    rosELO  = 3.0d+00;

return;

