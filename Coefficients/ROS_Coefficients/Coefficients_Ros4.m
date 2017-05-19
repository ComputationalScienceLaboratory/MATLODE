function [ rosMethod rosELO rosS rosName ] = Coefficients_Ros4( RS4 )
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Filename: Coefficients_Ros4.m
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
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Coefficients_Ros4:
%   Initializes coefficients for Ros4 Rosenbrock method.
% 
%       Stages  Order   Stability Property
%       4       4       L-stable
%
% Coefficients_Ros4: INPUT ARGUMENTS
%   RS4 (integer):
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
    rosMethod = RS4;
    
	% Method name
    rosName = 'Ros4';

    % Number of stages
    rosS = 4;

    %~~~> Y_stage_i ~ Y( T + H*Alpha_i )
    ros_Alpha(1) = 0.D0;
    ros_Alpha(2) = 0.114640000000000d+01;
    ros_Alpha(3) = 0.6552168638155900d+00;
    ros_Alpha(4) = ros_Alpha(3);
    
    %~~~> Gamma_i = \sum_j  gamma_{i,j}    
    ros_Gamma(1) = 0.5728200000000000d+00;
    ros_Gamma(2) =-0.1769193891319233d+01;
    ros_Gamma(3) = 0.7592633437920482d+00;
    ros_Gamma(4) =-0.1049021087100450d+00;
    
%~~~> The coefficient matrices A and C are strictly lower triangular.
%   The lower triangular (subdiagonal) elements are stored in row-wise order:
%   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
%   The general mapping formula is:
%       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )    
%       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )  
    ros_A(1) = 0.2000000000000000d+01;
    ros_A(2) = 0.1867943637803922d+01;
    ros_A(3) = 0.2344449711399156d+00;
    ros_A(4) = ros_A(2);
    ros_A(5) = ros_A(3);
    ros_A(6) = 0.0D0;

    ros_C(1) =-0.7137615036412310d+01;
    ros_C(2) = 0.2580708087951457d+01;
    ros_C(3) = 0.6515950076447975d+00;
    ros_C(4) =-0.2137148994382534d+01;
    ros_C(5) =-0.3214669691237626d+00;
    ros_C(6) =-0.6949742501781779d+00;
    
    %~~~> M_i = Coefficients for new step solution
    ros_M(1) = 0.2255570073418735d+01;
    ros_M(2) = 0.2870493262186792d+00;
    ros_M(3) = 0.4353179431840180d+00;
    ros_M(4) = 0.1093502252409163d+01;
    
    %~~~> E_i  = Coefficients for error estimator    
    ros_E(1) =-0.2815431932141155d+00;
    ros_E(2) =-0.7276199124938920d-01;
    ros_E(3) =-0.1082196201495311d+00;
    ros_E(4) =-0.1093502252409163d+01;
    
    %~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
    %   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
    ros_NewF(1)  = true;
    ros_NewF(2)  = true;
    ros_NewF(3)  = true;
    ros_NewF(4)  = false;

    %~~~> ros_ELO  = estimator of local order - the minimum between the
    %    main and the embedded scheme orders plus 1
    rosELO  = 4.0d0;

return;
