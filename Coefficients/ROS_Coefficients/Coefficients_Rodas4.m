function [ rosMethod, rosELO, rosS, rosName ] = Coefficients_Rodas4( RD4 )
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Filename: Coefficients_Rodas4.m
%
% Original Author:
%
% File Create Date:
%
% Input Arguments:
%   Name        Type
%   L3C         integer
%
% Output Arguments:
%   Name        Type
%   rkMethod    integer
%   rkELO       double
%   rkS         integer
%
% Modification History:
%   Date        Developer         Email             Action  
%   7/17/2012   Tony D'Augustine  adaug13@vt.edu    Translated Fortran90 Rodas4() method into MATLAB.
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Coefficients_Rodas4:
%   Initializes coefficients for Rodas4 Rosenbrock method.
%
% Coefficients_Rodas4: INPUT ARGUMENTS
%   RD4 (integer): Rosenbrock method number.
%
% Coefficients_Rodas4: OUTPUT ARGUMENTS
%   rosMethod (integer): Rosenbrock method number.
%   rosELO (doulbe): Estimator of local error.
%   rosS (integer): Nunmber of stages.
%
% Coefficients_Rodas4: GLOBAL VARIABLES
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

    % Rosenbrock method number
    rosMethod = RD4;

    % Number of stages
    rosS = 6;
    
    % rosELO  = estimator of local order - the minimum between the
    %           main and the embedded scheme orders plus 1
    rosELO = 4.0d0;

    % Method name
    rosName = 'Rodas4';

%~~~> Y_stage_i ~ Y( T + H*Alpha_i )
    ros_Alpha(1) = 0.000d0;
    ros_Alpha(2) = 0.386d0;
    ros_Alpha(3) = 0.210d0; 
    ros_Alpha(4) = 0.630d0;
    ros_Alpha(5) = 1.000d0;
    ros_Alpha(6) = 1.000d0;
        
%~~~> Gamma_i = \sum_j  gamma_{i,j}    
    ros_Gamma(1) = 0.2500000000000000d+00;
    ros_Gamma(2) =-0.1043000000000000d+00;
    ros_Gamma(3) = 0.1035000000000000d+00;
    ros_Gamma(4) =-0.3620000000000023d-01;
    ros_Gamma(5) = 0.0d0;
    ros_Gamma(6) = 0.0d0;

%~~~> The coefficient matrices A and C are strictly lower triangular.
%   The lower triangular (subdiagonal) elements are stored in row-wise order:
%   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
%   The general mapping formula is:  A(i,j) = ros_A( (i-1)*(i-2)/2 + j )    
%                  C(i,j) = ros_C( (i-1)*(i-2)/2 + j )  
     
    ros_A(1) = 0.1544000000000000d+01;
    ros_A(2) = 0.9466785280815826d+00;
    ros_A(3) = 0.2557011698983284d+00;
    ros_A(4) = 0.3314825187068521d+01;
    ros_A(5) = 0.2896124015972201d+01;
    ros_A(6) = 0.9986419139977817d+00;
    ros_A(7) = 0.1221224509226641d+01;
    ros_A(8) = 0.6019134481288629d+01;
    ros_A(9) = 0.1253708332932087d+02;
    ros_A(10) =-0.6878860361058950d+00;
    ros_A(11) = ros_A(7);
    ros_A(12) = ros_A(8);
    ros_A(13) = ros_A(9);
    ros_A(14) = ros_A(10);
    ros_A(15) = 1.0d+00;

    ros_C(1) =-0.5668800000000000d+01;
    ros_C(2) =-0.2430093356833875d+01;
    ros_C(3) =-0.2063599157091915d+00;
    ros_C(4) =-0.1073529058151375d+00;
    ros_C(5) =-0.9594562251023355d+01;
    ros_C(6) =-0.2047028614809616d+02;
    ros_C(7) = 0.7496443313967647d+01;
    ros_C(8) =-0.1024680431464352d+02;
    ros_C(9) =-0.3399990352819905d+02;
    ros_C(10) = 0.1170890893206160d+02;
    ros_C(11) = 0.8083246795921522d+01;
    ros_C(12) =-0.7981132988064893d+01;
    ros_C(13) =-0.3152159432874371d+02;
    ros_C(14) = 0.1631930543123136d+02;
    ros_C(15) =-0.6058818238834054d+01;

%~~~> M_i = Coefficients for new step solution
    ros_M(1) = ros_A(7);
    ros_M(2) = ros_A(8);
    ros_M(3) = ros_A(9);
    ros_M(4) = ros_A(10);
    ros_M(5) = 1.0d+00;
    ros_M(6) = 1.0d+00;

%~~~> E_i  = Coefficients for error estimator    
    ros_E(1) = 0.0d+00;
    ros_E(2) = 0.0d+00;
    ros_E(3) = 0.0d+00;
    ros_E(4) = 0.0d+00;
    ros_E(5) = 0.0d+00;
    ros_E(6) = 1.0d+00;

%~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
%     or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
    ros_NewF(1) = true;
    ros_NewF(2) = true;
    ros_NewF(3) = true;
    ros_NewF(4) = true;
    ros_NewF(5) = true;
    ros_NewF(6) = true;

return;

