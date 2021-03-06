function [ rosMethod rosELO rosS rosName ] = Coefficients_Ros3( RS3 )
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Filename: Coefficients_Ros3.m
%
% Original Author:
%
% File Create Date:
%
% Input Arguments:
%   Name        Type
%               integer
%
% Output Arguments:
%   Name        Type
%   rosMethod   integer
%   rosELO      double
%   rosS        integer
%
% Modification History:
%   Date        Developer         Email             Action  
%   7/17/2012   Tony D'Augustine  adaug13@vt.edu    Translated Fortran90 Ros3() method into MATLAB.
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Coefficients_Ros3:
%   Initializes coefficients for Ros3 Rosenbrock method.
% 
%       Stages  Order   Stability Property
%       2       2       L-stable
%
% Coefficients_Ros3: INPUT ARGUMENTS
%
% Coefficients_Ros3: OUTPUT ARGUMENTS
%   rosMethod (integer): Rosenbrock method number.
%   rosELO (double):
%   rosS (integer):
%
% Coefficients_Ros3: GLOBAL VARIABLES
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

    rosMethod = RS3;

    % Method name
    rosName = 'Ros3';

    % Number of stages
    rosS = 3;
   
%~~~> The coefficient matrices A and C are strictly lower triangular.
%   The lower triangular (subdiagonal) elements are stored in row-wise order:
%   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
%   The general mapping formula is:
%       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )    
%       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )  
   ros_A(1)= 1.d0;
   ros_A(2)= 1.d0;
   ros_A(3)= 0.d0;

   ros_C(1) = -0.10156171083877702091975600115545d+01;
   ros_C(2) =  0.40759956452537699824805835358067d+01;
   ros_C(3) =  0.92076794298330791242156818474003d+01;
   
%~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
%   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
   ros_NewF(1) = true;
   ros_NewF(2) = true;
   ros_NewF(3) = false;
%~~~> M_i = Coefficients for new step solution
   ros_M(1) =  0.1d+01;
   ros_M(2) =  0.61697947043828245592553615689730d+01;
   ros_M(3) = -0.42772256543218573326238373806514d+00;
% E_i = Coefficients for error estimator    
   ros_E(1) =  0.5d+00;
   ros_E(2) = -0.29079558716805469821718236208017d+01;
   ros_E(3) =  0.22354069897811569627360909276199d+00;
%~~~> ros_ELO = estimator of local order - the minimum between the
%    main and the embedded scheme orders plus 1
   rosELO = 3.0d0;   
%~~~> Y_stage_i ~ Y( T + H*Alpha_i )
   ros_Alpha(1)= 0.0d+00;
   ros_Alpha(2)= 0.43586652150845899941601945119356d+00;
   ros_Alpha(3)= 0.43586652150845899941601945119356d+00;
%~~~> Gamma_i = \sum_j  gamma_{i,j}    
   ros_Gamma(1)= 0.43586652150845899941601945119356d+00;
   ros_Gamma(2)= 0.24291996454816804366592249683314d+00;
   ros_Gamma(3)= 0.21851380027664058511513169485832d+01;

return;
