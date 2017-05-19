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
%     % Number of stages
    rosS = 4;
    ros_A =[2 1.41921731745576 -0.25923221167297 4.18476048231916 -0.285192017355496 2.29428036027904];ros_Alpha =[0 0.871733043016918 0.731579957788852 1];ros_C =[-4.58856072055809 -4.18476048231916 0.285192017355496 -6.36817920012836 -6.79562094446684 2.87009860433106];ros_E =[0.277749947647967;-1.403239895176;1.77263012766755;0.5];ros_Gamma =[0.435866521508459 -0.435866521508459 -0.413333376233886 -5.55111512312578e-17];ros_M =[4.18476048231916;-0.285192017355496;2.29428036027904;1];
   
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
