function [ erkMethod rkELO rkS erkName ] = Coefficients_Erk23( RK2 )
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Filename: Coefficients_Erk23.m
%
% Original Author:
%
% File Creation Date:
%
% Input Arguments:
%   Name            Type
%   RK2             integer
%
% Output Arguments:
%   Name            Type
%   erkMethod       integer
%   rkELO           double
%   rkS             integer
% 
% Modification History:
%   Date        Developer         Email             Action  
%   7/17/2012   Tony D'Augustine  adaug13@vt.edu    Translated Fortran90 Erk23() method into MATLAB.
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Coefficients_Erk23:
%   Initializes coefficients for the RK2(3) Runge-Kutta method. This is a
%   2(3) Fehlberg-type method: second order accurate with a third order
%   error control mechanism.
%
%       Stages  Order
%       3       2
%
% Coefficients_Erk23: INPUT ARGUEMENTS
%   RK2 (integer):
%
% Coefficients_Erk23: OUTPUT ARGUMENTS
%   erkMethod (integer):
%   rkELO (double):
%   rkS (integer):
%
% Coefficients_Erk23: GLOBAL VARIABLES
%   rkA (double): 3x3 matrix
%   rkB (double): 1x3 matrix
%   rkC (double): 1x3 matrix
%   rkE (double)" 1x3 matrix
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    global rkA rkB rkC rkE
    
    erkMethod = RK2;
    rkELO = 3.0;
    rkS = 3;
    erkName = 'Erk23';
    
    rkA(2,1) = 1.d0;
    rkA(3,1) = 0.25d0;
    rkA(3,2) = 0.25d0;
    
    rkB(1) = 0.5d0;
    rkB(2) = 0.5d0;
    rkB(3) = 0.0;
    
    rkC(1) = 0.0;
    rkC(2) = 1.d0;
    rkC(3) = 0.5d0;

    rkE(1) = 1.d0/3.d0;
    rkE(2) = 1.d0/3.d0;
    rkE(3) = -2.d0/3.d0;

return;

