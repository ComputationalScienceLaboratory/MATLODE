function [ erkMethod erkELO erkS erkName ] = Coefficients_Erk43( RK4 )
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Filename: Coefficients_Erk43.m
%
% Original Author:
%
% File Create Date:
%
% Input Arguments:
%   Name    Type 
%   RK4     integer
%
% Output Arguments:
%   Name        Type
%   erkMethod   integer
%   erkELO      double
%   erkS        integer
%   erkName     string
%
% Modification History:
%   Date        Developer         Email             Action  
%   7/17/2012   Tony D'Augustine  adaug13@vt.edu    Translated Fortran90 Dopri5() method into MATLAB.
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Coefficients_Erk43:
%   Initializes coefficients for the Erk43 Runge-Kutta method.
%
%       Stages  Order
%       5       4
%
% Coefficients_Erk43: INPUT ARGUMENTS
%   RK4 (integer):
%
% Coefficients_Erk43: OUTPUT ARGUMENTS
%   erkMethod (integer):
%       erkELO (double):
%        erkS (integer):
%      erkName (string):
%
% Coefficients_Erk43: GLOBAL VARIABLES
%   rkA (double):
%   rkB (double):
%   rkC (double):
%   rkD (double):
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    global rkA rkB rkC rkE

    erkMethod = RK4;
    erkELO = 5.0;
    erkS = 5;
    erkName = 'Erk43';

    rkA(2,1) = 1.d0/3.d0;
    rkA(3,1) = -1.d0/3.d0;
    rkA(3,2) = 1.d0;
    rkA(4,1) = 1.d0;
    rkA(4,2) =-1.d0;
    rkA(4,3) = 1.d0;
    rkA(5,1) = 1.d0/8.d0;
    rkA(5,2) = 3.d0/8.d0;
    rkA(5,3) = 3.d0/8.d0;
    rkA(5,4) = 1.d0/8.d0;

    rkB(1)   = 1.0/8.0;
    rkB(2)   = 3.0/8.0;
    rkB(3)   = 3.0/8.0;
    rkB(4)   = 1.0/8.0;
    rkB(5)   = 0;

    rkC(2)   = 1.d0/3.d0;
    rkC(3)   = 2.d0/3.d0;
    rkC(4)   = 1.d0;
    rkC(5)   = 1.d0;

    rkE(1)   =  1.d0/24.d0;
    rkE(2)   =  -1.d0/8.d0;
    rkE(3)   =  1.d0/8.d0;
    rkE(4)   =  1.d0/8.d0;
    rkE(5)   =  -1.d0/6.d0;
    
return;

