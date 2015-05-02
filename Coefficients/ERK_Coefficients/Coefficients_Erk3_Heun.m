function [ erkMethod, erkELO, erkS, erkName ] = Coefficients_Erk3_Heun( RK3 )
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Filename: Coefficients_Erk3_Heun.m
%
% Original Author:
%
% File Create Date:
%
% Input Arguments:
%
% Output Arguments:
%
% Modification History:
%   Date        Developer         Email             Action  
%   7/17/2012   Tony D'Augustine  adaug13@vt.edu    Translated Fortran90 Erk3_Heun() method into MATLAB.
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Coefficients_Erk3_Heun:
%   Initializes coefficients for the Erk3_Heun Runge-Kutta method.
%
%       Stages  Order
%       ?       ?
%
% Coefficients_Erk3_Heun: INPUT ARGUMENTS
%   RK3 (integer)
%
% Coefficients_Erk3_Heun: OUTPUT ARGUMENTS
%   erkMethod (integer):
%   rkELO (double):
%   rkS (integer):
%
% Coefficients_Erk3_Heun: GLOBAL VARIABLES
%   rkA (double):
%   rkB (double):
%   rkC (double):
%   rkD (double):
%   
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    global rkA rkB rkC rkE

    erkMethod = RK3;
    erkELO = 3.0;
    erkS = 4;
    erkName = 'Erk3 Heun';
    
    rkA(2,1) = 1.d0 / 3.d0;
    rkA(3,2) = 2.d0 / 3.d0;
    rkA(4,1) = 0.25d0;
    rkA(4,3) = 0.75d0;
    
    rkB(1) = 0.25d0;
    rkB(2) = 0.0;
    rkB(3) = 0.75d0;
    rkB(4) = 0.0;
    
    rkC(1) = 0.0;
    rkC(2) = 1.d0 / 3.d0;
    rkC(3) = 2.d0 / 3.d0;
    rkC(4) = 0.0;

    rkE(1) = 0.25d0;
    rkE(2) = -0.5d0;
    rkE(3) = 0.25d0;
    rkE(4) = 0.0;

return;

