function [ erkMethod, erkELO, erkS, erkName ] = Coefficients_Dopri5( RK5 )
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Filename: Coefficients_Dopri5.m
%
% Original Author:
%
% File Create Date:
%
% Input Arguments:
%   Name    Type
%   RK5     integer
%
% Output Arguments:
%   Name        Type 
%   erkMethod   integer
%   rkELO       double
%   rkS         integer
%   erkName     string
%
% Modification History:
%   Date        Developer         Email             Action  
%   7/17/2012   Tony D'Augustine  adaug13@vt.edu    Translated Fortran90 Dopri5() method into MATLAB.
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Coefficients_Dopri5:
%   Initializes coefficients for the Dopri5 Explict Runge-Kutta method.
%       
%       Stages  Order
%       7       5
%
% Coefficients_Dopri5: INPUT ARGUMENTS
%   RK5 (integer):
%
% Coefficients_Dopri5: OUTPUT ARGUMENTS
%   erkMethod (integer):
%        rkELO (double):
%         rkS (integer):
%      erkName (string):
%
% Coefficients_Dopri5: GLOBAL VARIABLES
%   rkA (double):
%   rkB (double):
%   rkC (double):
%   rkD (double):
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%
% Authored by Tony D'Augustine, Adrian Sandu, and Hong Zhang.
% Computational Science Laboratory, Virginia Tech.
% ©2015 Virginia Tech Intellectual Properties, Inc.
%%
    global rkA rkB rkC rkE
    
    % Clear all global variables to ensire correct matrix size and values.
    rkA = []; rkB = []; rkC = []; rkE =[];

    erkMethod = RK5;
    erkELO = 6.0;
    erkS = 7;
%     erkS = 6;
    erkName = 'Dopri5';

    rkA(1,1) = 0; 
    rkA(1,2) = 0; 
    rkA(1,3) = 0; 
    rkA(1,4) = 0; 
    rkA(1,5) = 0;
    rkA(1,6) = 0; 
    rkA(1,7) = 0;

    rkA(2,1) = .2;
    rkA(2,2) = 0; 
    rkA(2,3) = 0; 
    rkA(2,4) = 0;
    rkA(2,5) = 0;
    rkA(2,6) = 0;
    rkA(2,7) = 0;

    rkA(3,1) = 3.0/40.0;
    rkA(3,2) = 9.0/40.0;
    rkA(3,3) = 0;
    rkA(3,4) = 0;
    rkA(3,5) = 0;
    rkA(3,6) = 0;
    rkA(3,7) = 0;

    rkA(4,1) = 44.0/45.0;
    rkA(4,2) = -56.0/15.0;
    rkA(4,3) = 32.0/9.0;
    rkA(4,4) = 0;
    rkA(4,5) = 0;
    rkA(4,6) = 0;
    rkA(4,7) = 0;

    rkA(5,1) = 19372.0/6561.0;
    rkA(5,2) = -25360.0/2187.0;
    rkA(5,3) = 64448.0/6561.0;
    rkA(5,4) = -212.0/729.0;
    rkA(5,5) = 0;
    rkA(5,6) = 0;
    rkA(5,7) = 0;

    rkA(6,1) = 9017.0 / 3168.0;
    rkA(6,2) = -355.0 / 33.0;
    rkA(6,3) = 46732.0 / 5247.0;
    rkA(6,4) = 49.0 / 176.0;
    rkA(6,5) = -5103.0 / 18656.0;
    rkA(6,6) = 0;
    rkA(6,7) = 0;

    rkA(7,1) = 35.0 / 384.0;
    rkA(7,2) = 0;
    rkA(7,3) = 500.0 / 1113.0;
    rkA(7,4) = 125.0 / 192.0;
    rkA(7,5) = -2187.0 / 6784.0;
    rkA(7,6) = 11.0 / 84.0;
    rkA(7,7) = 0;
             
    rkB   = rkA(7,:);              

    rkC(1)   = 0.0;
    rkC(2)   = .20;
    rkC(3)   = .30;
    rkC(4)   = .80;
    rkC(5)   = 8.0/9.0;
    rkC(6)   = 1.0;
    rkC(7)   = 1.0;

    rkE(1)   = 71.0/57600.0;
    rkE(2)   = 0;
    rkE(3)   = -71.0/16695.0;
    rkE(4)   = 71.0/1920.0;
    rkE(5)   = -17253.0/339200.0;
    rkE(6)   = 22.0/525.0;
    rkE(7)   = -1.0/40.0;

return;

