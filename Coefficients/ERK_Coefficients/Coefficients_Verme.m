function [ erkMethod rkELO rkS erkName ] = Coefficients_Verme( RK6 )
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Filename: Coefficients_Verme.m
%
% Original Author:
%
% File Create Date:
%
% Input Arguments:
%   Name    Type 
%   RK6     integer
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
%   7/17/2012   Tony D'Augustine  adaug13@vt.edu    Translated Fortran90 Verme() method into MATLAB.
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Coefficients_Verme:
%   Initializes coefficients for the Verme Runge-Kutta method.
%
%       Stages  Order
%       8       6
%
% Coefficients_Verme: INPUT ARGUMENTS
%   RK6 (integer):
%
% Coefficients_Verme: OUTPUT ARGUMENTS
%   erkMethod (integer):
%   rkELO (double):
%   rkS (integer):
%   erkName (string):
%
% Coefficients_Verme: GLOBAL VARIABLES
%   rkA (double):
%   rkB (double):
%   rkC (double):
%   rkD (double):
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    global rkA rkB rkC rkE

    erkMethod = RK6;
    rkELO = 7.0;
    rkS = 8;
    erkName = 'Verme';

    rkA(1,1) = 0;
    rkA(1,2) = 0;
    rkA(1,3) = 0;
    rkA(1,4) = 0;
    rkA(1,5) = 0;
    rkA(1,6) = 0;
    rkA(1,7) = 0;
    rkA(1,8) = 0;

    rkA(2,1) = 1.0/6.0;
    rkA(2,2) = 0;
    rkA(2,3) = 0;
    rkA(2,4) = 0;
    rkA(2,5) = 0;
    rkA(2,6) = 0;
    rkA(2,7) = 0;
    rkA(2,8) = 0;

    rkA(3,1) = 4.0/75.0;
    rkA(3,2) = 16.d0/75.d0;
    rkA(4,1) = 5.d0/6.d0;
    rkA(4,2) = -8.d0/3.d0;
    rkA(4,3) = 5.d0/2.d0;
    rkA(5,1) = -165.d0/64.d0;
    rkA(5,2) = 55.d0/6.d0;
    rkA(5,3) = -425.d0/64.d0;
    rkA(5,4) = 85.d0/96.d0;
    rkA(6,1) = 12.d0/5.d0;
    rkA(6,2) = -8.d0;
    rkA(6,3) = 4015.d0/612.d0;
    rkA(6,4) = -11.d0/36.d0;
    rkA(6,5) = 88.d0/255.d0;
    rkA(7,1) = -8263.d0/15000.d0;
    rkA(7,2) = 124.d0/75.d0;
    rkA(7,3) = -643.d0/680.d0;
    rkA(7,4) = -81.d0/250.d0;
    rkA(7,5) = 2484.d0/10625.d0;
    rkA(7,6) = 0.d0;
    rkA(8,1) = 3501.d0/1720.d0;
    rkA(8,2) = -300.d0/43.d0;
    rkA(8,3) = 297275.d0/52632.d0;
    rkA(8,4) = -319.d0/2322.d0;
    rkA(8,5) = 24068.d0/84065.d0;
    rkA(8,6) = 0.d0;
    rkA(8,7) = 3850.d0/26703.d0;

    rkB(1) = 3.d0/40.d0;
    rkB(2) = 0.d0;
    rkB(3) = 875.d0/2244.d0;
    rkB(4) = 23.d0/72.d0;
    rkB(5) = 264.d0/1955.d0;
    rkB(6) = 0.d0;
    rkB(7) = 125.d0/11592.d0;
    rkB(8) = 43.d0/616.d0;

    rkC(2) = 1.d0/6.d0;
    rkC(3) = 4.d0/15.d0;
    rkC(4) = 2.d0/3.d0;
    rkC(5) = 5.d0/6.d0;
    rkC(6) = 1.0d0;
    rkC(7) = 1.d0/15.d0;
    rkC(8) = 1.d0;

    rkE(1) = -1.d0/160.d0;
    rkE(2) = 0.d0;
    rkE(3) = 875.d0/2244.d0 - 2375.d0/5984.d0;
    rkE(4) = 23.d0/72.d0 - 5.d0/16.d0;
    rkE(5) = 264.d0/1955.d0 - 12.d0/85.d0;
    rkE(6) = -3.d0/44.d0;
    rkE(7) = 125.d0/11592.d0;
    rkE(8) = 43.d0/616.d0;

return;

