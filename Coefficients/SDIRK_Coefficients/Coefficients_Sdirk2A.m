function [ sdirkMethod sdirkELO sdirkS sdirkName ] = Coefficients_Sdirk2A( S2A )
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Filename: Coefficients_Sdirk2A.m
%
% Original Author:
%
% File Create Date:
%
% Input Arguments:
%   Name            Type
%                   Integer
%
% Output Arguments:
%   Name            Type
%   sdirkMethod     integer
%   sdirkELO        double
%   sdirkS          integer
%
% Modification History:
%   Date        Developer         Email             Action  
%   7/17/2012   Tony D'Augustine  adaug13@vt.edu    Translated Fortran90 Sdirk2A() method into MATLAB.
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Coefficients_Sdirk2A:
%   Initializes coefficients for Sdirk2A method.
%
%       Stages  Order
%
%
% Coefficients_Sdirk2A: INPUT ARGUMENTS
%   S2A (integer):
%
% Coefficients_Sdirk2A: OUTPUT ARGUMENTS
%   sdirkMethod (integer):
%   sdirkELO (double):
%   sdirkS (integer):
%
% Coefficients_Sdirk2A: GLOBAL VARIABLES
%   rkA (double):
%   rkB (double):
%   rkC (double):
%   rkD (double):
%   rkE (double):
%   rkBhat (double):
%   rkAlpha (double):
%   rkGamma (double):
%   rkTheta (double):
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    global rkA rkB rkC rkD rkE
    global rkBhat
    global rkAlpha rkGamma rkTheta

    sdirkMethod = S2A;

    % Method name
    sdirkName = 'Sdirk2A';

    % Number of stages
    sdirkS = 2;

    % Method coefficients
    rkGamma = .2928932188134524755991556378951510d0;

    rkA(1,1) = .2928932188134524755991556378951510d0;
    rkA(2,1) = .7071067811865475244008443621048490d0;
    rkA(2,2) = .2928932188134524755991556378951510d0;

    rkB(1)   = .7071067811865475244008443621048490d0;
    rkB(2)   = .2928932188134524755991556378951510d0;

    rkBhat(1)= .6666666666666666666666666666666667d0;
    rkBhat(2)= .3333333333333333333333333333333333d0;

    rkC(1)   = 0.292893218813452475599155637895151d0;
    rkC(2)   = 1.0d0;

    % Ynew = Yold + h*Sum_i {rkB_i*k_i} = Yold + Sum_i {rkD_i*Z_i}
    rkD(1)   = 0.0d0;
    rkD(2)   = 1.0d0;

    % Err = h * Sum_i {(rkB_i-rkBhat_i)*k_i} = Sum_i {rkE_i*Z_i}
    rkE(1)   =  0.4714045207910316829338962414032326d0;
    rkE(2)   = -0.1380711874576983496005629080698993d0;

    % Local order of Err estimate
    sdirkELO    = 2;

    % h*Sum_j {rkA_ij*k_j} = Sum_j {rkTheta_ij*Z_j}
    rkTheta(2,1) = 2.414213562373095048801688724209698d0;

    % Starting value for Newton iterations: Z_i^0 = Sum_j {rkAlpha_ij*Z_j}
    rkAlpha(2,1) = 3.414213562373095048801688724209698d0;

return;
