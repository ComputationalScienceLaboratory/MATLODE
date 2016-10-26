function [ sdirkMethod sdirkELO sdirkS sdirkName ] = Coefficients_Sdirk2B( S2B )
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Filename: Coefficients_Sdirk2B.m
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
%   7/17/2012   Tony D'Augustine  adaug13@vt.edu    Translated Fortran90 Sdirk2B() method into MATLAB.
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Coefficients_Sdirk2B:
%   Initializes coefficients for Sdirk3A method.
%
%       Stages  Order
%
%
% Coefficients_Sdirk2B: INPUT ARGUMENTS
%   S2B (integer):
%
% Coefficients_Sdirk2B: OUTPUT ARGUMENTS
%   sdirkMethod (integer):
%   sdirkELO (double):
%   sdirkS (double):
%
% Coefficients_Sdirk2B: GLOBAL VARIABLES
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

    sdirkMethod = S2B;

    % Method name
    sdirkName = 'Sdirk2B';

    % Number of stages
    sdirkS      = 2;

    % Method coefficients
    rkGamma  = 1.707106781186547524400844362104849d0;

    rkA(1,1) = 1.707106781186547524400844362104849d0;
    rkA(2,1) = -.707106781186547524400844362104849d0;
    rkA(2,2) = 1.707106781186547524400844362104849d0;

    rkB(1)   = -.707106781186547524400844362104849d0;
    rkB(2)   = 1.707106781186547524400844362104849d0;

    rkBhat(1)= .6666666666666666666666666666666667d0;
    rkBhat(2)= .3333333333333333333333333333333333d0;

    rkC(1)   = 1.707106781186547524400844362104849d0;
    rkC(2)   = 1.0d0;

    % Ynew = Yold + h*Sum_i {rkB_i*k_i} = Yold + Sum_i {rkD_i*Z_i}
    rkD(1)   = 0.0d0;
    rkD(2)   = 1.0d0;

    % Err = h * Sum_i {(rkB_i-rkBhat_i)*k_i} = Sum_i {rkE_i*Z_i}
    rkE(1)   = -.4714045207910316829338962414032326d0;
    rkE(2)   =  .8047378541243650162672295747365659d0;

    % Local order of Err estimate
    sdirkELO    = 2;

    % h*Sum_j {rkA_ij*k_j} = Sum_j {rkTheta_ij*Z_j}
    rkTheta(2,1) = -.414213562373095048801688724209698d0;

    % Starting value for Newton iterations: Z_i^0 = Sum_j {rkAlpha_ij*Z_j}
    rkAlpha(2,1) = .5857864376269049511983112757903019d0;

return;
