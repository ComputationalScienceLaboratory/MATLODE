function [ sdirkMethod sdirkELO sdirkS sdirkName ] = Coefficients_Sdirk4B( S4B )
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Filename: Coefficients_Sdirk4B.m
%
% Original Author:
%
% File Create Date:
%
% Input Arguments:
%   Name            Type
%   S4B             Integer
%
% Output Arguments:
%   Name            Type
%   sdirkMethod     integer
%   sdirkELO        double
%   sdirkS          integer
%
% Modification History:
%   Date        Developer         Email             Action  
%   7/17/2012   Tony D'Augustine  adaug13@vt.edu    Translated Fortran90 Sdirk4B() method into MATLAB.
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Coefficients_Sdirk4B:
%   Initializes coefficients for Sdirk4B method.
%
%       Stages  Order
%
%
% Coefficients_Sdirk4B: INPUT ARGUMENTS
%   S4B (integer):
%
% Coefficients_Sdirk4B: OUTPUT ARGUMENTS
%   sdirkMethod (integer):
%   sdirkELO (double):
%   sdirkS (integer):
%
% Coefficients_Sdirk4B: GLOBAL VARIABLES
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

    sdirkMethod = S4B;

    % Method name
    sdirkName = 'Sdirk4B';

    % Number of stages
    sdirkS = 5;

    % Method coefficients
    rkGamma = .25d0;

    rkA(1,1) = 0.25d0;
    rkA(2,1) = 0.5d00;
    rkA(2,2) = 0.25d0;
    rkA(3,1) = 0.34d0;
    rkA(3,2) =-0.40d-1;
    rkA(3,3) = 0.25d0;
    rkA(4,1) = 0.2727941176470588235294117647058824d0;
    rkA(4,2) =-0.5036764705882352941176470588235294d-1;
    rkA(4,3) = 0.2757352941176470588235294117647059d-1;
    rkA(4,4) = 0.25d0;
    rkA(5,1) = 1.041666666666666666666666666666667d0;
    rkA(5,2) =-1.020833333333333333333333333333333d0;
    rkA(5,3) = 7.812500000000000000000000000000000d0;
    rkA(5,4) =-7.083333333333333333333333333333333d0;
    rkA(5,5) = 0.25d0;

    rkB(1)   =  1.041666666666666666666666666666667d0;
    rkB(2)   = -1.020833333333333333333333333333333d0;
    rkB(3)   =  7.812500000000000000000000000000000d0;
    rkB(4)   = -7.083333333333333333333333333333333d0;
    rkB(5)   =  0.250000000000000000000000000000000d0;

    rkBhat(1)=  1.069791666666666666666666666666667d0;
    rkBhat(2)= -0.894270833333333333333333333333333d0;
    rkBhat(3)=  7.695312500000000000000000000000000d0;
    rkBhat(4)= -7.083333333333333333333333333333333d0;
    rkBhat(5)=  0.212500000000000000000000000000000d0;

    rkC(1)   = 0.25d0;
    rkC(2)   = 0.75d0;
    rkC(3)   = 0.55d0;
    rkC(4)   = 0.50d0;
    rkC(5)   = 1.00d0;

    % Ynew = Yold + h*Sum_i {rkB_i*k_i} = Yold + Sum_i {rkD_i*Z_i}
    rkD(1)   = 0.0d0;
    rkD(2)   = 0.0d0;
    rkD(3)   = 0.0d0;
    rkD(4)   = 0.0d0;
    rkD(5)   = 1.0d0;

    % Err = h * Sum_i {(rkB_i-rkBhat_i)*k_i} = Sum_i {rkE_i*Z_i}
    rkE(1)   =  0.5750d0;
    rkE(2)   =  0.2125d0;
    rkE(3)   = -4.6875d0;
    rkE(4)   =  4.2500d0;
    rkE(5)   =  0.1500d0;

    % Local order of Err estimate
    sdirkELO    = 4;

    % h*Sum_j {rkA_ij*k_j} = Sum_j {rkTheta_ij*Z_j}
    rkTheta(2,1) = 2.d0;
    rkTheta(3,1) = 1.680000000000000000000000000000000d0;
    rkTheta(3,2) = -.1600000000000000000000000000000000d0;
    rkTheta(4,1) = 1.308823529411764705882352941176471d0;
    rkTheta(4,2) = -.1838235294117647058823529411764706d0;
    rkTheta(4,3) = 0.1102941176470588235294117647058824d0;
    rkTheta(5,1) = -3.083333333333333333333333333333333d0;
    rkTheta(5,2) = -4.291666666666666666666666666666667d0;
    rkTheta(5,3) =  34.37500000000000000000000000000000d0;
    rkTheta(5,4) = -28.33333333333333333333333333333333d0;

    % Starting value for Newton iterations: Z_i^0 = Sum_j {rkAlpha_ij*Z_j}
    rkAlpha(2,1) = 3.0;
    rkAlpha(3,1) = .8800000000000000000000000000000000d0;
    rkAlpha(3,2) = .4400000000000000000000000000000000d0;
    rkAlpha(4,1) = .1666666666666666666666666666666667d0;
    rkAlpha(4,2) = -.8333333333333333333333333333333333d-1;
    rkAlpha(4,3) = .9469696969696969696969696969696970d0;
    rkAlpha(5,1) = -6.d0;
    rkAlpha(5,2) = 9.d0;
    rkAlpha(5,3) = -56.81818181818181818181818181818182d0;
    rkAlpha(5,4) = 54.d0;

return;
