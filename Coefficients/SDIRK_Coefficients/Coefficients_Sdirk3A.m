function [ sdirkMethod sdirkELO sdirkS sdirkName ] = Coefficients_Sdirk3A( S3A )
% Filename: Coefficients_Sdirk3A.m
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
%   7/17/2012   Tony D'Augustine  adaug13@vt.edu    Translated Fortran90 Sdirk3A() method into MATLAB.
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Coefficients_Sdirk3A:
%   Initializes coefficients for Sdirk3A method.
%
%       Stages  Order
%
%
% Coefficients_Sdirk3A: INPUT ARGUMENTS
%   S3A (integer):
%
% Coefficients_Sdirk3A: OUTPUT ARGUMENTS
%   sdirkMethod (integer):
%   sdirkELO (double):
%   sdirkS (double):
%
% Coefficients_Sdirk3A: GLOBAL VARIABLES
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

    sdirkMethod = S3A;

    % Method name
    sdirkName = 'Sdirk3A';

    % Number of stages
    sdirkS = 3;

    % Method coefficients
    rkGamma = .2113248654051871177454256097490213d0;

    rkA(1,1) = .2113248654051871177454256097490213d0;
    rkA(2,1) = .2113248654051871177454256097490213d0;
    rkA(2,2) = .2113248654051871177454256097490213d0;
    rkA(3,1) = .2113248654051871177454256097490213d0;
    rkA(3,2) = .5773502691896257645091487805019573d0;
    rkA(3,3) = .2113248654051871177454256097490213d0;

    rkB(1)   = .2113248654051871177454256097490213d0;
    rkB(2)   = .5773502691896257645091487805019573d0;
    rkB(3)   = .2113248654051871177454256097490213d0;

    rkBhat(1)= .2113248654051871177454256097490213d0;
    rkBhat(2)= .6477918909913548037576239837516312d0;
    rkBhat(3)= .1408832436034580784969504064993475d0;

    rkC(1)   = .2113248654051871177454256097490213d0;
    rkC(2)   = .4226497308103742354908512194980427d0;
    rkC(3)   = 1.d0;

    % Ynew = Yold + h*Sum_i {rkB_i*k_i} = Yold + Sum_i {rkD_i*Z_i}
    rkD(1)   = 0.d0;
    rkD(2)   = 0.d0;
    rkD(3)   = 1.d0;

    % Err = h * Sum_i {(rkB_i-rkBhat_i)*k_i} = Sum_i {rkE_i*Z_i}
    rkE(1)   =  0.9106836025229590978424821138352906d0;
    rkE(2)   = -1.244016935856292431175815447168624d0;
    rkE(3)   =  0.3333333333333333333333333333333333d0;

    % Local order of Err estimate
    sdirkELO    = 2;

    % h*Sum_j {rkA_ij*k_j} = Sum_j {rkTheta_ij*Z_j}
    rkTheta(2,1) =  1.0d0;
    rkTheta(3,1) = -1.732050807568877293527446341505872d0;
    rkTheta(3,2) =  2.732050807568877293527446341505872d0;

    % Starting value for Newton iterations: Z_i^0 = Sum_j {rkAlpha_ij*Z_j}
    rkAlpha(2,1) =   2.0d0;
    rkAlpha(3,1) = -12.92820323027550917410978536602349d0;
    rkAlpha(3,2) =   8.83012701892219323381861585376468d0;

return;
