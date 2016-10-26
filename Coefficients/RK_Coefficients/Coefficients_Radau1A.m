function [ rkMethod rkELO rkS rkName ] = Coefficients_Radau1A( R1A )
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Filename: Coefficients_Radau1A.m
%
% Original Author:
%
% File Create Date:
%
% Input Arguments:
%   Name        Type
%   Gau         integer
%
% Output Arguments:
%   Name        Type
%   rkMethod    integer
%   rkELO       double
%   rkS         integer
%
% Modification History:
%   Date        Developer         Email             Action  
%   7/17/2012   Tony D'Augustine  adaug13@vt.edu    Translated Fortran90 Radau1A_Coefficients() method 
%                                                   into MATLAB.
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Coefficients_Radau1A:
%   Initializes coefficients for Radau1A Runge-Kutta method.
%
%       Stages  Order   Stability Properties
%       3       5       Stiffly-accurate
%
% Coefficients_Radau1A: INPUT ARGUMENTS
%   R1A (integer):
%
% Coefficients_Radau!A: OUTPUT ARGUMENTS
%   rkMethod (integer):
%   rkELO (integer):
%   rkS (integer):
%
% Coefficients_Radau1A: GLOBAL VARIABLES
%   rkA (double)
%   rkB (double)
%   rkC (double)
%   rkD (double)
%   rkE (double)
%   rkBgam (double)
%   rkBhat (double)
%   rkAlpha (double)
%   rkBeta (double)
%   rkGamma (double)
%   rkTheta (double)
%   rkT (double)
%   rkTinv (double)
%   rkTinvAinv (double)
%   rkAinvT (double)
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    global rkA rkB rkC rkD rkE
    global rkBgam rkBhat
    global rkAlpha rkBeta rkGamma rkTheta
    global rkT rkTinv rkTinvAinv rkAinvT

    % The coefficients of the Radau1A method
    
    rkMethod = R1A;
    rkS = 'irrelevant';
    b0 = 0.1;
    rkName = 'Radau1A';

    rkA(1,1) =  .1111111111111111111111111111111111d0;
    rkA(1,2) = -.1916383190435098943442935597058829d0;
    rkA(1,3) =  .8052720793239878323318244859477174d-1;
    rkA(2,1) =  .1111111111111111111111111111111111d0;
    rkA(2,2) =  .2920734116652284630205027458970589d0;
    rkA(2,3) = -.481334970546573839513422644787591d-1;
    rkA(3,1) =  .1111111111111111111111111111111111d0;
    rkA(3,2) =  .5370223859435462728402311533676479d0;
    rkA(3,3) =  .1968154772236604258683861429918299d0;

    rkB(1) = .1111111111111111111111111111111111d0;
    rkB(2) = .5124858261884216138388134465196080d0;
    rkB(3) = .3764030627004672750500754423692808d0;

    rkC(1) = 0.d0;
    rkC(2) = .3550510257216821901802715925294109d0;
    rkC(3) = .8449489742783178098197284074705891d0;

    % Classical error estimator, embedded solution: 
    rkBhat(1) = b0;
    rkBhat(2) = .11111111111111111111111111111111111d0-b0;
    rkBhat(3) = .51248582618842161383881344651960810d0;
    rkBhat(4) = .37640306270046727505007544236928079d0;

    % New solution: H* Sum B_j*f(Z_j) = Sum D_j*Z_j
    rkD(1) = .3333333333333333333333333333333333d0;
    rkD(2) = -.8914115380582557157653087040196127d0;
    rkD(3) = .1558078204724922382431975370686279d1;

    % Classical error estimator: 
    % H* Sum (b_j-bhat_j) f(Z_j) = H*E(0)*F(0) + Sum E_j Z_j
    rkE(1) =   .2748888295956773677478286035994148d0*b0;
    rkE(2) = -1.374444147978386838739143017997074d0*b0;
    rkE(3) = -1.335337922441686804550326197041126d0*b0;
    rkE(4) =   .235782604058977333559011782643466d0*b0;

    % Sdirk error estimator
    rkBgam(1) = 0.0d0;
    rkBgam(2) = .1948150124588532186183490991130616d-1;
    rkBgam(3) = .7575249005733381398986810981093584d0;
    rkBgam(4) = -.518952314149008295083446116200793d-1;
    rkBgam(5) = .2748888295956773677478286035994148d0;

    % H* Sum Bgam_j*f(Z_j) = H*Bgam(0)*f(0) + Sum Theta_j*Z_j
    rkTheta(1) = -1.224370034375505083904362087063351d0;
    rkTheta(2) = .9340045331532641409047527962010133d0;
    rkTheta(3) = .4656990124352088397561234800640929d0;

    % ELO = local order of classical error estimator 
    rkELO = 4.0d0;

    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %~~~> Diagonalize the RK matrix:               
    % rkTinv * inv(rkA) * rkT =          
    %           |  rkGamma      0           0     |
    %           |      0      rkAlpha   -rkBeta   |
    %           |      0      rkBeta     rkAlpha  |
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    rkGamma = 3.637834252744495732208418513577775d0;
    rkAlpha = 2.681082873627752133895790743211112d0;
    rkBeta  = 3.050430199247410569426377624787569d0;

    rkT(1,1) =  .424293819848497965354371036408369d0;
    rkT(1,2) = -.3235571519651980681202894497035503d0;
    rkT(1,3) = -.522137786846287839586599927945048d0;
    rkT(2,1) =  .57594609499806128896291585429339d-1;
    rkT(2,2) =  .3148663231849760131614374283783d-2;
    rkT(2,3) =  .452429247674359778577728510381731d0;
    rkT(3,1) = 1.d0;
    rkT(3,2) = 1.d0;
    rkT(3,3) = 0.d0;

    rkTinv(1,1) = 1.233523612685027760114769983066164d0;
    rkTinv(1,2) = 1.423580134265707095505388133369554d0;
    rkTinv(1,3) = .3946330125758354736049045150429624d0;
    rkTinv(2,1) = -1.233523612685027760114769983066164d0;
    rkTinv(2,2) = -1.423580134265707095505388133369554d0;
    rkTinv(2,3) = .6053669874241645263950954849570376d0;
    rkTinv(3,1) = -.1484438963257383124456490049673414d0;
    rkTinv(3,2) = 2.038974794939896109682070471785315d0;
    rkTinv(3,3) = -.544501292892686735299355831692542d-1;

    rkTinvAinv(1,1) =  4.487354449794728738538663081025420d0;
    rkTinvAinv(1,2) =  5.178748573958397475446442544234494d0;
    rkTinvAinv(1,3) =  1.435609490412123627047824222335563d0;
    rkTinvAinv(2,1) = -2.854361287939276673073807031221493d0;
    rkTinvAinv(2,2) = -1.003648660720543859000994063139137d+1;
    rkTinvAinv(2,3) =  1.789135380979465422050817815017383d0;
    rkTinvAinv(3,1) = -4.160768067752685525282947313530352d0;
    rkTinvAinv(3,2) =  1.124128569859216916690209918405860d0;
    rkTinvAinv(3,3) =  1.700644430961823796581896350418417d0;

    rkAinvT(1,1) = 1.543510591072668287198054583233180d0;
    rkAinvT(1,2) = -2.460228411937788329157493833295004d0;
    rkAinvT(1,3) = -.412906170450356277003910443520499d0;
    rkAinvT(2,1) = .209519643211838264029272585946993d0;
    rkAinvT(2,2) = 1.388545667194387164417459732995766d0;
    rkAinvT(2,3) = 1.20339553005832004974976023130002d0;
    rkAinvT(3,1) = 3.637834252744495732208418513577775d0;
    rkAinvT(3,2) = 2.681082873627752133895790743211112d0;
    rkAinvT(3,3) = -3.050430199247410569426377624787569d0;

return;

