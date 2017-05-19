function [ rkMethod rkELO rkS rkName ] = Coefficients_Radau2A( R2A, ICNTRL )
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Filename: Coefficients_Radau2A.m
%
% Original Author:
%
% File Create Date:
%
% Input Arguments:
%   Name        Type
%   R2A         integer
%
% Output Arguments:
%   Name        Type
%   rkMethod    integer
%   rkELO       double
%   rkS         integer
%
% Modification History:
%   Date        Developer         Email             Action  
%   7/17/2012   Tony D'Augustine  adaug13@vt.edu    Translated Fortran90 Radau2A_Coefficients() method 
%                                                   into MATLAB.
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Coefficients_Radau2A:
%   Initializes coefficients for Radau2A Runge-Kutta method.
%
%       Stages  Order   Stability Property
%       3       5       Stiffly-accurate
%
% Coefficients_Radau2A:
%   R2A (integer):
%
% Coefficients_Radau2A: INPUT ARGUMENTS
%   rkMethod (integer):
%   rkELO (integer):
%   rkS (integer):
%
% Coefficients_Radau2A: GLOBAL VARIABLES
%   rkA (double)
%   rkB (double)
%   rkC (double)
%   rkD (double)
%   rkE (double)
%   rkBgam (double)
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
    global rkBgam
    global rkAlpha rkBeta rkGamma rkTheta
    global rkT rkTinv rkTinvAinv rkAinvT

    % b0 = 1.0d0
    if ( ICNTRL.SdirkError )
        b0 = 0.2d-1;
    else
        b0 = 0.5d-1;
    end

    % The coefficients of the Radau2A method
    rkMethod = R2A;
    rkS = 'Irrelevent';
    rkName = 'Radau2A';

    rkA(1,1) =  1.968154772236604258683861429918299d-1;
    rkA(1,2) = -6.55354258501983881085227825696087d-2;
    rkA(1,3) =  2.377097434822015242040823210718965d-2;
    rkA(2,1) =  3.944243147390872769974116714584975d-1;
    rkA(2,2) =  2.920734116652284630205027458970589d-1;
    rkA(2,3) = -4.154875212599793019818600988496743d-2;
    rkA(3,1) =  3.764030627004672750500754423692808d-1;
    rkA(3,2) =  5.124858261884216138388134465196080d-1;
    rkA(3,3) =  1.111111111111111111111111111111111d-1;

    rkB(1) = 3.764030627004672750500754423692808d-1;
    rkB(2) = 5.124858261884216138388134465196080d-1;
    rkB(3) = 1.111111111111111111111111111111111d-1;

    rkC(1) = 1.550510257216821901802715925294109d-1;
    rkC(2) = 6.449489742783178098197284074705891d-1;
    rkC(3) = 1.0d0;

    % New solution: H* Sum B_j*f(Z_j) = Sum D_j*Z_j
    rkD(1) = 0.0d0;
    rkD(2) = 0.0d0;
    rkD(3) = 1.0d0;

    % Classical error estimator: 
    % H* Sum (B_j-Bhat_j)*f(Z_j) = H*E(0)*f(0) + Sum E_j*Z_j
    rkE(1) = 1.0d0*b0;
    rkE(2) = -10.04880939982741556246032950764708d0*b0;
    rkE(3) = 1.382142733160748895793662840980412d0*b0;
    rkE(4) = -.3333333333333333333333333333333333d0*b0;

    % Sdirk error estimator
    rkBgam(1) = b0;
    rkBgam(2) = .3764030627004672750500754423692807d0-1.558078204724922382431975370686279d0*b0;
    rkBgam(3) = .8914115380582557157653087040196118d0*b0+.5124858261884216138388134465196077d0;
    rkBgam(4) = -.1637777184845662566367174924883037d0-.3333333333333333333333333333333333d0*b0;
    rkBgam(5) = .2748888295956773677478286035994148d0;

    % H* Sum Bgam_j*f(Z_j) = H*Bgam(0)*f(0) + Sum Theta_j*Z_j
    rkTheta(1) = -1.520677486405081647234271944611547d0-10.04880939982741556246032950764708d0*b0;
    rkTheta(2) = 2.070455145596436382729929151810376d0+1.382142733160748895793662840980413d0*b0;
    rkTheta(3) = -.3333333333333333333333333333333333d0*b0-.3744441479783868387391430179970741d0;

    % Local order of error estimator 
    if ( b0 == 0.0 )
        rkELO  = 6.0d0;
    else	
        rkELO  = 4.0d0;
    end	

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

    rkT(1,1) =  9.443876248897524148749007950641664d-2;
    rkT(1,2) = -1.412552950209542084279903838077973d-1;
    rkT(1,3) = -3.00291941051474244918611170890539d-2;
    rkT(2,1) =  2.502131229653333113765090675125018d-1;
    rkT(2,2) =  2.041293522937999319959908102983381d-1;
    rkT(2,3) =  3.829421127572619377954382335998733d-1;
    rkT(3,1) =  1.0d0;
    rkT(3,2) =  1.0d0;
    rkT(3,3) =  0.0d0;

    rkTinv(1,1) =  4.178718591551904727346462658512057d0;
    rkTinv(1,2) =  3.27682820761062387082533272429617d-1;
    rkTinv(1,3) =  5.233764454994495480399309159089876d-1;
    rkTinv(2,1) = -4.178718591551904727346462658512057d0;
    rkTinv(2,2) = -3.27682820761062387082533272429617d-1;
    rkTinv(2,3) =  4.766235545005504519600690840910124d-1;
    rkTinv(3,1) = -5.02872634945786875951247343139544d-1;
    rkTinv(3,2) =  2.571926949855605429186785353601676d0;
    rkTinv(3,3) = -5.960392048282249249688219110993024d-1;

    rkTinvAinv(1,1) =  1.520148562492775501049204957366528d+1;
    rkTinvAinv(1,2) =  1.192055789400527921212348994770778d0;
    rkTinvAinv(1,3) =  1.903956760517560343018332287285119d0;
    rkTinvAinv(2,1) = -9.669512977505946748632625374449567d0;
    rkTinvAinv(2,2) = -8.724028436822336183071773193986487d0;
    rkTinvAinv(2,3) =  3.096043239482439656981667712714881d0;
    rkTinvAinv(3,1) = -1.409513259499574544876303981551774d+1;
    rkTinvAinv(3,2) =  5.895975725255405108079130152868952d0;
    rkTinvAinv(3,3) = -1.441236197545344702389881889085515d-1;

    rkAinvT(1,1) = .3435525649691961614912493915818282d0;
    rkAinvT(1,2) = -.4703191128473198422370558694426832d0;
    rkAinvT(1,3) = .3503786597113668965366406634269080d0;
    rkAinvT(2,1) = .9102338692094599309122768354288852d0;
    rkAinvT(2,2) = 1.715425895757991796035292755937326d0;
    rkAinvT(2,3) = .4040171993145015239277111187301784d0;
    rkAinvT(3,1) = 3.637834252744495732208418513577775d0;
    rkAinvT(3,2) = 2.681082873627752133895790743211112d0;
    rkAinvT(3,3) = -3.050430199247410569426377624787569d0;

return;

