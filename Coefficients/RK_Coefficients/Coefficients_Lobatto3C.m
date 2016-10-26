function [ rkMethod, rkELO, rkS, rkName ] = Coefficients_Lobatto3C( L3C, ICNTRL )
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Filename: Coefficients_Lobatto3C.m
%
% Original Author:
%
% File Create Date:
%
% Input Arguments:
%   Name        Type
%   L3C         integer
%
% Output Arguments:
%   Name        Type
%   rkMethod    integer
%   rkELO       double
%   rkS         integer
%
% Modification History:
%   Date        Developer         Email             Action  
%   7/17/2012   Tony D'Augustine  adaug13@vt.edu    Translated Fortran90 Lobatto3C_Coefficients() method into MATLAB.
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Coefficients_Lobatto3C:
%   Initializes coefficients for Lobatto3C Runge-Kutta method.
%
%       Stages  Order   Stability Properties
%       3       4       Stiffly-accurrate
%
% Coefficients_Lobatto3C: INPUT ARGUMENTS
%   L3C (integer):
%
% Coefficients_Lobatto3C: OUTPUT ARGUMENTS
%   rkMethod (integer):
%   rkELO (integer):
%   rkS (integer):
%
% Coefficients_Lobatto3C: GLOBAL VARIABLES
%   rkA (double)
%   rkB (double)
%   rkC (double)
%   rkD (double)
%   rkE (double)
%   rkbgam (double)
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

    rkMethod = L3C;
    rkS = 'irrelevent';
    rkName = 'Lobatto3C';

    % b0 = 1.0d0
    if ( ICNTRL.SdirkError )
        b0 = 0.2d0;
    else
        b0 = 0.5d0;
    end
    
    % The coefficients of the Lobatto3C method
    rkA(1,1) =  .1666666666666666666666666666666667d0;
    rkA(1,2) = -.3333333333333333333333333333333333d0;
    rkA(1,3) =  .1666666666666666666666666666666667d0;
    rkA(2,1) =  .1666666666666666666666666666666667d0;
    rkA(2,2) =  .4166666666666666666666666666666667d0;
    rkA(2,3) = -.8333333333333333333333333333333333d-1;
    rkA(3,1) =  .1666666666666666666666666666666667d0;
    rkA(3,2) =  .6666666666666666666666666666666667d0;
    rkA(3,3) =  .1666666666666666666666666666666667d0;

    rkB(1) = .1666666666666666666666666666666667d0;
    rkB(2) = .6666666666666666666666666666666667d0;
    rkB(3) = .1666666666666666666666666666666667d0;

    rkC(1) = 0.0d0;
    rkC(2) = 0.5d0;
    rkC(3) = 1.0d0;

    % Classical error estimator, embedded solution: 
    rkBhat(1) = b0;
    rkBhat(2) = .16666666666666666666666666666666667d0-b0;
    rkBhat(3) = .66666666666666666666666666666666667d0;
    rkBhat(4) = .16666666666666666666666666666666667d0;

    % New solution: h Sum_j b_j f(Z_j) = sum d_j Z_j
    rkD(1) = 0.0d0;
    rkD(2) = 0.0d0;
    rkD(3) = 1.0d0;

    % Classical error estimator: 
    %   H* Sum (B_j-Bhat_j)*f(Z_j) = H*E(0)*f(0) + Sum E_j*Z_j
    rkE(1) =   .3808338772072650364017425226487022*b0;
    rkE(2) = -1.142501631621795109205227567946107*b0;
    rkE(3) = -1.523335508829060145606970090594809*b0;
    rkE(4) =   .3808338772072650364017425226487022*b0;

    % Sdirk error estimator
    rkBgam(1) = b0;
    rkBgam(2) = .1666666666666666666666666666666667d0-1.d0*b0;
    rkBgam(3) = .6666666666666666666666666666666667d0;
    rkBgam(4) = -.2141672105405983697350758559820354d0;
    rkBgam(5) = .3808338772072650364017425226487021d0;

    % H* Sum Bgam_j*f(Z_j) = H*Bgam(0)*f(0) + Sum Theta_j*Z_j
    rkTheta(1) = -3.d0*b0-.3808338772072650364017425226487021d0;
    rkTheta(2) = -4.d0*b0+1.523335508829060145606970090594808d0;
    rkTheta(3) = -.142501631621795109205227567946106d0+b0;

    % Local order of error estimator 
    if ( b0==0.0d0 )
        rkELO  = 5.0d0;
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

    rkGamma = 2.625816818958466716011888933765284d0;
    rkAlpha = 1.687091590520766641994055533117359d0;
    rkBeta  = 2.508731754924880510838743672432351d0;

    rkT(1,1) = 1.d0;
    rkT(1,2) = 1.d0;
    rkT(1,3) = 0.d0;
    rkT(2,1) = .4554100411010284672111720348287483d0;
    rkT(2,2) = -.6027050205505142336055860174143743d0;
    rkT(2,3) = -.4309321229203225731070721341350346d0;
    rkT(3,1) = 2.195823345445647152832799205549709d0;
    rkT(3,2) = -1.097911672722823576416399602774855d0;
    rkT(3,3) = .7850032632435902184104551358922130d0;

    rkTinv(1,1) = .4205559181381766909344950150991349d0;
    rkTinv(1,2) = .3488903392193734304046467270632057d0;
    rkTinv(1,3) = .1915253879645878102698098373933487d0;
    rkTinv(2,1) = .5794440818618233090655049849008650d0;
    rkTinv(2,2) = -.3488903392193734304046467270632057d0;
    rkTinv(2,3) = -.1915253879645878102698098373933487d0;
    rkTinv(3,1) = -.3659705575742745254721332009249516d0;
    rkTinv(3,2) = -1.463882230297098101888532803699806d0;
    rkTinv(3,3) = .4702733607340189781407813565524989d0;

    rkTinvAinv(1,1) = 1.104302803159744452668648155627548d0;
    rkTinvAinv(1,2) = .916122120694355522658740710823143d0;
    rkTinvAinv(1,3) = .5029105849749601702795812241441172d0;
    rkTinvAinv(2,1) = 1.895697196840255547331351844372453d0;
    rkTinvAinv(2,2) = 3.083877879305644477341259289176857d0;
    rkTinvAinv(2,3) = -1.502910584974960170279581224144117d0;
    rkTinvAinv(3,1) = .8362439183082935036129145574774502d0;
    rkTinvAinv(3,2) = -3.344975673233174014451658229909802d0;
    rkTinvAinv(3,3) = .312908409479233358005944466882642d0;

    rkAinvT(1,1) = 2.625816818958466716011888933765282d0;
    rkAinvT(1,2) = 1.687091590520766641994055533117358d0;
    rkAinvT(1,3) = -2.508731754924880510838743672432351d0;
    rkAinvT(2,1) = 1.195823345445647152832799205549710d0;
    rkAinvT(2,2) = -2.097911672722823576416399602774855d0;
    rkAinvT(2,3) = .7850032632435902184104551358922130d0;
    rkAinvT(3,1) = 5.765829871932827589653709477334136d0;
    rkAinvT(3,2) = .1170850640335862051731452613329320d0;
    rkAinvT(3,3) = 4.078738281412060947659653944216779d0;

return;

