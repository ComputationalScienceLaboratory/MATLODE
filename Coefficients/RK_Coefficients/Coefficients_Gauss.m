function [ rkMethod rkELO rkS rkName ] = Coefficients_Gauss( GAU )
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Filename: Coefficients_Gauss.m
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
%   7/17/2012   Tony D'Augustine  adaug13@vt.edu    Translated Fortran90 Gauss_Coefficients() method into MATLAB.
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Coefficients_Gauss:
%   Initializes coefficients for Gauss Runge-Kutta method.
%
%       Stages  Order   Stability Property
%       3       6       Stiffly-accurrate
%
% Coefficients_Gauss: INPUT ARGUMENTS
%   GAU (integer):
%
% Coefficients_Gauss: OUTPUT ARGUEMENTS
%   rkMethod (integer):
%   rkELO (integer):
%   rkS (integer):
%
% Coefficients_Gauss: GLOBAL VARIABLES
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
    

    % The coefficients of the Gauss method
    rkMethod = GAU;
    rkS = 'irrelevant';
    rkName = 'Gauss';

    % b0 = 4.0d0
    b0 = 0.1d0;

    % The coefficients of the Gauss method

    rkA(1,1) =  .1388888888888888888888888888888889d0;
    rkA(1,2) = -.359766675249389034563954710966045d-1;
    rkA(1,3) =  .97894440153083260495800422294756d-2;
    rkA(2,1) =  .3002631949808645924380249472131556d0;
    rkA(2,2) =  .2222222222222222222222222222222222d0;
    rkA(2,3) = -.224854172030868146602471694353778d-1;
    rkA(3,1) =  .2679883337624694517281977355483022d0;
    rkA(3,2) =  .4804211119693833479008399155410489d0;
    rkA(3,3) =  .1388888888888888888888888888888889d0;

    rkB(1) = .2777777777777777777777777777777778d0;
    rkB(2) = .4444444444444444444444444444444444d0;
    rkB(3) = .2777777777777777777777777777777778d0;

    rkC(1) = .1127016653792583114820734600217600d0;
    rkC(2) = .5000000000000000000000000000000000d0;
    rkC(3) = .8872983346207416885179265399782400d0;

    % Classical error estimator, embedded solution: 
    rkBhat(1) = b0;
    rkBhat(2) =-1.4788305577012361475298775666303999d0*b0 + .27777777777777777777777777777777778d0;
    rkBhat(3) =  .44444444444444444444444444444444444d0 + .66666666666666666666666666666666667d0*b0;
    rkBhat(4) = -.18783610896543051913678910003626672d0*b0 + .27777777777777777777777777777777778d0;

    % New solution: h Sum_j b_j f(Z_j) = sum d_j Z_j
    rkD(1) = .1666666666666666666666666666666667d1;
    rkD(2) = -.1333333333333333333333333333333333d1;
    rkD(3) = .1666666666666666666666666666666667d1;

    % Classical error estimator: 
    %   H* Sum (B_j-Bhat_j)*f(Z_j) = H*E(0)*f(0) + Sum E_j*Z_j
    rkE(1) = .2153144231161121782447335303806954d0*b0;
    rkE(2) = -2.825278112319014084275808340593191d0*b0;
    rkE(3) = .2870858974881495709929780405075939d0*b0;
    rkE(4) = -.4558086256248162565397206448274867d-1*b0;

    % Sdirk error estimator
    rkBgam(1) = 0.d0;
    rkBgam(2) = .2373339543355109188382583162660537d0;
    rkBgam(3) = .5879873931885192299409334646982414d0;
    rkBgam(4) = -.4063577064014232702392531134499046d-1;
    rkBgam(5) = .2153144231161121782447335303806955d0;

    % H* Sum Bgam_j*f(Z_j) = H*Bgam(0)*f(0) + Sum Theta_j*Z_j
    rkTheta(1) = -2.594040933093095272574031876464493d0;
    rkTheta(2) = 1.824611539036311947589425112250199d0;
    rkTheta(3) = .1856563166634371860478043996459493d0;

    % ELO = local order of classical error estimator 
    rkELO = 4.0d0;

    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %~~~> Diagonalize the RK matrix:               
    % rkTinv * inv(rkA) * rkT =          
    %           |  rkGamma      0           0     |
    %           |      0      rkAlpha   -rkBeta   |
    %           |      0      rkBeta     rkAlpha  |
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    rkGamma = 4.644370709252171185822941421408064d0;
    rkAlpha = 3.677814645373914407088529289295970d0;
    rkBeta  = 3.508761919567443321903661209182446d0;

    rkT(1,1) =  .7215185205520017032081769924397664d-1;
    rkT(1,2) = -.8224123057363067064866206597516454d-1;
    rkT(1,3) = -.6012073861930850173085948921439054d-1;
    rkT(2,1) =  .1188325787412778070708888193730294d0;
    rkT(2,2) =  .5306509074206139504614411373957448d-1;
    rkT(2,3) =  .3162050511322915732224862926182701d0;
    rkT(3,1) = 1.d0;
    rkT(3,2) = 1.d0;
    rkT(3,3) = 0.d0;

    rkTinv(1,1) =  5.991698084937800775649580743981285d0;
    rkTinv(1,2) =  1.139214295155735444567002236934009d0;
    rkTinv(1,3) =   .4323121137838583855696375901180497d0;
    rkTinv(2,1) = -5.991698084937800775649580743981285d0;
    rkTinv(2,2) = -1.139214295155735444567002236934009d0;
    rkTinv(2,3) =   .5676878862161416144303624098819503d0;
    rkTinv(3,1) = -1.246213273586231410815571640493082d0;
    rkTinv(3,2) =  2.925559646192313662599230367054972d0;
    rkTinv(3,3) =  -.2577352012734324923468722836888244d0;

    rkTinvAinv(1,1) =  27.82766708436744962047620566703329d0;
    rkTinvAinv(1,2) =   5.290933503982655311815946575100597d0;
    rkTinvAinv(1,3) =   2.007817718512643701322151051660114d0;
    rkTinvAinv(2,1) = -17.66368928942422710690385180065675d0;
    rkTinvAinv(2,2) = -14.45491129892587782538830044147713d0;
    rkTinvAinv(2,3) =   2.992182281487356298677848948339886d0;
    rkTinvAinv(3,1) = -25.60678350282974256072419392007303d0;
    rkTinvAinv(3,2) =   6.762434375611708328910623303779923d0;
    rkTinvAinv(3,3) =   1.043979339483109825041215970036771d0;

    rkAinvT(1,1) = .3350999483034677402618981153470483d0;
    rkAinvT(1,2) = -.5134173605009692329246186488441294d0;
    rkAinvT(1,3) = .6745196507033116204327635673208923d-1;
    rkAinvT(2,1) = .5519025480108928886873752035738885d0;
    rkAinvT(2,2) = 1.304651810077110066076640761092008d0;
    rkAinvT(2,3) = .9767507983414134987545585703726984d0;
    rkAinvT(3,1) = 4.644370709252171185822941421408064d0;
    rkAinvT(3,2) = 3.677814645373914407088529289295970d0;
    rkAinvT(3,3) = -3.508761919567443321903661209182446d0;

return;

