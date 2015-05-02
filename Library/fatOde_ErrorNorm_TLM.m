function [ FWD_Err ] = fatOde_ErrorNorm_TLM( NVAR, FWD_Err, Yerr_TLM, OPTIONS )
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Filename: fatOde_ErrorNorm_TLM.m
%
% Original Author: 
%
% File Creation Date: 
%
% Input Arguments:
%
% Output Arguments:
% 
% Modification History:
%   Date        Developer         Email             Action  
%   7/17/2012   Tony D'Augustine  adaug13@vt.edu    Translated Fortran90
%                                                   ErrorNorm_TLM() to MATLAB.
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% fatOde_ErrorNorm_TLM:
%   Calculate the error norm for the tangent linear model implementation.
%
% fatOde_ErrorNorm_TLM: INPUT ARGUMENTS
%   NVAR (integer):
%   NTLM (integer):
%   Y_TLM (integer):
%
% fatOde_ErrorNorm_TLM: OUTPUT ARGUMENTS
%   FWD_Err (double):
%
% fatOde_ErrorNorm_TLM: SYNTAX
%
% fatOde_ErrorNorm_TLM: EXAMPLE
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%     SCAL = 1.0 ./ OPTIONS.AbsTol_TLM+OPTIONS.RelTol_TLM .* abs(Yerr_TLM);
%     FWD_Err = max([FWD_Err, sqrt(sum((OPTIONS.Y_TLM.*SCAL).^2,1)./NVAR), 1.0d-10] );
    
    
   NTLM = max(size(OPTIONS.Y_TLM));

    for itlm=1:NTLM
        SCAL_TLM = fatOde_ErrorScale( NVAR, OPTIONS.ITOL, OPTIONS.AbsTol_TLM(:,itlm), ...
            OPTIONS.RelTol_TLM(:,itlm), Yerr_TLM(:,itlm) );
        
        % Perform error norm
        Err = sum((OPTIONS.Y_TLM(:,itlm).*SCAL_TLM).^2,1);
        Err = max( sqrt( Err/double(NVAR) ), 1.0d-10 );
        
        FWD_Err = max( FWD_Err, Err );
    end

return;

