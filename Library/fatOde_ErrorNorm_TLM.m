%% fatOde_ErrorNorm_TLM
%
% <html>
%   <div>
%       <img style="float: right" src="../../MATLODE_LOGO.png" height="150px"></img>
%   </div>
% </html>
%
%% Syntax
%
%% Input Parameters
%
%% Output Parameters
%
%% Description
%
%% Reference
% [1] Tony D'Augustine, Adrian Sandu. MATLODE: A MATLAB ODE Solver and
%     Sensitivity Analysis Toolbox. Submitted to ACM TOMS.
%
%%
%  Authored by Tony D'Augustine, Adrian Sandu, and Hong Zhang.
%  Computational Science Laboratory, Virginia Tech.
%  ©2015 Virginia Tech Intellectual Properties, Inc.
%
function [ FWD_Err ] = fatOde_ErrorNorm_TLM( NVAR, FWD_Err, Yerr_TLM, OPTIONS )

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

%% Major Modification History
% <html>
% <table border=1>
%   <tr>
%       <td><b>Date</b></td>
%       <td>Developer</td>
%       <td>Email</td>
%       <td>Action</td>
%   </tr>
%   <tr>
%       <td>1/1/2014</td>
%       <td>Tony D'Augustine</td>
%       <td>adaug13@vt.edu</td>
%       <td>Release MATLODE_v2.0.00</td>
%   </tr>
% </table>
% </html>
% 
%%
% <html>
%   <div>
%       <img style="float: right" src="../../CSL_LogoWithName_1.png" height="50px"></img>
%   </div>
% </html>
