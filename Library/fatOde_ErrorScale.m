%% fatOde_ErrorScale
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
function SCAL = fatOde_ErrorScale( NVAR, ITOL, AbsTol, RelTol, Y0 )
    SCAL = zeros(NVAR,1);

    if ( ITOL == 0 )
        for i = 1:NVAR
            SCAL(i) = 1.0 / ( AbsTol(1) + RelTol(1) * abs( Y0(i) ) );
        end
    else
        for i = 1:NVAR
            SCAL(i) = 1.0 / ( AbsTol(i) + RelTol(i) * abs( Y0(i) ) );
        end
    end
%     SCAL = 1.0 ./ ( AbsTol + RelTol .* abs(Y0) );

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
