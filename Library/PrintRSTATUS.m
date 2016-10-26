%% PrintRSTATUS
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
function PrintRSTATUS( RSTATUS )

    Ntexit_formatSpec = '    Ntexit: %6.4f \n';
    Nhacc_formatSpec  = '     Nhacc: %6.4f \n';
    Nhnew_formatSpec  = '     Nhnew: %6.4f \n';
    Nhexit_formatSpec = '    Nhexit: %6.4f \n';
    Etime_formatSpec  = '     Etime: %6.4f \n';
    
    fprintf( '\n\nRSTATUS = \n\n' );
    fprintf( Ntexit_formatSpec, RSTATUS.Ntexit );
    fprintf( Nhacc_formatSpec, RSTATUS.Nhacc );
    fprintf( Nhnew_formatSpec, RSTATUS.Nhnew );
    fprintf( Etime_formatSpec, RSTATUS.Etime );
    
    if ( ~isempty(RSTATUS.Nhexit) )
        fprintf( Nhexit_formatSpec, RSTATUS.Nhexit );
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