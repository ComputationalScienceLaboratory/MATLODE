%% PrintISTATUS
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
function PrintISTATUS( ISTATUS )

    Nfun_formatSpec = '    Nfun: %6.0f \n';
    Njac_formatSpec = '    Njac: %6.0f \n';
    Nstp_formatSpec = '    Nstp: %6.0f \n';
    Nacc_formatSpec = '    Nacc: %6.0f \n';
    Nrej_formatSpec = '    Nrej: %6.0f \n';
    Ndec_formatSpec = '    Ndec: %6.0f \n';
    Nsol_formatSpec = '    Nsol: %6.0f \n';
    Nsng_formatSpec = '    Nsng: %6.0f \n';
%    Nchk_formatSpec = '    Nchk: %6.0f \n';
    
    fprintf( '\n\nISTATUS = \n\n' );
    fprintf( Nfun_formatSpec, ISTATUS.Nfun );
    fprintf( Njac_formatSpec, ISTATUS.Njac );
    fprintf( Nstp_formatSpec, ISTATUS.Nstp );
    fprintf( Nacc_formatSpec, ISTATUS.Nacc );
    fprintf( Nrej_formatSpec, ISTATUS.Nrej );
    fprintf( Ndec_formatSpec, ISTATUS.Ndec );
    fprintf( Nsol_formatSpec, ISTATUS.Nsol );
    fprintf( Nsng_formatSpec, ISTATUS.Nsng );
%    fprintf( Nchk_formatSpec, ISTATUS.Nchk );

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