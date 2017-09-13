%% ISTATUS_Add
%
% <html>
%   <div>
%       <img style="float: right" src="../../../MATLODE_LOGO.png" height="150px"></img>
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
function [ ISTATUS ] = ISTATUS_Add( ISTATUS1, ISTATUS2 )

    ISTATUS.Nfun = ISTATUS1.Nfun + ISTATUS2.Nfun;
    ISTATUS.Njac = ISTATUS1.Njac + ISTATUS2.Njac;
    ISTATUS.Nstp = ISTATUS1.Nstp + ISTATUS2.Nstp;
    ISTATUS.Nacc = ISTATUS1.Nacc + ISTATUS2.Nacc;
    ISTATUS.Nrej = ISTATUS1.Nrej + ISTATUS2.Nrej;
    ISTATUS.Ndec = ISTATUS1.Ndec + ISTATUS2.Ndec;
    ISTATUS.Nsol = ISTATUS1.Nsol + ISTATUS2.Nsol;
    ISTATUS.Nsng = ISTATUS1.Nsng + ISTATUS2.Nsng;
    ISTATUS.Nchk = ISTATUS1.Nchk + ISTATUS2.Nchk;
    ISTATUS.Nkdim = ISTATUS1.Nkdim + ISTATUS2.Nkdim;
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
%       <img style="float: right" src="../../../CSL_LogoWithName_1.png" height="50px"></img>
%   </div>
% </html>
