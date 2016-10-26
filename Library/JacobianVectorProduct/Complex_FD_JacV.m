%% Complex_FD_JacV
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
function [ JacV, Stats ] = Complex_FD_JacV( T, Y, V, OdeFunction, Stats )
    yDimension = length(Y);

    % Complex finite difference
    JacV = zeros(yDimension,1);
    canonical = eye(yDimension);
    %epsilon = (1/(yDimension*norm(V)))*sum(sqrt(eps).*Y+sqrt(eps));    
    epsilon = sqrt((1+norm(Y))*eps)/norm(V);
    for i=1:yDimension
        JacV = JacV + ((imag(OdeFunction(T,Y+1i*canonical(:,i)*epsilon))/epsilon))*V(i);
    end
    Stats.Nfun = Stats.Nfun + yDimension;
    
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
