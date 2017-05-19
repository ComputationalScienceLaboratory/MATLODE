%% fatOde_ROS_PrepareMatrix
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
function [ H, ISING, e, ISTATUS ] = fatOde_ROS_PrepareMatrix( NVAR, H, Direction, gam, fjac, ISTATUS )

    Nconsecutive = 0;
    ISING = true;

    while ( ISING )
        % Construct Ghimj = 1/(H*ham) - Jac0
        ghinv = 1.0/(Direction*H*gam);
        
        % Compute LU decomposition
        [ ISING, e ] = lss_decomp( NVAR, ghinv, fjac );
        ISTATUS.Ndec = ISTATUS.Ndec + 1;
        
        if ( ISING == 0 ) % If successful, done.
            ISING = false;
        else % If unsuccessful half the step size; If 5 consecutive fails then return.
            ISTATUS.Nsng = ISTATUS.Nsng + 1;
            Nconsecutive = Nconsecutive + 1;
            ISING = true;
            str = [ 'Warning: LU Decomposition returned ISING = ', num2str(ISING) ];
            disp(str);
            if ( Nconsecutive <= 5 ) % Less than 5 consecutive failed decompositions
                H = H*0.5;
            else
                % More than 5 consecutive failed decompositions
                return;
            end % Nconsecutive
        end % ising
    end % while singular

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
