function [ H, ISING, e, ISTATUS ] = fatOde_ROS_PrepareMatrix( NVAR, H, Direction, gam, fjac, ISTATUS )
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Filename: fatOde_ROS_PrepareMatrix.m
%
% Original Author:
%
% File Creation Date:
%
% Input Arguments:
%   Name        Type
%   NVAR        integer
%   H           double
%   Direction   integer
%   gam         double
%   fjac        double
%   ISTATUS     struct
%   
% Output Arguments:
%   Name        Type
%   H           double
%   ISING       integer
%   e           double
%   ISTATUS     struct
%
% Modification History:
%   Date        Developer         Email             Action  
%   7/17/2012   Tony D'Augustine  adaug13@vt.edu 
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% fatOde_ROS_PrepareMatrix:
%
% fatOde_PrepareMatrix: INPUT ARGUMENTS
%        NVAR (integer):
%            H (double):
%   Direction (integer):
%          gam (double):
%         fjac (double):
%      ISTATUS (struct):
%
% fatOde_ROS_PrepareMatrix: OUTPUT ARGUMENTS
%         H (double):
%    ISING (integer):
%         e (double):
%   ISTATUS (struct):
%
% fatOde_ROS_PrepareMatrix: SYNTAX
%   [ H ISING e ISTATUS ] = fatOde_PrepareMatrix( NVAR, H, Direction, gam, fjac, ISTATUS );
%
% fatOde_ROS_PrepareMatrix: EXAMPLE
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

