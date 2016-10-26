function stack_ptr = rosPush( NVAR, T, H, Ystage, K, ICNTRL, Coefficient, stack_ptr )
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Filename: rosPush.m
%
% Original Author:
% 
% File Create Date:
%
% Input Arguments:
%   Name    Type
%
% Output Arguments:
%   Name    Type
%
% Modification History:
%   Date    Developer   Email   Action
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% rosPush:
%   Pop Explicit Runge-Kutta stack variables.
%
% rosPush: INPUT ARGUMENTS
%
% rosPush: OUTPUT ARGUMENTS
%
% rosPush: GLOBAL VARIABLES
%
% rosPush: SYNTAX
%
% rosPush: EXAMPLE
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    global chk_H chk_T chk_Y chk_K
    
    stack_ptr = stack_ptr + 1;
    if ( stack_ptr > ICNTRL.Max_no_steps )
        disp('Push failed: buffer overflow');
    end
    
    chk_H(stack_ptr) = H;
    chk_T(stack_ptr) = T;
    chk_Y(1:NVAR*Coefficient.NStage,stack_ptr) = Ystage(1:NVAR*Coefficient.NStage);
    chk_K(1:NVAR*Coefficient.NStage,stack_ptr) = K(1:NVAR*Coefficient.NStage);

return;

