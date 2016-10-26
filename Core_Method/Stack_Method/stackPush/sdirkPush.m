function stack_ptr = sdirkPush( NVAR, rkS, T, H, Y, Z, Max_no_steps, stack_ptr )
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Filename: sdirkPush.m
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
% sdirkPush:
%   Pop Explicit Runge-Kutta stack variables.
%
% sdirkPush: INPUT ARGUMENTS
%
% sdirkPush: OUTPUT ARGUMENTS
%
% sdirkPush: GLOBAL VARIABLES
%
% sdirkPush: SYNTAX
%
% sdirkPush: EXAMPLE
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    global chk_H chk_T chk_Y chk_Z

    stack_ptr = stack_ptr + 1;
    if ( stack_ptr > Max_no_steps )
        error( 'Push failed: buffer overflow' );
    end
    chk_H( stack_ptr ) = H;
    chk_T( stack_ptr ) = T;
    chk_Y( 1:NVAR, stack_ptr ) = Y(1:NVAR);
    chk_Z( 1:NVAR, 1:rkS, stack_ptr ) = Z(1:NVAR, 1:rkS);

return;

