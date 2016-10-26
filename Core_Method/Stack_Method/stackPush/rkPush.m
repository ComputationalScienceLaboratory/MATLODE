function [ stack_ptr ] = rkPush( NVAR, T, H, Y, Zstage, NewIt, ICNTRL, stack_ptr )
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Filename: rkPush.m
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
% rkPush:
%   Pop Explicit Runge-Kutta stack variables.
%
% rkPush: INPUT ARGUMENTS
%
% rkPush: OUTPUT ARGUMENTS
%
% rkPush: GLOBAL VARIABLES
%
% rkPush: SYNTAX
%
% rkPush: EXAMPLE
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    global chk_H chk_T chk_Y chk_Z chk_NiT
    
    stack_ptr = stack_ptr + 1;

    if ( stack_ptr > ICNTRL.Max_no_steps )
        error('Push failed: buffer overflow');
    end
    
    chk_H( stack_ptr ) = H;
    chk_T( stack_ptr ) = T;
    chk_Y( 1:NVAR, stack_ptr ) = Y;
    chk_Z(1:3*NVAR,stack_ptr ) = Zstage(1:3*NVAR);
    chk_NiT( stack_ptr ) = NewIt;
    
return;

