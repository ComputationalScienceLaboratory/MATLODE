function [ T, H, Y, Zstage, NewIt, stack_ptr ] = rkPop( NVAR, stack_ptr )
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Filename: rkPop.m
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
% rkPop:
%   Pop Explicit Runge-Kutta stack variables.
%
% rkPop: INPUT ARGUMENTS
%
% rkPop: OUTPUT ARGUMENTS
%
% rkPop: GLOBAL VARIABLES
%
% rkPop: SYNTAX
%
% rkPop: EXAMPLE
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    global chk_H chk_T chk_Y chk_Z chk_NiT

    if ( stack_ptr <= 0 )
        error( 'Pop failed: empty buffer' );
    end
    
    H = chk_H( stack_ptr );
    T = chk_T( stack_ptr );
    Y = chk_Y( 1:NVAR,stack_ptr );
    Zstage = chk_Z(1:3*NVAR,stack_ptr);
    NewIt = chk_NiT( stack_ptr );
    
    stack_ptr = stack_ptr - 1;

return;

