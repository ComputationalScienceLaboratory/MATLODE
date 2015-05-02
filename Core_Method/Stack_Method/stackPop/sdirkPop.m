function [ T, H, Y, Z, stack_ptr ] = sdirkPop( NVAR, rkS, stack_ptr )
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Filename: sdirkkPop.m
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
% sdirkPop:
%   Pop Explicit Runge-Kutta stack variables.
%
% erkPop: INPUT ARGUMENTS
%
% erkPop: OUTPUT ARGUMENTS
%
% erkPop: GLOBAL VARIABLES
%
% erkPop: SYNTAX
%
% erkPop: EXAMPLE
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    global chk_H chk_T chk_Y chk_Z

    if ( stack_ptr <= 0 )
        error( 'Pop failed: empty buffer' );
    end
    H = chk_H( stack_ptr );
    T = chk_T( stack_ptr );
    Y(1:NVAR) = chk_Y( 1:NVAR, stack_ptr );
    Z(1:NVAR, 1:rkS) = chk_Z(1:NVAR, 1:rkS, stack_ptr );
    
    Y = Y'; % temp fix
    
    stack_ptr = stack_ptr - 1;

return;
