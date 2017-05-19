function [ T, H, Y, K, stack_ptr ] = rosPop( NVAR, rkS, stack_ptr )
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Filename: rosPop.m
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
% rosPop:
%   Pop Explicit Runge-Kutta stack variables.
%
% rosPop: INPUT ARGUMENTS
%
% rosPop: OUTPUT ARGUMENTS
%
% rosPop: GLOBAL VARIABLES
%
% rosPop: SYNTAX
%
% rosPop: EXAMPLE
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    global chk_H chk_T chk_Y chk_K

    if ( stack_ptr <= 0 )
        disp('Pop failed: empty buffer');
    end
    
    H = chk_H(stack_ptr);
    T = chk_T(stack_ptr);
    Y = chk_Y(1:NVAR*rkS,stack_ptr);
    K = chk_K(1:NVAR*rkS,stack_ptr);
    
    stack_ptr = stack_ptr - 1;

return;

