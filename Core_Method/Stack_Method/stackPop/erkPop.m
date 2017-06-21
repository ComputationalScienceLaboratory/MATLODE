function [ T, H, Y, Z, stack_ptr ] = erkPop( NVAR, rkS, stack_ptr )
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Filename: erkPop.m
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
% erkPop:
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
        % for debugging
        %TODO WE MAY WANT TO ADD THIS. The same function was inside
        %ERK_ADJ2_Discretintegrator.
%      if ( OPTIONS.displaySteps == true )
%         str = ['Backtracking. Time = ', num2str(T), '; Stepsize = ', num2str(H)];
%         disp(str);
%      end    
    
    Y(1:NVAR) = chk_Y( 1:NVAR, stack_ptr );
    Z(1:NVAR, 1:rkS) = chk_Z(1:NVAR, 1:rkS, stack_ptr );
    
    stack_ptr = stack_ptr - 1;

return;

