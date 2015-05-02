function SCAL = fatOde_ErrorScale( NVAR, ITOL, AbsTol, RelTol, Y0 )
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Filename: fatOde_ErrorScale.m
%
% Original Author: 
%
% File Creation Date: 
%
% Input Arguments:
%   Name    Type
%   NVAR    integer
%   ITOL    
%   AbsTol  double
%   RelTol  double
%   Y0      double
%
% Output Arguments:
%   Name    Type
%   SCAL    double
% 
% Modification History:
%   Date        Developer         Email             Action  
%   7/17/2012   Tony D'Augustine  adaug13@vt.edu    Translated Fortran90
%                                                   ErrorScale() to MATLAB.
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% fatOde_ErrorNorm:
%   Calculate the error norm.
%
% fatOde_ErrorNorm: INPUT ARGUMENTS
%    NVAR (integer):
%          ITOL ( ):
%   AbsTol (double):
%   RelTol (double):
%       Y0 (double):
%   
% fatOde_ErrorNorm: OUTPUT ARGUMENT
%   SCAL (double):
%
% fatOde_ErrorNorm: SYNTAX
%   [ SCAL ] = fatOde_ErrorScale( NVAR, ITOL, AbsTol, RelTol, Y0 )
%
% fatOde_ErrorNorm: EXAMPLE
%   SCAL = fatOde)ErrorScale( NVAR, OPTIONS.ITOL, OPTIONS.AbsTol, OPTIONS.RelTol, Y0 );
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SCAL = zeros(NVAR,1);

    if ( ITOL == 0 )
        for i = 1:NVAR
            SCAL(i) = 1.0 / ( AbsTol(1) + RelTol(1) * abs( Y0(i) ) );
        end
    else
        for i = 1:NVAR
            SCAL(i) = 1.0 / ( AbsTol(i) + RelTol(i) * abs( Y0(i) ) );
        end
    end
%     SCAL = 1.0 ./ ( AbsTol + RelTol .* abs(Y0) );

return; 

