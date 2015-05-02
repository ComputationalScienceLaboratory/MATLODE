function [ dFdT ISTATUS ] = fatOde_FunctionTimeDerivative( T, Roundoff, Y, Fcn0, OdeFunction, ISTATUS )
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Filename: fatOde_FunctionTimeDerivative.m
%
% Original Author: 
%
% File Creation Date: 
%
% Input Arguments:
%   Name        Type
%   T           double
%   Roundoff    double
%   Y           double
%   Fcn0        double
%   ISTATUS     struct
%
% Output Arguments:
%   Name        Type
%   dFdT        double
%   ISTATUS     struct
% 
% Modification History:
%   Date        Developer         Email             Action  
%   7/17/2012   Tony D'Augustine  adaug13@vt.edu    Translated Fortran90
%                                                   ros_FunTimeDerivative() to MATLAB.
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% fatOde_FunctionTimeDerivative:
%   Calculate the time partial derivative of the function by finite
%   differences.
%
% fatOde_FunctionTimeDerivative: INPUT ARGUMENTS
%          T (double):
%   Roundoff (double):
%          Y (double):
%       Fcn0 (double):
%    ISTATUS (struct):
%
% fatOde_FunctionTimeDerivative: OUTPUT ARGUMENTS
%      dFdT (double):
%   ISTATUS (struct):
%
% fatOde_FunctionTimeDerivative: SYNTAX
%   [ dFdT ISTATUS ] = fatOde_FunctionTimeDerivative( T, Roundoff, Y, Fcn0, ISTATUS )
%
% fatOde_FunctionTimeDerivative: EXAMPLE
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    % local variables
    DeltaMin = 1.0d-6;

    Delta = sqrt(Roundoff)*max(DeltaMin,abs(T));
    dFdT = OdeFunction( T+Delta, Y );
    ISTATUS.Nfun = ISTATUS.Nfun + 1;
    dFdT = dFdT - Fcn0;
    dFdT = dFdT/Delta;

return;