% Testing Suite: General

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   General Setup
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% analysis: -1: solution
%           -2: error
analysis = -1;

% integrator: -1: ERK_FWD_DRIVER_Integrator
%             -2: ROS_FWD_DRIVER_Integrator
%             -3: RK_FWD_DRIVER_Integrator
%             -4: SDIRK_FWD_DRIVER_Integrator
%
%             -5: ERK_TLM_DRIVER_Integrator
%             -6: ROS_TLM_DRIVER_Integrator
%             -7: RK_TLM_DRIVER_Integrator
%             -8: SDIRK_TLM_DRIVER_Integrator
%
%              -9: ERK_ADJ_DRIVER_Integrator
%             -10: ROS_ADJ_DRIVER_Integrator
%             -11: RK_ADJ_DRIVER_Integrator
%             -12: SDIRK_ADJ_DRIVER_Integrator              
integrator = -10;

% General Requirements:
Absolute_Tolerance = PLACE_HOLDER;
Relative_Tolerance = PLACE_HOLDER;

% General Options:

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   Forward Setup
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   Tangent Linear Setup
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Tangent Linear Requirements
Absolute_Tolerance_TLM = PLACE_HOLDER;
Relative_Tolerance_TLM = PLACE_HOLDER;



%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   Adjoint Setup
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Adjoint Requirements:
Absolute_Tolerance_ADJ = PLACE_HOLDER;
Relative_Tolerance_ADJ = PLACE_HOLDER;

% ADJ1: cost function sensitivities w.r.t. initial conditions
% ADJ2: cost function senstiivities w.r.t. initial conditions and parameters
%
%                 | ADJ1 w/ Quad | ADJ1 no Quad | ADJ2 w/ Quad | ADJ2 no Quad
% -----------------------------------------------------------------------------
% Jacp            |   REQUIRED   |   REQUIRED   |              |
% Mu              |   REQUIRED   |   REQUIRED   |              |
% QFun            |   REQUIRED   |              |   REQUIRED   |
% Quadrature      |   REQUIRED   |              |   REQUIRED   |
% DRDY            |   REQUIRED   |              |   REQUIRED   |
% DRDP            |   REQUIRED   |              |              |
%
% Additional Adjoint Requirements: (Rosenbrock Only)
%                 | ADJ1 w/ Quad | ADJ1 no Quad | ADJ2 w/ Quad | ADJ2 no Quad
% -----------------------------------------------------------------------------
% Hesstr_vec      |   REQUIRED   |   REQUIRED   |   REQUIRED   |   REQUIRED
% Hesstr_vec_f_py |   REQUIRED   |   REQUIRED   |              |
% Hesstr_vec_r    |   REQUIRED   |              |   REQUIRED   |
% Hesstr_vec_r_py |   REQUIRED   |              |              |

Jacp       = PLACE_HOLDER;
Mu         = PLACE_HOLDER;
QFun       = PLACE_HOLDER;
Quadrature = PLACE_HOLDER;
DRDY       = PLACE_HOLDER;
DRDP       = PLACE_HOLDER;

Hesstr_vec      = PLACE_HOLDER;
Hesstr_vec_f_py = PLACE_HOLDER;
Hesstr_vec_r    = PLACE_HOLDER;
Hesstr_vec_r_py = PLACE_HOLDER;



