%% OPTIONS_UserSupplied_RK_ADJ
%
% <html>
%   <div>
%       <img style="float: right" src="../../../../../MATLODE_LOGO.png" height="150px"></img>
%   </div>
% </html>
%
%% Syntax
%    [Options] = OPTIONS_UserSupplied_RK_ADJ(Options_U)
%
%% Input Parameters
% |Options_U|: User supplied option struct
%
%% Output Parameters
% |Options|: Option struct with only necessary option parameters
%
%% Description
% Warns user if option parameter is not used.
%
%% Reference
% [1] Tony D'Augustine, Adrian Sandu. MATLODE: A MATLAB ODE Solver and
%     Sensitivity Analysis Toolbox. Submitted to ACM TOMS.
%
%%
%  Authored by Tony D'Augustine, Adrian Sandu, and Hong Zhang.
%  Computational Science Laboratory, Virginia Tech.
%  ©2015 Virginia Tech Intellectual Properties, Inc.
%
function [ OPTIONS_U ] = OPTIONS_UserSupplied_RK_ADJ( OPTIONS_U )

    if ( ~isempty(OPTIONS_U.AbsTol_TLM) == true )
        warning('MatlODE:configuration','AbsTol_TLM is not used in integrator.');
        OPTIONS_U.AbsTol_TLM = [];
    end
    if ( ~isempty(OPTIONS_U.Autonomous) == true )
        warning('MatlODE:configuration','Autonomous is not used in integrator.');
        OPTIONS_U.Autonomous = [];
    end
    if ( ~isempty(OPTIONS_U.DirectTLM) == true )
        warning('MatlODE:configuration','DirectTLM is not used in integrator.');
        OPTIONS_U.DirectTLM = [];
    end
    if ( ~isempty(OPTIONS_U.Hess_vec) == true )
        warning('MatlODE:configuration','Hess_vec is not used in integrator.');
        OPTIONS_U.Hess_vec = [];
    end
    if ( ~isempty(OPTIONS_U.Hesstr_vec) == true )
        warning('MatlODE:configuration','Hesstr_vec is not used in integrator.');
        OPTIONS_U.Hesstr_vec = [];
    end
    if ( ~isempty(OPTIONS_U.Hesstr_vec_f_py) == true );
        warning('MatlODE:configuration','Hesstr_vec_f_py is not used in integrator.');
        OPTIONS_U.Hesstr_vec_f_py = [];
    end
    if ( ~isempty(OPTIONS_U.Hesstr_vec_r) == true )
        warning('MatlODE:configuration','Hesstr_vec_r is not used in integrator.');
        OPTIONS_U.Hesstr_vec_r = [];
    end
    if ( ~isempty(OPTIONS_U.Hesstr_vec_r_py) == true )
        warning('MatlODE:configuration','Hesstr_vec_r_py is not used in integrator.');
        OPTIONS_U.Hesstr_vec_r_py = [];
    end
    if ( ~isempty(OPTIONS_U.MatrixFree) == true )
        warning('MatlODE:configuration','MatrixFree is not used in integrator.');
        OPTIONS_U.MatrixFree = [];
    end                    
    if ( ~isempty(OPTIONS_U.NewtonTol) == true )
        warning('MatlODE:configuration','NewtonTol is not used in integrator.');
        OPTIONS_U.NewtonTol = [];
    end
    if ( ~isempty(OPTIONS_U.Qmax) == true )
        warning('MatlODE:configuration','Qmax is not used in integrator.');
        OPTIONS_U.Qmax = [];
    end


end

%%
% <html>
%   <div>
%       <img style="float: right" src="../../../../../CSL_LogoWithName_1.png" height="50px"></img>
%   </div>
% </html>