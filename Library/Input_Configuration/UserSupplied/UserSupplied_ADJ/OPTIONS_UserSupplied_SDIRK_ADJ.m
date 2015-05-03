%% OPTIONS_UserSupplied_SDIRK_ADJ
%
% <html>
%   <div>
%       <img style="float: right" src="../../../../../MATLODE_LOGO.png" height="150px"></img>
%   </div>
% </html>
%
%% Syntax
%    [Options] = OPTIONS_UserSupplied_SDIRK_ADJ(Options_U)
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
% Copyright 2015 Computational Science Laboratory
function [ OPTIONS_U ] = OPTIONS_UserSupplied_SDIRK_ADJ( OPTIONS_U )

    if ( ~isempty(OPTIONS_U.AbsTol_TLM) == true )
        warning('MatlODE:configuration','AbsTol_TLM is not used in integrator.');
        OPTIONS_U.AbsTol_TLM = [];
    end
    if ( ~isempty(OPTIONS_U.Desired_Mode) == true );
        warning('MatlODE:configuration','Desired_Mode is not used in integrator.');
        OPTIONS_U.Desired_Mode = [];
    end
    if ( ~isempty(OPTIONS_U.DirectTLM) == true )
        warning('MatlODE:configuration','DirectTLM is not used in integrator.');
        OPTIONS_U.DirectTLM = [];
    end
    if ( ~isempty(OPTIONS_U.FDAprox) == true )
        warning('MatlODE:configuration','FDAprox is not used in integrator.');
        OPTIONS_U.FDAprox = [];
    end
    if ( ~isempty(OPTIONS_U.FDIncrement) == true )
        warning('MatlODE:configuration','FDIncrement is not used in integrator.');
        OPTIONS_U.FDIncrement = [];
    end                    
    if ( ~isempty(OPTIONS_U.Gustafsson) == true )
        warning('MatlODE:configuration','Gustafsson is not used in integrator.');
        OPTIONS_U.Gustafsson = [];
    end
    if ( ~isempty(OPTIONS_U.Hess_vec) == true )
        warning('MatlODE:configuration','Hess_vec is not used in integrator.');
        OPTIONS_U.Hess_vec = [];
    end
    if ( ~isempty(OPTIONS_U.Hesstr_vec) == true )
        warning('MatlODE:configuration','Hesstr_vec is not used in integrator.');
        OPTIONS_U.Hesstr_vec = [];
    end
    if ( ~isempty(OPTIONS_U.Hesstr_vec_f_py) == true )
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
    if ( ~isempty(OPTIONS_U.TLMNewtonEst) == true )
        warning('MatlODE:configuration','TLMNewtonEst is not used in integrator.');
        OPTIONS_U.TLMNewtonEst = [];
    end
    if ( ~isempty(OPTIONS_U.TLMTruncErr) == true )
        warning('MatlODE:configuration','TLMTruncErr is not used in integrator.');
        OPTIONS_U.TLMTruncErr = [];
    end
    if ( ~isempty(OPTIONS_U.Qmax) == true )
        warning('MatlODE:configuration','Qmax is not used in integrator.');
        OPTIONS_U.Qmax = [];
    end
    if ( ~isempty(OPTIONS_U.Qmin) == true )
        warning('MatlODE:configuration','Qmin is not used in integrator.');
        OPTIONS_U.Qmin = [];
    end
    if ( ~isempty(OPTIONS_U.RelTol_TLM) == true )
        warning('MatlODE:configuration','RelTol_TLM is not used in integrator.');
        OPTIONS_U.RelTol_TLM = [];
    end
    if ( ~isempty(OPTIONS_U.SdirkError) == true )
        warning('MatlODE:configuration','SdirkError is not used in integrator.');
        OPTIONS_U.SdirkError = [];
    end
    if ( ~isempty(OPTIONS_U.TLMTruncErr) == true )
        warning('MatlODE:configuration','TLMTruncErr is not used in integrator.');
        OPTIONS_U.TLMTruncErr = [];
    end
    if ( ~isempty(OPTIONS_U.Y_TLM) == true )
        warning('MatlODE:configuration','Y_TLM is not used in integrator.');
        OPTIONS_U.Y_TLM = [];
    end       

end

%%
% <html>
%   <div>
%       <img style="float: right" src="../../../../../CSL_LogoWithName_1.png" height="50px"></img>
%   </div>
% </html>