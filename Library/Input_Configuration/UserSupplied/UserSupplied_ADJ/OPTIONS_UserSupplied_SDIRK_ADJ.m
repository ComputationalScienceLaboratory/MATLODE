function [ OPTIONS_U ] = OPTIONS_UserSupplied_SDIRK_ADJ( OPTIONS_U )

    if ( ~isempty(OPTIONS_U.AbsTol_TLM) == true )
        warning('MatlODE:configuration','AbsTol_TLM is not used in integrator.');
        OPTIONS_U.AbsTol_TLM = [];
    end
    if ( ~isempty(OPTIONS_U.AdjointSolve) == true )
        warning('MatlODE:configuration','AdjointSolve is not used in integrator.');
        OPTIONS_U.AdjointSolve = [];
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

