function [ OPTIONS_U ] = OPTIONS_UserSupplied_ERK_ADJ( OPTIONS_U )

    if ( ~isempty(OPTIONS_U.AbsTol_ADJ) == true )
        warning('MatlODE:configuration','AbsTol_ADJ is not used in integrator.');
        OPTIONS_U.AbsTol_ADJ = [];
    end
    if ( ~isempty(OPTIONS_U.AbsTol_TLM) == true )
        warning('MatlODE:configuration','Abs_TLM is not used in integrator.');
        OPTIONS_U.AbsTol_TLM = [];
    end                    
    if ( ~isempty(OPTIONS_U.AdjointSolve) == true )
        warning('MatlODE:configuration','AdjointSolve is not used in integrator.');
        OPTIONS_U.AdjointSolve = [];
    end
    if ( ~isempty(OPTIONS_U.FDAprox) == true )
        warning('MatlODE:configuration','FDAprox is not used in integrator.');
        OPTIONS_U.FDAprox = [];
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
    if ( ~isempty(OPTIONS_U.MatrixFree) == true )
        warning('MatlODE:configuration','MatrixFree is not used in integrator.');
        OPTIONS_U.MatrixFree = [];
    end                    
    if ( ~isempty(OPTIONS_U.NewtonMaxIt) == true )
        warning('MatlODE:configuration','NewtonMaxIt is not used in integrator.');
        OPTIONS_U.NewtonMaxIt = [];
    end
    if ( ~isempty(OPTIONS_U.NewtonTol) == true )
        warning('MatlODE:configuration','NewtonTol is not used in integrator.');
        OPTIONS_U.NewtonTol = [];
    end
    if ( ~isempty(OPTIONS_U.Qmax) == true )
        warning('MatlODE:configuration','Qmax is not used in integrator.');
        OPTIONS_U.Qmax = [];
    end
    if ( ~isempty(OPTIONS_U.Qmin) == true )
        warning('MatlODE:configuration','Qmin is not used in integrator.');
        OPTIONS_U.Qmin = [];
    end
    if ( ~isempty(OPTIONS_U.RelTol_ADJ) == true )
        warning('MatlODE:configuration','RelTol_ADJ is not used in integrator.');
        OPTIONS_U.RelTol_ADJ = [];
    end
    if ( ~isempty(OPTIONS_U.RelTol_TLM) == true )
        warning('MatlODE:configuration','RelTol_TLM is not used in integrator.');
        OPTIONS_U.RelTol_TLM = [];
    end
    if ( ~isempty(OPTIONS_U.SaveLU) == true )
        warning('MatlODE:configuration','SaveLU is not used in integrator.');
        OPTIONS_U.SaveLU = [];
    end
    if ( ~isempty(OPTIONS_U.SdirkError) == true )
        warning('MatlODE:configuration','SdirkError is not used in integrator.');
        OPTIONS_U.SdirkError = [];
    end
    if ( ~isempty(OPTIONS_U.StartNewton) == true )
        warning('MatlODE:configuration','StartNewtion is not used in integrator.');
        OPTIONS_U.StartNewton = [];
    end
    if ( ~isempty(OPTIONS_U.ThetaMin) == true )
        warning('MatlODE:configuration','ThetaMin is not used in integrator.');
        OPTIONS_U.ThetaMin = [];
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

