function [ OPTIONS_U ] = OPTIONS_UserSupplied_RK_ADJ( OPTIONS_U )

    if ( ~isempty(OPTIONS_U.AbsTol_TLM) == true )
        warning('MatlODE:configuration','AbsTol_TLM is not used in integrator.');
        OPTIONS_U.AbsTol_TLM = [];
    end
    if ( ~isempty(OPTIONS_U.AdjointSolve) == true )
        warning('MatlODE:configuration','AdjointSolve is not used in integrator.');
        OPTIONS_U.AdjointSolve = [];
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

