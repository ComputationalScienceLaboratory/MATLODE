function [ OPTIONS_U ] = OPTIONS_UserSupplied_RK_FWD( OPTIONS_U )

    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %   On/Off Warnings
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
    if ( OPTIONS_U.WarningConfig == false )
        warning('off','MatlODE:configuration');
    else 
        warning('on','MatlODE:configuration');
    end

    if ( ~isempty(OPTIONS_U.AbsTol_ADJ) == true )
        warning('MatlODE:configuration','AbsTol_ADJ is not used in integrator.');
        OPTIONS_U.AbsTol_ADJ = [];
    end
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
    if ( ~isempty(OPTIONS_U.DirectADJ) == true )
        warning('MatlODE:configuration','DirectADJ is not used in integrator.');
        OPTIONS_U.DirectADJ = [];
    end
    if ( ~isempty(OPTIONS_U.DirectTLM) == true )
        warning('MatlODE:configuration','DirectTLM is not used in integrator.');
        OPTIONS_U.DirectTLM = [];
    end
    if ( ~isempty(OPTIONS_U.DRDP) == true )
        warning('MatlODE:configuration','DRDP is not used in integrator.');
        OPTIONS_U.DRDP = [];
    end
    if ( ~isempty(OPTIONS_U.DRDY) == true )
        warning('MatlODE:configuration','DRDY is not used in integrator.');
        OPTIONS_U.DRDY = [];
    end
    if ( ~isempty(OPTIONS_U.FDAprox) == true )
        warning('MatlODE:configuration','FDAprox is not used in integrator.');
        OPTIONS_U.FDAprox = [];
    end
    if ( ~isempty(OPTIONS_U.FDIncrement) == true )
        warning('MatlODE:configuration','FDIncrement is not used in integrator.');
        OPTIONS_U.FDIncrement = [];
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
    if ( isempty(OPTIONS_U.Jacobian) == true )
        error('Jacobian is required.');
    end
    if ( ~isempty(OPTIONS_U.TLMNewtonEst) == true )
        warning('MatlODE:configuration','TLMNewtonEst is not used in integrator.');
        OPTIONS_U.TLMNewtonEst = [];
    end
    if ( ~isempty(OPTIONS_U.TLMTruncErr) == true )
        warning('MatlODE:configuration','TLMTruncErr is not used in integrator.');
        OPTIONS_U.TLMTruncErr = [];
    end
    if ( ~isempty(OPTIONS_U.MatrixFree) == true )
        warning('MatlODE:configuration','MatrixFree is not used in integrator.');
        OPTIONS_U.MatrixFree = [];
    end
    if ( ~isempty(OPTIONS_U.Mu) == true )
        warning('MatlODE:configuration','Mu is not used in integrator.');
        OPTIONS_U.Mu = [];
    end
    if ( ~isempty(OPTIONS_U.QFun) == true )
        warning('MatlODE:configuration','QFun is not used in integrator.');
        OPTIONS_U.QFun = [];
    end
    if ( ~isempty(OPTIONS_U.Quadrature) == true )
        warning('MatlODE:configuration','Quadrature is not used in integrator.');
        OPTIONS_U.Quadrature = [];
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
    if ( ~isempty(OPTIONS_U.TLMTruncErr) == true )
        warning('MatlODE:configuration','TLMTruncErr is not used in integrator.');
        OPTIONS_U.TLMTruncErr = [];
    end
    if ( ~isempty(OPTIONS_U.Y_TLM) == true )
        warning('MatlODE:configuration','Y_TLM is not used in integrator.');
        OPTIONS_U.Y_TLM = [];
    end


end

