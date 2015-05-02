function [ OPTIONS_U ] = OPTIONS_UserSupplied_ERK_TLM( OPTIONS_U )

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
    if ( ~isempty(OPTIONS_U.AdjointSolve) == true )
        warning('MatlODE:configuration','AbjointSolve is not used in integrator.');
        OPTIONS_U.AdjointSolve = [];
    end
    if ( ~isempty(OPTIONS_U.Autonomous) == true ) 
        warning('MatlODE:configuration','Autonomous is not used in integrator.');
        OPTIONS_U.Autonomous = [];
    end
    if ( ~isempty(OPTIONS_U.Desired_Mode) == true )
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
        warning('MatlODE:configuration','Hesstr_vec_py is not used in integrator.');
        OPTIONS_U.Hesstr_vec_py = [];
    end
    if ( ~isempty(OPTIONS_U.Jacp) == true )
        warning('MatlODE:configuration','Jacp is not used in integrator.');
        OPTIONS_U.Jacp = [];
    end
    if ( ~isempty(OPTIONS_U.Lambda) == true )
        warning('MatlODE:configuration','Lambda is not used in integrator.');
        OPTIONS_U.Lambda = [];
    end
    if ( ~isempty(OPTIONS_U.MatrixFree) == true )
        warning('MatlODE:configuration','MatrixFree is not used in integrator.');
        OPTIONS_U.MatrixFree = [];
    end                    
    if ( ~isempty(OPTIONS_U.Mu) == true )
        warning('MatlODE:configuration','Mu is not used in integrator.');
        OPTIONS_U.Mu = [];
    end
    if ( ~isempty(OPTIONS_U.NewtonMaxIt) == true )
        warning('MatlODE:configuration','NewtonMaxIt not used in integrator.');
        OPTIONS_U.NewtonMaxIt = [];
    end
    if ( ~isempty(OPTIONS_U.NewtonTol) == true )
        warning('MatlODE:configuration','NewtonTol is not used in integrator.');
        OPTIONS_U.NewtonTol = [];
    end
    if ( ~isempty(OPTIONS_U.QFun) == true )
        warning('MatlODE:configuration','QFun is not used in integrator.');
        OPTIONS_U.QFun = [];
    end
    if ( ~isempty(OPTIONS_U.Qmax) == true )
        warning('MatlODE:configuration','Qmax is not used in integrator.');
        OPTIONS_U.Qmax = [];
    end
    if ( ~isempty(OPTIONS_U.Qmin) == true )
        warning('MatlODE:configuration','Qmin is not used in integrator.');
        OPTIONS_U.Qmin = [];
    end
    if ( ~isempty(OPTIONS_U.Quadrature) == true )
        warning('MatlODE:configuration','Quadrature is not used in integrator.');
        OPTIONS_U.Quadrature = [];
    end
    if ( ~isempty(OPTIONS_U.RelTol_ADJ) == true )
        warning('MatlODE:configuration','RelTol_ADJ is not used in integrator.');
        OPTIONS_U.RelTol_ADJ = [];
    end
    if ( ~isempty(OPTIONS_U.SaveLU) == true )
        warning('MatlODE:configuration','SaveLU is not used in integrator.');
        OPTIONS_U.SaveLU = [];
    end
    if ( ~isempty(OPTIONS_U.SdirkError) == true )
        warning('MatlODE:configuration','SdirkError is not used in integrator.');
        OPTIONS_U.SdirkErr = [];
    end
    if ( ~isempty(OPTIONS_U.StartNewton) == true )
        warning('MatlODE:configuration','StartNewton is not used in integrator.');
        OPTIONS_U.StartNewton = [];
    end
    if ( ~isempty(OPTIONS_U.ThetaMin) == true )
        warning('MatlODE:configuration','ThetaMin is not used in integrator.');
        OPTIONS_U.ThetaMin = [];
    end
    if ( ~isempty(OPTIONS_U.TLMNewtonEst) == true )
        warning('MatlODE:configuration','TLMNewtonEst is not used in integrator.');
        OPTIONS_U.TLMNewtonEst = [];
    end


end

