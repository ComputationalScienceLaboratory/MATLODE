%% OPTIONS_Merge
%
% <html>
%   <div>
%       <img style="float: right" src="../../../MATLODE_LOGO.png" height="150px"></img>
%   </div>
% </html>
%
%% Syntax
%    [Options] = OPTIONS_Merge(Options_U,Options)
%
%% Input Parameters
% |Options_U|: User supplied option struct
%
% |Options|: Fine tuned option struct
%
%% Output Parameters
% |Options|: The merged output option struct user and fine tuned option struct
% paramters.
%
%% Description
% Merges two option structs together.
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
function [ OPTIONS ] = OPTIONS_Merge( OPTIONS_U, OPTIONS )

    if ( ~isempty(OPTIONS_U.AbsTol) == true )
        OPTIONS.AbsTol = OPTIONS_U.AbsTol;
    end
    if ( ~isempty(OPTIONS_U.AbsTol_ADJ) == true )
        OPTIONS.AbsTol_ADJ = OPTIONS_U.AbsTol_ADJ;
    end
    if ( ~isempty(OPTIONS_U.AbsTol_TLM) == true )
        OPTIONS.AbsTol_TLM = OPTIONS_U.AbsTol_TLM;
    end
    if ( ~isempty(OPTIONS_U.Adaptive_Krylov) == true )
        OPTIONS.Adaptive_Krylov = OPTIONS_U.Adaptive_Krylov;
    end
    if ( ~isempty(OPTIONS_U.Autonomous) == true )
        OPTIONS.Autonomous = OPTIONS_U.Autonomous;
    end
    if ( ~isempty(OPTIONS_U.ChunkSize) == true )
        OPTIONS.ChunkSize = OPTIONS_U.ChunkSize;
    end
    if ( ~isempty(OPTIONS_U.Desired_Mode) == true )
        OPTIONS.Desired_Mode = OPTIONS_U.Desired_Mode;
    end
    if ( ~isempty(OPTIONS_U.DirectADJ) == true )
        OPTIONS.DirectADJ = OPTIONS_U.DirectADJ;
    end
    if ( ~isempty(OPTIONS_U.DirectTLM) == true )
        OPTIONS.DirectTLM = OPTIONS_U.DirectTLM;
    end
    if ( ~isempty(OPTIONS_U.displayStats) == true )
        OPTIONS.displayStats = OPTIONS_U.displayStats;
    end
    if ( ~isempty(OPTIONS_U.displayStats) == true )
        OPTIONS.displayStats = OPTIONS_U.displayStats;
    end
    if ( ~isempty(OPTIONS_U.displaySteps) == true )
        OPTIONS.displaySteps = OPTIONS_U.displaySteps;      
    end
    if ( ~isempty(OPTIONS_U.DRDP) == true )
        OPTIONS.DRDP = OPTIONS_U.DRDP;
    end
    if ( ~isempty(OPTIONS_U.DRDY) == true )
        OPTIONS.DRDY = OPTIONS_U.DRDY;
    end
    if ( ~isempty(OPTIONS_U.FacMax) == true )
        OPTIONS.FacMax = OPTIONS_U.FacMax;
    end
    if ( ~isempty(OPTIONS_U.FacMin) == true )
        OPTIONS.FacMin = OPTIONS_U.FacMin;
    end
    if ( ~isempty(OPTIONS_U.FacRej) == true )
        OPTIONS.FacRej = OPTIONS_U.FacRej;
    end
    if ( ~isempty(OPTIONS_U.FacSafe) == true )
        OPTIONS.FacSafe = OPTIONS_U.FacSafe;
    end
    if ( ~isempty(OPTIONS_U.FDAprox) == true )
        OPTIONS.FDAprox = OPTIONS_U.FDAprox;
    end
    if ( ~isempty(OPTIONS_U.FDIncrement) == true )
        OPTIONS.FDIncrement = OPTIONS_U.FDIncrement;
    end
    if ( ~isempty(OPTIONS_U.Gustafsson) == true )
        OPTIONS.Gustafsson = OPTIONS_U.Gustafsson;
    end
    if ( ~isempty(OPTIONS_U.GMRES_TOL) == true )
        OPTIONS.GMRES_TOL = OPTIONS_U.GMRES_TOL;
    end
    
    if ( ~isempty(OPTIONS_U.GMRES_MaxIt) == true )
        OPTIONS.GMRES_MaxIt = OPTIONS_U.GMRES_MaxIt;
    end
    if ( ~isempty(OPTIONS_U.GMRES_Restart) == true )
        OPTIONS.GMRES_Restart = OPTIONS_U.GMRES_Restart;
    end
    if ( ~isempty(OPTIONS_U.GMRES_P) == true )
        OPTIONS.GMRES_P = OPTIONS_U.GMRES_P;
    end
    if ( ~isempty(OPTIONS_U.Hess_vec) == true )
        OPTIONS.Hess_vec = OPTIONS_U.Hess_vec;
    end
    if ( ~isempty(OPTIONS_U.Hesstr_vec) == true )
        OPTIONS.Hesstr_vec = OPTIONS_U.Hesstr_vec;
    end
    if ( ~isempty(OPTIONS_U.Hesstr_vec_f_py) == true )
        OPTIONS.Hesstr_vec_f_py = OPTIONS_U.Hesstr_vec_f_py;
    end
    if ( ~isempty(OPTIONS_U.Hesstr_vec_r) == true )
        OPTIONS.Hesstr_vec_r = OPTIONS_U.Hesstr_vec_r;
    end
    if ( ~isempty(OPTIONS_U.Hesstr_vec_r_py) == true )
        OPTIONS.Hesstr_vec_r_py = OPTIONS_U.Hesstr_vec_r_py;
    end
    if ( ~isempty(OPTIONS_U.Hmax) == true )
        OPTIONS.Hmax = OPTIONS_U.Hmax;
    end
    if ( ~isempty(OPTIONS_U.Hmin) == true )
        OPTIONS.Hmin = OPTIONS_U.Hmin;
    end
    if ( ~isempty(OPTIONS_U.Hstart) == true )
        OPTIONS.Hstart = OPTIONS_U.Hstart;
    end
    if ( ~isempty(OPTIONS_U.Jacp) == true )
        OPTIONS.Jacp = OPTIONS_U.Jacp;
    end
    if ( ~isempty(OPTIONS_U.Lambda) == true )
        OPTIONS.Lambda = OPTIONS_U.Lambda;
    end
    if ( ~isempty(OPTIONS_U.NTLM) == true )
        OPTIONS.NTLM = OPTIONS_U.NTLM;
    end
    if ( ~isempty(OPTIONS_U.TLMNewtonEst) == true )
        OPTIONS.TLMNewtonEst = OPTIONS_U.TLMNewtonEst;
    end
    if ( ~isempty(OPTIONS_U.ITOL) == true )
        OPTIONS.ITOL = OPTIONS_U.ITOL;
    end
    if ( ~isempty(OPTIONS_U.IsJACSymm) == true )
        OPTIONS.IsJACSymm = OPTIONS_U.IsJACSymm;
    end
    if ( ~isempty(OPTIONS_U.Jacobian) == true )
        OPTIONS.Jacobian = OPTIONS_U.Jacobian;
    end
    if ( ~isempty(OPTIONS_U.Max_no_steps) == true )
        OPTIONS.Max_no_steps = OPTIONS_U.Max_no_steps;
    end
    if ( ~isempty(OPTIONS_U.MatrixFree) == true )
        OPTIONS.MatrixFree = OPTIONS_U.MatrixFree;
    end
    if ( ~isempty(OPTIONS_U.Method) == true )
        OPTIONS.Method = OPTIONS_U.Method;
    end
    if ( ~isempty(OPTIONS_U.Mu) == true )
        OPTIONS.Mu = OPTIONS_U.Mu;
    end
    if ( ~isempty(OPTIONS_U.NBasisVectors) == true )
        OPTIONS.NBasisVectors = OPTIONS_U.NBasisVectors;
    end
    if ( ~isempty(OPTIONS_U.NADJ) == true )
        OPTIONS.NADJ = OPTIONS_U.NADJ;
    end
    if ( ~isempty(OPTIONS_U.NewtonMaxIt) == true )
        OPTIONS.NewtonMaxIt = OPTIONS_U.NewtonMaxIt;
    end
    if ( ~isempty(OPTIONS_U.NewtonTol) == true )
        OPTIONS.NewtonTol = OPTIONS_U.NewtonTol;
    end
    if ( ~isempty(OPTIONS_U.NP) == true )
        OPTIONS.NP = OPTIONS_U.NP;
    end    
    if ( ~isempty(OPTIONS_U.QFun) == true )
        OPTIONS.QFun = OPTIONS_U.QFun;
    end
    if ( ~isempty(OPTIONS_U.Qmax) == true )
        OPTIONS.Qmax = OPTIONS_U.Qmax;
    end
    if ( ~isempty(OPTIONS_U.Qmin) == true )
        OPTIONS.Qmin = OPTIONS_U.Qmin;
    end
    if ( ~isempty(OPTIONS_U.Quadrature) == true )
        OPTIONS.Quadrature = OPTIONS_U.Quadrature;
    end
    if ( ~isempty(OPTIONS_U.RelTol) == true )
        OPTIONS.RelTol = OPTIONS_U.RelTol;
    end
    if ( ~isempty(OPTIONS_U.RelTol_ADJ) == true )
        OPTIONS.RelTol_ADJ = OPTIONS_U.RelTol_ADJ;
    end
    if ( ~isempty(OPTIONS_U.RelTol_TLM) == true )
        OPTIONS.RelTol_TLM = OPTIONS_U.RelTol_TLM;
    end
    if ( ~isempty(OPTIONS_U.SaveLU) == true )
        OPTIONS.SaveLU = OPTIONS_U.SaveLU;
    end
    if ( ~isempty(OPTIONS_U.SdirkError) == true )
        OPTIONS.SdirkError = OPTIONS_U.SdirkError;
    end
    if ( ~isempty(OPTIONS_U.storeCheckpoint) == true )
        OPTIONS.storeCheckpoint = OPTIONS_U.storeCheckpoint;
    end
    if ( ~isempty(OPTIONS_U.StartNewton) == true )
        OPTIONS.StartNewton = OPTIONS_U.StartNewton;
    end
    if ( ~isempty(OPTIONS_U.ThetaMin) == true )
        OPTIONS.ThetaMin = OPTIONS_U.ThetaMin;
    end
    if ( ~isempty(OPTIONS_U.TLMTruncErr) == true )
        OPTIONS.TLMTruncErr = OPTIONS_U.TLMTruncErr;
    end
    if ( ~isempty(OPTIONS_U.Y_TLM) == true )
        OPTIONS.Y_TLM = OPTIONS_U.Y_TLM;
    end
    if ( ~isempty(OPTIONS_U.LU) == true )
        OPTIONS.LU = OPTIONS_U.LU;
    end

end

%%
% <html>
%   <div>
%       <img style="float: right" src="../../../CSL_LogoWithName_1.png" height="50px"></img>
%   </div>
% </html>