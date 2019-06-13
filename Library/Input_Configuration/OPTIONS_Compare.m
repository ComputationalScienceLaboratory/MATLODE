%% OPTIONS_Compare
%
% <html>
%   <div>
%       <img style="float: right" src="../../../MATLODE_LOGO.png" height="150px"></img>
%   </div>
% </html>
%
%% Syntax
%    OPTIONS_Compare(Options1,Options2)
%
%% Input Parameters
% |Options1|: Option struct to be compared to Options2
%
% |Options2|: Option struct to be compared to Options1
%
%% Output Parameters
% No output parameters. Function prints to command window.
%
%% Description
% Compared two MATLODE option structs and prints results to command window.
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
function OPTIONS_Compare( option1, option2 )

    % AbsTol
    if ( isequal(option1.AbsTol,option2.AbsTol) )
        fprintf(  ['         AbsTol : ',mat2str(option1.AbsTol),' | ',mat2str(option2.AbsTol),'\n']);
    else
        fprintf(2,['         AbsTol : ',mat2str(option1.AbsTol),' | ',mat2str(option2.AbsTol),'\n']);
    end
    
    % AbsTol_ADJ
    if ( isequal(option1.AbsTol_ADJ,option2.AbsTol_ADJ) )
        fprintf(  ['     AbsTol_ADJ : ',mat2str(option1.AbsTol_ADJ),' | ',mat2str(option2.AbsTol_ADJ),'\n']);
    else
        fprintf(2,['     AbsTol_ADJ : ',mat2str(option1.AbsTol_ADJ),' | ',mat2str(option2.AbsTol_ADJ),'\n']);
    end    
    
    % AbsTol_TLM
    if ( isequal(option1.AbsTol_TLM,option2.AbsTol_TLM) )
        fprintf(  ['     AbsTol_TLM : ',mat2str(option1.AbsTol_TLM),' | ',mat2str(option2.AbsTol_TLM),'\n']);
    else
        fprintf(2,['     AbsTol_TLM : ',mat2str(option1.AbsTol_TLM),' | ',mat2str(option2.AbsTol_TLM),'\n']);
    end        
    
    % AdjointSolve
    if ( isequal(option1.AdjointSolve,option2.AdjointSolve) )
        fprintf(  ['   AdjointSolve : ',num2str(option1.AdjointSolve),' | ',num2str(option2.AdjointSolve),'\n']);
    else
        fprintf(2,['   AdjointSolve : ',num2str(option1.AdjointSolve),' | ',num2str(option2.AdjointSolve),'\n']);
    end        
    
    % Autonomous
    if ( isequal(option1.Autonomous,option2.Autonomous) )
        fprintf(  ['     Autonomous : ',num2str(option1.Autonomous),' | ',num2str(option2.Autonomous),'\n']);
    else
        fprintf(2,['     Autonomous : ',num2str(option1.Autonomous),' | ',num2str(option2.Autonomous),'\n']);
    end       
    
    % ChunkSize
    if ( isequal(option1.ChunkSize,option2.ChunkSize) )
        fprintf(  ['      ChunkSize : ',num2str(option1.ChunkSize),' | ',num2str(option2.ChunkSize),'\n']);
    else
        fprintf(2,['      ChunkSize : ',num2str(option1.ChunkSize),' | ',num2str(option2.ChunkSize),'\n']);
    end      
    
    % Desired_Mode
    if ( isequal(option1.Desired_Mode,option2.Desired_Mode) )
        fprintf(  ['   Desired_Mode : ',num2str(option1.Desired_Mode),' | ',num2str(option2.Desired_Mode),'\n']);
    else
        fprintf(2,['   Desired_Mode : ',num2str(option1.Desired_Mode),' | ',num2str(option2.Desired_Mode),'\n']);
    end  
    
    % DirectADJ
    if ( isequal(option1.DirectADJ,option2.DirectADJ) )
        fprintf(  ['      DirectADJ : ',num2str(option1.DirectADJ),' | ',num2str(option2.DirectADJ),'\n']);
    else
        fprintf(2,['      DirectADJ : ',num2str(option1.DirectADJ),' | ',num2str(option2.DirectADJ),'\n']);
    end      
    
    % DirectTLM
    if ( isequal(option1.DirectTLM,option2.DirectTLM) )
        fprintf(  ['      DirectTLM : ',num2str(option1.DirectTLM),' | ',num2str(option2.DirectTLM),'\n']);
    else
        fprintf(2,['      DirectTLM : ',num2str(option1.DirectTLM),' | ',num2str(option2.DirectTLM),'\n']);
    end      
    
    % displayStats
    if ( isequal(option1.displayStats,option2.displayStats) )
        fprintf(  ['   displayStats : ',num2str(option1.displayStats),' | ',num2str(option2.displayStats),'\n']);
    else
        fprintf(2,['   displayStats : ',num2str(option1.displayStats),' | ',num2str(option2.displayStats),'\n']);
    end          
     
    % displaySteps
    if ( isequal(option1.displaySteps,option2.displaySteps) )
        fprintf(  ['   displaySteps : ',num2str(option1.displaySteps),' | ',num2str(option2.displaySteps),'\n']);
    else
        fprintf(2,['   displaySteps : ',num2str(option1.displaySteps),' | ',num2str(option2.displaySteps),'\n']);
    end      
    
    % DRDP
    if ( isequal(option1.DRDP,option2.DRDP) )
        if ( isa(option1.DRDP,'function_handle') && isa(option2.DRDP,'function_handle') )
            fprintf(  ['           DRDP : ',func2str(option1.DRDP),' | ',func2str(option2.DRDP),'\n']);
        else
            fprintf(  ['           DRDP : ',num2str(option1.DRDP),' | ',num2str(option2.DRDP),'\n']);
        end
    else
        if ( isa(option1.DRDP,'function_handle') && isa(option2.DRDP,'function_handle') )
            fprintf(2,['           DRDP : ',func2str(option1.DRDP),' | ',func2str(option2.DRDP),'\n']);
        elseif ( isa(option1.DRDP,'function_handle') && ~isa(option2.DRDP,'function_handle') )
            fprintf(2,['           DRDP : ',func2str(option1.DRDP),' | ',num2str(option2.DRDP),'\n']);
        elseif ( ~isa(option1.DRDP,'function_handle') && isa(option2.DRDP,'function_handle') )
            fprintf(2,['           DRDP : ',num2str(option1.DRDP),' | ',func2str(option2.DRDP),'\n']);
        elseif ( ~isa(option1.DRDP,'function_handle') && ~isa(option2.DRDP,'function_handle') )
            fprintf(2,['           DRDP : ',num2str(option1.DRDP),' | ',num2str(option2.DRDP),'\n']);
        end
    end      
    
    % DRDY
    if ( isequal(option1.DRDY,option2.DRDY) )
        if ( isa(option1.DRDY,'function_handle') && isa(option2.DRDY,'function_handle') )
            fprintf(  ['           DRDY : ',func2str(option1.DRDY),' | ',func2str(option2.DRDY),'\n']);
        else
            fprintf(  ['           DRDY : ',num2str(option1.DRDY),' | ',num2str(option2.DRDY),'\n']);
        end
    else
        if ( isa(option1.DRDY,'function_handle') && isa(option2.DRDY,'function_handle') )
            fprintf(2,['           DRDY : ',func2str(option1.DRDY),' | ',func2str(option2.DRDY),'\n']);
        elseif ( isa(option1.DRDY,'function_handle') && ~isa(option2.DRDY,'function_handle') )
            fprintf(2,['           DRDY : ',func2str(option1.DRDY),' | ',num2str(option2.DRDY),'\n']);
        elseif ( ~isa(option1.DRDY,'function_handle') && isa(option2.DRDY,'function_handle') )
            fprintf(2,['           DRDY : ',num2str(option1.DRDY),' | ',func2str(option2.DRDY),'\n']);
        elseif ( ~isa(option1.DRDY,'function_handle') && ~isa(option2.DRDY,'function_handle') )
            fprintf(2,['           DRDY : ',num2str(option1.DRDY),' | ',num2str(option2.DRDY),'\n']);
        end
    end      
    
    % FacMax
    if ( isequal(option1.FacMax,option2.FacMax) )
        fprintf(  ['         FacMax : ',num2str(option1.FacMax),' | ',num2str(option2.FacMax),'\n']);
    else
        fprintf(2,['         FacMax : ',num2str(option1.FacMax),' | ',num2str(option2.FacMax),'\n']);
    end      
    
    % FacMin
    if ( isequal(option1.FacMin,option2.FacMin) )
        fprintf(  ['         FacMin : ',num2str(option1.FacMin),' | ',num2str(option2.FacMin),'\n']);
    else
        fprintf(2,['         FacMin : ',num2str(option1.FacMin),' | ',num2str(option2.FacMin),'\n']);
    end          
    
    % FacRej
    if ( isequal(option1.FacRej,option2.FacRej) )
        fprintf(  ['         FacRej : ',num2str(option1.FacRej),' | ',num2str(option2.FacRej),'\n']);
    else
        fprintf(2,['         FacRej : ',num2str(option1.FacRej),' | ',num2str(option2.FacRej),'\n']);
    end    
    
    % FacSafe
    if ( isequal(option1.FacSafe,option2.FacSafe) )
        fprintf(  ['        FacSafe : ',num2str(option1.FacSafe),' | ',num2str(option2.FacSafe),'\n']);
    else
        fprintf(2,['        FacSafe : ',num2str(option1.FacSafe),' | ',num2str(option2.FacSafe),'\n']);
    end        
    
    % FacSafeHigh
    if ( isequal(option1.FacSafeHigh,option2.FacSafeHigh) )
        fprintf(  ['    FacSafeHigh : ',num2str(option1.FacSafeHigh),' | ',num2str(option2.FacSafeHigh),'\n']);
    else
        fprintf(2,['    FacSafeHigh : ',num2str(option1.FacSafeHigh),' | ',num2str(option2.FacSafeHigh),'\n']);
    end        
    
    % FacSafeLow
    if ( isequal(option1.FacSafeLow,option2.FacSafeLow) )
        fprintf(  ['     FacSafeLow : ',num2str(option1.FacSafeLow),' | ',num2str(option2.FacSafeLow),'\n']);
    else
        fprintf(2,['     FacSafeLow : ',num2str(option1.FacSafeLow),' | ',num2str(option2.FacSafeLow),'\n']);
    end        
    
    % FDAprox
    if ( isequal(option1.FDAprox,option2.FDAprox) )
        fprintf(  ['        FDAprox : ',num2str(option1.FDAprox),' | ',num2str(option2.FDAprox),'\n']);
    else
        fprintf(2,['        FDAprox : ',num2str(option1.FDAprox),' | ',num2str(option2.FDAprox),'\n']);
    end       
    
    % FDIncrement
    if ( isequal(option1.FDIncrement,option2.FDIncrement) )
        fprintf(  ['    FDIncrement : ',num2str(option1.FDIncrement),' | ',num2str(option2.FDIncrement),'\n']);
    else
        fprintf(2,['    FDIncrement : ',num2str(option1.FDIncrement),' | ',num2str(option2.FDIncrement),'\n']);
    end       
    
    % Gustafsson
    if ( isequal(option1.Gustafsson,option2.Gustafsson) )
        fprintf(  ['     Gustafsson : ',num2str(option1.Gustafsson),' | ',num2str(option2.Gustafsson),'\n']);
    else
        fprintf(2,['     Gustafsson : ',num2str(option1.Gustafsson),' | ',num2str(option2.Gustafsson),'\n']);
    end        
    
    % Hess_vec
    
    % Hesstr_vec
    
    % Hesstr_vec_f_py
    if ( isequal(option1.Hesstr_vec_f_py,option2.Hesstr_vec_f_py) )
        if ( isa(option1.Hesstr_vec_f_py,'function_handle') && isa(option2.Hesstr_vec_f_py,'function_handle') )
            fprintf(  ['Hesstr_vec_f_py : ',func2str(option1.Hesstr_vec_f_py),' | ',func2str(option2.Hesstr_vec_f_py),'\n']);
        else
            fprintf(  ['Hesstr_vec_f_py : ',num2str(option1.Hesstr_vec_f_py),' | ',num2str(option2.Hesstr_vec_f_py),'\n']);
        end
    else
        if ( isa(option1.Hesstr_vec_f_py,'function_handle') && isa(option2.Hesstr_vec_f_py,'function_handle') )
            fprintf(2,['Hesstr_vec_f_py : ',func2str(option1.Hesstr_vec_f_py),' | ',func2str(option2.Hesstr_vec_f_py),'\n']);
        elseif ( isa(option1.Hesstr_vec_f_py,'function_handle') && ~isa(option2.Hesstr_vec_f_py,'function_handle') )
            fprintf(2,['Hesstr_vec_f_py : ',func2str(option1.Hesstr_vec_f_py),' | ',num2str(option2.Hesstr_vec_f_py),'\n']);
        elseif ( ~isa(option1.Hesstr_vec_f_py,'function_handle') && isa(option2.Hesstr_vec_f_py,'function_handle') )
            fprintf(2,['Hesstr_vec_f_py : ',num2str(option1.Hesstr_vec_f_py),' | ',func2str(option2.Hesstr_vec_f_py),'\n']);
        elseif ( ~isa(option1.Hesstr_vec_f_py,'function_handle') && ~isa(option2.Hesstr_vec_f_py,'function_handle') )
            fprintf(2,['Hesstr_vec_f_py : ',num2str(option1.Hesstr_vec_f_py),' | ',num2str(option2.Hesstr_vec_f_py),'\n']);
        end
    end       
 
    % Hesstr_vec_r
    if ( isequal(option1.Hesstr_vec_r,option2.Hesstr_vec_r) )
        if ( isa(option1.Hesstr_vec_r,'function_handle') && isa(option2.Hesstr_vec_r,'function_handle') )
            fprintf(  ['   Hesstr_vec_r : ',func2str(option1.Hesstr_vec_r),' | ',func2str(option2.Hesstr_vec_r),'\n']);
        else
            fprintf(  ['   Hesstr_vec_r : ',num2str(option1.Hesstr_vec_r),' | ',num2str(option2.Hesstr_vec_r),'\n']);
        end
    else
        if ( isa(option1.Hesstr_vec_r,'function_handle') && isa(option2.Hesstr_vec_r,'function_handle') )
            fprintf(2,['   Hesstr_vec_r : ',func2str(option1.Hesstr_vec_r),' | ',func2str(option2.Hesstr_vec_r),'\n']);
        elseif ( isa(option1.Hesstr_vec_r,'function_handle') && ~isa(option2.Hesstr_vec_r,'function_handle') )
            fprintf(2,['   Hesstr_vec_r : ',func2str(option1.Hesstr_vec_r),' | ',num2str(option2.Hesstr_vec_r),'\n']);
        elseif ( ~isa(option1.Hesstr_vec_r,'function_handle') && isa(option2.Hesstr_vec_r,'function_handle') )
            fprintf(2,['   Hesstr_vec_r : ',num2str(option1.Hesstr_vec_r),' | ',func2str(option2.Hesstr_vec_r),'\n']);
        elseif ( ~isa(option1.Hesstr_vec_r,'function_handle') && ~isa(option2.Hesstr_vec_r,'function_handle') )
            fprintf(2,['   Hesstr_vec_r : ',num2str(option1.Hesstr_vec_r),' | ',num2str(option2.Hesstr_vec_r),'\n']);
        end
    end      
    
    % Hesstr_vec_r_py
    if ( isequal(option1.Hesstr_vec_r_py,option2.Hesstr_vec_r_py) )
        if ( isa(option1.Hesstr_vec_r_py,'function_handle') && isa(option2.Hesstr_vec_r_py,'function_handle') )
            fprintf(  ['Hesstr_vec_r_py : ',func2str(option1.Hesstr_vec_r_py),' | ',func2str(option2.Hesstr_vec_r_py),'\n']);
        else
            fprintf(  ['Hesstr_vec_r_py : ',num2str(option1.Hesstr_vec_r_py),' | ',num2str(option2.Hesstr_vec_r_py),'\n']);
        end
    else
        if ( isa(option1.Hesstr_vec_r_py,'function_handle') && isa(option2.Hesstr_vec_r_py,'function_handle') )
            fprintf(2,['Hesstr_vec_r_py : ',func2str(option1.Hesstr_vec_r_py),' | ',func2str(option2.Hesstr_vec_r_py),'\n']);
        elseif ( isa(option1.Hesstr_vec_r_py,'function_handle') && ~isa(option2.Hesstr_vec_r_py,'function_handle') )
            fprintf(2,['Hesstr_vec_r_py : ',func2str(option1.Hesstr_vec_r_py),' | ',num2str(option2.Hesstr_vec_r_py),'\n']);
        elseif ( ~isa(option1.Hesstr_vec_r_py,'function_handle') && isa(option2.Hesstr_vec_r_py,'function_handle') )
            fprintf(2,['Hesstr_vec_r_py : ',num2str(option1.Hesstr_vec_r_py),' | ',func2str(option2.Hesstr_vec_r_py),'\n']);
        elseif ( ~isa(option1.Hesstr_vec_r_py,'function_handle') && ~isa(option2.Hesstr_vec_r_py,'function_handle') )
            fprintf(2,['Hesstr_vec_r_py : ',num2str(option1.Hesstr_vec_r_py),' | ',num2str(option2.Hesstr_vec_r_py),'\n']);
        end
    end      
    
   
    % Hmax
    if ( isequal(option1.Hmax,option2.Hmax) )
        fprintf(  ['           Hmax : ',num2str(option1.Hmax),' | ',num2str(option2.Hmax),'\n']);
    else
        fprintf(2,['           Hmax : ',num2str(option1.Hmax),' | ',num2str(option2.Hmax),'\n']);
    end        
    
    % Hmin
    if ( isequal(option1.Hmin,option2.Hmin) )
        fprintf(  ['           Hmin : ',num2str(option1.Hmin),' | ',num2str(option2.Hmin),'\n']);
    else
        fprintf(2,['           Hmin : ',num2str(option1.Hmin),' | ',num2str(option2.Hmin),'\n']);
    end    
    
    % Hstart
    if ( isequal(option1.Hstart,option2.Hstart) )
        fprintf(  ['         Hstart : ',num2str(option1.Hstart),' | ',num2str(option2.Hstart),'\n']);
    else
        fprintf(2,['         Hstart : ',num2str(option1.Hstart),' | ',num2str(option2.Hstart),'\n']);
    end        
    
    % Jacp
    %
    if ( isequal(option1.Jacp,option2.Jacp) )
        if ( isa(option1.Jacp,'function_handle') && isa(option2.Jacp,'function_handle') )
            fprintf(  ['           Jacp : ',func2str(option1.Jacp),' | ',func2str(option2.Jacp),'\n']);
        else
            fprintf(  ['           Jacp : ',num2str(option1.Jacp),' | ',num2str(option2.Jacp),'\n']);
        end
    else
        if ( isa(option1.Jacp,'function_handle') && isa(option2.Jacp,'function_handle') )
            fprintf(2,['          Jacp : ',func2str(option1.Jacp),' | ',func2str(option2.Jacp),'\n']);
        elseif ( isa(option1.Jacp,'function_handle') && ~isa(option2.Jacp,'function_handle') )
            fprintf(2,['          Jacp : ',func2str(option1.Jacp),' | ',num2str(option2.Jacp),'\n']);
        elseif ( ~isa(option1.Jacp,'function_handle') && isa(option2.Jacp,'function_handle') )
            fprintf(2,['          Jacp : ',num2str(option1.Jacp),' | ',func2str(option2.Jacp),'\n']);
        elseif ( ~isa(option1.Jacp,'function_handle') && ~isa(option2.Jacp,'function_handle') )
            fprintf(2,['          Jacp : ',num2str(option1.Jacp),' | ',num2str(option2.Jacp),'\n']);
        end
    end      
    
    % TLMNewtonEst
    if ( isequal(option1.TLMNewtonEst,option2.TLMNewtonEst) )
        fprintf(  ['   TLMNewtonEst : ',num2str(option1.TLMNewtonEst),' | ',num2str(option2.TLMNewtonEst),'\n']);
    else
        fprintf(2,['   TLMNewtonEst : ',num2str(option1.TLMNewtonEst),' | ',num2str(option2.TLMNewtonEst),'\n']);
    end        
    
    % ITOL
    if ( isequal(option1.ITOL,option2.ITOL) )
        fprintf(  ['           ITOL : ',num2str(option1.ITOL),' | ',num2str(option2.ITOL),'\n']);
    else
        fprintf(2,['           ITOL : ',num2str(option1.ITOL),' | ',num2str(option2.ITOL),'\n']);
    end        
    
    
    % Jacobian
    if ( isequal(option1.Jacobian,option2.Jacobian) )
        if ( isa(option1.Jacobian,'function_handle') && isa(option2.Jacobian,'function_handle') )
            fprintf(  ['       Jacobian : ',func2str(option1.Jacobian),' | ',func2str(option2.Jacobian),'\n']);
        else
            fprintf(  ['       Jacobian : ',num2str(option1.Jacobian),' | ',num2str(option2.Jacobian),'\n']);
        end
    else
        if ( isa(option1.Jacobian,'function_handle') && isa(option2.Jacobian,'function_handle') )
            fprintf(2,['       Jacobian : ',func2str(option1.Jacobian),' | ',func2str(option2.Jacobian),'\n']);
        elseif ( isa(option1.Jacobian,'function_handle') && ~isa(option2.Jacobian,'function_handle') )
            fprintf(2,['       Jacobian : ',func2str(option1.Jacobian),' | ',num2str(option2.Jacobian),'\n']);
        elseif ( ~isa(option1.Jacobian,'function_handle') && isa(option2.Jacobian,'function_handle') )
            fprintf(2,['       Jacobian : ',num2str(option1.Jacobian),' | ',func2str(option2.Jacobian),'\n']);
        elseif ( ~isa(option1.Jacobian,'function_handle') && ~isa(option2.Jacobian,'function_handle') )
            fprintf(2,['       Jacobian : ',num2str(option1.Jacobian),' | ',num2str(option2.Jacobian),'\n']);
        end
    end    
    
    % Lambda
    if ( isequal(option1.Lambda,option2.Lambda) )
        if ( isa(option1.Lambda,'function_handle') && isa(option2.Lambda,'function_handle') )
            fprintf(  ['         Lambda : ',func2str(option1.Lambda),' | ',func2str(option2.Lambda),'\n']);
        else
            fprintf(  ['         Lambda : ',mat2str(option1.Lambda),' | ',mat2str(option2.Lambda),'\n']);
        end
    else
        if ( isa(option1.Lambda,'function_handle') && isa(option2.Lambda,'function_handle') )
            fprintf(2,['         Lambda : ',func2str(option1.Lambda),' | ',func2str(option2.Lambda),'\n']);
        elseif ( isa(option1.Lambda,'function_handle') && ~isa(option2.Lambda,'function_handle') )
            fprintf(2,['         Lambda : ',func2str(option1.Lambda),' | ',num2str(option2.Lambda),'\n']);
        elseif ( ~isa(option1.Lambda,'function_handle') && isa(option2.Lambda,'function_handle') )
            fprintf(2,['         Lambda : ',num2str(option1.Lambda),' | ',func2str(option2.Lambda),'\n']);
        elseif ( ~isa(option1.Lambda,'function_handle') && ~isa(option2.Lambda,'function_handle') )
            fprintf(2,['         Lambda : ',num2str(option1.Lambda),' | ',num2str(option2.Lambda),'\n']);
        end
    end      
    
    % MatrixFree
    if ( isequal(option1.MatrixFree,option2.MatrixFree) )
        fprintf(  ['     MatrixFree : ',num2str(option1.MatrixFree),' | ',num2str(option2.MatrixFree),'\n']);
    else
        fprintf(2,['     MatrixFree : ',num2str(option1.MatrixFree),' | ',num2str(option2.MatrixFree),'\n']);
    end        
    
    % Max_no_steps
    if ( isequal(option1.Max_no_steps,option2.Max_no_steps) )
        fprintf(  ['   Max_no_steps : ',num2str(option1.Max_no_steps),' | ',num2str(option2.Max_no_steps),'\n']);
    else
        fprintf(2,['   Max_no_steps : ',num2str(option1.Max_no_steps),' | ',num2str(option2.Max_no_steps),'\n']);
    end        
    
    % MaxOrder
    if ( isequal(option1.MaxOrder,option2.MaxOrder) )
        fprintf(  ['       MaxOrder : ',num2str(option1.MaxOrder),' | ',num2str(option2.MaxOrder),'\n']);
    else
        fprintf(2,['       MaxOrder : ',num2str(option1.MaxOrder),' | ',num2str(option2.MaxOrder),'\n']);
    end        
    
    % Method
    if ( isequal(option1.Method,option2.Method) )
        fprintf(  ['         Method : ',num2str(option1.Method),' | ',num2str(option2.Method),'\n']);
    else
        fprintf(2,['         Method : ',num2str(option1.Method),' | ',num2str(option2.Method),'\n']);
    end        
    
    % Mu
    if ( isequal(option1.Mu,option2.Mu) )
        if ( isa(option1.Mu,'function_handle') && isa(option2.Mu,'function_handle') )
            fprintf(  ['             Mu : ',func2str(option1.Mu),' | ',func2str(option2.Mu),'\n']);
        else
            fprintf(  ['             Mu : ',num2str(option1.Mu),' | ',num2str(option2.Mu),'\n']);
        end
    else
        if ( isa(option1.Mu,'function_handle') && isa(option2.Mu,'function_handle') )
            fprintf(2,['             Mu : ',func2str(option1.Mu),' | ',func2str(option2.Mu),'\n']);
        elseif ( isa(option1.Mu,'function_handle') && ~isa(option2.Mu,'function_handle') )
            fprintf(2,['             Mu : ',func2str(option1.Mu),' | ',num2str(option2.Mu),'\n']);
        elseif ( ~isa(option1.Mu,'function_handle') && isa(option2.Mu,'function_handle') )
            fprintf(2,['             Mu : ',num2str(option1.Mu),' | ',func2str(option2.Mu),'\n']);
        elseif ( ~isa(option1.Mu,'function_handle') && ~isa(option2.Mu,'function_handle') )
            fprintf(2,['             Mu : ',num2str(option1.Mu),' | ',num2str(option2.Mu),'\n']);
        end
    end      
    
    % NADJ
    if ( isequal(option1.NADJ,option2.NADJ) )
        fprintf(  ['           NADJ : ',num2str(option1.NADJ),' | ',num2str(option2.NADJ),'\n']);
    else
        fprintf(2,['           NADJ : ',num2str(option1.NADJ),' | ',num2str(option2.NADJ),'\n']);
    end        
    
    % NewtonMaxIt
    if ( isequal(option1.NewtonMaxIt,option2.NewtonMaxIt) )
        fprintf(  ['    NewtonMaxIt : ',num2str(option1.NewtonMaxIt),' | ',num2str(option2.NewtonMaxIt),'\n']);
    else
        fprintf(2,['    NewtonMaxIt : ',num2str(option1.NewtonMaxIt),' | ',num2str(option2.NewtonMaxIt),'\n']);
    end        
    
    % NewtonTol
    if ( isequal(option1.NewtonTol,option2.NewtonTol) )
        fprintf(  ['      NewtonTol : ',num2str(option1.NewtonTol),' | ',num2str(option2.NewtonTol),'\n']);
    else
        fprintf(2,['      NewtonTol : ',num2str(option1.NewtonTol),' | ',num2str(option2.NewtonTol),'\n']);
    end        
    
    % NP
    if ( isequal(option1.NP,option2.NP) )
        fprintf(  ['             NP : ',num2str(option1.NP),' | ',num2str(option2.NP),'\n']);
    else
        fprintf(2,['             NP : ',num2str(option1.NP),' | ',num2str(option2.NP),'\n']);
    end        
    
    % NTLM
    if ( isequal(option1.NTLM,option2.NTLM) )
        fprintf(  ['           NTLM : ',num2str(option1.NTLM),' | ',num2str(option2.NTLM),'\n']);
    else
        fprintf(2,['           NTLM : ',num2str(option1.NTLM),' | ',num2str(option2.NTLM),'\n']);
    end        
    
    % QFun
    if ( isequal(option1.QFun,option2.QFun) )
        if ( isa(option1.QFun,'function_handle') && isa(option2.QFun,'function_handle') )
            fprintf(  ['           QFun : ',func2str(option1.QFun),' | ',func2str(option2.QFun),'\n']);
        else 
            fprintf(  ['           QFun : ',num2str(option1.QFun),' | ',num2str(option2.QFun),'\n']);
        end
    else
        if ( isa(option1.QFun,'function_handle') && isa(option2.QFun,'function_handle') )
            fprintf(2,['          QFun : ',func2str(option1.QFun),' | ',func2str(option2.QFun),'\n']);
        elseif ( isa(option1.QFun,'function_handle') && ~isa(option2.QFun,'function_handle') )
            fprintf(2,['          QFun : ',func2str(option1.QFun),' | ',num2str(option2.QFun),'\n']);
        elseif ( ~isa(option1.QFun,'function_handle') && isa(option2.QFun,'function_handle') )
            fprintf(2,['          QFun : ',num2str(option1.QFun),' | ',func2str(option2.QFun),'\n']);
        elseif ( ~isa(option1.QFun,'function_handle') && ~isa(option2.QFun,'function_handle') )
            fprintf(2,['          QFun : ',num2str(option1.QFun),' | ',num2str(option2.QFun),'\n']);
        end
    end      
    
    % Qmax
    if ( isequal(option1.Qmax,option2.Qmax) )
        fprintf(  ['           Qmax : ',num2str(option1.Qmax),' | ',num2str(option2.Qmax),'\n']);
    else 
        fprintf(2,['           Qmax : ',num2str(option1.Qmax),' | ',num2str(option2.Qmax),'\n']);
    end        
    
    % Qmin
    if ( isequal(option1.Qmin,option2.Qmin) )
        fprintf(  ['           Qmin : ',num2str(option1.Qmin),' | ',num2str(option2.Qmin),'\n']);
    else
        fprintf(2,['           Qmin : ',num2str(option1.Qmin),' | ',num2str(option2.Qmin),'\n']);
    end        
    
    % Quadrature
    if ( isequal(option1.Quadrature,option2.Quadrature) )
        if ( isa(option1.Quadrature,'function_handle') && isa(option2.Quadrature,'function_handle') )
            fprintf(  ['     Quadrature : ',func2str(option1.Quadrature),' | ',func2str(option2.Quadrature),'\n']);
        else
            fprintf(  ['     Quadrature : ',num2str(option1.Quadrature),' | ',num2str(option2.Quadrature),'\n']);
        end
    else
        if ( isa(option1.Quadrature,'function_handle') && isa(option2.Quadrature,'function_handle') )
            fprintf(2,['     Quadrature : ',func2str(option1.Quadrature),' | ',func2str(option2.Quadrature),'\n']);
        elseif ( isa(option1.Quadrature,'function_handle') && ~isa(option2.Quadrature,'function_handle') )
            fprintf(2,['     Quadrature : ',func2str(option1.Quadrature),' | ',num2str(option2.Quadrature),'\n']);
        elseif ( ~isa(option1.Quadrature,'function_handle') && isa(option2.Quadrature,'function_handle') )
            fprintf(2,['     Quadrature : ',num2str(option1.Quadrature),' | ',func2str(option2.Quadrature),'\n']);
        elseif ( ~isa(option1.Quadrature,'function_handle') && ~isa(option2.Quadrature,'function_handle') )
            fprintf(2,['     Quadrature : ',num2str(option1.Quadrature),' | ',num2str(option2.Quadrature),'\n']);
        end
    end      
    
    % RelTol
    if ( isequal(option1.RelTol,option2.RelTol) )
        fprintf(  ['         RelTol : ',mat2str(option1.RelTol),' | ',mat2str(option2.RelTol),'\n']);
    else
        fprintf(2,['         RelTol : ',mat2str(option1.RelTol),' | ',mat2str(option2.RelTol),'\n']);
    end        
    
    % RelTol_ADJ
    if ( isequal(option1.RelTol_ADJ,option2.RelTol_ADJ) )
        fprintf(  ['     RelTol_ADJ : ',mat2str(option1.RelTol_ADJ),' | ',mat2str(option2.RelTol_ADJ),'\n']);
    else
        fprintf(2,['     RelTol_ADJ : ',mat2str(option1.RelTol_ADJ),' | ',mat2str(option2.RelTol_ADJ),'\n']);
    end
    
    % RelTol_TLM
    if ( isequal(option1.RelTol_TLM,option2.RelTol_TLM) )
        fprintf(  ['     RelTol_TLM : ',mat2str(option1.RelTol_TLM),' | ',mat2str(option2.RelTol_TLM),'\n']);
    else
        fprintf(2,['     RelTol_TLM : ',mat2str(option1.RelTol_TLM),' | ',mat2str(option2.RelTol_TLM),'\n']);
    end        
    
    % SaveLU
    if ( isequal(option1.SaveLU,option2.SaveLU) )
        fprintf(  ['         SaveLU : ',num2str(option1.SaveLU),' | ',num2str(option2.SaveLU),'\n']);
    else
        fprintf(2,['         SaveLU : ',num2str(option1.SaveLU),' | ',num2str(option2.SaveLU),'\n']);
    end        
    
    % SdirkError
    if ( isequal(option1.SdirkError,option2.SdirkError) )
        fprintf(  ['     SdirkError : ',num2str(option1.SdirkError),' | ',num2str(option2.SdirkError),'\n']);
    else
        fprintf(2,['     SdirkError : ',num2str(option1.SdirkError),' | ',num2str(option2.SdirkError),'\n']);
    end        
    
    % storeCheckpoint
    if ( isequal(option1.storeCheckpoint,option2.storeCheckpoint) )
        fprintf(  ['storeCheckpoint : ',num2str(option1.storeCheckpoint),' | ',num2str(option2.storeCheckpoint),'\n']);
    else
        fprintf(2,['storeCheckpoint : ',num2str(option1.storeCheckpoint),' | ',num2str(option2.storeCheckpoint),'\n']);
    end        
    
    % StartNewton
    if ( isequal(option1.StartNewton,option2.StartNewton) )
        fprintf(  ['    StartNewton : ',num2str(option1.StartNewton),' | ',num2str(option2.StartNewton),'\n']);
    else
        fprintf(2,['    StartNewton : ',num2str(option1.StartNewton),' | ',num2str(option2.StartNewton),'\n']);
    end        
    
    % ThetaMin
    if ( isequal(option1.ThetaMin,option2.ThetaMin) )
        fprintf(  ['       ThetaMin : ',num2str(option1.ThetaMin),' | ',num2str(option2.ThetaMin),'\n']);
    else
        fprintf(2,['       ThetaMin : ',num2str(option1.ThetaMin),' | ',num2str(option2.ThetaMin),'\n']);
    end        
    
    % TLMTruncErr
    if ( isequal(option1.TLMTruncErr,option2.TLMTruncErr) )
        fprintf(  ['    TLMTruncErr : ',num2str(option1.TLMTruncErr),' | ',num2str(option2.TLMTruncErr),'\n']);
    else
        fprintf(2,['    TLMTruncErr : ',num2str(option1.TLMTruncErr),' | ',num2str(option2.TLMTruncErr),'\n']);
    end        
    
    % WarningConfig
    if ( isequal(option1.WarningConfig,option2.WarningConfig) )
        fprintf(  ['  WarningConfig : ',num2str(option1.WarningConfig),' | ',num2str(option2.WarningConfig),'\n']);
    else
        fprintf(2,['  WarningConfig : ',num2str(option1.WarningConfig),' | ',num2str(option2.WarningConfig),'\n']);
    end        
    
    % Y_TLM
    if ( isequal(option1.Y_TLM,option2.Y_TLM) )
        if ( isa(option1.Y_TLM,'function_handle') && isa(option2.Y_TLM,'function_handle') )
            fprintf(  ['          Y_TLM : ',func2str(option1.Y_TLM),' | ',func2str(option2.Y_TLM),'\n']);
        else
            fprintf(  ['          Y_TLM : ',mat2str(option1.Y_TLM),' | ',mat2str(option2.Y_TLM),'\n']);
        end
    else
        if ( isa(option1.Y_TLM,'function_handle') && isa(option2.Y_TLM,'function_handle') )
            fprintf(2,['          Y_TLM : ',func2str(option1.Y_TLM),' | ',func2str(option2.Y_TLM),'\n']);
        elseif ( isa(option1.Y_TLM,'function_handle') && ~isa(option2.Y_TLM,'function_handle') )
            fprintf(2,['          Y_TLM : ',func2str(option1.Y_TLM),' | ',mat2str(option2.Y_TLM),'\n']);
        elseif ( ~isa(option1.Y_TLM,'function_handle') && isa(option2.Y_TLM,'function_handle') )
            fprintf(2,['          Y_TLM : ',mat2str(option1.Y_TLM),' | ',func2str(option2.Y_TLM),'\n']);
        elseif ( ~isa(option1.Y_TLM,'function_handle') && ~isa(option2.Y_TLM,'function_handle') )
            fprintf(2,['          Y_TLM : ',mat2str(option1.Y_TLM),' | ',mat2str(option2.Y_TLM),'\n']);
        end
    end      
    
end

%%
% <html>
%   <div>
%       <img style="float: right" src="../../../CSL_LogoWithName_1.png" height="50px"></img>
%   </div>
% </html>

