%% MATLODE_OPTIONS
%
%% Syntax
%            MATLODE_OPTIONS
%            MATLODE_OPTIONS('family','implementation')
%  Options = MATLODE_OPTIONS('key1',value1,'key2',value2,...,'keyN',valueN)
%  Options = MATLODE_OPTIONS(Options,'key1',value1,'key2',value2,...,'keyN',valueN)
%
%% Input Parameters
% |'family'| specifies which family of integrator configuration parameters
% to print to the command window. (i.e. 'ERK', 'EXP', 'EXPK', 'RK', 'ROK',
% 'ROS', 'SDIRK')
%
% |'implementation'| specifies which implementation of integrator
% configuration parameters to print to the command window. (i.e. 'FWD",
% 'TLM', 'ADJ')
%
% |'key'| specifies which configuration parameter to assign |value| to.
%
% |value| is the value assigned to |'key'| in the _('key',value)_ option
% struc.
%
% |Options| is the option configuration _('key',value)_ struct.
%
%% Output Parameters
% |Options| is the option configuration _('key',value)_ struct.
%
%% Description
% A _('key',value)_ pair option structure for fine tuning MATLODE's library
% of ordinary differential equation integrators.
%
% |MATLODE_OPTIONS| displays ALL the option configuration parameter
% descriptions.
%
% |MATLODE_OPTIONS('family','implementation')| displays the optional and
% required configuration parameters associated with 'family' and
% 'implementation'.
%
% |Options = MATLODE_OPTIONS('key1',value1,'key2',value2,...,'keyN',valueN)|
%
% |Options = MATLODE_OPTIONS(Options,'key1',value1,'key2',value2,...,'keyN',valueN)|
%
%% Example
% To create a new option struct specifying the absolute and relative
% tolerance, pass the appropriate _('key',value)_ to |MATLODE_OPTIONS|.
%
%   Options = MATLODE_OPTIONS('AbsTol',1e-5,'RelTol',1e-5);
%
% To alter the option struct, pass the struct as the first parameter. 
%
%  Options = MATLODE_OPTIONS(Options,'NBasisVectors',4);
%
%% Reference
%
function options = MATLODE_OPTIONS(varargin)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   Add option descriptions here so you don't have to change all cases
%   individually
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    AbsTol_str          = '         AbsTol: [ Absolute error tolerance for Forward integrators. ]\n';
    AbsTol_ADJ_str      = '     AbsTol_ADJ: [ Absolute Newton iteration tolerance for solving adjoint system. ]\n';
    AbsTol_TLM_str      = '     AbsTol_TLM: [ Absolute error estimation for TLM at Newton stages. ]\n';
    Autonomous_str      = '     Autonomous: [ ]\n';
    ChunkSize_str       = '      ChunkSize: [ Appended memory block size ]\n';
    DirectADJ_str       = '      DirectADJ: [ Determines whether direct adjoint sensitivity analysis is performed. ]\n';
    DirectTLM_str       = '      DirectTLM: [ Determines whether direct tangent linear sensitivity analysis is performed. ]\n';
    displayStats_str    = '   displayStats: [ Determines whether statistics are displayed ]\n';
    displaySteps_str    = '   displaySteps: [ Determines whether steps are displayed ]\n';
    DRDP_str            = '           DRDP: [ User defined function: Derivative of r w.r.t. parameters. ]\n';
    DRDY_str            = '           DRDY: [ User defined function: Derivative of r w.r.t. y vector. ]\n';
    FacMax_str          = '         FacMax: [ Step size upper bound change ratio. ]\n';
    FacMin_str          = '         FacMin: [ Step size lower bound change ratio. ]\n';
    FacRej_str          = '         FacRej: [ Decrease step size factor after two successive rejections. ]\n';
    FacSafe_str         = '        FacSafe: [ Fact by which the new step is slightly smaller than the predicted value. ]\n';
    FDAprox_str         = '        FDAprox: [ Determines whether Jacobian Vector products by finite difference is used. ]\n';
    FDIncrement_str     = '    FDIncrement: [ ]\n';
    Gustafsson_str      = '     Gustafsson: [ An alternative error controller approach which may be advantageous depending on the model characteristics. ]\n';
    Hess_vec_str        = '       Hess_vec: [ ]\n';
    Hesstr_vec_str      = '     Hesstr_vec: [ ]\n';
    Hesstr_vec_f_py_str = 'Hesstr_vec_f_py: [ ]\n';
    Hesstr_vec_r_str    = '   Hesstr_vec_r: [ ]\n';
    Hesstr_vec_r_py_str = 'Hesstr_vec_r_py: [ ]\n';
    Hmax_str            = '           Hmax: [ Step size upper bound. ]\n';
    Hmin_str            = '           Hmin: [ Step size lower bound. ]\n';
    Hstart_str          = '         Hstart: [ Initial step size. ]\n';
    ITOL_str            = '           ITOL: [ ]\n';
    Jacobian_str        = '       Jacobian: [ Jacobian function (user supplied) ]\n';
    Jacp_str            = '           Jacp: [ ]\n';
    Lambda_str          = '         Lambda: [ ]\n';
    MatrixFree_str      = '     MatrixFree: [ ]\n';
    Max_no_steps_str    = '   Max_no_steps: [ Maximum nunber of steps upper bound ]\n';
    Method_str          = '         Method: [ Determines which coefficients to use ]\n';
    Mu_str              = '             Mu: [ ]\n';
    NADJ_str            = '           NADJ: [ The number of cost functionals for which adjoints are evaluated simultaneously ]\n';
    NBasisVectors       = '  NBasisVectors: [ ] \n';
    NewtonMaxit_str     = '    NewtonMaxIt: [ Maximum number of newton iterations performed ]\n';
    NewtonTol_str       = '      NewtonTol: [ Newtons method stopping criterion. ]\n';
    NP_str              = '             NP: [ ]\n';
    NTLM_str            = '           NTLM: [ The number of tangent linear sensitivity coefficients ]\n';
    QFun_str            = '           QFun: [ ]\n';
    Qmax_str            = '           Qmax: [ Predicted step size to current step size upper bound ratio. ]\n';
    Qmin_str            = '           Qmin: [ Predicted step size to current step size lower bound ratio. ]\n';
    Quadrature_str      = '     Quadrature: [ ]\n';
    RelTol_str          = '         RelTol: [ Relative error tolerance ]\n';
    RelTol_ADJ_str      = '     RelTol_ADJ: [ Relative Newton iteration tolerance for solving adjoint system. ]\n';
    RelTol_TLM_str      = '     RelTol_TLM: [ Relative error estimation for TLM at Newton stages ]\n';
    SaveLU_str          = '         SaveLU: [ ]\n';
    SdirkError_str      = '     SdirkError: [ ]\n';
    StartNewton_str     = '    StartNewton: [ ]\n';    
    storeCheckpoint_str = 'storeCheckpoint: [ Determines whether intermediate values are stored ]\n';
    ThetaMin_str        = '       ThetaMin: [ Factir deciding whether the Jacobian should be recomputed. ]\n';
    TLMNewtonEst_str    = '   TLMNewtonEst: [ ]\n';
    TLMtruncErr_str     = '    TLMtruncErr: [ ]\n';
    WarningConfig_str   = '  WarningConfig: [ Determines whether warning messages are displayed during option configuration. ]\n';
    Y_TLM_str           = '          Y_TLM: [ Contains the sensitivities of Y with respect to the specified coefficients ]\n';

% Print out all possible values of properties.
if (nargin == 0) && (nargout == 0)
    fprintf( AbsTol_str );
    fprintf( AbsTol_ADJ_str );
    fprintf( AbsTol_TLM_str ); 
    fprintf( Autonomous_str );
    fprintf( ChunkSize_str );
        
    fprintf( DirectADJ_str );
    fprintf( DirectTLM_str );
    fprintf( displayStats_str );
    fprintf( displaySteps_str );   
    fprintf( DRDP_str );
    fprintf( DRDY_str );     
    fprintf( FacMax_str );
    fprintf( FacMin_str );
    fprintf( FacRej_str );
    fprintf( FacSafe_str );
    fprintf( FDAprox_str );
    fprintf( FDIncrement_str );
    fprintf( Gustafsson_str );
    fprintf( Hess_vec_str );  
    fprintf( Hesstr_vec_str );    
    fprintf( Hesstr_vec_f_py_str );
    fprintf( Hesstr_vec_r_str );    
    fprintf( Hesstr_vec_r_py_str );
    fprintf( Hmax_str );
    fprintf( Hmin_str );
    fprintf( Hstart_str );
    fprintf( ITOL_str );    
    fprintf( Jacobian_str );    
    fprintf( Jacp_str );    
    fprintf( Lambda_str );
    fprintf( MatrixFree_str );
    fprintf( Max_no_steps_str );
    fprintf( Method_str );
    fprintf( Mu_str );    
    fprintf( NADJ_str );  
    fprintf( NewtonMaxit_str );
    fprintf( NewtonTol_str );
    fprintf( NP_str );
    fprintf( NTLM_str );      
    fprintf( QFun_str );    
    fprintf( Qmax_str );
    fprintf( Qmin_str );
    fprintf( Quadrature_str );    
    fprintf( RelTol_str );
    fprintf( RelTol_ADJ_str );
    fprintf( RelTol_TLM_str );    
    fprintf( SaveLU_str );
    fprintf( SdirkError_str );
    fprintf( StartNewton_str );
    fprintf( storeCheckpoint_str );    
    fprintf( ThetaMin_str );        
    fprintf( TLMNewtonEst_str );
    fprintf( TLMtruncErr_str );
    fprintf( WarningConfig_str );
    fprintf( Y_TLM_str );    
    fprintf( '\n' );
  return;
end

% print out family or implementation specific properties
if (nargin == 1) && (nargout == 0)
    firstArg = varargin{1};
    switch ( firstArg )
        case 'ERK'
            fprintf( '\n' );
            fprintf( AbsTol_str );
            fprintf( AbsTol_ADJ_str );
            fprintf( AbsTol_TLM_str );
            fprintf( ChunkSize_str );
            
            fprintf( displayStats_str );
            fprintf( displaySteps_str );
            fprintf( DRDP_str );
            fprintf( DRDY_str );
            fprintf( FacMax_str );
            fprintf( FacMin_str );
            fprintf( FacRej_str );
            fprintf( FacSafe_str);
            fprintf( Hmax_str );
            fprintf( Hmin_str );
            fprintf( Hstart_str );
            fprintf( ITOL_str );
            fprintf( Jacobian_str );
            fprintf( Jacp_str );
            fprintf( Lambda_str );
            fprintf( Max_no_steps_str );
            fprintf( Method_str );
            fprintf( Mu_str );
            fprintf( QFun_str );
            fprintf( Quadrature_str );
            fprintf( RelTol_str );
            fprintf( RelTol_TLM_str );
            fprintf( storeCheckpoint_str );
            fprintf( TLMtruncErr_str );
            fprintf( WarningConfig_str );
            fprintf( Y_TLM_str );
            fprintf( '\n' );
        case 'RK'
        case 'ROS'
        case 'SDIRK'
        case 'FWD'
        case 'TLM'
        case 'ADJ'
        otherwise
            disp('ERROR_1: see fatOde_OPTIONS_Set');
    end
    return;
end

if (nargin == 2) && (nargout == 0)
    firstArg = varargin{1};
    secondArg = varargin{2};
    if ( strcmp(firstArg,'ERK') ||  strcmp(firstArg,'RK') ||  strcmp(firstArg,'ROS') ||  strcmp(firstArg,'SDIRK') )
        if ( strcmp(secondArg,'FWD') || strcmp(secondArg,'TLM') || strcmp(secondArg,'ADJ') )
            family = firstArg;
            implementation = secondArg;
        else
            error('Internal Error: Input parameters must be (family,implementation) or (implementation,family)');
        end
    elseif ( strcmp(secondArg,'ERK') ||  strcmp(secondArg,'RK') ||  strcmp(secondArg,'ROS') ||  strcmp(secondArg,'SDIRK') )
        if ( strcmp(firstArg,'FWD') || strcmp(firstArg,'TLM') || strcmp(firstArg,'ADJ') ) 
            family = secondArg;
            implementation = firstArg;
        else
            error('Internal Error: Input parameters must be (family,implementation) or (implementation,family)');
        end
    end
    
    switch ( family )
        case 'ERK'
            switch ( implementation )
                case 'FWD'
                    fprintf( '\n' );
                    fprintf( AbsTol_str );
                    fprintf( ChunkSize_str );
                    fprintf( displayStats_str );
                    fprintf( displaySteps_str );
                    fprintf( FacMax_str );
                    fprintf( FacMin_str );
                    fprintf( FacRej_str );
                    fprintf( FacSafe_str );
                    fprintf( Hmax_str );
                    fprintf( Hmin_str );
                    fprintf( Hstart_str );
                    fprintf( ITOL_str );
                    fprintf( Max_no_steps_str );
                    fprintf( Method_str );
                    fprintf( RelTol_str );
                    fprintf( storeCheckpoint_str );
                    fprintf( WarningConfig_str );
                    fprintf( '\n' );
                case 'TLM'
                    fprintf( '\n' );
                    fprintf( AbsTol_str );
                    fprintf( AbsTol_TLM_str );
                    fprintf( ChunkSize_str );
                    fprintf( displayStats_str );
                    fprintf( displaySteps_str );
                    fprintf( FacMax_str );
                    fprintf( FacMin_str );
                    fprintf( FacRej_str );
                    fprintf( FacSafe_str );
                    fprintf( FDAprox_str );
                    fprintf( FDIncrement_str );
                    fprintf( Hmax_str );
                    fprintf( Hmin_str );
                    fprintf( Hstart_str );
                    fprintf( ITOL_str );
                    fprintf( Jacobian_str );
                    fprintf( Max_no_steps_str );
                    fprintf( Method_str );
                    fprintf( RelTol_str );
                    fprintf( RelTol_TLM_str );
                    fprintf( storeCheckpoint_str );
                    fprintf( TLMtruncErr_str );
                    fprintf( WarningConfig_str );
                    fprintf( Y_TLM_str );
                    fprintf( '\n' );
                case 'ADJ'
                    fprintf( '\n' );
                    fprintf( AbsTol_str );
                    fprintf( AbsTol_ADJ_str );
                    fprintf( ChunkSize_str );
                    
                    fprintf( displayStats_str );
                    fprintf( displaySteps_str );
                    fprintf( DRDP_str );
                    fprintf( DRDY_str );
                    fprintf( FacMax_str );
                    fprintf( FacMin_str );
                    fprintf( FacRej_str );
                    fprintf( FacSafe_str );
                    fprintf( Hmax_str );
                    fprintf( Hmin_str );
                    fprintf( Hstart_str );
                    fprintf( ITOL_str );
                    fprintf( Jacobian_str );
                    fprintf( Jacp_str );
                    fprintf( Lambda_str );
                    fprintf( Max_no_steps_str );
                    fprintf( Method_str );
                    fprintf( Mu_str );
                    fprintf( QFun_str );
                    fprintf( Quadrature_str );
                    fprintf( RelTol_str );
                    fprintf( RelTol_ADJ_str );
                    fprintf( storeCheckpoint_str );
                    fprintf( WarningConfig_str );
                    fprintf( '\n' );
                otherwise
                    error('Input parameter implementation is invalid. This error should never occur.');
            end
        case 'RK'
            switch ( implementation )
                case 'FWD'
                    fprintf( '\n' );
                    fprintf( AbsTol_str );
                    fprintf( ChunkSize_str );
                    fprintf( displayStats_str );
                    fprintf( displaySteps_str );                    
                    fprintf( FacMax_str );
                    fprintf( FacMin_str );
                    fprintf( FacRej_str );
                    fprintf( FacSafe_str );                    
                    fprintf( Gustafsson_str );                    
                    fprintf( Hmax_str );
                    fprintf( Hmin_str );
                    fprintf( Hstart_str );                    
                    fprintf( ITOL_str );
                    fprintf( Jacobian_str );
                    fprintf( Max_no_steps_str );
                    fprintf( Method_str );
                    fprintf( NewtonMaxit_str );
                    fprintf( NewtonTol_str );
                    fprintf( Qmax_str );
                    fprintf( Qmin_str );                     
                    fprintf( RelTol_str );
                    fprintf( RelTol_TLM_str );                   
                    fprintf( SdirkError_str );
                    fprintf( StartNewton_str );
                    fprintf( storeCheckpoint_str );
                    fprintf( WarningConfig_str );
                    fprintf( ThetaMin_str );
                    fprintf( '\n' );
                case 'TLM'
                    fprintf( '\n' );
                    fprintf( AbsTol_str );
                    fprintf( AbsTol_TLM_str );
                    fprintf( ChunkSize_str );
                    fprintf( DirectTLM_str );
                    fprintf( displayStats_str );
                    fprintf( displaySteps_str );
                    fprintf( FacMax_str );
                    fprintf( FacMin_str );
                    fprintf( FacRej_str );
                    fprintf( FacSafe_str );
                    fprintf( Gustafsson_str );
                    fprintf( Hmax_str );
                    fprintf( Hmin_str );
                    fprintf( Hstart_str );
                    fprintf( ITOL_str );
                    fprintf( Jacobian_str );
                    fprintf( Max_no_steps_str );
                    fprintf( Method_str );
                    fprintf( NewtonMaxit_str );
                    fprintf( NewtonTol_str );
                    fprintf( Qmax_str );
                    fprintf( Qmin_str );
                    fprintf( RelTol_str );
                    fprintf( RelTol_TLM_str );
                    fprintf( SdirkError_str );
                    fprintf( StartNewton_str );
                    fprintf( storeCheckpoint_str );
                    fprintf( ThetaMin_str );
                    fprintf( TLMNewtonEst_str );
                    fprintf( TLMtruncErr_str );
                    fprintf( WarningConfig_str );
                    fprintf( Y_TLM_str );
                    fprintf( '\n' );
                case 'ADJ'
                    fprintf( '\n' );
                    fprintf( AbsTol_str );
                    fprintf( AbsTol_ADJ_str );
                    fprintf( ChunkSize_str );
                    
                    fprintf( DirectADJ_str );
                    fprintf( displayStats_str );
                    fprintf( displaySteps_str );  
                    fprintf( DRDP_str );
                    fprintf( DRDY_str );
                    fprintf( FacMax_str );
                    fprintf( FacMin_str );
                    fprintf( FacRej_str );
                    fprintf( FacSafe_str );
                    fprintf( Gustafsson_str );
                    fprintf( Hmax_str );
                    fprintf( Hmin_str );
                    fprintf( Hstart_str );                    
                    fprintf( ITOL_str );
                    fprintf( Jacobian_str );
                    fprintf( Jacp_str );
                    fprintf( Lambda_str );
                    fprintf( Max_no_steps_str );
                    fprintf( Method_str );
                    fprintf( Mu_str );
                    fprintf( NewtonMaxit_str );
                    fprintf( QFun_str );
                    fprintf( Qmax_str );
                    fprintf( Qmin_str );
                    fprintf( Quadrature_str );                    
                    fprintf( RelTol_str );
                    fprintf( SaveLU_str );
                    fprintf( SdirkError_str );                    
                    fprintf( StartNewton_str );
                    fprintf( storeCheckpoint_str );
                    fprintf( storeCheckpoint_str );
                    fprintf( WarningConfig_str );
                    fprintf( '\n' );
                otherwise
                    error('Input parameter implementation is invalid. This error should never occur.');
            end            
        case 'ROS'
            switch ( implementation )
                case 'FWD'
                    fprintf( '\n' );
                    fprintf( AbsTol_str );                    
                    fprintf( ChunkSize_str );
                    fprintf( displayStats_str );
                    fprintf( displaySteps_str );
                    fprintf( FacMax_str );
                    fprintf( FacMin_str );
                    fprintf( FacRej_str );
                    fprintf( FacSafe_str );
                    fprintf( Hmax_str );
                    fprintf( Hmin_str );
                    fprintf( Hstart_str );
                    fprintf( ITOL_str );
                    fprintf( Jacobian_str );
                    fprintf( Max_no_steps_str );
                    fprintf( Method_str );
                    fprintf( RelTol_str );
                    fprintf( storeCheckpoint_str );
                    fprintf( WarningConfig_str );
                    fprintf( '\n' );
                case 'TLM'
                    fprintf( '\n' );
                    fprintf( AbsTol_str );
                    fprintf( AbsTol_TLM_str );
                    fprintf( ChunkSize_str );
                    fprintf( displayStats_str );
                    fprintf( displaySteps_str );
                    fprintf( FacMax_str );
                    fprintf( FacMin_str );
                    fprintf( FacRej_str );
                    fprintf( FacSafe_str );
                    fprintf( Hmax_str );
                    fprintf( Hmin_str );
                    fprintf( Hstart_str );
                    fprintf( ITOL_str );
                    fprintf( Jacobian_str );
                    fprintf( Max_no_steps_str );
                    fprintf( Method_str );
                    fprintf( RelTol_str );
                    fprintf( RelTol_TLM_str );
                    fprintf( storeCheckpoint_str );
                    fprintf( TLMtruncErr_str );
                    fprintf( WarningConfig_str );
                    fprintf( Y_TLM_str );
                    fprintf( '\n' );
                case 'ADJ'
                    fprintf( '\n' );
                    fprintf( AbsTol_str );
                    fprintf( AbsTol_ADJ_str );
                    fprintf( ChunkSize_str );
                    
                    fprintf( displayStats_str );
                    fprintf( displaySteps_str );
                    fprintf( DRDP_str );
                    fprintf( DRDY_str );
                    fprintf( FacMax_str );
                    fprintf( FacMin_str );
                    fprintf( FacRej_str );
                    fprintf( FacSafe_str );
                    fprintf( Hess_vec_str );
                    fprintf( Hesstr_vec_str );
                    fprintf( Hesstr_vec_f_py_str );
                    fprintf( Hesstr_vec_r_str );
                    fprintf( Hesstr_vec_r_py_str );
                    fprintf( Hmax_str );
                    fprintf( Hmin_str );
                    fprintf( Hstart_str );                    
                    fprintf( ITOL_str );
                    fprintf( Jacobian_str );
                    fprintf( Jacp_str );
                    fprintf( Lambda_str );                    
                    fprintf( Max_no_steps_str );
                    fprintf( Method_str );
                    fprintf( Mu_str );
                    fprintf( QFun_str );
                    fprintf( Quadrature_str );                    
                    fprintf( RelTol_str );
                    fprintf( storeCheckpoint_str );
                    fprintf( WarningConfig_str );
                    fprintf( '\n' );
                otherwise
                    error('Input parameter implementation is invalid. This error should never occur.');
            end            
        case 'SDIRK'
            switch ( implementation )
                case 'FWD'
                    fprintf( '\n' );
                    fprintf( AbsTol_str );                    
                    fprintf( ChunkSize_str );
                    fprintf( displayStats_str );
                    fprintf( displaySteps_str );                    
                    fprintf( FacMax_str );
                    fprintf( FacMin_str );
                    fprintf( FacRej_str );
                    fprintf( FacSafe_str );                    
                    fprintf( Hmax_str );
                    fprintf( Hmin_str );
                    fprintf( Hstart_str );
                    fprintf( ITOL_str );
                    fprintf( Jacobian_str );
                    fprintf( MatrixFree_str );
                    fprintf( Max_no_steps_str );
                    fprintf( Method_str );
                    fprintf( NewtonMaxit_str );
                    fprintf( Qmax_str );
                    fprintf( Qmin_str );
                    fprintf( RelTol_str );
                    fprintf( StartNewton_str );
                    fprintf( storeCheckpoint_str );
                    fprintf( ThetaMin_str );
                    fprintf( WarningConfig_str );
                    fprintf( '\n' );
                case 'TLM'
                    fprintf( '\n' );
                    fprintf( AbsTol_str );
                    fprintf( AbsTol_TLM_str );
                    fprintf( ChunkSize_str );
                    fprintf( DirectTLM_str );
                    fprintf( displayStats_str );
                    fprintf( displaySteps_str );
                    fprintf( FacMax_str );
                    fprintf( FacMin_str );
                    fprintf( FacRej_str );
                    fprintf( FacSafe_str );                    
                    fprintf( Hmax_str );
                    fprintf( Hmin_str );
                    fprintf( Hstart_str );                    
                    fprintf( ITOL_str );
                    fprintf( Jacobian_str );
                    fprintf( Max_no_steps_str );
                    fprintf( Method_str );
                    fprintf( NewtonMaxit_str );
                    fprintf( RelTol_str );
                    fprintf( RelTol_TLM_str );
                    fprintf( StartNewton_str );
                    fprintf( storeCheckpoint_str );
                    fprintf( TLMNewtonEst_str );
                    fprintf( TLMtruncErr_str ); 
                    fprintf( WarningConfig_str );
                    fprintf( Y_TLM_str );
                    fprintf( '\n' );
                case 'ADJ'
                    fprintf( '\n' );
                    fprintf( AbsTol_str );
                    fprintf( AbsTol_ADJ_str );
                    fprintf( ChunkSize_str );
                    
                    fprintf( DirectADJ_str );
                    fprintf( displayStats_str );
                    fprintf( displaySteps_str );
                    fprintf( DRDP_str );
                    fprintf( DRDY_str );
                    fprintf( FacMax_str );
                    fprintf( FacMin_str );
                    fprintf( FacRej_str );
                    fprintf( FacSafe_str );                    
                    fprintf( Hmax_str );
                    fprintf( Hmin_str );
                    fprintf( Hstart_str );                    
                    fprintf( ITOL_str );
                    fprintf( Lambda_str );                    
                    fprintf( Jacobian_str );
                    fprintf( Jacp_str );
                    fprintf( MatrixFree_str );
                    fprintf( Max_no_steps_str );
                    fprintf( Method_str );
                    fprintf( Mu_str );
                    fprintf( NewtonMaxit_str );
                    fprintf( QFun_str );
                    fprintf( Quadrature_str );                    
                    fprintf( RelTol_str );
                    fprintf( SaveLU_str );
                    fprintf( StartNewton_str );
                    fprintf( storeCheckpoint_str );
                    fprintf( WarningConfig_str );
                    fprintf( '\n' );
                otherwise
                    error('Input parameter implementation is invalid. This error should never occur.');
            end            
        otherwise
            error('Input parameter family is invalid. This error should never occur.');
    end
        
    return;
end

Names = [
    'AbsTol             '
    'AbsTol_ADJ         '
    'AbsTol_TLM         '
    'AdaptiveArnoldiTol '
    'Autonomous         '
    'BasisExtended      '
    'BiOrthogonalLanczos'
    'BlockSize          '
    'ChunkSize          '
    'Desired_Mode       '
    'DirectADJ          '
    'DirectTLM          '
    'displayStats       '
    'displaySteps       '
    'DRDP               '
    'DRDY               '
    'FacMax             '
    'FacMin             '
    'FacRej             '
    'FacSafe            '
    'FDAprox            '
    'FDIncrement        '
    'GMRES_TOL          '
    'GMRES_MaxIt        '
    'GMRES_Restart      '
    'GMRES_P            '
    'Gustafsson         '
    'Hess_vec           '
    'Hesstr_vec         '
    'Hesstr_vec_f_py    '
    'Hesstr_vec_r       '     
    'Hesstr_vec_r_py    '
    'Hmax               '
    'Hmin               '
    'Hstart             '
    'IOArnoldi          '
    'JacobianAdjointVec '
    'Jacp               '
    'TLMNewtonEst       '
    'ITOL               '
    'Jacobian           '
    'Lambda             '
    'MatrixFree         '
    'Max_no_steps       '
    'Method             '
    'Mu                 '
    'NADJ               '
    'NBasisVectors      '
    'NewtonMaxIt        '
    'NewtonTol          '
    'NP                 '
    'NRecycledVectors   '
    'NTLM               '
    'QFun               '
    'Qmax               '
    'Qmin               '
    'Quadrature         '
    'RecycleBasisSize   '
    'RelTol             '
    'RelTol_ADJ         '
    'RelTol_TLM         '
    'SaveLU             '
    'SdirkError         '
    'OneStepIntegrator  '
    'storeCheckpoint    '
    'StartNewton        '
    'ThetaMin           '
    'TLMTruncErr        '
    'WarningConfig      '
    'Y_TLM              '
    ];
m = size(Names,1);
names = lower(Names);

% Combine all leading options structures o1, o2, ... in odeset(o1,o2,...).
options = [];
for j = 1:m
  options.(deblank(Names(j,:))) = [];
end
i = 1;
while i <= nargin
  arg = varargin{i};
  if ischar(arg)                         % arg is an option name
    break;
  end
  if ~isempty(arg)                      % [] is a valid options argument
    if ~isa(arg,'struct')
      error(message('MATLAB:odeset:NoPropNameOrStruct', i));
    end
    for j = 1:m
      if any(strcmp(fieldnames(arg),deblank(Names(j,:))))
        val = arg.(deblank(Names(j,:)));
      else
        val = [];
      end
      if ~isempty(val)
        options.(deblank(Names(j,:))) = val;
      end
    end
  end
  i = i + 1;
end

% A finite state machine to parse name-value pairs.
if rem(nargin-i+1,2) ~= 0
  error(message('MATLAB:odeset:ArgNameValueMismatch'));
end
expectval = 0;                          % start expecting a name, not a value
while i <= nargin
  arg = varargin{i};
    
  if ~expectval
    if ~ischar(arg)
      error(message('MATLAB:odeset:NoPropName', i));
    end
    
    lowArg = lower(arg);
    j = strmatch(lowArg,names);
    if isempty(j)                       % if no matches
      error(message('MATLAB:odeset:InvalidPropName', arg));
    elseif length(j) > 1                % if more than one match
      % Check for any exact matches (in case any names are subsets of others)
      k = strmatch(lowArg,names,'exact');
      if length(k) == 1
        j = k;
      else
            matches = deblank(Names(j(1),:));
        for k = j(2:length(j))'
                matches = [matches ', ' deblank(Names(k,:))]; %#ok<AGROW>
        end
            error(message('MATLAB:odeset:AmbiguousPropName',arg,matches));
      end
    end
    expectval = 1;                      % we expect a value next
    
  else
    options.(deblank(Names(j,:))) = arg;
    expectval = 0;
      
  end
  i = i + 1;
end

if expectval
  error(message('MATLAB:odeset:NoValueForProp', arg));
end
%% Major Modification History
% <html>
% <table border=1>
%   <tr>
%       <td><b>Date</b></td>
%       <td>Developer</td>
%       <td>Email</td>
%       <td>Action</td>
%   </tr>
%   <tr>
%       <td>1/1/2014</td>
%       <td>Tony D'Augustine</td>
%       <td>adaug13@vt,edu</td>
%       <td>Release MATLODE</td>
%   </tr>
% </table>
% </html>
% 
