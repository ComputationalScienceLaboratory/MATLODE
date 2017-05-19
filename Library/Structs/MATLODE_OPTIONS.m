%% MATLODE_OPTIONS
%
% <html>
%   <div>
%       <img style="float: right" src="../../../MATLODE_LOGO.png" height="150px"></img>
%   </div>
% </html>
%
%% Syntax
%            MATLODE_OPTIONS
%  Options = MATLODE_OPTIONS('key1',value1,'key2',value2,...,'keyN',valueN)
%  Options = MATLODE_OPTIONS(Options,'key1',value1,'key2',value2,...,'keyN',valueN)
%
%% Input Parameters
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
%% Contact Information
%%
% Dr. Adrian Sandu                 | Phone: (540) 231-2193 | Email: sandu@cs.vt.edu
%%
% Tony D'Augustine                 | Phone: (540) 231-6186 | Email: adaug13@vt.edu 
%%
% Computational Science Laboratory | Phone: (540) 231-6186
%% Reference
% [1] Tony D'Augustine, Adrian Sandu. MATLODE: A MATLAB ODE Solver and
%     Sensitivity Analysis Toolbox. Submitted to ACM TOMS.
%
%%
%  Authored by Tony D'Augustine, Adrian Sandu, and Hong Zhang.
%  Computational Science Laboratory, Virginia Tech.
%  ©2015 Virginia Tech Intellectual Properties, Inc.
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
    FacSafe_str         = '        FacSafe: [ Factor by which the new step is slightly smaller than the predicted value. ]\n';
    FDAprox_str         = '        FDAprox: [ Determines whether Jacobian Vector products by finite difference is used. ]\n';
    Gustafsson_str      = '     Gustafsson: [ An alternative error controller approach which may be advantageous depending on the model characteristics. ]\n';
    Hess_vec_str        = '       Hess_vec: [ User defined function: H * v ]\n';
    Hesstr_vec_str      = '     Hesstr_vec: [ User defined function: H^T * v ]\n';
    Hesstr_vec_f_py_str = 'Hesstr_vec_f_py: [ User defined function: (d(f_p^T * u)/dy) * k ]\n';
    Hesstr_vec_r_str    = '   Hesstr_vec_r: [ User defined function: (d(r_y^T * u)/dy) * k ]\n';
    Hesstr_vec_r_py_str = 'Hesstr_vec_r_py: [ User defined function: (d(r_p^T *u)/dy) *k ]\n';
    Hmax_str            = '           Hmax: [ Step size upper bound. ]\n';
    Hmin_str            = '           Hmin: [ Step size lower bound. ]\n';
    Hstart_str          = '         Hstart: [ Initial step size. ]\n';
    ITOL_str            = '           ITOL: [ Depreciated: Determines whether tolerances are vector or scalar ]\n';
    Jacobian_str        = '       Jacobian: [ User defined function: Jacobian ]\n';
    Jacp_str            = '           Jacp: [ User defined function: df/dp ]\n';
    Lambda_str          = '         Lambda: [ Adjoint sensitivity matrix ]\n';
    MatrixFree_str      = '     MatrixFree: [ Determines whether Jacobian is approximated ]\n';
    Max_no_steps_str    = '   Max_no_steps: [ Maximum number of steps upper bound ]\n';
    Method_str          = '         Method: [ Determines which coefficients to use ]\n';
    Mu_str              = '             Mu: [ Mu vector for sensitivity analysis ]\n';
    NADJ_str            = '           NADJ: [ Depreciated: The number of cost functionals for which adjoints are evaluated simultaneously ]\n';
    NBasisVectors       = '  NBasisVectors: [ Number of basis vectors ] \n';
    NewtonMaxit_str     = '    NewtonMaxIt: [ Maximum number of newton iterations performed ]\n';
    NewtonTol_str       = '      NewtonTol: [ Newtons method stopping criterion. ]\n';
    NP_str              = '             NP: [ Number of parameters ]\n';
    NTLM_str            = '           NTLM: [ Depreciated: The number of tangent linear sensitivity coefficients ]\n';
    QFun_str            = '           QFun: [ User defined function: r function ]\n';
    Qmax_str            = '           Qmax: [ Predicted step size to current step size upper bound ratio. ]\n';
    Qmin_str            = '           Qmin: [ Predicted step size to current step size lower bound ratio. ]\n';
    Quadrature_str      = '     Quadrature: [ Initial quadrature ]\n';
    RelTol_str          = '         RelTol: [ Relative error tolerance ]\n';
    RelTol_ADJ_str      = '     RelTol_ADJ: [ Relative Newton iteration tolerance for solving adjoint system. ]\n';
    RelTol_TLM_str      = '     RelTol_TLM: [ Relative error estimation for TLM at Newton stages ]\n';
    SaveLU_str          = '         SaveLU: [ Determines whether to save during LU factorization ]\n';
    SdirkError_str      = '     SdirkError: [ Alternative error criterion ]\n';
    StartNewton_str     = '    StartNewton: [ Determines whether newton iterations are performed ]\n';    
    storeCheckpoint_str = 'storeCheckpoint: [ Determines whether intermediate values are stored ]\n';
    ThetaMin_str        = '       ThetaMin: [ Factor deciding whether the Jacobian should be recomputed. ]\n';
    TLMNewtonEst_str    = '   TLMNewtonEst: [ Determines whether to user a tangent linear scaling factor in newton interation ]\n';
    TLMtruncErr_str     = '    TLMtruncErr: [ Determiens whether to incorpate sensitivity truncation error ]\n';
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
    fprintf( NBasisVectors );
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

Names = [
    'AbsTol           '
    'AbsTol_ADJ       '
    'AbsTol_TLM       '
    'Autonomous       '
    'ChunkSize        '
    'Desired_Mode     '
    'DirectADJ        '
    'DirectTLM        '
    'displayStats     '
    'displaySteps     '
    'DRDP             '
    'DRDY             '
    'FacMax           '
    'FacMin           '
    'FacRej           '
    'FacSafe          '
    'FDAprox          '
    'FDIncrement      '
    'GMRES_TOL        '
    'Gustafsson       '
    'Hess_vec         '
    'Hesstr_vec       '
    'Hesstr_vec_f_py  '
    'Hesstr_vec_r     '     
    'Hesstr_vec_r_py  '
    'Hmax             '
    'Hmin             '
    'Hstart           '
    'Jacp             '
    'TLMNewtonEst     '
    'ITOL             '
    'Jacobian         '
    'Lambda           '
    'MatrixFree       '
    'Max_no_steps     '
    'Method           '
    'Mu               '
    'NADJ             '
    'NBasisVectors    '
    'NewtonMaxIt      '
    'NewtonTol        '
    'NP               '
    'NTLM             '
    'QFun             '
    'Qmax             '
    'Qmin             '
    'Quadrature       '
    'RelTol           '
    'RelTol_ADJ       '
    'RelTol_TLM       '
    'SaveLU           '
    'SdirkError       '
    'OneStepIntegrator'
    'storeCheckpoint  '
    'StartNewton      '
    'ThetaMin         '
    'TLMTruncErr      '
    'WarningConfig    '
    'Y_TLM            '
    'LU               '
    ];
m = size(Names,1);
names = lower(Names);

% Combine all leading options structures o1, o2, ... in odeset(o1,o2,...).
options = [];
for j = 1:m
  options.(deblank(Names(j,:))) = [];
end

if ( ~isempty(varargin) )
    if ( isa(varargin{1},'struct') )
        if ( rem(nargin,2) ~= 1 )
            error('Invalid input arguments. See help MATLODE_OPTIONS.');
        end

        options = varargin{1};

        for i=2:2:nargin
            arg = varargin{i};
            val = varargin{i+1};

            lowArg = lower(arg);
            j = strmatch(lowArg,names,'exact');
            if ~isempty(j)
                options.(deblank(Names(j,:))) = val;
            end

        end    

    else
        if ( rem(nargin,2) ~= 0 )
            error('Invalid input arguments. See help MATLODE_OPTIONS.');
        end

        for i=1:2:nargin
            arg = varargin{i};
            val = varargin{i+1};

            lowArg = lower(arg);
            j = strmatch(lowArg,names,'exact');
            if ~isempty(j)
                options.(deblank(Names(j,:))) = val;
            end
        end


    end
end

options = orderfields(options);


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
%       <td>adaug13@vt.edu</td>
%       <td>Release MATLODE_v2.0.00</td>
%   </tr>
% </table>
% </html>
% 
%%
% <html>
%   <div>
%       <img style="float: right" src="../../../CSL_LogoWithName_1.png" height="50px"></img>
%   </div>
% </html>

