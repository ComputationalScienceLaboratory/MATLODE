%% OPTIONS_Configuration
%
% <html>
%   <div>
%       <img style="float: right" src="../../../MATLODE_LOGO.png" height="150px"></img>
%   </div>
% </html>
%%
%% Reference
% [1] Tony D'Augustine, Adrian Sandu. MATLODE: A MATLAB ODE Solver and
%     Sensitivity Analysis Toolbox. Submitted to ACM TOMS.
%
%%
%  Authored by Tony D'Augustine, Adrian Sandu, and Hong Zhang.
%  Computational Science Laboratory, Virginia Tech.
%  ©2015 Virginia Tech Intellectual Properties, Inc.
%
function [ OPTIONS, Coefficient ] = OPTIONS_Configuration( OPTIONS_U, family, implementation, y0, tspan  )

    %% Toggle Warning Configuration
    % To toggle user supplied warning system, toggle option parameter
    % 'WarningConfig'.
    %
    if ( OPTIONS_U.WarningConfig == false )
        warning('off','MatlODE:configuration');
    else 
        warning('on','MatlODE:configuration');
    end

    %% User Supplied Options and Fine Tuning
    % Each integrator is combed to guarentee applicable option settings are
    % used and available by the user. Then required settings are fined
    % tuned in prepartion for merging to the official option setting. 
    switch ( implementation )
        case 'FWD'
            switch ( family )
                case 'ERK'
                    [ OPTIONS_U ] = OPTIONS_UserSupplied_ERK_FWD( OPTIONS_U );
                    [ OPTIONS ]   = OPTIONS_FineTuned_ERK_FWD;
                case 'EXPK'
                    [ OPTIONS_U ] = OPTIONS_UserSupplied_EXPK_FWD( OPTIONS_U );
                    [ OPTIONS ]   = OPTIONS_FineTuned_EXPK_FWD;
                case 'EXP'
                    [ OPTIONS_U ] = OPTIONS_UserSupplied_EXP_FWD( OPTIONS_U );
                    [ OPTIONS ]   = OPTIONS_FineTuned_EXP_FWD;
                case 'RK'
                    [ OPTIONS_U ] = OPTIONS_UserSupplied_RK_FWD( OPTIONS_U );
                    [ OPTIONS ]   = OPTIONS_FineTuned_RK_FWD;
                case 'ROK'
                    [ OPTIONS_U ] = OPTIONS_UserSupplied_ROK_FWD( OPTIONS_U );
                    [ OPTIONS ]   = OPTIONS_FineTuned_ROK_FWD;
                case 'ROS'
                    [ OPTIONS_U ] = OPTIONS_UserSupplied_ROS_FWD( OPTIONS_U );
                    [ OPTIONS ]   = OPTIONS_FineTuned_ROS_FWD;
                case 'SDIRK'
                    [ OPTIONS_U ] = OPTIONS_UserSupplied_SDIRK_FWD( OPTIONS_U );
                    [ OPTIONS ]   = OPTIONS_FineTuned_SDIRK_FWD;
                otherwise
                    error('Internal Error: Family parameter is invalid.');
            end
        case 'TLM'
            switch ( family )
                case 'ERK'
                    [ OPTIONS_U ] = OPTIONS_UserSupplied_ERK_TLM( OPTIONS_U );
                    [ OPTIONS ]   = OPTIONS_FineTuned_ERK_TLM;
                case 'RK'
                    [ OPTIONS_U ] = OPTIONS_UserSupplied_RK_TLM( OPTIONS_U );
                    [ OPTIONS ]   = OPTIONS_FineTuned_RK_TLM;                    
                case 'ROS'
                    [ OPTIONS_U ] = OPTIONS_UserSupplied_ROS_TLM( OPTIONS_U );
                    [ OPTIONS ]   = OPTIONS_FineTuned_ROS_TLM;                    
                case 'SDIRK'
                    [ OPTIONS_U ] = OPTIONS_UserSupplied_SDIRK_TLM( OPTIONS_U );
                    [ OPTIONS ]   = OPTIONS_FineTuned_SDIRK_TLM;                    
                otherwise
                    error('Internal Error: Family parameter is invalid.');                    
            end 
        case 'ADJ'
            switch ( family )
                case 'ERK'
                    [ OPTIONS_U ] = OPTIONS_UserSupplied_ERK_ADJ( OPTIONS_U );
                    [ OPTIONS ]   = OPTIONS_FineTuned_ERK_ADJ;                    
                case 'RK'
                    [ OPTIONS_U ] = OPTIONS_UserSupplied_RK_ADJ( OPTIONS_U );
                    [ OPTIONS ]   = OPTIONS_FineTuned_RK_ADJ;                                        
                case 'ROS'
                    [ OPTIONS_U ] = OPTIONS_UserSupplied_ROS_ADJ( OPTIONS_U );
                    [ OPTIONS ]   = OPTIONS_FineTuned_ROS_ADJ;                                        
                case 'SDIRK'
                    [ OPTIONS_U ] = OPTIONS_UserSupplied_SDIRK_ADJ( OPTIONS_U );
                    [ OPTIONS ]   = OPTIONS_FineTuned_SDIRK_ADJ;                                        
                otherwise
                    error('Internal Error: Family parameter is invalid.');                    
            end 
        otherwise
            error('Internal Error: Implementation parameter is invalid.');            
    end
    
    %% Merge and General Configuration
    % All integrator user option settings must merge with fine tuning. Once
    % merged, default settings configuration will run to guarentee uniform
    % default configurations across all integrators.
    [ OPTIONS ]               = OPTIONS_Merge(OPTIONS_U,OPTIONS);
    [ OPTIONS, Coefficient ]  = OPTIONS_GeneralConfiguration(OPTIONS,family,implementation,y0,tspan);    
    %% Flags 
   
    OPTIONS.Family         = family;
    OPTIONS.Implementation = implementation;
    OPTIONS.Events         = OPTIONS_U.Events;
    
    
    % TODO: LATER OPTION.Method should be a string not a number
end

%%
% <html>
%   <div>
%       <img style="float: right" src="../../../CSL_LogoWithName_1.png" height="50px"></img>
%   </div>
% </html>