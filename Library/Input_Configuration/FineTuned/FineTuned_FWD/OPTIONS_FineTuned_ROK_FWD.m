%% OPTIONS_FineTuned_ROK_FWD
%
% <html>
% Up: <a href="../../../../html/Library.html">Library</a>
% </html>
%
%% Rosenbrock Krylov (ROK) Forward (FWD) Fine Tuning
% The following is the default fine tuning for the ROK FWD integrator.
% All option settings fine tuned to 0 will use the default
% configuration settings. 
%
% *IMPORTANT:* If an option setting is required by integrator, the
% corresponding setting must be fined to 0 or an appropriate value as
% desired. 
%
% Copyright 2015 Computational Science Laboratory
function [ OPTIONS ] = OPTIONS_FineTuned_ROK_FWD
    OPTIONS = MATLODE_OPTIONS( ...
        'AbsTol',          0, ...
        'Autonomous',      0, ...
        'ChunkSize',       0, ...
        'displaySteps',    0, ...        
        'FacMax',          0, ...
        'FacMin',          0, ...
        'FacRej',          0, ...
        'FacSafe',         0, ...        
        'Hmin',            0, ...
        'Hstart',          0, ...
        'Hmax',            0, ...
        'NBasisVectors',   4, ...        
        'Max_no_steps',    0, ...
        'MatrixFree',      0, ...
        'Method',          0, ...
        'RelTol',          0, ...
        'ITOL',            0, ...
        'storeCheckpoint', 0 );
        
end

