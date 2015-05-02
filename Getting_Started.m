%% Getting Started
%
% <html>
% Up: <a href="../index.html">MATLODE Toolbox</a>
% </html>
%% Installation (Recommended)
% # Download MATLODE.mltbx
% # Navigate to MATLODE.mltbx in directory
% # Right click on MATLODE.mltbx
% # Click install
%
%% Alternative Installation
% # Download MATLODE.tar.gz
% # Navigate to MATLODE.tar.gz
% # Extract all files
% # Right click on root folder MATLODE
% # Click Add to Path -> Selected Folders and Subfolders
% 
%% Uninstall (Recommended)
% # Click Add-Ons -> Manage Custom Toolboxes -> Uninstall MATLODE
%
% *NOTE:* Only applies if MATLODE.tltbx is installed as a custom MATLAB
% toolbox.
%
%% Alternative Uninstall Method 1
% # Naviage to MATLAB/Toolboxes
% # Delete MATLODE
% # Delete all MATLODE/* paths from MATLAB/Toolboxes/toolboxFolders.txt 
%
% *NOTE:* Only applies if MATLODE.tltbx is installed as a custom MATLAB
% toolbox.
%
%% Alternative Uninstall Method 2
% # Naviage to MATLODE root folder
% # Right click MATLODE root folder
% # Click Remove from Path -> Selected Folders and Subfolders
%
% *NOTE:* Only applies if MATLODE is installed by adding MATLODE to the
% current matlabpath.
%
%% Troubleshooting
% Included in MATLODE are example scripts for every integrator. To
% troubleshoot a specific integrator, execute the appropriate script listed in 
% the command window. To edit the script, execut 'edit scriptFileName' in
% the command window.
%
% <html>
%   <h3>All Forward Integrators</h3>
%   <ul>
%      <li><a href="../Testing_Material/html/Testing_Simple_FWD.html">Testing_Simple_FWD</a></li>
%   </ul>
%   <h3>All Tangent Linear Integrators</h3>
%   <ul>
%      <li><a href="../Testing_Material/html/Testing_Simple_TLM.html">Testing_Simple_TLM</a></li>
%   </ul>
%   <h3>All Adjoint Integrators</h3>
%   <ul>
%      <li><a href="../Testing_Material/html/Testing_Simple_ADJ.html">Testing_Simple_ADJ</a></li>
%   </ul>
% </html>
%
%% Simple Example Problem
% Execute the following commands into the command window to run Van Der Pol using
% a forward explicit Runge-Kutta intregration scheme.
Ode_Function        = @arenstorfOrbit_Function;
Time_Interval       = [ 0 17.0652166 ];
Y0                  = [0.994; 0; 0; -2.00158510637908252240537862224];

[ ~, Y ] = MATLODE_ERK_FWD_Integrator(Ode_Function,Time_Interval,Y0)

%%
% Copyright 2015 Computational Science Laboratory