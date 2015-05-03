%% User Guide
%
% <html>
%   <div>
%       <img style="float: right" src="../MATLODE_LOGO.png" height="150px"></img>
%   </div>
% </html>
%
% <html>
% Up: <a href="../index.html">MATLODE Toolbox</a>
% </html>
%
% <html>
% <table cellspacing="0" style="max-width:90%; min-width=80%">
%   <tr>
%     <th style="text-align: center; background-color:#D3D3D3">ODE Solver</th>
%     <th style="text-align: center; background-color:#D3D3D3">Problem Characteristics</th>
%     <th style="text-align: center; background-color:#D3D3D3">Forward Solution</th>
%     <th style="text-align: center; background-color:#D3D3D3">Tangent Linear Sensitivity</th>
%     <th style="text-align: center; background-color:#D3D3D3">Adjoint Sensitivity</th>
%     <th style="text-align: center; background-color:#D3D3D3">Method</th>
%   </tr>
%   <tr>
%     <td><a href="../Driver/Driver_ADJ/html/MATLODE_ERK_ADJ_Integrator.html">MATLODE_ERK_ADJ_Integrator</a></td>
%     <td>Nonstiff differential equations</td>
%     <td>&#x2713</td>
%     <td>&#x2717</td>
%     <td>&#x2713</td>
%     <td>Runge-Kutta</td>
%   </tr>
%   <tr>
%     <td><a href="../Driver/Driver_ADJ/html/MATLODE_RK_ADJ_Integrator.html">MATLODE_RK_ADJ_Integrator</a></td>
%     <td>Stiff differential equations</td>
%     <td>&#x2713</td>
%     <td>&#x2717</td>
%     <td>&#x2713</td>
%     <td>Runge-Kutta</td>
%   </tr>
%   <tr>
%     <td><a href="../Driver/Driver_ADJ/html/MATLODE_ROS_ADJ_Integrator.html">MATLODE_ROS_ADJ_Integrator</a></td>
%     <td>Stiff differential equations</td>
%     <td>&#x2713</td>
%     <td>&#x2717</td>
%     <td>&#x2713</td>
%     <td>Rosenbrock</td>
%   </tr>
%   <tr>
%     <td><a href="../Driver/Driver_ADJ/html/MATLODE_SDIRK_ADJ_Integrator.html">MATLODE_SDIRK_ADJ_Integrator</a></td>
%     <td>Stiff differential equations</td>
%     <td>&#x2713</td>
%     <td>&#x2717</td>
%     <td>&#x2713</td>
%     <td>Runge-Kutta</td>
%   </tr>
%   <tr>
%     <td><a href="../Driver/Driver_FWD/html/MATLODE_ERK_FWD_Integrator.html">MATLODE_ERK_FWD_Integrator</a></td>
%     <td>Nonstiff differential equations</td>
%     <td>&#x2713</td>
%     <td>&#x2717</td>
%     <td>&#x2717</td>
%     <td>Runge-Kutta</td>
%   </tr>
%   <tr>
%     <td><a href="../Driver/Driver_FWD/html/MATLODE_RK_FWD_Integrator.html">MATLODE_RK_FWD_Integrator</a></td>
%     <td>Stiff differential equations</td>
%     <td>&#x2713</td>
%     <td>&#x2717</td>
%     <td>&#x2717</td>
%     <td>Runge-Kutta</td>
%   </tr>
%   <tr>
%     <td><a href="../Driver/Driver_FWD/html/MATLODE_ROS_FWD_Integrator.html">MATLODE_ROS_FWD_Integrator</a></td>
%     <td>Stiff differential equations</td>
%     <td>&#x2713</td>
%     <td>&#x2717</td>
%     <td>&#x2717</td>
%     <td>Rosenbrock</td>
%   </tr>
%   <tr>
%     <td><a href="../Driver/Driver_FWD/html/MATLODE_SDIRK_FWD_Integrator.html">MATLODE_SDIRK_FWD_Integrator</a></td>
%     <td>Stiff differential equations</td>
%     <td>&#x2713</td>
%     <td>&#x2717</td>
%     <td>&#x2717</td>
%     <td>Runge-Kutta</td>
%   </tr>
%   <tr>
%     <td><a href="../Driver/Driver_TLM/html/MATLODE_ERK_TLM_Integrator.html">MATLODE_ERK_TLM_Integrator</a></td>
%     <td>Nonstiff differential equations</td>
%     <td>&#x2713</td>
%     <td>&#x2713</td>
%     <td>&#x2717</td>
%     <td>Runge-Kutta</td>
%   </tr>
%   <tr>
%     <td><a href="../Driver/Driver_TLM/html/MATLODE_RK_TLM_Integrator.html">MATLODE_RK_TLM_Integrator</a></td>
%     <td>Stiff differential equations</td>
%     <td>&#x2713</td>
%     <td>&#x2713</td>
%     <td>&#x2717</td>
%     <td>Runge-Kutta</td>
%   </tr>
%   <tr>
%     <td><a href="../Driver/Driver_TLM/html/MATLODE_ROS_TLM_Integrator.html">MATLODE_ROS_TLM_Integrator</a></td>
%     <td>Stiff differential equations</td>
%     <td>&#x2713</td>
%     <td>&#x2713</td>
%     <td>&#x2717</td>
%     <td>Rosenbrock</td>
%   </tr>
%   <tr>
%     <td><a href="../Driver/Driver_TLM/html/MATLODE_SDIRK_TLM_Integrator.html">MATLODE_SDIRK_TLM_Integrator</a></td>
%     <td>Stiff differential equations</td>
%     <td>&#x2713</td>
%     <td>&#x2713</td>
%     <td>&#x2717</td>
%     <td>Runge-Kutta</td>
%   </tr>
% </table> 
% </html>
%
%% Methods
% <html>
% <table cellspacing="50">
%   <tr>
%   <td style="border:none">
%   <h3>Method: Explicit Runge-Kutta (ERK)</h3>
%   <table cellspacing="0">
%   <tr>
%      <th style="text-align: center; background-color:#D3D3D3">Value</th>
%      <th style="text-align: center; background-color:#D3D3D3">Configuration</th>
%      <th style="text-align: center; background-color:#D3D3D3">Stages</th>
%      <th style="text-align: center; background-color:#D3D3D3">Order</th>
%      <th style="text-align: center; background-color:#D3D3D3">Stability Properties</th>
%   </tr>
%   <tr>
%      <td>0 (default)</td>
%      <td>Dopri5</td>
%      <td>7</td>
%      <td>5</td>
%      <td>conditionally-stable</td>
%   </tr>
%   <tr>
%      <td>1</td>
%      <td>Erk23</td>
%      <td>3</td>
%      <td>2</td>
%      <td>conditionally-stable</td>
%   </tr>
%   <tr>
%      <td>2</td>
%      <td>Erk3 Heun</td>
%      <td>4</td>
%      <td>3</td>
%      <td>conditionally-stable</td>
%   </tr>
%   <tr>
%      <td>3</td>
%      <td>Erk43</td>
%      <td>5</td>
%      <td>4</td>
%      <td>conditionally-stable</td>
%   </tr>
%   <tr>
%      <td>4</td>
%      <td>Dopri5</td>
%      <td>7</td>
%      <td>5</td>
%      <td>conditionally-stable</td>
%   </tr>
%   <tr>
%      <td>5</td>
%      <td>Dopri853</td>
%      <td>12</td>
%      <td>8</td>
%      <td>conditionally-stable</td>
%   </tr>
%   </table>
%   </td>
%   <td style="border:none">
%   <h3>Method: Explicit Runge-Kutta (ERK)</h3>
%   <table cellspacing="0">
%   <tr>
%      <th style="text-align: center; background-color:#D3D3D3">Value</th>
%      <th style="text-align: center; background-color:#D3D3D3">Configuration</th>
%      <th style="text-align: center; background-color:#D3D3D3">Stages</th>
%      <th style="text-align: center; background-color:#D3D3D3">Order</th>
%      <th style="text-align: center; background-color:#D3D3D3">Stability Properties</th>
%   </tr>
%   <tr>
%      <td>0 (default)</td>
%      <td>Dopri5</td>
%      <td>7</td>
%      <td>5</td>
%      <td>conditionally-stable</td>
%   </tr>
%   <tr>
%      <td>1</td>
%      <td>Erk23</td>
%      <td>3</td>
%      <td>2</td>
%      <td>conditionally-stable</td>
%   </tr>
%   <tr>
%      <td>2</td>
%      <td>Erk3 Heun</td>
%      <td>4</td>
%      <td>3</td>
%      <td>conditionally-stable</td>
%   </tr>
%   <tr>
%      <td>3</td>
%      <td>Erk43</td>
%      <td>5</td>
%      <td>4</td>
%      <td>conditionally-stable</td>
%   </tr>
%   <tr>
%      <td>4</td>
%      <td>Dopri5</td>
%      <td>7</td>
%      <td>5</td>
%      <td>conditionally-stable</td>
%   </tr>
%   <tr>
%      <td>5</td>
%      <td>Dopri853</td>
%      <td>12</td>
%      <td>8</td>
%      <td>conditionally-stable</td>
%   </tr>
%   </table>
%   </td>
% <tr>
% </table>
% </html>
%% Forward Integrators
% <html>
%   <h3> MATLODE_OPTIONS: Forward Integrator
%   <table cellspacing="0">
%      <tr>
%         <th style="text-align: center; background-color:#D3D3D3">Key</th>
%         <th style="text-align: center; background-color:#D3D3D3">Value</th>
%         <th style="text-align: center; background-color:#D3D3D3" >Description</th>
%         <th style="text-align: center; background-color:#D3D3D3">MATLODE_ERK_FWD_Integrator</th>
%         <th style="text-align: center; background-color:#D3D3D3">MATLODE_RK_FWD_Integrator</th>
%         <th style="text-align: center; background-color:#D3D3D3">MATLODE_ROS_FWD_Integrator</th>
%         <th style="text-align: center; background-color:#D3D3D3">MATLODE_SDIRK_FWD_Integrator</th>
%      </tr>
%      <tr><td>AbsTol</td><td>double, double[]</td><td>Absolute error tolerance for Forward integrators</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td></tr>
%      <tr><td>AbsTol_ADJ</td><td>double, double[]</td><td>Absolute Newton iteration tolerance for solving adjoint system</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td></tr>
%      <tr><td>AbsTol_TLM</td><td>double, double[]</td><td>Absolute error estimation for TLM at Newton stages</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td></tr>
%      <tr><td>ChunkSize</td><td>integer<td>Appended memory block size</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td></tr>
%      <tr><td>DirectADJ</td><td>boolean<td>Determines whether direct adjoint sensitivity analysis is performed</td><td></td><td></td><td></td><td></td></tr>
%      <tr><td>DirectTLM</td><td>boolean<td>Determines whether direct tangent linear sensitivity analysis is performed</td><td></td><td></td><td></td><td></td></tr>
%      <tr><td>DisplayStats</td><td>boolean<td>Determines whether direct tangent linear sensitivity analysis is performed</td><td></td><td></td><td></td><td></td></tr>
%      <tr><td>DisplaySteps</td><td>boolean<td>Determines whether direct tangent linear sensitivity analysis is performed</td><td></td><td></td><td></td><td></td></tr>
%      <tr><td>DRDP</td><td>function handle<td>Derivitive of r w.r.t. parameters</td><td></td><td></td><td></td><td></td></tr>
%      <tr><td>DRDY</td><td>function handle<td>Derivitive of r w.r.t. y vector</td><td></td><td></td><td></td><td></td></tr>
%   </table>
% </html>
%
%% Tangent Linear Integrators
%
%% Adjoint Integrators
%

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
%       <td>Release MATLODE_v2.0.00</td>
%   </tr>
% </table>
% </html>
%
%%
% <html>
%   <div>
%       <img style="float: right" src="../CSL_LogoWithName_1.png" height="50px"></img>
%   </div>
% </html>
%
%%
% Copyright 2015 Computational Science Laboratory
