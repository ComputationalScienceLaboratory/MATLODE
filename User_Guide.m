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
%%
% <html>
% <table cellspacing="0" style="max-width:90%; min-width=80%">
%   <tr>
%     <th style="text-align: center; background-color:#D3D3D3">ODE Solver</th>
%     <th style="text-align: center; background-color:#D3D3D3">Example</th>
%     <th style="text-align: center; background-color:#D3D3D3">Problem Characteristics</th>
%     <th style="text-align: center; background-color:#D3D3D3">Forward Solution</th>
%     <th style="text-align: center; background-color:#D3D3D3">Tangent Linear Sensitivity</th>
%     <th style="text-align: center; background-color:#D3D3D3">Adjoint Sensitivity</th>
%     <th style="text-align: center; background-color:#D3D3D3">Method</th>
%   </tr>
%   <tr>
%     <td><a href="../Driver/Driver_ADJ/html/MATLODE_ERK_ADJ_Integrator.html">MATLODE_ERK_ADJ_Integrator</a></td>
%     <td><a href="MATLODE_Example_ERK_ADJ_Integrator.html">MATLODE_Example_ERK_ADJ_Integrator</a></td>
%     <td>Nonstiff differential equations</td>
%     <td>&#x2713</td>
%     <td>&#x2717</td>
%     <td>&#x2713</td>
%     <td>Runge-Kutta</td>
%   </tr>
%   <tr>
%     <td><a href="../Driver/Driver_ADJ/html/MATLODE_RK_ADJ_Integrator.html">MATLODE_RK_ADJ_Integrator</a></td>
%     <td><a href="MATLODE_Example_RK_ADJ_Integrator.html">MATLODE_Example_RK_ADJ_Integrator</a></td>
%     <td>Stiff differential equations</td>
%     <td>&#x2713</td>
%     <td>&#x2717</td>
%     <td>&#x2713</td>
%     <td>Runge-Kutta</td>
%   </tr>
%   <tr>
%     <td><a href="../Driver/Driver_ADJ/html/MATLODE_ROS_ADJ_Integrator.html">MATLODE_ROS_ADJ_Integrator</a></td>
%     <td><a href="MATLODE_Example_ROS_ADJ_Integrator.html">MATLODE_Example_ROS_ADJ_Integrator</a></td>
%     <td>Stiff differential equations</td>
%     <td>&#x2713</td>
%     <td>&#x2717</td>
%     <td>&#x2713</td>
%     <td>Rosenbrock</td>
%   </tr>
%   <tr>
%     <td><a href="../Driver/Driver_ADJ/html/MATLODE_SDIRK_ADJ_Integrator.html">MATLODE_SDIRK_ADJ_Integrator</a></td>
%     <td><a href="MATLODE_Example_SDIRK_ADJ_Integrator.html">MATLODE_Example_SDIRK_ADJ_Integrator</a></td>
%     <td>Stiff differential equations</td>
%     <td>&#x2713</td>
%     <td>&#x2717</td>
%     <td>&#x2713</td>
%     <td>Runge-Kutta</td>
%   </tr>
%   <tr>
%     <td><a href="../Driver/Driver_FWD/html/MATLODE_ERK_FWD_Integrator.html">MATLODE_ERK_FWD_Integrator</a></td>
%     <td><a href="MATLODE_Example_ERK_FWD_Integrator.html">MATLODE_Example_ERK_FWD_Integrator</a></td>
%     <td>Nonstiff differential equations</td>
%     <td>&#x2713</td>
%     <td>&#x2717</td>
%     <td>&#x2717</td>
%     <td>Runge-Kutta</td>
%   </tr>
%   <tr>
%     <td><a href="../Driver/Driver_FWD/html/MATLODE_RK_FWD_Integrator.html">MATLODE_RK_FWD_Integrator</a></td>
%     <td><a href="MATLODE_Example_RK_FWD_Integrator.html">MATLODE_Example_RK_FWD_Integrator</a></td>
%     <td>Stiff differential equations</td>
%     <td>&#x2713</td>
%     <td>&#x2717</td>
%     <td>&#x2717</td>
%     <td>Runge-Kutta</td>
%   </tr>
%   <tr>
%     <td><a href="../Driver/Driver_FWD/html/MATLODE_ROS_FWD_Integrator.html">MATLODE_ROS_FWD_Integrator</a></td>
%     <td><a href="MATLODE_Example_ROS_FWD_Integrator.html">MATLODE_Example_ROS_FWD_Integrator</a></td>
%     <td>Stiff differential equations</td>
%     <td>&#x2713</td>
%     <td>&#x2717</td>
%     <td>&#x2717</td>
%     <td>Rosenbrock</td>
%   </tr>
%   <tr>
%     <td><a
%     href="../Driver/Driver_FWD/html/MATLODE_SDIRK_FWD_Integrator.html">MATLODE_SDIRK_FWD_Integrator</a></td>
%     <td><a href="MATLODE_Example_SDIRK_FWD_Integrator.html">MATLODE_Example_SDIRK_FWD_Integrator</a></td>
%     <td>Stiff differential equations</td>
%     <td>&#x2713</td>
%     <td>&#x2717</td>
%     <td>&#x2717</td>
%     <td>Runge-Kutta</td>
%   </tr>
%   <tr>
%     <td><a href="../Driver/Driver_TLM/html/MATLODE_ERK_TLM_Integrator.html">MATLODE_ERK_TLM_Integrator</a></td>
%     <td><a href="MATLODE_Example_ERK_TLM_Integrator.html">MATLODE_Example_ERK_TLM_Integrator</a></td>
%     <td>Nonstiff differential equations</td>
%     <td>&#x2713</td>
%     <td>&#x2713</td>
%     <td>&#x2717</td>
%     <td>Runge-Kutta</td>
%   </tr>
%   <tr>
%     <td><a href="../Driver/Driver_TLM/html/MATLODE_RK_TLM_Integrator.html">MATLODE_RK_TLM_Integrator</a></td>
%     <td><a href="MATLODE_Example_RK_TLM_Integrator.html">MATLODE_Example_RK_TLM_Integrator</a></td>
%     <td>Stiff differential equations</td>
%     <td>&#x2713</td>
%     <td>&#x2713</td>
%     <td>&#x2717</td>
%     <td>Runge-Kutta</td>
%   </tr>
%   <tr>
%     <td><a href="../Driver/Driver_TLM/html/MATLODE_ROS_TLM_Integrator.html">MATLODE_ROS_TLM_Integrator</a></td>
%     <td><a href="MATLODE_Example_ROS_TLM_Integrator.html">MATLODE_Example_ROS_TLM_Integrator</a></td>
%     <td>Stiff differential equations</td>
%     <td>&#x2713</td>
%     <td>&#x2713</td>
%     <td>&#x2717</td>
%     <td>Rosenbrock</td>
%   </tr>
%   <tr>
%     <td><a href="../Driver/Driver_TLM/html/MATLODE_SDIRK_TLM_Integrator.html">MATLODE_SDIRK_TLM_Integrator</a></td>
%     <td><a href="MATLODE_Example_SDIRK_TLM_Integrator.html">MATLODE_Example_SDIRK_TLM_Integrator</a></td>
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
%   <h3>Method: Implicit Runge-Kutta (RK)</h3>
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
%      <td>Lobatto3C</td>
%      <td>3</td>
%      <td>4</td>
%      <td>L-stable</td>
%   </tr>
%   <tr>
%      <td>1</td>
%      <td>Radua2A</td>
%      <td>3</td>
%      <td>5</td>
%      <td>L-stable</td>
%   </tr>
%   <tr>
%      <td>2</td>
%      <td>Lobatto3C</td>
%      <td>3</td>
%      <td>4</td>
%      <td>L-stable</td>
%   </tr>
%   <tr>
%      <td>3</td>
%      <td>Gauss</td>
%      <td>3</td>
%      <td>6</td>
%      <td>weakly L-stable</td>
%   </tr>
%   <tr>
%      <td>4</td>
%      <td>Radau1A</td>
%      <td>3</td>
%      <td>5</td>
%      <td>L-stable</td>
%   </tr>
%   </table>
%   </td>
% </tr>
% <tr>
%   <td style="border:none">
%   <h3>Method: Rosenbrock (ROS)</h3>
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
%      <td>Ros4</td>
%      <td>4</td>
%      <td>4</td>
%      <td>L-stable</td>
%   </tr>
%   <tr>
%      <td>1</td>
%      <td>Ros2</td>
%      <td>2</td>
%      <td>2</td>
%      <td>L-stable</td>
%   </tr>
%   <tr>
%      <td>2</td>
%      <td>Ros3</td>
%      <td>3</td>
%      <td>3</td>
%      <td>L-stable</td>
%   </tr>
%   <tr>
%      <td>3</td>
%      <td>Ros4</td>
%      <td>4</td>
%      <td>4</td>
%      <td>L-stable</td>
%   </tr>
%   <tr>
%      <td>4</td>
%      <td>Rodas3</td>
%      <td>4</td>
%      <td>3</td>
%      <td>L-stable</td>
%   </tr>
%   <tr>
%      <td>5</td>
%      <td>Rodas4</td>
%      <td>6</td>
%      <td>4</td>
%      <td>L-stable</td>
%   </tr>
%   </table>
%   </td>
%   <td style="border:none">
%   <h3>Method: Singly Diagonally Implicit Runge-Kutta (SDIRK)</h3>
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
%      <td>Sdirk4A</td>
%      <td>5</td>
%      <td>4</td>
%      <td>L-stable</td>
%   </tr>
%   <tr>
%      <td>1</td>
%      <td>Sdrik2A</td>
%      <td>2</td>
%      <td>2</td>
%      <td>L-stable</td>
%   </tr>
%   <tr>
%      <td>2</td>
%      <td>Sdirk2B</td>
%      <td>2</td>
%      <td>2</td>
%      <td>L-stable</td>
%   </tr>
%   <tr>
%      <td>3</td>
%      <td>Sdirk3A</td>
%      <td>3</td>
%      <td>2</td>
%      <td>L-stable</td>
%   </tr>
%   <tr>
%      <td>4</td>
%      <td>Sdirk4A</td>
%      <td>5</td>
%      <td>4</td>
%      <td>L-stable</td>
%   </tr>
%   <tr>
%      <td>5</td>
%      <td>Sdirk4B</td>
%      <td>5</td>
%      <td>4</td>
%      <td>L-stable</td>
%   </tr>
%   </table>
%   </td>
% </tr>
% </table>
% <p>Return to <a href="User_Guide.html">Top</a></p>
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
%      <tr><td>DirectADJ</td><td>boolean<td>Determines whether direct adjoint sensitivity analysis is performed</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td></tr>
%      <tr><td>DirectTLM</td><td>boolean<td>Determines whether direct tangent linear sensitivity analysis is performed</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td></tr>
%      <tr><td>DisplayStats</td><td>boolean<td>Determines whether statistics are displayed</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td></tr>
%      <tr><td>DisplaySteps</td><td>boolean<td>Determines whether steps are displayedd</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td></tr>
%      <tr><td>DRDP</td><td>function handle<td>Derivitive of r w.r.t. parameters</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td></tr>
%      <tr><td>DRDY</td><td>function handle<td>Derivitive of r w.r.t. y vector</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td></tr>
%      <tr><td>FacMax</td><td>double<td>Step size upper bound change ratio</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td></tr>
%      <tr><td>FacMin</td><td>double<td>Step size lower bound change ratio</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td></tr>
%      <tr><td>FacSafe</td><td>double<td>Factor by which the new step is slightly smaller than the predicted value</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td></tr>
%      <tr><td>FDAprox</td><td>boolean<td>Determines whether Jacobian vector products by finite difference is used</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td></tr>
%      <tr><td>Gustafsson</td><td>boolean<td>An alternative error controller approach which may be advantageous depending on the model characteristics</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td></tr>
%      <tr><td>Hess_vec</td><td>function handle<td>H * v</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td></tr>
%      <tr><td>Hesstr_vec</td><td>function handle<td>H^T * v</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td></tr>
%      <tr><td>Hesstr_vec_f_py</td><td>function handle<td>(d(f_p^T * u)/dy) * k</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td></tr>
%      <tr><td>Hesstr_vec_r</td><td>function handle<td>(d(r_y^T * u)/dy) * k</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td></tr>
%      <tr><td>Hesstr_vec_r_py</td><td>function handle<td>(d(r_p^T *u)/dy) *k</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td></tr>
%      <tr><td>Hmax</td><td>double<td>Step size upper bound</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td></tr>
%      <tr><td>Hmin</td><td>double<td>Step size lower bound</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td></tr>
%      <tr><td>Hstart</td><td>double<td>Initial step size</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td></tr>
%      <tr><td>ITOL</td><td>boolean<td>Depreciated: Tolerances are scalor or vector</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td></tr>
%      <tr><td>Jacobian</td><td>function handle<td>User defined function: Jacobian</td><td>&#x2717</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td></tr>
%      <tr><td>Jacp</td><td>function handle<td>User defined function: df/dp </td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td></tr>
%      <tr><td>Lambda</td><td>double[]<td>Adjoint sensitivity matrix</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td></tr>
%      <tr><td>MatrixFree</td><td>boolean<td>Determines whether Jacobian is approximated</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td><td>&#x2713</td></tr>
%      <tr><td>Max_no_steps</td><td>integer<td>Maximum number of steps upper bound</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td></tr>
%      <tr><td>Method</td><td>integer<td>Determines which coefficients to use</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td></tr>
%      <tr><td>Mu</td><td>double[]<td>Mu vector for sensitivity analysis</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td></tr>
%      <tr><td>NBasisVectors</td><td>integer<td>Number of basis vectors</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td></tr>
%      <tr><td>NewtonMaxIt</td><td>integer<td>Maximum number of newton iterations performed</td><td>&#x2717</td><td>&#x2713</td><td>&#x2717</td><td>&#x2713</td></tr>
%      <tr><td>NewtonTol</td><td>double<td>Newton method stopping criterion</td><td>&#x2717</td><td>&#x2713</td><td>&#x2717</td><td>&#x2713</td></tr>
%      <tr><td>NP</td><td>integer<td>Number of parameters</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td></tr>
%      <tr><td>QFun</td><td>function handle<td>r function</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td></tr>
%      <tr><td>Qmax</td><td>double<td>Predicted step size to current step size upper bound ratio</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td></tr>
%      <tr><td>Qmin</td><td>double<td>Predicted step size to current step size lower bound ratio</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td></tr>
%      <tr><td>Quadrature</td><td>double[]<td>Initial Quadrature</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td></tr>
%      <tr><td>RelTol</td><td>double, double[]<td>Relative error tolerance</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td></tr>
%      <tr><td>RelTol_ADJ</td><td>double, double[]<td>Relative Newton iteration tolerance for solving adjoint system</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td></tr>
%      <tr><td>RelTol_TLM</td><td>double, double[]<td>Relative error estimation for TLM at Newton stages</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td></tr>
%      <tr><td>SaveLU</td><td>boolean<td>Determines whether to save during LU factorization</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td></tr>
%      <tr><td>SdirkError</td><td>boolean<td>Alternative error criterionn</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td></tr>
%      <tr><td>StartNewton</td><td>boolean<td>Determines whether newton iterations are performed</td><td>&#x2717</td><td>&#x2713</td><td>&#x2717</td><td>&#x2713</td></tr>
%      <tr><td>StoreCheckpoint</td><td>boolean<td>Determines whether intermediate values are stored</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td></tr>
%      <tr><td>ThetaMin</td><td>double<td>Factor deciding whether the Jacobian should be recomputed</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td></tr>
%      <tr><td>TLMNewtonEst</td><td>boolean<td>Determines whether to user a tangent linear scaling factor in newton interation</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td></tr>
%      <tr><td>TLMtruncErr</td><td>boolean<td>Determiens whether to incorpate sensitivity truncation error</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td></tr>
%      <tr><td>WarningConfig</td><td>boolean<td>Determines whether warning messages are displayed during option configuration</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td></tr>
%      <tr><td>Y_TLM</td><td>double[]<td>Contains the sensitivities of Y with respect to the specified coefficients</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td></tr>
%   </table>
% <p>Return to <a href="User_Guide.html">Top</a></p>
% </html>
%
%% Tangent Linear Integrators
% <html>
%   <h3> MATLODE_OPTIONS: Tangent Linear Integrator
%   <table cellspacing="0">
%      <tr>
%         <th style="text-align: center; background-color:#D3D3D3">Key</th>
%         <th style="text-align: center; background-color:#D3D3D3">Value</th>
%         <th style="text-align: center; background-color:#D3D3D3" >Description</th>
%         <th style="text-align: center; background-color:#D3D3D3">MATLODE_ERK_TLM_Integrator</th>
%         <th style="text-align: center; background-color:#D3D3D3">MATLODE_RK_TLM_Integrator</th>
%         <th style="text-align: center; background-color:#D3D3D3">MATLODE_ROS_TLM_Integrator</th>
%         <th style="text-align: center; background-color:#D3D3D3">MATLODE_SDIRK_TLM_Integrator</th>
%      </tr>
%      <tr><td>AbsTol</td><td>double, double[]</td><td>Absolute error tolerance for Forward integrators</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td></tr>
%      <tr><td>AbsTol_ADJ</td><td>double, double[]</td><td>Absolute Newton iteration tolerance for solving adjoint system</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td></tr>
%      <tr><td>AbsTol_TLM</td><td>double, double[]</td><td>Absolute error estimation for TLM at Newton stages</td><td>&#x2717</td><td>&#x2717</td><td>&#x2713</td><td>&#x2717</td></tr>
%      <tr><td>ChunkSize</td><td>integer<td>Appended memory block size</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td></tr>
%      <tr><td>DirectADJ</td><td>boolean<td>Determines whether direct adjoint sensitivity analysis is performed</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td></tr>
%      <tr><td>DirectTLM</td><td>boolean<td>Determines whether direct tangent linear sensitivity analysis is performed</td><td>&#x2717</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td></tr>
%      <tr><td>DisplayStats</td><td>boolean<td>Determines whether statistics are displayed</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td></tr>
%      <tr><td>DisplaySteps</td><td>boolean<td>Determines whether steps are displayedd</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td></tr>
%      <tr><td>DRDP</td><td>function handle<td>Derivitive of r w.r.t. parameters</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td></tr>
%      <tr><td>DRDY</td><td>function handle<td>Derivitive of r w.r.t. y vector</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td></tr>
%      <tr><td>FacMax</td><td>double<td>Step size upper bound change ratio</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td></tr>
%      <tr><td>FacMin</td><td>double<td>Step size lower bound change ratio</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td></tr>
%      <tr><td>FacSafe</td><td>double<td>Factor by which the new step is slightly smaller than the predicted value</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td></tr>
%      <tr><td>FDAprox</td><td>boolean<td>Determines whether Jacobian vector products by finite difference is used</td><td>&#x2713</td><td>&#x2713</td><td>&#x2717</td><td>&#x2717</td></tr>
%      <tr><td>Gustafsson</td><td>boolean<td>An alternative error controller approach which may be advantageous depending on the model characteristics</td><td>&#x2717</td><td>&#x2713</td><td>&#x2717</td><td>&#x2717</td></tr>
%      <tr><td>Hess_vec</td><td>function handle<td>H * v</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td></tr>
%      <tr><td>Hesstr_vec</td><td>function handle<td>H^T * v</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td></tr>
%      <tr><td>Hesstr_vec_f_py</td><td>function handle<td>(d(f_p^T * u)/dy) * k</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td></tr>
%      <tr><td>Hesstr_vec_r</td><td>function handle<td>(d(r_y^T * u)/dy) * k</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td></tr>
%      <tr><td>Hesstr_vec_r_py</td><td>function handle<td>(d(r_p^T *u)/dy) *k</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td></tr>
%      <tr><td>Hmax</td><td>double<td>Step size upper bound</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td></tr>
%      <tr><td>Hmin</td><td>double<td>Step size lower bound</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td></tr>
%      <tr><td>Hstart</td><td>double<td>Initial step size</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td></tr>
%      <tr><td>ITOL</td><td>boolean<td>Depreciated: Tolerances are scalor or vector</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td></tr>
%      <tr><td>Jacobian</td><td>function handle<td>User defined function: Jacobian</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td></tr>
%      <tr><td>Jacp</td><td>function handle<td>User defined function: df/dp </td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td></tr>
%      <tr><td>Lambda</td><td>double[]<td>Adjoint sensitivity matrix</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td></tr>
%      <tr><td>MatrixFree</td><td>boolean<td>Determines whether Jacobian is approximated</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td><td>&#x2713</td></tr>
%      <tr><td>Max_no_steps</td><td>integer<td>Maximum number of steps upper bound</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td></tr>
%      <tr><td>Method</td><td>integer<td>Determines which coefficients to use</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td></tr>
%      <tr><td>Mu</td><td>double[]<td>Mu vector for sensitivity analysis</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td></tr>
%      <tr><td>NBasisVectors</td><td>integer<td>Number of basis vectors</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td></tr>
%      <tr><td>NewtonMaxIt</td><td>integer<td>Maximum number of newton iterations performed</td><td>&#x2717</td><td>&#x2713</td><td>&#x2717</td><td>&#x2713</td></tr>
%      <tr><td>NewtonTol</td><td>double<td>Newton method stopping criterion</td><td>&#x2717</td><td>&#x2713</td><td>&#x2717</td><td>&#x2713</td></tr>
%      <tr><td>NP</td><td>integer<td>Number of parameters</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td></tr>
%      <tr><td>QFun</td><td>function handle<td>r function</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td></tr>
%      <tr><td>Qmax</td><td>double<td>Predicted step size to current step size upper bound ratio</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td></tr>
%      <tr><td>Qmin</td><td>double<td>Predicted step size to current step size lower bound ratio</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td></tr>
%      <tr><td>Quadrature</td><td>double[]<td>Initial Quadrature</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td></tr>
%      <tr><td>RelTol</td><td>double, double[]<td>Relative error tolerance</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td></tr>
%      <tr><td>RelTol_ADJ</td><td>double, double[]<td>Relative Newton iteration tolerance for solving adjoint system</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td></tr>
%      <tr><td>RelTol_TLM</td><td>double, double[]<td>Relative error estimation for TLM at Newton stages</td><td>&#x2717</td><td>&#x2717</td><td>&#x2713</td><td>&#x2717</td></tr>
%      <tr><td>SaveLU</td><td>boolean<td>Determines whether to save during LU factorization</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td></tr>
%      <tr><td>SdirkError</td><td>boolean<td>Alternative error criterionn</td><td>&#x2717</td><td>&#x2713</td><td>&#x2717</td><td>&#x2717</td></tr>
%      <tr><td>StartNewton</td><td>boolean<td>Determines whether newton iterations are performed</td><td>&#x2717</td><td>&#x2713</td><td>&#x2717</td><td>&#x2713</td></tr>
%      <tr><td>StoreCheckpoint</td><td>boolean<td>Determines whether intermediate values are stored</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td></tr>
%      <tr><td>ThetaMin</td><td>double<td>Factor deciding whether the Jacobian should be recomputed</td><td>&#x2717</td><td>&#x2713</td><td>&#x2717</td><td>&#x2713</td></tr>
%      <tr><td>TLMNewtonEst</td><td>boolean<td>Determines whether to user a tangent linear scaling factor in newton interation</td><td>&#x2717</td><td>&#x2713</td><td>&#x2717</td><td>&#x2713</td></tr>
%      <tr><td>TLMtruncErr</td><td>boolean<td>Determiens whether to incorpate sensitivity truncation error</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td></tr>
%      <tr><td>WarningConfig</td><td>boolean<td>Determines whether warning messages are displayed during option configuration</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td></tr>
%      <tr><td>Y_TLM</td><td>double[]<td>Contains the sensitivities of Y with respect to the specified coefficients</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td></tr>
%   </table>
% <p>Return to <a href="User_Guide.html">Top</a></p>
% </html>
%
%% Adjoint Integrators
% <html>
%   <h3> MATLODE_OPTIONS: Adjoint Integrator
%   <table cellspacing="0">
%      <tr>
%         <th style="text-align: center; background-color:#D3D3D3">Key</th>
%         <th style="text-align: center; background-color:#D3D3D3">Value</th>
%         <th style="text-align: center; background-color:#D3D3D3" >Description</th>
%         <th style="text-align: center; background-color:#D3D3D3">MATLODE_ERK_ADJ_Integrator</th>
%         <th style="text-align: center; background-color:#D3D3D3">MATLODE_RK_ADJ_Integrator</th>
%         <th style="text-align: center; background-color:#D3D3D3">MATLODE_ROS_ADJ_Integrator</th>
%         <th style="text-align: center; background-color:#D3D3D3">MATLODE_SDIRK_ADJ_Integrator</th>
%      </tr>
%      <tr><td>AbsTol</td><td>double, double[]</td><td>Absolute error tolerance for Forward integrators</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td></tr>
%      <tr><td>AbsTol_ADJ</td><td>double, double[]</td><td>Absolute Newton iteration tolerance for solving adjoint system</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td><td>&#x2713</td></tr>
%      <tr><td>AbsTol_TLM</td><td>double, double[]</td><td>Absolute error estimation for TLM at Newton stages</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td></tr>
%      <tr><td>ChunkSize</td><td>integer<td>Appended memory block size</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td></tr>
%      <tr><td>DirectADJ</td><td>boolean<td>Determines whether direct adjoint sensitivity analysis is performed</td><td>&#x2717</td><td>&#x2713</td><td>&#x2717</td><td>&#x2713</td></tr>
%      <tr><td>DirectTLM</td><td>boolean<td>Determines whether direct tangent linear sensitivity analysis is performed</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td></tr>
%      <tr><td>DisplayStats</td><td>boolean<td>Determines whether statistics are displayed</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td></tr>
%      <tr><td>DisplaySteps</td><td>boolean<td>Determines whether steps are displayedd</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td></tr>
%      <tr><td>DRDP</td><td>function handle<td>Derivitive of r w.r.t. parameters</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td></tr>
%      <tr><td>DRDY</td><td>function handle<td>Derivitive of r w.r.t. y vector</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td></tr>
%      <tr><td>FacMax</td><td>double<td>Step size upper bound change ratio</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td></tr>
%      <tr><td>FacMin</td><td>double<td>Step size lower bound change ratio</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td></tr>
%      <tr><td>FacSafe</td><td>double<td>Factor by which the new step is slightly smaller than the predicted value</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td></tr>
%      <tr><td>FDAprox</td><td>boolean<td>Determines whether Jacobian vector products by finite difference is used</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td></tr>
%      <tr><td>Gustafsson</td><td>boolean<td>An alternative error controller approach which may be advantageous depending on the model characteristics</td><td>&#x2717</td><td>&#x2713</td><td>&#x2717</td><td>&#x2717</td></tr>
%      <tr><td>Hess_vec</td><td>function handle<td>H * v</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td></tr>
%      <tr><td>Hesstr_vec</td><td>function handle<td>H^T * v</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td></tr>
%      <tr><td>Hesstr_vec_f_py</td><td>function handle<td>(d(f_p^T * u)/dy) * k</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td></tr>
%      <tr><td>Hesstr_vec_r</td><td>function handle<td>(d(r_y^T * u)/dy) * k</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td></tr>
%      <tr><td>Hesstr_vec_r_py</td><td>function handle<td>(d(r_p^T *u)/dy) *k</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td></tr>
%      <tr><td>Hmax</td><td>double<td>Step size upper bound</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td></tr>
%      <tr><td>Hmin</td><td>double<td>Step size lower bound</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td></tr>
%      <tr><td>Hstart</td><td>double<td>Initial step size</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td></tr>
%      <tr><td>ITOL</td><td>boolean<td>Depreciated: Tolerances are scalor or vector</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td></tr>
%      <tr><td>Jacobian</td><td>function handle<td>User defined function: Jacobian</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td></tr>
%      <tr><td>Jacp</td><td>function handle<td>User defined function: df/dp </td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td></tr>
%      <tr><td>Lambda</td><td>double[]<td>Adjoint sensitivity matrix</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td></tr>
%      <tr><td>MatrixFree</td><td>boolean<td>Determines whether Jacobian is approximated</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td><td>&#x2713</td></tr>
%      <tr><td>Max_no_steps</td><td>integer<td>Maximum number of steps upper bound</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td></tr>
%      <tr><td>Method</td><td>integer<td>Determines which coefficients to use</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td></tr>
%      <tr><td>Mu</td><td>double[]<td>Mu vector for sensitivity analysis</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td></tr>
%      <tr><td>NBasisVectors</td><td>integer<td>Number of basis vectors</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td></tr>
%      <tr><td>NewtonMaxIt</td><td>integer<td>Maximum number of newton iterations performed</td><td>&#x2717</td><td>&#x2713</td><td>&#x2717</td><td>&#x2713</td></tr>
%      <tr><td>NewtonTol</td><td>double<td>Newton method stopping criterion</td><td>&#x2717</td><td>&#x2713</td><td>&#x2717</td><td>&#x2713</td></tr>
%      <tr><td>NP</td><td>integer<td>Number of parameters</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td></tr>
%      <tr><td>QFun</td><td>function handle<td>r function</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td></tr>
%      <tr><td>Qmax</td><td>double<td>Predicted step size to current step size upper bound ratio</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td></tr>
%      <tr><td>Qmin</td><td>double<td>Predicted step size to current step size lower bound ratio</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td></tr>
%      <tr><td>Quadrature</td><td>double[]<td>Initial Quadrature</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td></tr>
%      <tr><td>RelTol</td><td>double, double[]<td>Relative error tolerance</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td></tr>
%      <tr><td>RelTol_ADJ</td><td>double, double[]<td>Relative Newton iteration tolerance for solving adjoint system</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td><td>&#x2713</td></tr>
%      <tr><td>RelTol_TLM</td><td>double, double[]<td>Relative error estimation for TLM at Newton stages</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td></tr>
%      <tr><td>SaveLU</td><td>boolean<td>Determines whether to save during LU factorization</td><td>&#x2717</td><td>&#x2713</td><td>&#x2717</td><td>&#x2713</td></tr>
%      <tr><td>SdirkError</td><td>boolean<td>Alternative error criterionn</td><td>&#x2717</td><td>&#x2713</td><td>&#x2717</td><td>&#x2717</td></tr>
%      <tr><td>StartNewton</td><td>boolean<td>Determines whether newton iterations are performed</td><td>&#x2717</td><td>&#x2713</td><td>&#x2717</td><td>&#x2713</td></tr>
%      <tr><td>StoreCheckpoint</td><td>boolean<td>Determines whether intermediate values are stored</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td></tr>
%      <tr><td>ThetaMin</td><td>double<td>Factor deciding whether the Jacobian should be recomputed</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td></tr>
%      <tr><td>TLMNewtonEst</td><td>boolean<td>Determines whether to user a tangent linear scaling factor in newton interation</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td></tr>
%      <tr><td>TLMtruncErr</td><td>boolean<td>Determiens whether to incorpate sensitivity truncation error</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td></tr>
%      <tr><td>WarningConfig</td><td>boolean<td>Determines whether warning messages are displayed during option configuration</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td><td>&#x2713</td></tr>
%      <tr><td>Y_TLM</td><td>double[]<td>Contains the sensitivities of Y with respect to the specified coefficients</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td><td>&#x2717</td></tr>
%   </table>
% <p>Return to <a href="User_Guide.html">Top</a></p>
% </html>
%

%% Contact Information
%%
% Dr. Adrian Sandu                 | Phone: (540) 231-2193 | Email: sandu@cs.vt.edu
%%
% Tony D'Augustine                 | Phone: (540) 231-6186 | Email: adaug13@vt.edu 
%%
% Computational Science Laboratory | Phone: (540) 231-6186

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
