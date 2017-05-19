% MATLODE
%: A MATLAB ODE Solver and Sensitivity Analysis Toolbox
% Commands for forward integration and sensitivity analysis
%
% Forward Integrators
%   MATLODE_ERK_FWD_Integrator - Explicit Runge-Kutta forward integrator
%   MATLODE_RK_FWD_Integrator - Implicit Runge-Kutta forward integrator
%   MATLODE_ROS_FWD_Integrator - Rosenbrock forward integrator
%   MATLODE_SDIRK_FWD_Integrator - Singly Diagonally Implicit Runge-Kutta integrator
%
% Tangent Linear Sensitivity Analysis
%   MATLODE_ERK_TLM_Integrator - Explicit Runge-Kutta tangent linear integrator
%   MATLODE_RK_TLM_Integrator - Implicit Runge-Kutta tangent linear integrator
%   MATLODE_ROS_TLM_Integrator - Rosenbrock tangent linear integrator
%   MATLODE_SDIRK_TLM_Integrator - Singly Diagonally Implicit Runge-Kutta integrator
%
% Adjoint Sensitivity Analysis 
%   MATLODE_ERK_ADJ_Integrator - Explicit Runge-Kutta adjoint integrator
%   MATLODE_RK_ADJ_Integrator - Implicit Runge-Kutta adjoint integrator
%   MATLODE_ROS_ADJ_Integrator - Rosenbrock adjoint integrator
%   MATLODE_SDIRK_ADJ_Integrator - Singly Diagonally Implicit Runge-Kutta adjoint integrator
%
% Option Struct
%   MATLODE_OPTIONS - Option struct for fine tuning all integrators
%
% Execute 'doc -classic MATLODE' for additional documentation.
%
% Copyright 2015 Computational Science Laboratory
