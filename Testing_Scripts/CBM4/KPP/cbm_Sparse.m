% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%                                                                  
% Sparse Data Definition File                                      
%                                                                  
% Generated by KPP-2.2.3 symbolic chemistry Kinetics PreProcessor  
%       (http://www.cs.vt.edu/~asandu/Software/KPP)                
% KPP is distributed under GPL, the general public licence         
%       (http://www.gnu.org/copyleft/gpl.html)                     
% (C) 1995-1997, V. Damian & A. Sandu, CGRER, Univ. Iowa           
% (C) 1997-2005, A. Sandu, Michigan Tech, Virginia Tech            
%     With important contributions from:                           
%        M. Damian, Villanova University, USA                      
%        R. Sander, Max-Planck Institute for Chemistry, Mainz, Germany
%                                                                  
% File                 : cbm_Sparse.m                              
% Time                 : Fri Mar 15 14:06:05 2013                  
% Working directory    : /home/sandu/kpp-2.2.3/examples/Cbm_matlab 
% Equation file        : cbm.kpp                                   
% Output root filename : cbm                                       
%                                                                  
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




%  ----------> Sparse Jacobian Data                                

% LU_IROW - Row indexes of the LU Jacobian of variables            
 global LU_IROW;
% LU_ICOL - Column indexes of the LU Jacobian of variables         
 global LU_ICOL;
% LU_CROW - Compressed row indexes of the LU Jacobian of variables 
 global LU_CROW;
% LU_DIAG - Diagonal indexes of the LU Jacobian of variables       
 global LU_DIAG;


%  ----------> Sparse Hessian Data                                 

% IHESS_I - Index i of Hessian element d^2 f_i/dv_j.dv_k           
 global IHESS_I;
% IHESS_J - Index j of Hessian element d^2 f_i/dv_j.dv_k           
 global IHESS_J;
% IHESS_K - Index k of Hessian element d^2 f_i/dv_j.dv_k           
 global IHESS_K;


%  ----------> Sparse Stoichiometric Matrix                        

% STOICM - Stoichiometric Matrix in compressed column format       
 global STOICM;
% IROW_STOICM - Row indices in STOICM                              
 global IROW_STOICM;
% CCOL_STOICM - Beginning of columns in STOICM                     
 global CCOL_STOICM;
% ICOL_STOICM - Column indices in STOICM                           
 global ICOL_STOICM;


%  ----------> Sparse Data for Jacobian of Reactant Products       

% ICOL_JVRP - Column indices in JVRP                               
 global ICOL_JVRP;
% IROW_JVRP - Row indices in JVRP                                  
 global IROW_JVRP;
% CROW_JVRP - Beginning of rows in JVRP                            
 global CROW_JVRP;

