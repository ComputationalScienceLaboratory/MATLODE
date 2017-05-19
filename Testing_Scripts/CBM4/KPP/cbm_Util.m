% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%                                                                  
% Auxiliary Routines File                                          
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
% File                 : cbm_Util.m                                
% Time                 : Fri Mar 15 14:06:05 2013                  
% Working directory    : /home/sandu/kpp-2.2.3/examples/Cbm_matlab 
% Equation file        : cbm.kpp                                   
% Output root filename : cbm                                       
%                                                                  
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




% User INLINED Utility Functions                                   

% End INLINED Utility Functions                                    

% Utility Functions from KPP_HOME/util/util                        
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%                                                                  
% UTIL - Utility functions                                         
%   Arguments :                                                    
%                                                                  
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% ****************************************************************
%                            
% InitSaveData - Opens the data file for writing
%
% ****************************************************************

function InitSaveData ()

global cbm_FID

      cbm_FID = fopen('cbm.dat','w');

return %  InitSaveData

% End of InitSaveData function
% ****************************************************************

% ****************************************************************
%                            
% SaveData - Write LOOKAT species in the data file 
%
% ****************************************************************

function SaveData ()

global VAR FIX CFACTOR LOOKAT NLOOKAT cbm_FID

      C(1:32) = VAR(1:32);
      C(32+1:33) = FIX(1:1);
      
      fprintf(cbm_FID,'%12.5e,',C(LOOKAT(1:NLOOKAT)));

return %  SaveData

% End of SaveData function
% ****************************************************************

% ****************************************************************
%                            
% CloseSaveData - Close the data file 
%
% ****************************************************************

function CloseSaveData ()
global cbm_FID

      fclose( cbm_FID );

return %  CloseSaveData

% End of CloseSaveData function
% ****************************************************************


% End Utility Functions from KPP_HOME/util/util                    
% End of UTIL function                                             
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


