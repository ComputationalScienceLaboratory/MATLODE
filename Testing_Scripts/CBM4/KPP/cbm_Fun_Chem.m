
% Wrapper for calling the ODE function routine
% in a format required by Matlab's ODE integrators

function P = cbm_Fun_Chem(T, Y) 
     
  global TIME FIX RCONST ind_NO ind_NO2 
 
  Told = TIME;
  TIME = T;
  cbm_Update_SUN;
  cbm_Update_RCONST;
  
%  This line calls the Matlab ODE function routine  
  P = cbm_Fun( Y, FIX, RCONST );
  
%  To call the mex routine instead, comment the line above and uncomment the following line:
%  P = cbm_mex_Fun( Y, FIX, RCONST );

% Add emissions
  P(ind_NO)  = P(ind_NO)  + 5.0e6;
  P(ind_NO2) = P(ind_NO2) + 5.0e6;
  
  TIME = Told;

return
