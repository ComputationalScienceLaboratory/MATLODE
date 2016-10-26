% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%                                                                  
% ReactantProd - Reactant Products in each equation                
%   Arguments :                                                    
%      V         - Concentrations of variable species (local)      
%      F         - Concentrations of fixed species (local)         
%      ARP       - Reactant product in each equation               
%                                                                  
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%                                                                  
% Generated by KPP - symbolic chemistry Kinetics PreProcessor      
%     KPP is developed at CGRER labs University of Iowa by         
%     Valeriu Damian & Adrian Sandu                                
%                                                                  
% File                 : cbm_ReactantProd.m                        
% Time                 : Wed Dec 31 19:00:00 1969                  
% Working directory    : /home/sandu/kpp-2.2.3/examples/Cbm_matlab 
% Equation file        : cbm.kpp                                   
% Output root filename : cbm                                       
%                                                                  
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function  [ ARP ] =  cbm_ReactantProd ( V , F )


% Reactant Products in each equation are useful in the             
%     stoichiometric formulation of mass action law                
   ARP(1) = V(26) ;
   ARP(2) = V(29) ;
   ARP(3) = V(25)*V(31) ;
   ARP(4) = V(26)*V(29) ;
   ARP(5) = V(26)*V(29) ;
   ARP(6) = V(29)*V(31) ;
   ARP(7) = V(25)*V(26) ;
   ARP(8) = V(25) ;
   ARP(9) = V(25) ;
   ARP(10) = V(1) ;
   ARP(11) = V(1)*F(1) ;
   ARP(12) = V(25)*V(27) ;
   ARP(13) = V(25)*V(28) ;
   ARP(14) = V(30) ;
   ARP(15) = V(30)*V(31) ;
   ARP(16) = V(26)*V(30) ;
   ARP(17) = V(26)*V(30) ;
   ARP(18) = V(6)*F(1) ;
   ARP(19) = V(6) ;
   ARP(20) = V(31)*V(31) ;
   ARP(21) = V(26)*V(31)*F(1) ;
   ARP(22) = V(27)*V(31) ;
   ARP(23) = V(9) ;
   ARP(24) = V(9)*V(27) ;
   ARP(25) = V(9)*V(9) ;
   ARP(26) = V(26)*V(27) ;
   ARP(27) = V(12)*V(27) ;
   ARP(28) = V(28)*V(31) ;
   ARP(29) = V(26)*V(28) ;
   ARP(30) = V(10) ;
   ARP(31) = V(10)*V(27) ;
   ARP(32) = V(28)*V(28) ;
   ARP(33) = V(28)*V(28)*F(1) ;
   ARP(34) = V(2) ;
   ARP(35) = V(2)*V(27) ;
   ARP(36) = V(16)*V(27) ;
   ARP(37) = V(21)*V(27) ;
   ARP(38) = V(21) ;
   ARP(39) = V(21) ;
   ARP(40) = V(21)*V(29) ;
   ARP(41) = V(21)*V(30) ;
   ARP(42) = V(24)*V(29) ;
   ARP(43) = V(24)*V(27) ;
   ARP(44) = V(24)*V(30) ;
   ARP(45) = V(24) ;
   ARP(46) = V(31)*V(32) ;
   ARP(47) = V(26)*V(32) ;
   ARP(48) = V(3) ;
   ARP(49) = V(32)*V(32) ;
   ARP(50) = V(28)*V(32) ;
   ARP(51) = V(27) ;
   ARP(52) = V(20)*V(27) ;
   ARP(53) = V(13) ;
   ARP(54) = V(13) ;
   ARP(55) = V(13)*V(26) ;
   ARP(56) = V(23)*V(29) ;
   ARP(57) = V(23)*V(27) ;
   ARP(58) = V(23)*V(25) ;
   ARP(59) = V(23)*V(30) ;
   ARP(60) = V(17)*V(29) ;
   ARP(61) = V(17)*V(27) ;
   ARP(62) = V(17)*V(25) ;
   ARP(63) = V(5)*V(27) ;
   ARP(64) = V(11)*V(31) ;
   ARP(65) = V(11) ;
   ARP(66) = V(14)*V(27) ;
   ARP(67) = V(14)*V(30) ;
   ARP(68) = V(4)*V(26) ;
   ARP(69) = V(7)*V(27) ;
   ARP(70) = V(19)*V(27) ;
   ARP(71) = V(19) ;
   ARP(72) = V(19)*V(25) ;
   ARP(73) = V(15)*V(27) ;
   ARP(74) = V(15) ;
   ARP(75) = V(22)*V(29) ;
   ARP(76) = V(22)*V(27) ;
   ARP(77) = V(22)*V(25) ;
   ARP(78) = V(22)*V(30) ;
   ARP(79) = V(18)*V(31) ;
   ARP(80) = V(18)*V(18) ;
   ARP(81) = V(8)*V(31) ;
      
return

% End of ReactantProd function                                     
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


