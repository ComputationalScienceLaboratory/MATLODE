%% fatOde_ROS_ErrorNorm
%
% <html>
%   <div>
%       <img style="float: right" src="../../MATLODE_LOGO.png" height="150px"></img>
%   </div>
% </html>
%
%% Syntax
%
%% Input Parameters
%
%% Output Parameters
%
%% Description
%
%% Reference
% [1] Tony D'Augustine, Adrian Sandu. MATLODE: A MATLAB ODE Solver and
%     Sensitivity Analysis Toolbox. Submitted to ACM TOMS.
%
%%
%  Authored by Tony D'Augustine, Adrian Sandu, and Hong Zhang.
%  Computational Science Laboratory, Virginia Tech.
%  ©2015 Virginia Tech Intellectual Properties, Inc.
%
function [ ros_ErrorNorm ] = fatOde_ROS_ErrorNorm( NVAR, Y, Ynew, Yerr, VectorTol, OPTIONS )

    Err = 0.0;
    for i=1:NVAR
        Ymax = max( abs(Y(i)), abs(Ynew(i)) );
        if ( VectorTol )
            SCAL = OPTIONS.AbsTol(i) + OPTIONS.RelTol(i)*Ymax;
        else
            SCAL = OPTIONS.AbsTol(1) + OPTIONS.RelTol(1)*Ymax;
        end
        Err = Err + ( Yerr(i)/SCAL )^2;
    end
    Err = sqrt( Err/NVAR );
    
    %SCAL_1 = OPTIONS.AbsTol + OPTIONS.RelTol.*max(abs(Y),abs(Ynew));
    %Err_1 = sqrt(sum((Yerr./SCAL_1).^2)/NVAR);
    
%     if(~VectorTol) 
%         at = OPTIONS.AbsTol(1)*ones(1,NVAR);
%         rt = OPTIONS.RelTol(1)*ones(1,NVAR);
%     else
%         at = OPTIONS.AbsTol;
%         rt = OPTIONS.RelTol;
%     end
%     SCAL2 = at + rt.*(max(Y, Ynew)');
%     Err2 = norm(Yerr./SCAL2')/sqrt(NVAR);
%     if( VectorTol )
%         SCAL2 = OPTIONS.AbsTol + OPTIONS.RelTol.*max(Y,Ynew)';
%         Err2 = norm(Yerr./SCAL2')/sqrt(NVAR);
%     else
%         SCAL2 = OPTIONS.AbsTol(1)*ones(1,NVAR) + OPTIONS.RelTol(1)*max(Y, Ynew)';
%         Err2 = norm(Yerr./SCAL2')/sqrt(NVAR);
%     end
%     
%     diff = abs(Err - Err2)
%     
%     if ( diff > 1 )
%         keyboard
%     end
    
    ros_ErrorNorm = max( Err, 1d-10 );

return;

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
%       <td>adaug13@vt.edu</td>
%       <td>Release MATLODE_v2.0.00</td>
%   </tr>
% </table>
% </html>
% 
%%
% <html>
%   <div>
%       <img style="float: right" src="../../CSL_LogoWithName_1.png" height="50px"></img>
%   </div>
% </html>