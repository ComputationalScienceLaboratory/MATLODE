function [ ros_ErrorNorm ] = fatOde_ROS_ErrorNorm( NVAR, Y, Ynew, Yerr, VectorTol, OPTIONS )
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Filenanme: fatODe_ROS_ErrorNorm.m
%
% Original Author:
%
% File Creation Date:
%
% Input Arguments:
%
% Output Arguments:
%
% Modification History:
%   Date        Developer         Email             Action  
%   7/17/2012   Tony D'Augustine  adaug13@vt.edu 
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% fatOde_ROS_ErrorNorm:
%   Compute the "scaled norm" of the error vector Yerr.
%
% fatOde_ROS_ErrorNorm: INPUT ARGUMENTS
%
% fatOde_ROS_ErrorNorm: OUTPUT ARGUMENTS
%
% fatOde_ROS_ErrorNorm: SYNTAX
%
% fatOde_ROS_ErrorNorm: EXAMPLE
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

