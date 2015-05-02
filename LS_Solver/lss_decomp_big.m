function [ e_big ] = lss_decomp_big( NVAR, H, fjac1, fjac2, fjac3 )
%LSS_DECOMP_BIG Summary of this function goes here
%   Detailed explanation goes here

    global rkA

%     idx = 1;
%     for j=1:3*NVAR
%         for i=1:3*NVAR
%             row = mod(i-1,NVAR)+1;
%             col = mod(j-1,NVAR)+1;
%             if ( j <= NVAR )
%                 current = -H*rkA((i-row)/NVAR+1,1)*fjac1(row,col);
%             elseif ( j <= 2*NVAR )
%                 current = -H*rkA((i-row)/NVAR+1,2)*fjac2(row,col);
%             else
%                 current = -H*rkA((i-row)/NVAR+1,3)*fjac3(row,col);
%             end
%             if ( current ~= 0.0 )
%                 if ( i == j )
%                     ax_big(idx) = current + 1;
%                 else
%                     ax_big(idx) = current;
%                 end
%                 idx = idx + 1;
%             end
%         end
%     end

%     for j=1:NVAR
%         for i=1:NVAR
%           ax_big(i,j)               = -H*rkA(1,1)*fjac1(i,j);
%           ax_big(NVAR+i,j)          = -H*rkA(2,1)*fjac1(i,j);
%           ax_big(2*NVAR+i,j)        = -H*rkA(3,1)*fjac1(i,j);
%           ax_big(i,NVAR+j)          = -H*rkA(1,2)*fjac2(i,j);
%           ax_big(NVAR+i,NVAR+j)     = -H*rkA(2,2)*fjac2(i,j);
%           ax_big(2*NVAR+i,NVAR+j)   = -H*rkA(3,2)*fjac2(i,j);
%           ax_big(i,2*NVAR+j)        = -H*rkA(1,3)*fjac3(i,j);
%           ax_big(NVAR+i,2*NVAR+j)   = -H*rkA(2,3)*fjac3(i,j);
%           ax_big(2*NVAR+i,2*NVAR+j) = -H*rkA(3,3)*fjac3(i,j);
%         end
%     end
    
    HrkA = -H.*rkA;
    e_big = [ HrkA(1,1).*fjac1 HrkA(1,2).*fjac2 HrkA(1,3).*fjac3; ...
                    HrkA(2,1).*fjac1 HrkA(2,2).*fjac2 HrkA(2,3).*fjac3; ...
                    HrkA(3,1).*fjac1 HrkA(3,2).*fjac2 HrkA(3,3).*fjac3 ];
    
                
%     rRMS(ax_big,ax_big_temp)
    
    
%     for i=1:3*NVAR
%         e_big(i,i) = e_big(i,i) + 1;
%     end
    
    e_big = e_big + eye(3*NVAR,3*NVAR);
               
return;