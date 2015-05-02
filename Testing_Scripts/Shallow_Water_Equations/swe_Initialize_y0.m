function [ y0, H ] = swe_Initialize_y0()

    [ NX, NY, ~, ~, ~, midgrid] = swe_Constant_Parameters();
    
    H = ones(NX+2,NY+2); 
    U = zeros(NX+2,NY+2);  
    V = zeros(NX+2,NY+2);
    
    U(NX/2-midgrid:NX/2+midgrid,NX/2-midgrid:NX/2+midgrid)=1;
    V(NX/2-midgrid:NX/2+midgrid,NX/2-midgrid:NX/2+midgrid)=1;
    
    %X Momentum, Y Momentum
    for i=1:NX
        for j=1:NY
            H(i,j) = 1 + exp(-((i-NX/3)^2+(j-2*NY/3)^2)/36);
        end
    end
  
   NVAR=(NX+2)*(NY+2);
   y0=zeros(NVAR,1);
  
   y0(1:NVAR)=reshape(transpose(H),NVAR,1);
   y0(NVAR+1:2*NVAR)=reshape(transpose(U),NVAR,1);
   y0(2*NVAR+1:3*NVAR)=reshape(transpose(V),NVAR,1);

return;

