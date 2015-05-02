function swe_simpleplot(tout,yout)
[ NX, NY] = swe_Constant_Parameters();


[surfplot,top,stop]=initgraphics();
bb=size(yout);
NVAR=(NX+2)*(NY+2);
movieFlag = 1;
while movieFlag==1
    for k=1:bb(1)
        H=reshape(yout(k,1:NVAR),NX+2,NY+2);H=H';
        U=reshape(yout(k,NVAR+1:2*NVAR),NX+2,NY+2);U=U';
        V=reshape(yout(k,2*NVAR+1:3*NVAR),NX+2,NY+2);V=V';
        i=2:NX+1;
        j=2:NY+1;
        C = abs(U(i,j)) + abs(V(i,j));  % Color shows momemtum
        tv = norm(C,'fro');
        set(surfplot,'zdata',H(i,j),'cdata',C);
        set(top,'string',sprintf('t = %6.2f,  tv = %6.2f',tout(k),tv))
        drawnow
        if (get(stop,'value')==1)
            movieFlag = 0;
            break;
        end
    end
end

function [surfplot,top,stop] = initgraphics()
% INITGRAPHICS  Initialize graphics for waterwave.
% [surfplot,top,start,stop] = initgraphics(n)
% returns handles to a surface plot, its title, and two uicontrol toggles.
[ NX, NY] = swe_Constant_Parameters();

   clf;
   shg;
   set(gcf,'numbertitle','off','name','Shallow_water')
   x = (0:NX-1)/(NX-1);
   surfplot = surf(x,x,ones(NX,NY),zeros(NX,NY));
   grid off
   axis([0 1 0 1 -1 3])
   caxis([-1 1])
   shading faceted
   c = (1:NX)'/NX;
   cyan = [0*c c c];
   colormap(cyan)
   top = title('Click start');
   %start = uicontrol('position',[20 20 80 20],'style','toggle','string','start');
   stop = uicontrol('position',[120 20 80 20],'style','toggle','string','stop');
  

