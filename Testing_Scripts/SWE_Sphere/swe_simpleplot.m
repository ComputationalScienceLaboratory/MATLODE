function swe_simpleplot(yout,~)
global imax jmax

%It can be called by swe_simpleplot(yout, 0)
% Input the final solution to visualize. This code is written specifically for 
%Shallow water equations, with the first third being U, second third being V and 
%the last third being H. But it can be changed easily by adding/deleting components 
%in the for loop. It should also be noted that if we change NX1, NY1, NX and NY then 
%the NVAR also changes. Here 3*NVAR represents the number of variables in yout.
%The code was written mainly to visualize a series of solution snapshots.
%So yout should of size K X 3NVAR. K should be atleast 1.

NX1=jmax;NY1=imax;NX=jmax;NY=imax;

NVAR=NX1*NY1;

bb=size(yout);
%Creates the map
[x, y] =meshgrid(linspace(-85.58,83.65,NX),linspace(-180,190.17,NY));

for k=1:bb(1);
    u=reshape(yout(k,1:NVAR),NY1,NX1);
    v=reshape(yout(k,NVAR+1:2*NVAR),NY1,NX1);
    p=reshape(yout(k,2*NVAR+1:3*NVAR),NY1,NX1);
    
    
    %plots H
    
    figure
    coast =load('coast');
    axesm eckert4; framem; gridm;
    patchesm(coast.lat,coast.long,'FaceColor',[1 1 1])
    contourfm(x',y',p','LineStyle', 'none');
    title('Height');
    colormap(jet(11));
    colorbar
    geoshow(coast.lat, coast.long, 'DisplayType', 'line', 'Color', [0 0 0]);
    
    %plots U
    figure
    coast = load('coast');
    axesm eckert4; framem; gridm;
    patchesm(coast.lat,coast.long,'FaceColor',[1 1 1])
    contourfm(x',y',u', 'LineStyle', 'none' );
    title('U');
    colormap(jet(11));
    colorbar;
    geoshow(coast.lat, coast.long, 'DisplayType', 'line', 'Color', [0 0 0]);

    %plots V
    figure
    coast = load('coast');
    axesm eckert4; framem; gridm;
    patchesm(coast.lat,coast.long,'FaceColor',[1 1 1])
    contourfm(x',y',v', 'LineStyle', 'none' );
    title('V');
    colormap(jet(11));
    colorbar;
    geoshow(coast.lat, coast.long, 'DisplayType', 'line', 'Color', [0 0 0]);
    
   
end

