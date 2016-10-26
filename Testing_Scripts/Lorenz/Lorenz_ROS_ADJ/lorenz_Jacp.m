function fpjac = lorenz_Jacp( ~, y )

    fpjac(1,1) = 0;
    fpjac(1,2) = 0;
    fpjac(1,3) = -y(1)+y(2);
    
    fpjac(2,1) = 0;
    fpjac(2,2) = y(1);
    fpjac(2,3) = 0;
    
    fpjac(3,1) = -y(3);
    fpjac(3,2) = 0;
    fpjac(3,3) = 0;

return;

