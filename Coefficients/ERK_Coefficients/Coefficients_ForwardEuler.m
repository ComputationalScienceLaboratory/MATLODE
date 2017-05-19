function [ erkMethod, erkELO, erkS, erkName ] = Coefficients_ForwardEuler( RK1 )

    global rkA rkB rkC

    erkMethod = RK1;
    erkS = 1;
    erkName = 'Euler';
    erkELO = 0; % temp
    
    rkA(1,1) = 0;
    
    rkB(1,1) = 1;
    
    rkC(1) = 0;
    

end

