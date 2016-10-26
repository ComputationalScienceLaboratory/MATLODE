function [ erkMethod, erkELO, erkS, erkName ] = Coefficients_BackwardEuler( RK1 )

    global rkA rkB rkC

    erkMethod = RK1;
    erkS = 1;
    erkName = 'BackwardEuler';
    erkELO = 1; % tmep fix
    
    rkA(1,1) = 1;
    
    rkB(1,1) = 1;
    
    rkC(1) = 1;
    

end

