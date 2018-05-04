function [AF, AS, BTF, BTS, BThatF, BThatS, Gamma, methodOrder, embeddedOrder, Afs, Asf] = EXIM_MRGARK_FWD_Coupling(n)

switch(n)
        
    case 2
        
        %Two Stage Second order
        Gamma = 1 - 1/sqrt(2);
        
        AF = [0, 0;
            2/3, 0];
        AS = [Gamma, 0;
            1/sqrt(2), Gamma];
        
        BTF = [1/4; 3/4];
        BTS = [1 - Gamma; Gamma];
        
        BThatF = [1; 0];
        BThatS = [3/5; 2/5];
        
        methodOrder = 2;
        embeddedOrder = 1;
        
        Afs = @AfastslowDepenTwoStageSecondOrder;
        Asf = @AslowfastDepenTwoStageSecondOrder;
        
    case 3
        Gamma = (1/2).*(2+(-1).*2.^(1/2).*cos((1/3).*acot(2.*2.^(1/2)))+6.^(1/2).* ...
            sin((1/3).*acot(2.*2.^(1/2))));
        
        AF = [0, 0, 0;
            1/2, 0, 0;
            0, 3/4, 0];
        
        AS = [...
            Gamma,0,0;
            ...
            (3+(-12).*Gamma+6.*Gamma.^2).^(-1).*(2+(-12).*Gamma+ ...
            18.*Gamma.^2+(-6).*Gamma.^3),Gamma,0;
            ...
            (1+(-4).*Gamma).*(4+(-24).* ...
            Gamma+36.*Gamma.^2+(-12).*Gamma.^3).^(-1),(-3/4).*(1+(-4).*Gamma+ ...
            2.*Gamma.^2).^2.*((-1)+6.*Gamma+(-9).*Gamma.^2+3.*Gamma.^3).^(-1), ...
            Gamma];
        
        BTF = [2/9; 1/3; 4/9];
        BTS = [...
            (1+(-4).*Gamma).*(4+(-24).*Gamma+36.*Gamma.^2+(-12).*Gamma.^3).^(-1);
            ...
            (-3/4).*(1+(-4).*Gamma+2.*Gamma.^2).^2.*((-1)+6.*Gamma+(-9).* ...
            Gamma.^2+3.*Gamma.^3).^(-1);
            ...
            Gamma];
        
        BThatF = [1/40; 37/40; 1/20];
        BThatS = [...
            (1+(-6).*Gamma+6.*Gamma.^2).*(4+(-24).*Gamma+36.*Gamma.^2 + ...
            (-12).*Gamma.^3).^(-1);
            ...
            (3/4).*((-1)+6.*Gamma+(-9).*Gamma.^2+3.*Gamma.^3) ...
            .^(-1).*((-1)+6.*Gamma+(-10).*Gamma.^2+4.*Gamma.^3);
            ...
            0];
        
        methodOrder = 3;
        embeddedOrder = 2;
        
        Afs = @AfastslowDepenThreeStageThirdOrder;
        Asf = @AslowfastDepenThreeStageThirdOrder;
    case 4
        Gamma = 1/4;
        
        AF = [... 
            0,0,0,0,0,0;
            (1/4),0,0,0,0,0;
            (3/32),(9/32),0,0,0,0;
            (1932/2197),( ...
            -7200/2197),(7296/2197),0,0,0;(439/216),(-8),(3680/513),( ...
            -845/4104),0,0;
            ...
            (-8/27),2,(-3544/2565),(1859/4104),(-11/40),0];
        
        AS = [(1/4),0,0,0,0,;...
            (13/20),(1/4),0,0,0,;
            (580/1287),(-175/5148),(1/4),0,0,;
            (12698/37375),(-201/2990),(891/11500),(1/4),0;
            (944/1365),(-400/819),(99/35),(-575/252),(1/4)];
        
        BTF = [25/216; 0; 1408/2565; 2197/4104; -1/5; 0];
        BTS = [ 944/1365; -400/819; 99/35; -575/252; 1/4];
        
        BThatF = [16/135; 0; 6656/12825; 28561/56430; -9/50; 2/55];
        BThatS = [41911/60060; -83975/144144; 3393/1120; -27025/11088; 103/352];
        
        methodOrder = 4;
        embeddedOrder = 3;
        
        Afs = @AfastslowDepenSixStageFourthOrder;
        Asf = @AslowfastDepenSixStageFourthOrder;
        
        
    otherwise
        [AF, AS, BTF, BTS, BThatF, BThatS, Gamma, methodOrder, embeddedOrder, Afs, Asf] = EXIM_MRGARK_FWD_Coupling(4);
end

end

function Afslambda = AfastslowDepenTestStageSecondOrder(lambda, ~, M)

     Afslambda = [...
         (lambda - 1)/M;
         (3*lambda - 1)/(3*M) ...
         ];
end

function Asflambda = AslowfastDepenTestStageSecondOrder(lambda, ~, M)

    if lambda == 1
        
        Asflambda = [M/2, 0];
        
    else
        
        Asflambda = [0, 0];
        
    end

end

function Afslambda = AfastslowDepenTwoStageSecondOrder(lambda, ~, M)

     Afslambda = [...
         (lambda - 1)/M, 0;
         (3*lambda - 1)/(3*M), 0 ...
         ];
end

function Asflambda = AslowfastDepenTwoStageSecondOrder(lambda, ~, M)

    if lambda == 1
        
        Asflambda = [ ...
            M - (M/(sqrt(2))), 0;
            1/4, 3/4 ...
            ];
        
    else
        
        Asflambda = [0, 0; 1/4, 3/4];
        
    end

end

function Afslambda = AfastslowDepenThreeStageThirdOrder(lambda, gamma, M)

    Afslambda = [ ...
        ((-1)+lambda).*M.^(-1), 0, 0;
        ...
        (1/2).*((-1)+2.*lambda).*M.^(-1), 0, 0;
        ...
        (1/16).*((-1)+6.*gamma+(-9).*gamma.^2+3.*gamma.^3).^(-1).*M.^(-1).* ...
        (4+gamma.^3.*(42+(-60).*lambda)+72.*gamma.^2.*((-1)+lambda)+(-16) ...
        .*lambda+gamma.*(3+42.*lambda)+9.*(1+(-4).*gamma+2.*gamma.^2).*M), ...
        (-9/16).*(1+(-4).*gamma+2.*gamma.^2).*((-1)+6.*gamma+(-9).* ...
        gamma.^2+3.*gamma.^3).^(-1).*M.^(-1).*(gamma.*(3+(-6).*lambda)+M), ...
        0];
end

function Asflambda = AslowfastDepenThreeStageThirdOrder(lambda, gamma, M)

if lambda == 1
    
    Asflambda = [...
        gamma.*M, 0, 0;
        ...
        (1/9).*(1+(-4).*gamma+2.*gamma.^2).^(-2).*M.*(3.*(2 + ...
        (-17).*gamma+46.*gamma.^2+(-42).*gamma.^3+12.*gamma.^4)+(-4).*(1+( ...
        -9).*gamma+27.*gamma.^2+(-30).*gamma.^3+9.*gamma.^4).*M), ...
        (4/9).*(1+(-4).*gamma+2.*gamma.^2).^(-2).*(1+(-9).*gamma+27.*gamma.^2+( ...
        -30).*gamma.^3+9.*gamma.^4).*M.^2, 0;
        ...
        (2/9),(1/3),(4/9)];
else
    
    Asflambda = [
        0,0,0;
        0,0,0;
        (2/9),(1/3),(4/9)
        ];
end

end

function Afslambda = AfastslowDepenSixStageFourthOrder(lambda, ~, M)

    Afslambda = [ ...
        ((-1)+lambda).*M.^(-1),0,0,0,0,;
        (1/4).*((-3)+4.*lambda).*M.^(-1), ...
        0,0,0,0,;
        (1/416).*M.^(-2).*(90.*((-1)+lambda)+((-335)+551.*lambda) ...
        .*M+(-90).*M.^2+45.*M.^3),(-15/416).*M.^(-2).*((-6)+6.*lambda+(-5) ...
        .*M+9.*lambda.*M+(-6).*M.^2+3.*M.^3),0,0,0,;...
        (1/2197).*M.^(-2).*(( ...
        -3960).*((-1)+lambda)+((-3709)+6517.*lambda).*M+(-2880).*M.^2+ ...
        1440.*M.^3),(-60/2197).*M.^(-2).*(66+(-66).*lambda+(-59).*M+72.* ...
        lambda.*M+(-48).*M.^2+24.*M.^3),0,0,0,; ...
        (1/273).*M.^(-2).*((-1155).* ...
        ((-1)+lambda)+((-529)+386.*lambda).*M+(-362).*M.^2+560.*M.^3),( ...
        5/1638).*(2439+1386.*((-1)+lambda).*M.^(-2)+(2590+(-4046).*lambda) ...
        .*M.^(-1)+(-672).*M),(-33/28).*M.^(-1).*(7+(-14).*lambda+11.*M),( ...
        575/252).*M.^(-1).*(1+(-2).*lambda+3.*M),0,; ...
        0,0,(-109/32)+(165/32) ...
        .*M.^(-2)+(-75/8).*M.^(-1)+5.*M,(1/32).*(109+(-165).*M.^(-2)+4.*( ...
        71+8.*lambda).*M.^(-1)+(-160).*M),0,;];
end

function Asflambda = AslowfastDepenSixStageFourthOrder(lambda, ~, M)

if lambda == 1
    
    Asflambda = [...
        (1/4).*M,0,0,0,0,0;
        (1/100).*(90+(-169).*M).*M,(169/100).*M.^2,0, 0,0,0;
        (1/198).*(132+(-155).*M).*M,(155/198).*M.^2,0,0,0,0;
        (1/920).*(552+(-497).*M).*M,(14/23).*M.^2,(-896/10925).*M.^2,(1183/87400) ...
        .*M.^2,0,0;(25/216),0,(1408/2565),(2197/4104),(-1/5),0];
else
    
    Asflambda = [0,0,0,0,0,0;
        0,0,0,0,0,0;
        0,0,0,0,0,0;
        0,0,0,0,0,0;
        (25/216),0,(1408/2565),(2197/4104),(-1/5),0];
end

end




