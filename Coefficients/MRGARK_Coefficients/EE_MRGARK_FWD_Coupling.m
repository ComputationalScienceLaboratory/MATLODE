function [A, BT, BThat, methodOrder, embeddedOrder, Afs, Asf] = EE_MRGARK_FWD_Coupling(n)

switch(n)
    case 1
        
        %Two Stage parameterized
        A = [0, 0; 2/3, 0];
        BT = [1/4; 3/4];
        BThat = [1; 0];
        
        methodOrder = 2;
        embeddedOrder = 1;
        
        Afs = @AfastslowDepenTwoStageSecondOrderP;
        Asf = @AslowfastDepenTwoStageSecondOrderP;
        
    case 2
        
        %Two Stage Second order
        A = [0, 0;
            2/3, 0];
        
        BT = [1/4; 3/4];  
        BThat = [1; 0];
        
        methodOrder = 2;
        embeddedOrder = 1;
        
        Afs = @AfastslowDepenTwoStageSecondOrder;
        Asf = @AslowfastDepenTwoStageSecondOrder;
        
    case 3
      
        %Three Stage Third order method
        A = [0, 0, 0; 
            1/2, 0, 0;
            0, 3/4, 0];
        
        BT = [2/9; 1/3; 4/9];
        BThat = [1/40; 37/40; 1/20];
        
        methodOrder = 3;
        embeddedOrder = 2;
        
        Afs = @AfastslowDepenThreeStageThirdOrder;
        Asf = @AslowfastDepenThreeStageThirdOrder;
        
    case 4
        
        %Five Stage Fourth Order method
        A = [0, 0, 0, 0, 0;
            (2/5), 0, 0, 0, 0;
            (-3/20), (3/4), 0, 0, 0;
            (19/44), (-15/44), (10/11), 0, 0;
            (11/72), (25/72), (25/72), (11/72), 0];
        
        BT = [11/72; 25/72; 25/72; 11/72; 0];
        BThat = [1251515/8970912; 3710105/8970912; 2519695/8970912;
            61105/8970912; 119041/747576];
        
        methodOrder = 4;
        embeddedOrder = 3;
        
        Afs = @AfastslowDepenFiveStageFourthOrder;
        Asf = @AslowfastDepenFiveStageFourthOrder;
        
    otherwise
        %Three Stage Third order method
        A = [0, 0, 0; 
            1/2, 0, 0;
            0, 3/4, 0];
        
        BT = [2/9; 1/3; 4/9];
        BThat = [1/40; 37/40; 1/20];
        
        methodOrder = 3;
        embeddedOrder = 2;
        
        Afs = @AfastslowDepenThreeStageThirdOrder;
        Asf = @AslowfastDepenThreeStageThirdOrder;
              
end

end

function Afslambda = AfastslowDepenTwoStageSecondOrderP(lambda, M)

if M == 1
    Afslambda = AfastslowDepenTwoStageSecondOrder(lambda, M);
else
    l2 = floor(2*M/3);
    
    if lambda <= l2;
        
        Afslambda = [(lambda - 1)/M, 0; (lambda - 1/3)/M, 0];
        
    else
        
        Afslambda = [ ...
            (1/4).*((-1)+lambda).*M.^(-1),(3/4).*((-1)+lambda).*M.^(-1);( ...
            1/12).*M.^(-1).*((-7)+15.*lambda+4.*M.^2.*((-1).*M+floor((2/3).*M) ...
            ).^(-1)),(1/4).*M.^(-1).*(1+(-1).*lambda+4.*M.^2.*(3.*M+(-3).* ...
            floor((2/3).*M)).^(-1))];
    end
end

end


function Asflambda = AslowfastDepenTwoStageSecondOrderP(lambda, M)

if M == 1
    Asflambda = AslowfastDepenTwoStageSecondOrder(lambda, M);
else
    l2 = floor(2*M/3);
    
    if lambda <= l2
        Asflambda = [0, 0;
            (M * (1 - 2*M + 3*l2)) / (6*l2), ...
            (M * (3 + 2*M - 3*l2))/(6*l2)];
    else
        Asflambda = [0, 0; 0, 0];
    end
end

end

function Afslambda = AfastslowDepenTwoStageSecondOrder(lambda, M)

     if lambda == 1
         
         Afslambda = [ 0, 0; 2/(3*M), 0]; 
       
     else
         
        Afslambda = [ ... 
            (3*M^3 - 11*M^2 + 20*lambda*M - 20*M - 20*lambda + 20) / (20 * (M-1)*M) ...
            , - (M * (3*M - 11))/ (20*(M-1));
            ...
            (-3*M^3 - 9*M^2 + 60*lambda*M - 20*M - 60*lambda + 20) / (60 * (M-1) * M) ...
            , M*(M+3)/(20 * (M-1))
            ];
         
     end
end

function Asflambda = AslowfastDepenTwoStageSecondOrder(lambda, M)

if M == 1
    
    Asflambda = [0, 0; 2/3, 0];
else
    if lambda == 1
        
        Asflambda = [0, 0; -(1/3)*(M - 2)*M, M^2/3];
        
    else
        
        Asflambda = [0, 0; 0, 0];
        
    end
end

end
 function Afslambda = AfastslowDepenThreeStageThirdOrder(lambda, M)
       
 if M == 1
     Afslambda = [0, 0, 0; 
            1/2, 0, 0;
            0, 3/4, 0];
 else
        if lambda == 1
            Afslambda = [0, 0, 0; 1/(2*M), 0, 0; 0, 3/(4*M), 0];
            
        else
            
            Afslambda = [ ...
                (3*M^3 - 8*M^2 + 6*lambda*M - 6*lambda + 6) / (6 * (M-1)*M) ...
                , (-3*M^2 + 8*M - 6)/ (6 * (M -1)) ...
                , 0; ...
                ...
                (-2*M^2 + 6*lambda*M - 3*M - 6*lambda + 3) / (6 * (M-1)*M), ...
                M / (3 * (M-1)) ...
                , 0; ...
                ...
                (-3*M^3 + 2*M^2 + 12*lambda*M - 9*M - 12*lambda + 12) / (12 * (M-1)*M), ...
                (3*M^3 - 2*M^2 + 6*M - 9) / (12 * (M-1)*M) ...
                , 0
            ];
        end
 end
            
 end
    
 function Asflambda = AslowfastDepenThreeStageThirdOrder(lambda, M)
 
 if M == 1
          Asflambda = [0, 0, 0; 
            1/2, 0, 0;
            0, 3/4, 0];
 else
        if lambda == 1
            Asflambda = [ 0, 0, 0; ...
                (-1/66)*M*(16*M - 33), 8 * M^2 / 33,  0; ...
                ...
                (1/264 * (11*M^4 - 22*M^3 + 26*M^2 + 11*M + 44)) ...
                ,(1/88 * (-11*M^4 + 22*M^3 - 16*M^2 - 11*M + 22)) ...
                ,(1/12 * (M^4 - 2*M^3 + M^2 + M + 4))
                ];
        else
            
            Asflambda = [0, 0, 0; ...
                0, 0, 0; ...
                ...
                (-M^4 + 2*M^3 + 2*M^2 + 3*M - 4)/ (24*(M-1)) ...
                , (1/8) * (M^3 - M^2 -M + 2) ...
                , (-M^4 + 2 * M^3 - M^2 + 3 * M - 4) / (12 * (M - 1))
                
                ];
        end
 end
 end
 
 %Five stages, Method Order = 4, Embedded Orer = 3
 function Asflambda = AslowfastDepenFiveStageFourthOrder(lambda, M)
 
 if M == 1
     Asflambda = [0, 0, 0, 0, 0;
            (2/5), 0, 0, 0, 0;
            (-3/20), (3/4), 0, 0, 0;
            (19/44), (-15/44), (10/11), 0, 0;
            (11/72), (25/72), (25/72), (11/72), 0];
        
 elseif lambda == 1
     Asflambda = [...
         0,0,0,0,0; ...
         ...
         (2/5).*M,0,0,0,0; ...
         ...
         (3/20).* (4 + (-5) .* M) .* M, (3/4) .* M.^2, 0, 0, 0; ...
         ...
         (1/44).*M.*(44+(-81).*M+56.*M.^2),(-5/44).*M.^2.*((-13)+16.*M) ...
         ,(-5/11).*((-3)+M).*M.^2,((-1)+M).*M.^2,0; ...
         ...
         (11/72), (25/72), (25/72), (11/72), 0 ...
         ];
     
 else
     
     Asflambda = [ ...
         0,0,0,0,0; ...
         0,0,0,0,0; ...
         0,0,0,0,0; ...
         0,0,0,0,0; ...
         (11/72), (25/72), (25/72),(11/72),0 ...
         ];
     
 end
 
end

function Afslambda = AfastslowDepenFiveStageFourthOrder(lambda, M)

    if M == 1
        Afslambda = [0, 0, 0, 0, 0;
            (2/5), 0, 0, 0, 0;
            (-3/20), (3/4), 0, 0, 0;
            (19/44), (-15/44), (10/11), 0, 0;
            (11/72), (25/72), (25/72), (11/72), 0];
        
    elseif lambda == 1
        Afslambda = [ ...
            0,0,0,0,0; ...
            ...
            (2/5).*M.^(-1),0,0,0,0; ...
            ...
            (80 .* M +(-60) .* M.^2) .^ (-1) .* (3 + (-66) .* M + 90 .* M.^2 +(-30).*M.^3) ...
            ,(16 .* M + (-12) .* M.^2) .^ (-1) .* (9 + 6 .* M +(-18) .* M.^2 + 6 .* M.^3),0,0,0; ...
            ...
            0,(88 .* M + (-66) .* M.^2).^(-1).*(249 + (-348) .* M + 150 .* M.^2 + (-30).*M.^3), ...
            (-1) .* (88 .* M + (-66) .* M.^2).^(-1) .* (161 + (-282) .* M + 150 .* M.^2 + ...
            (-30).*M.^3),0,0; ...
            ...
            (11/72).*M.^(-1),(25/72).*M.^(-1),(25/72).*M.^(-1),(11/72).*M.^(-1),0];
    
    else
        %Super monster. Copied from Mathematica. pls no touch
        Afslambda = [...
            (11/72).*M.^(-1).*((-1)+lambda),(25/72).*M.^(-1).*((-1)+lambda) ...
            ,(25/72).* M.^(-1).*((-1)+lambda),(11/72).*M.^(-1).*((-1) + lambda),0;...
            ...
            (1/450).*((-1)+M).^(-1).*M.^(-1).*(776+(-450).*M.^2 ...
            +(-956).*lambda+M.*((-497)+956.*lambda)),...
            (1/450).*((-1)+M).^(-1).*M.^(-1).*(450.*M.^2+M.*(227+(-506).*lambda)...
            + 506.*((-1)+lambda)),0,0,0; ...
            ...
            (1/600).*M.^(-1).*(4+(-7).*M+3.* M.^2).^(-1).* ((-900).*M.^3 ...
            + 3.*M.^2.*(739+413.*lambda)+2.*((-781)+826.*lambda) ...
            +(-1).*M.*(97+2891.*lambda)),(1/600).*M.^(-1).*(4+(-7).*M + ...
            3.*M.^2).^(-1).*(602 + 900.*M.^3+M.*(1777+(-1309).*lambda) ...
            +748.*lambda+33.*M.^2.*((-89)+17.*lambda)),0,0,0;...
            ...
            0,(1/22).*M.^(-1).*(4+(-7).*M+3.*M.^2).^(-1).*((-90).*M.^3+3.* ...
            (39+44.*lambda)+M.^2.*(197+99.*lambda)+(-1).*M.*(205+231.*...
            lambda)),(1/792).*M.^(-1).*(4 +(-7).*M+3.*M.^2).^(-1).*(3240.*M.^3 ...
            +(-15).*M.^2.*(497+55.*lambda)+(-4).*(1174+275.*lambda) ...
            +M.*(8227+1925.*lambda)),(-11/72).*M.^(-1).*((-1)+lambda),0;...
            ...
            (11/72).*M.^(-1).*lambda,(25/72).*M.^(-1).*lambda,...
            (25/72).*M.^(-1).*lambda, (11/72).*M.^(-1).*lambda,0....
            ];
    end
end