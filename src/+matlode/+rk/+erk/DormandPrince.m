classdef DormandPrince < matlode.rk.erk.ERK
%Method: Dormand PRince
% p = 5 s = 7 pe = 4
%Reference: 
    
    methods

        function obj = DormandPrince(datatype)
            arguments
				datatype(1,1) string = 'double';
			end
            
            caster = @(x) matlode.util.CoefficentTransformers.transform(x,datatype);
            
            
            a = caster(join(['[[0,          0,           0,          0,        0,           0,     0];',...
                              '[1/5,        0,           0,          0,        0,           0,     0];',...
                              '[3/40,       9/40,        0,          0,        0,           0,     0];',...
                              '[44/45,      -56/15,      32/9,       0,        0,           0,     0];',...
                              '[19372/6561, -25360/2187, 64448/6561, -212/729, 0,           0,     0];',...
                              '[9017/3168,  -355/33,     46732/5247, 49/176,   -5103/18656, 0,     0];',...
                              '[35/384,     0,           500/1113,   125/192,  -2187/6784,  11/84, 0]]']));
             
            b =  caster('[35/384, 0, 500/1113, 125/192, -2187/6784, 11/84, 0]');
            
            bHat = caster('[5179/57600, 0, 7571/16695, 393/640, -92097/339200, 187/2100, 1/40]');
            
            e = caster('[71/57600, 0, -71/16695, 71/1920, -17253/339200, 22/525, -1/40]');
            
            c = caster('[0, 1/5, 3/10, 4/5, 8/9, 1, 1]');
            
            order = 5;
            
            embbededOrder = 4;
            
            obj = obj@matlode.rk.erk.ERK(a, b, bHat, c, e, order, embbededOrder);
            
            %Shampien formula for dense output given in the Solving ODE
            %book on pg192 (6.12)
            %the Dense output has been transformed from a function form to
            %a matrix form
            %p* == 4 s* == 5
            denseMatrix = caster(join(['[[1, -(4034104133/1410260304), 105330401/33982176,      -(13107642775/11282082432), 6542295/470086768];',...
                                       '[ 0, 0,                        0,                       0,                          0];',...
                                       '[ 0, 132343189600/32700410799, -(833316000/131326951),  91412856700/32700410799,    -(523383600/10900136933)];',...
                                       '[ 0, -(115792950/29380423),    185270875/16991088,      -(12653452475/1880347072),  98134425/235043384];',...
                                       '[ 0, 70805911779/24914598704,  -(4531260609/600351776), 988140236175/199316789632,  -(14307999165/24914598704)];',...
                                       '[ 0, -(331320693/205662961),   31361737/7433601,        -(2426908385/822651844),    97305120/205662961];',...
                                       '[ 0, 44764047/29380423,        -(1532549/353981),       90730570/29380423,          -(8293050/29380423)]]']));
            obj.DenseOut = matlode.denseoutput.RKDense(denseMatrix);
           
        end
        
        function sol = integrate(obj, f, tspan, y0, varargin)
            sol = integrate@matlode.rk.erk.ERK(obj, f, tspan, y0, 'StepSizeController', matlode.stepsizecontroller.StandardController, 'Dense', obj.DenseOut, varargin{:});
        end
    end
end

