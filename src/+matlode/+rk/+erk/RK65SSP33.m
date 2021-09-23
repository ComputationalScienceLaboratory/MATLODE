classdef RK65SSP33 < matlode.rk.erk.ERK
    properties (SetAccess = immutable)
        DenseOut
    end
    
    methods
        function obj = RK65SSP33(datatype)
            
            % A strong stability preserving RK method
            % Reference : Page 5.6, Table 4.12 a) https://core.ac.uk/download/pdf/56372752.pdf
            % Macdonald, Colin Barr. "Constructing high-order Runge-Kutta methods with embedded strong-stability-preserving pairs." 
            % PhD diss., Theses (Dept. of Mathematics)/Simon Fraser University, 2003.
            
            if nargin == 0
                datatype = 'double';
            end
            
            caster = @(x) matlode.util.CoefficentTransformers.transform(x,datatype);
            
            
            a = caster(join(['[[0,            0,            0,            0,            0,            0];',...
                             '[1,             0,            0,            0,            0,            0];',...
                             '[1/4,           1/4,          0,            0,            0,            0];',...
                             '[2046/15625,    -454/15625,   1533/15625,   0,            0,            0];',...
                             '[-739/5625,     511/5625,     -566/16875,   20/27,        0,            0];',...
                             '[11822/21875,   -6928/21875,  -4269/21875,  -4/7,         54/35,        0]]']));
             
            b =  caster('[1/6, 1/6, 2/3, 0, 0, 0]');
            
            bHat = caster('[1/24, 0, 0, 125/336, 27/56, 5/48]');
            
            e = caster('[1/8, 1/6, 2/3, -125/336, -27/56, -5/48]');
            
            c = caster('[0, 1, 1/2, 1/5, 2/3, 1]');
            
            order = 3;
            
            embbededOrder = 5;
            
            obj = obj@matlode.rk.erk.ERK(a, b, bHat, c, e, order, embbededOrder);
            
            % Fix the Dense thing.
            
%             %Shampien formula for dense output given in the Solving ODE
%             %book on pg192 (6.12)
%             %the Dense output has been transformed from a function form to
%             %a matrix form
%             %p* == 4 s* == 5
%             denseMatrix = caster(join(['[[1, -(4034104133/1410260304), 105330401/33982176,      -(13107642775/11282082432), 6542295/470086768];',...
%                                        '[ 0, 0,                        0,                       0,                          0];',...
%                                        '[ 0, 132343189600/32700410799, -(833316000/131326951),  91412856700/32700410799,    -(523383600/10900136933)];',...
%                                        '[ 0, -(115792950/29380423),    185270875/16991088,      -(12653452475/1880347072),  98134425/235043384];',...
%                                        '[ 0, 70805911779/24914598704,  -(4531260609/600351776), 988140236175/199316789632,  -(14307999165/24914598704)];',...
%                                        '[ 0, -(331320693/205662961),   31361737/7433601,        -(2426908385/822651844),    97305120/205662961];',...
%                                        '[ 0, 44764047/29380423,        -(1532549/353981),       90730570/29380423,          -(8293050/29380423)]]']));
%             obj.DenseOut = matlode.denseoutput.RKDense(denseMatrix);
           
        end
        
        function sol = integrate(obj, f, tspan, y0, varargin)
            sol = integrate@matlode.rk.erk.ERK(obj, f, tspan, y0, 'StepSizeController', matlode.stepsizecontroller.StandardController, 'Dense', obj.DenseOut, varargin{:});
        end
    end
end

