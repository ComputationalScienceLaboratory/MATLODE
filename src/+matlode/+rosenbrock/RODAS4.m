classdef RODAS4 < matlode.rosenbrock.Rosenbrock
%Method: RODAS4
% p = 4 s = 6 pe = 3
%Reference: 
    
    methods
        function obj = RODAS4(datatype)
            arguments
				datatype(1,1) string = 'double';
			end
            
            caster = @(x) matlode.util.CoefficentTransformers.transform(x,datatype);
			
            gammadia = caster('[1/4, 1/4, 1/4, 1/4, 1/4, 1/4]');
            alphasum = caster('[0, 0.386, 0.210, 0.630, 1, 1]');
            gammasum = caster('[1/4, -0.1043, 0.1035, -0.0362, 0, 0]');
            
            a = caster(join(['[[0,0,0,0,0,0];',...
							 '[1.544,0,0,0,0,0];',...
                             '[0.9466785280815826, 0.2557011698983284, 0, 0, 0, 0];',...
                             '[0.3314825187068521e1, 0.2896124015972201e1, 0.9986419139977817, 0, 0, 0];',...
                             '[0.1221224509226641e1, 0.6019134481288629e1, 0.1253708332932087e2, -0.6878860361058950, 0, 0];',...
							 '[0.1221224509226641d+01, 0.6019134481288629e1, 0.1253708332932087e2,  -0.6878860361058950, 1, 0]]']));
            
            c = caster(join(['[[0, 0, 0, 0, 0, 0];',...
                             '[-0.5668800000000000e1, 0, 0, 0, 0, 0];',...
                             '[-0.2430093356833875e1, -0.2063599157091915, 0, 0, 0, 0];',...
                             '[-0.1073529058151375, -0.9594562251023355e1, -0.2047028614809616e2, 0, 0, 0];',...
							 '[0.7496443313967647e1, -0.1024680431464352e2, -0.3399990352819905e2, 0.1170890893206160e2, 0, 0];',...
							 '[0.8083246795921522e1, -0.7981132988064893e1, -0.3152159432874371e2, 0.1631930543123136e2, -0.6058818238834054e1, 0]]']));
			
            m = caster('[0.1221224509226641d+01, 0.6019134481288629e1, 0.1253708332932087e2,  -0.6878860361058950, 1, 1]');
            e = caster('[0, 0, 0, 0, 0, 1]');
            
            order = 4;
            
            embbededOrder = 3;
            
            obj = obj@matlode.rosenbrock.Rosenbrock(gammadia, gammasum, alphasum, a, c, m, e, order, embbededOrder);
            
        end
        
        function sol = integrate(obj, f, tspan, y0, varargin)
            sol = integrate@matlode.rosenbrock.Rosenbrock(obj, f, tspan, y0, 'StepSizeController', matlode.stepsizecontroller.StandardController, 'Dense', obj.DenseOut, varargin{:});
        end
    end
end

