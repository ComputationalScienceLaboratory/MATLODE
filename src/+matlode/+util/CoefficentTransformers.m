classdef CoefficentTransformers
    methods (Static)
        
        function valuetype = transform(str, datatype)
            
            if isa(datatype, 'sym')
                valuetype = str2sym(str);
            else
                %Better than str2double as it can deal with arrays and sqrt
                valuetype = cast(str2num(str), datatype);
            end
        end
    end
end

