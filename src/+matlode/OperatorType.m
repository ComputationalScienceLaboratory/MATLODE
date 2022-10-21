classdef OperatorType
	%Defines operator type from R x F^n -> F^n
	% Or in the implicit case for jets R x F^n x F^n -> F^n
    enumeration
        Zero, Constant, TimeDependent, StateDependent, JetDependent, TSDep, SJDep, TJDep, TSJDep
    end

end

