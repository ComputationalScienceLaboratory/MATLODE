classdef LinearOperatorType
	%Defines operators R x F^n -> F^{n x n}
    %There are multiple kinds of linear operators
    %This class simplifies check what type the linear operator is
    %Does not identify if singular or not
    %TODO: Add Jet Depedency
    enumeration
        Empty, Zero, Identity, Constant, TimeDependent, StateDependent, TSDepedent
	end
end

