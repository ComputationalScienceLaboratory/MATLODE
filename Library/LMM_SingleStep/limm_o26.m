function [alpha, beta, gamma] = limm_o26(k)
%
switch(k)
  case 1
    coeffs = load('limm_o26_o1.mat');
    alpha.denominator_i = coeffs.alphaDenominatorPos;
    alpha.denominator_v = coeffs.alphaDenominatorVal;
    alpha.denominator_d = coeffs.alphaDenominatorDim;
    alpha.numerator_i = {coeffs.alphaNumeratorPos1};
    alpha.numerator_v = {coeffs.alphaNumeratorVal1};
    alpha.numerator_d = {coeffs.alphaNumeratorDim1};
    beta.denominator_i = coeffs.betaDenominatorPos;
    beta.denominator_v = coeffs.betaDenominatorVal;
    beta.denominator_d = coeffs.betaDenominatorDim;
    beta.numerator_i = {coeffs.betaNumeratorPos1};
    beta.numerator_v = {coeffs.betaNumeratorVal1};
    beta.numerator_d = {coeffs.betaNumeratorDim1};
    gamma.denominator_i = coeffs.gammaDenominatorPos;
    gamma.denominator_v = coeffs.gammaDenominatorVal;
    gamma.denominator_d = coeffs.gammaDenominatorDim;
    gamma.numerator_i = {coeffs.gammaNumeratorPos1, coeffs.gammaNumeratorPos2};
    gamma.numerator_v = {coeffs.gammaNumeratorVal1, coeffs.gammaNumeratorVal2};
    gamma.numerator_d = {coeffs.gammaNumeratorDim1, coeffs.gammaNumeratorDim2};
  case 2
    coeffs = load('limm_o26_o2.mat');
    alpha.denominator_i = coeffs.alphaDenominatorPos;
    alpha.denominator_v = coeffs.alphaDenominatorVal;
    alpha.denominator_d = coeffs.alphaDenominatorDim;
    alpha.numerator_i = {coeffs.alphaNumeratorPos1, coeffs.alphaNumeratorPos2};
    alpha.numerator_v = {coeffs.alphaNumeratorVal1, coeffs.alphaNumeratorVal2};
    alpha.numerator_d = {coeffs.alphaNumeratorDim1, coeffs.alphaNumeratorDim2};
    beta.denominator_i = coeffs.betaDenominatorPos;
    beta.denominator_v = coeffs.betaDenominatorVal;
    beta.denominator_d = coeffs.betaDenominatorDim;
    beta.numerator_i = {coeffs.betaNumeratorPos1, coeffs.betaNumeratorPos2};
    beta.numerator_v = {coeffs.betaNumeratorVal1, coeffs.betaNumeratorVal2};
    beta.numerator_d = {coeffs.betaNumeratorDim1, coeffs.betaNumeratorDim2};
    gamma.denominator_i = coeffs.gammaDenominatorPos;
    gamma.denominator_v = coeffs.gammaDenominatorVal;
    gamma.denominator_d = coeffs.gammaDenominatorDim;
    gamma.numerator_i = {coeffs.gammaNumeratorPos1, coeffs.gammaNumeratorPos2, coeffs.gammaNumeratorPos3};
    gamma.numerator_v = {coeffs.gammaNumeratorVal1, coeffs.gammaNumeratorVal2, coeffs.gammaNumeratorVal3};
    gamma.numerator_d = {coeffs.gammaNumeratorDim1, coeffs.gammaNumeratorDim2, coeffs.gammaNumeratorDim3};
  case 3
    coeffs = load('limm_o26_o3.mat');
    alpha.denominator_i = coeffs.alphaDenominatorPos;
    alpha.denominator_v = coeffs.alphaDenominatorVal;
    alpha.denominator_d = coeffs.alphaDenominatorDim;
    alpha.numerator_i = {coeffs.alphaNumeratorPos1, coeffs.alphaNumeratorPos2, coeffs.alphaNumeratorPos3};
    alpha.numerator_v = {coeffs.alphaNumeratorVal1, coeffs.alphaNumeratorVal2, coeffs.alphaNumeratorVal3};
    alpha.numerator_d = {coeffs.alphaNumeratorDim1, coeffs.alphaNumeratorDim2, coeffs.alphaNumeratorDim3};
    beta.denominator_i = coeffs.betaDenominatorPos;
    beta.denominator_v = coeffs.betaDenominatorVal;
    beta.denominator_d = coeffs.betaDenominatorDim;
    beta.numerator_i = {coeffs.betaNumeratorPos1, coeffs.betaNumeratorPos2, coeffs.betaNumeratorPos3};
    beta.numerator_v = {coeffs.betaNumeratorVal1, coeffs.betaNumeratorVal2, coeffs.betaNumeratorVal3};
    beta.numerator_d = {coeffs.betaNumeratorDim1, coeffs.betaNumeratorDim2, coeffs.betaNumeratorDim3};
    gamma.denominator_i = coeffs.gammaDenominatorPos;
    gamma.denominator_v = coeffs.gammaDenominatorVal;
    gamma.denominator_d = coeffs.gammaDenominatorDim;
    gamma.numerator_i = {coeffs.gammaNumeratorPos1, coeffs.gammaNumeratorPos2, coeffs.gammaNumeratorPos3, coeffs.gammaNumeratorPos4};
    gamma.numerator_v = {coeffs.gammaNumeratorVal1, coeffs.gammaNumeratorVal2, coeffs.gammaNumeratorVal3, coeffs.gammaNumeratorVal4};
    gamma.numerator_d = {coeffs.gammaNumeratorDim1, coeffs.gammaNumeratorDim2, coeffs.gammaNumeratorDim3, coeffs.gammaNumeratorDim4};
  case 4
    coeffs = load('limm_o26_o4.mat');
    alpha.denominator_i = coeffs.alphaDenominatorPos;
    alpha.denominator_v = coeffs.alphaDenominatorVal;
    alpha.denominator_d = coeffs.alphaDenominatorDim;
    alpha.numerator_i = {coeffs.alphaNumeratorPos1, coeffs.alphaNumeratorPos2, coeffs.alphaNumeratorPos3, coeffs.alphaNumeratorPos4};
    alpha.numerator_v = {coeffs.alphaNumeratorVal1, coeffs.alphaNumeratorVal2, coeffs.alphaNumeratorVal3, coeffs.alphaNumeratorVal4};
    alpha.numerator_d = {coeffs.alphaNumeratorDim1, coeffs.alphaNumeratorDim2, coeffs.alphaNumeratorDim3, coeffs.alphaNumeratorDim4};
    beta.denominator_i = coeffs.betaDenominatorPos;
    beta.denominator_v = coeffs.betaDenominatorVal;
    beta.denominator_d = coeffs.betaDenominatorDim;
    beta.numerator_i = {coeffs.betaNumeratorPos1, coeffs.betaNumeratorPos2, coeffs.betaNumeratorPos3, coeffs.betaNumeratorPos4};
    beta.numerator_v = {coeffs.betaNumeratorVal1, coeffs.betaNumeratorVal2, coeffs.betaNumeratorVal3, coeffs.betaNumeratorVal4};
    beta.numerator_d = {coeffs.betaNumeratorDim1, coeffs.betaNumeratorDim2, coeffs.betaNumeratorDim3, coeffs.betaNumeratorDim4};
    gamma.denominator_i = coeffs.gammaDenominatorPos;
    gamma.denominator_v = coeffs.gammaDenominatorVal;
    gamma.denominator_d = coeffs.gammaDenominatorDim;
    gamma.numerator_i = {coeffs.gammaNumeratorPos1, coeffs.gammaNumeratorPos2, coeffs.gammaNumeratorPos3, coeffs.gammaNumeratorPos4, coeffs.gammaNumeratorPos5};
    gamma.numerator_v = {coeffs.gammaNumeratorVal1, coeffs.gammaNumeratorVal2, coeffs.gammaNumeratorVal3, coeffs.gammaNumeratorVal4, coeffs.gammaNumeratorVal5};
    gamma.numerator_d = {coeffs.gammaNumeratorDim1, coeffs.gammaNumeratorDim2, coeffs.gammaNumeratorDim3, coeffs.gammaNumeratorDim4, coeffs.gammaNumeratorDim5};
  case 5
    coeffs = load('limm_o26_o5.mat');
    alpha.denominator_i = coeffs.alphaDenominatorPos;
    alpha.denominator_v = coeffs.alphaDenominatorVal;
    alpha.denominator_d = coeffs.alphaDenominatorDim;
    alpha.numerator_i = {coeffs.alphaNumeratorPos1, coeffs.alphaNumeratorPos2, coeffs.alphaNumeratorPos3, coeffs.alphaNumeratorPos4, coeffs.alphaNumeratorPos5};
    alpha.numerator_v = {coeffs.alphaNumeratorVal1, coeffs.alphaNumeratorVal2, coeffs.alphaNumeratorVal3, coeffs.alphaNumeratorVal4, coeffs.alphaNumeratorVal5};
    alpha.numerator_d = {coeffs.alphaNumeratorDim1, coeffs.alphaNumeratorDim2, coeffs.alphaNumeratorDim3, coeffs.alphaNumeratorDim4, coeffs.alphaNumeratorDim5};
    beta.denominator_i = coeffs.betaDenominatorPos;
    beta.denominator_v = coeffs.betaDenominatorVal;
    beta.denominator_d = coeffs.betaDenominatorDim;
    beta.numerator_i = {coeffs.betaNumeratorPos1, coeffs.betaNumeratorPos2, coeffs.betaNumeratorPos3, coeffs.betaNumeratorPos4, coeffs.betaNumeratorPos5};
    beta.numerator_v = {coeffs.betaNumeratorVal1, coeffs.betaNumeratorVal2, coeffs.betaNumeratorVal3, coeffs.betaNumeratorVal4, coeffs.betaNumeratorVal5};
    beta.numerator_d = {coeffs.betaNumeratorDim1, coeffs.betaNumeratorDim2, coeffs.betaNumeratorDim3, coeffs.betaNumeratorDim4, coeffs.betaNumeratorDim5};
    gamma.denominator_i = coeffs.gammaDenominatorPos;
    gamma.denominator_v = coeffs.gammaDenominatorVal;
    gamma.denominator_d = coeffs.gammaDenominatorDim;
    gamma.numerator_i = {coeffs.gammaNumeratorPos1, coeffs.gammaNumeratorPos2, coeffs.gammaNumeratorPos3, coeffs.gammaNumeratorPos4, coeffs.gammaNumeratorPos5, coeffs.gammaNumeratorPos6};
    gamma.numerator_v = {coeffs.gammaNumeratorVal1, coeffs.gammaNumeratorVal2, coeffs.gammaNumeratorVal3, coeffs.gammaNumeratorVal4, coeffs.gammaNumeratorVal5, coeffs.gammaNumeratorVal6};
    gamma.numerator_d = {coeffs.gammaNumeratorDim1, coeffs.gammaNumeratorDim2, coeffs.gammaNumeratorDim3, coeffs.gammaNumeratorDim4, coeffs.gammaNumeratorDim5, coeffs.gammaNumeratorDim6};
  case 6
    coeffs = load('limm_o26_o6.mat');
    alpha.denominator_i = coeffs.alphaDenominatorPos;
    alpha.denominator_v = coeffs.alphaDenominatorVal;
    alpha.denominator_d = coeffs.alphaDenominatorDim;
    alpha.numerator_i = {coeffs.alphaNumeratorPos1, coeffs.alphaNumeratorPos2, coeffs.alphaNumeratorPos3, coeffs.alphaNumeratorPos4, coeffs.alphaNumeratorPos5, coeffs.alphaNumeratorPos6};
    alpha.numerator_v = {coeffs.alphaNumeratorVal1, coeffs.alphaNumeratorVal2, coeffs.alphaNumeratorVal3, coeffs.alphaNumeratorVal4, coeffs.alphaNumeratorVal5, coeffs.alphaNumeratorVal6};
    alpha.numerator_d = {coeffs.alphaNumeratorDim1, coeffs.alphaNumeratorDim2, coeffs.alphaNumeratorDim3, coeffs.alphaNumeratorDim4, coeffs.alphaNumeratorDim5, coeffs.alphaNumeratorDim6};
    beta.denominator_i = coeffs.betaDenominatorPos;
    beta.denominator_v = coeffs.betaDenominatorVal;
    beta.denominator_d = coeffs.betaDenominatorDim;
    beta.numerator_i = {coeffs.betaNumeratorPos1, coeffs.betaNumeratorPos2, coeffs.betaNumeratorPos3, coeffs.betaNumeratorPos4, coeffs.betaNumeratorPos5, coeffs.betaNumeratorPos6};
    beta.numerator_v = {coeffs.betaNumeratorVal1, coeffs.betaNumeratorVal2, coeffs.betaNumeratorVal3, coeffs.betaNumeratorVal4, coeffs.betaNumeratorVal5, coeffs.betaNumeratorVal6};
    beta.numerator_d = {coeffs.betaNumeratorDim1, coeffs.betaNumeratorDim2, coeffs.betaNumeratorDim3, coeffs.betaNumeratorDim4, coeffs.betaNumeratorDim5, coeffs.betaNumeratorDim6};
    gamma.denominator_i = coeffs.gammaDenominatorPos;
    gamma.denominator_v = coeffs.gammaDenominatorVal;
    gamma.denominator_d = coeffs.gammaDenominatorDim;
    gamma.numerator_i = {coeffs.gammaNumeratorPos1, coeffs.gammaNumeratorPos2, coeffs.gammaNumeratorPos3, coeffs.gammaNumeratorPos4, coeffs.gammaNumeratorPos5, coeffs.gammaNumeratorPos6, coeffs.gammaNumeratorPos7};
    gamma.numerator_v = {coeffs.gammaNumeratorVal1, coeffs.gammaNumeratorVal2, coeffs.gammaNumeratorVal3, coeffs.gammaNumeratorVal4, coeffs.gammaNumeratorVal5, coeffs.gammaNumeratorVal6, coeffs.gammaNumeratorVal7};
    gamma.numerator_d = {coeffs.gammaNumeratorDim1, coeffs.gammaNumeratorDim2, coeffs.gammaNumeratorDim3, coeffs.gammaNumeratorDim4, coeffs.gammaNumeratorDim5, coeffs.gammaNumeratorDim6, coeffs.gammaNumeratorDim7};
  otherwise
    error('bad input k');
end
end
