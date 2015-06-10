%%
%  Authored by Tony D'Augustine, Adrian Sandu, and Hong Zhang.
%  Computational Science Laboratory, Virginia Tech.
%  ©2015 Virginia Tech Intellectual Properties, Inc.
%
function [ fjac1 ] = lss_jac1( t, y, Jacobian )
%LSS_JAC1: DO NOT USE THIS FUNCTION. ONLY HERE TO MAKE TRANSLATING EASIER

    fjac1 = Jacobian(t,y);

end

