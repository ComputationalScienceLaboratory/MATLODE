function [alpha, beta, gamma] = limm_tensor(a_struct, b_struct, g_struct, nosptensor)
%
if (exist('nosptensor', 'var') ~= 1)
    nosptensor = 0;
end

tlist = limm_tensor_helper({a_struct, b_struct, g_struct}, nosptensor);

alpha = tlist{1};
beta = tlist{2};
gamma = tlist{3};

end

