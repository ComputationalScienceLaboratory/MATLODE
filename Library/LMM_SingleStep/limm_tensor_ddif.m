function [c, d, e, gamma, err] = limm_tensor_ddif(c_struct, d_struct, e_struct, g_struct, r_struct, nosptensor)
%
if (exist('nosptensor', 'var') ~= 1)
    nosptensor = 0;
end

tlist = limm_tensor_helper({c_struct, d_struct, e_struct, g_struct, r_struct}, nosptensor);

c = tlist{1};
d = tlist{2};
e = tlist{3};
gamma = tlist{4};
err = tlist{5};

end

