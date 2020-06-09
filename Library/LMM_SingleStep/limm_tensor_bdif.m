function [c, d, e, gamma, r, err] = limm_tensor_bdif(c_struct, d_struct, e_struct, g_struct, r_struct, err_struct, nosptensor)
%
if (exist('nosptensor', 'var') ~= 1)
    nosptensor = 0;
end

tlist = limm_tensor_helper({c_struct, d_struct, e_struct, g_struct, r_struct, err_struct}, nosptensor);

c = tlist{1};
d = tlist{2};
e = tlist{3};
gamma = tlist{4};
r = tlist{5};
err = tlist{6};

end

