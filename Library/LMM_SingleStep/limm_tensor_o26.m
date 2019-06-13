function [alpha, beta, gamma] = limm_tensor_o26(a_struct, b_struct, g_struct, nosptensor)
%

if ((exist('nosptensor', 'var') == 1 && nosptensor ~= 0) || exist('sptensor', 'file') ~= 2)
    ad_ten.idx = a_struct.denominator_i; ad_ten.vals = a_struct.denominator_v'; ad_ten.dim = a_struct.denominator_d;
    bd_ten.idx = b_struct.denominator_i; bd_ten.vals = b_struct.denominator_v'; bd_ten.dim = b_struct.denominator_d;
    gd_ten.idx = g_struct.denominator_i; gd_ten.vals = g_struct.denominator_v'; gd_ten.dim = g_struct.denominator_d;
    
    an_ten = {};
    for i = 1:length(a_struct.numerator_v)
      an_ten{i}.idx = a_struct.numerator_i{i};
      an_ten{i}.vals = a_struct.numerator_v{i}';
      an_ten{i}.dim = a_struct.numerator_d{i};
    end
    bn_ten = {};
    for i = 1:length(b_struct.numerator_v)
      bn_ten{i}.idx = b_struct.numerator_i{i};
      bn_ten{i}.vals = b_struct.numerator_v{i}';
      bn_ten{i}.dim = b_struct.numerator_d{i};
    end
    gn_ten = {};
    for i = 1:length(g_struct.numerator_v)
      gn_ten{i}.idx = g_struct.numerator_i{i};
      gn_ten{i}.vals = g_struct.numerator_v{i}';
      gn_ten{i}.dim = g_struct.numerator_d{i};
    end
else
    ad_ten = sptensor(a_struct.denominator_i, a_struct.denominator_v', a_struct.denominator_d);
    bd_ten = sptensor(b_struct.denominator_i, b_struct.denominator_v', b_struct.denominator_d);
    gd_ten = sptensor(g_struct.denominator_i, g_struct.denominator_v', g_struct.denominator_d);

    an_ten = {};
    for i = 1:length(a_struct.numerator_v)
      an_ten{i} = sptensor(a_struct.numerator_i{i}, a_struct.numerator_v{i}', a_struct.numerator_d{i});
    end
    bn_ten = {};
    for i = 1:length(b_struct.numerator_v)
      bn_ten{i} = sptensor(b_struct.numerator_i{i}, b_struct.numerator_v{i}', b_struct.numerator_d{i});
    end
    gn_ten = {};
    for i = 1:length(g_struct.numerator_v)
      gn_ten{i} = sptensor(g_struct.numerator_i{i}, g_struct.numerator_v{i}', g_struct.numerator_d{i});
    end
end

alpha.d_ten = ad_ten;
alpha.n_ten = an_ten;
beta.d_ten = bd_ten;
beta.n_ten = bn_ten;
gamma.d_ten = gd_ten;
gamma.n_ten = gn_ten;

end

