function [ tlist ] = limm_tensor_helper( t_struct_list, nosptensor )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

if (exist('nosptensor', 'var') ~= 1)
    nosptensor = 0;
end

if ~iscell(t_struct_list)
    tlist = limm_tensor_helper({t_struct_list}, nosptensor);
    return;
end

tlist = {};
if (nosptensor ~= 0 || exist('sptensor', 'file') ~= 2)
    % homemade tensors
    for i = 1:length(t_struct_list)
        dten.idx = t_struct_list{i}.denominator_i;
        dten.vals = t_struct_list{i}.denominator_v';
        dten.dim = t_struct_list{i}.denominator_d;
        
        nten = {};
        for j = 1:length(t_struct_list{i}.numerator_v)
            nten{j}.idx = t_struct_list{i}.numerator_i{j};
            nten{j}.vals = t_struct_list{i}.numerator_v{j}';
            nten{j}.dim = t_struct_list{i}.numerator_d{j};
        end
        
        tlist{i} = struct;
        tlist{i}.d_ten = dten;
        tlist{i}.n_ten = nten;
    end
else
    % tensor toolbox
    for i = 1:length(t_struct_list)
        dten = sptensor(t_struct_list{i}.denominator_i, t_struct_list{i}.denominator_v', t_struct_list{i}.denominator_d);
        
        nten = {};
        for j = 1:length(t_struct_list{i}.numerator_v)
            nten{j} = sptensor(t_struct_list{i}.numerator_i{j}, t_struct_list{i}.numerator_v{j}', t_struct_list{i}.numerator_d{j});
        end
        
        tlist{i} = struct;
        tlist{i}.d_ten = dten;
        tlist{i}.n_ten = nten;
    end
end

end

