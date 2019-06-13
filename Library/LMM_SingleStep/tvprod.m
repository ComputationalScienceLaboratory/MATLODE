function [ out_ten ] = tvprod(tensor, vector, dim)
%tvprod 
%   
idxmask = true(size(tensor.dim)); idxmask(dim) = false;
vals = zeros(size(tensor.vals));
idx  = zeros(size(tensor.idx,1),size(tensor.idx,2)-1);
dims  = tensor.dim(idxmask);

% product
for i = 1:length(tensor.vals)
    vals(i) = vals(i) + tensor.vals(i) * vector(tensor.idx(i,dim));
    idx(i,:) = tensor.idx(i,idxmask);
end

% collapse repeated indices
[nidx, ~, loc] = unique(idx, 'rows');
nvals = accumarray(loc, vals, [size(nidx,1) 1]);

% remove zeros
zmask = nvals ~= 0.0;
nvals = nvals(zmask);
nidx  = nidx(zmask,:);

if isempty(nidx) || (size(nidx,1) == 1 && size(nidx,2) == 1 && nidx == 1)
    out_ten = nvals;
else
    out_ten.idx = nidx;
    out_ten.vals = nvals;
    out_ten.dim = dims;
end

end

