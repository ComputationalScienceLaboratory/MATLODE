function [ vals ] = limm_tensor_eval(w, coef_tensor_struct)
% Evaluates a set of multivariable rational functions with a common denominator

 denominator = coef_tensor_struct.d_ten;
 numerators = coef_tensor_struct.n_ten;

 n = zeros(1,length(numerators));
 
 if isempty(w) && ~isobject(denominator)
     denominator = denominator.vals;
     for j = 1:length(numerators)
         numerators{j} = numerators{j}.vals;
     end
 else
     if (isobject(denominator))
         dprod = @ttv;
         siz = @size;
     else
         dprod = @tvprod;
         siz = @(t)t.dim;
     end

     d_size = siz(denominator);
     if length(d_size) == 1 && d_size(1) == 1
         denominator = dprod(denominator, [1], 1);
     else
         for i = length(d_size):-1:1
             mon = w(end-i+1).^[0:d_size(i)-1];
             denominator = dprod(denominator, mon', i);
         end
     end
     denominator = double(denominator);

     for j = 1:length(numerators)
         n_size = siz(numerators{j});
         if length(n_size) == 1 && n_size(1) == 1
             numerators{j} = dprod(numerators{j}, [1], 1);
         else
             for i = length(n_size):-1:1
                 mon = w(end-i+1).^[0:n_size(i)-1];
                 numerators{j} = dprod(numerators{j}, mon', i);
             end
         end
         numerators{j} = double(numerators{j});
         if isempty(numerators{j})
             numerators{j} = 0.0;
         end
     end
 end

 vals = cell2mat(numerators)./denominator;
 
 %%%%% DEBUG %%%%%%
 assert(length(vals) == length(coef_tensor_struct.n_ten));
 
end

