function validFunc = scalarValidationFunc()
msg = 'Value must be a scalar';
validFunc = @(x) assert(isscalar(x), msg);
end

