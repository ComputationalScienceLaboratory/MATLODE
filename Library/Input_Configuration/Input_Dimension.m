function [ OPTIONS ] = Input_Dimension(Tspan, Y0, OdeFunction, OPTIONS)
%% The goal of this function is to update
% OPTIONS.NVAR (number of state variables),
% OPTIONS.NP (number of parameters used in the adjoint algorithm)
% OPTIONS.NADJ (number of cost function in the adjoint algorithm)
% OPTIONS.NTLM (number of parameters used in the TLM algorithm)
% and make sure of the dimension of the following functions or states:
% Y0 should be NVAR*1
% DRDP should return NADJ*NP
% DRDY should return NADJ*NVAR
% Jacp should return NVAR*NP
% Jacobian should return NVAR*NVAR
% OdeFunction should return NVAR*1
% Lambda_T should return NADJ*NVAR (d(terminal costfn)/d(state))
% Mu_T should return max(NADJ,1)*NP (d(terminal costfn)/d(parameters))

    if isempty(Tspan)
        error('The time span should be provided.');
    elseif ~(iscolumn(Tspan) || isrow(Tspan)) || numel(Tspan) < 2 
        error('The time span should be a vector with at least two elements.');
    end
    
    t0 = Tspan(1);

    %% For any integrator, we check that Y0 is a NVAR*1 vector
    if isempty(Y0)
        error('The Initial vector Y0 for the ODE system should be provided.')
    end

    if ~iscolumn(Y0)
        error('The Initial vector Y0 should be a column vector and cannot be a matrix')
    end

    OPTIONS.NVAR = size(Y0, 1);

    % For any integrator, we check that the ODE function returns a NVAR*1 vector
    if ~isa(OdeFunction, 'function_handle')
        error('OdeFunction function should be a function handle like: OdeFunction=(t,Y)function_name(t,Y,p1,p2,..).')
    end

    % Check the number of arguments the function accepts
    if nargin(OdeFunction) ~= 2
        error('ODE function should have exactly two arguments like: OdeFunction=(t,Y)function_name(t,Y,p1,p2,..).')
    end

    try
        Yp0 = OdeFunction(t0, Y0);
    catch ME
        error(strcat('Calling ODE function passing (t0, Y0) in order returned the following error:', ME.message));
    end

    if ~iscolumn(Yp0) || size(Yp0, 1) ~= OPTIONS.NVAR
        error('The ODE system should return a column vector NVAR*1')
    end

    %% For any integrator that requires the Jacobian of the ODE system, we check
    % it returns a NVAR*NVAR matrix
    if  ~isempty(OPTIONS.Jacobian)
        if ~isa(OPTIONS.Jacobian, 'function_handle')
            error('Jacobian function should be a function handle like: OPTIONS.Jacobian=(t,Y)function_name(t,Y,p1,p2,..).')
        end

        if nargin(OPTIONS.Jacobian) == 2
            try
                Jac = OPTIONS.Jacobian(t0, Y0);
            catch ME
                error(strcat('Calling Jacobian function passing (t0, Y0) in order returned the following error:', ME.message));
            end

            if size(Jac, 1) ~= OPTIONS.NVAR || size(Jac, 2) ~= OPTIONS.NVAR
                error('Jacobian function should return a NVAR*NVAR matrix')
            end
        elseif  nargin(OPTIONS.Jacobian) == 3
            if ~OPTIONS.MatrixFree
                error('Three argument Jacobian vector product functions require MatrixFree option to be set.')             
            end

            try
                randvec = rand(OPTIONS.NVAR, 1);
                Jv = OPTIONS.Jacobian(t0, Y0, randvec);
            catch ME
                error(strcat('Calling Jacobian function passing (t0, Y0, v) in order returned the following error:', ME.message));
            end

            if ~iscolumn(Jv) || size(Jv, 1) ~= OPTIONS.NVAR
                error('Jacobian-vector product function should return a column vector of size NVAR*1')
            end
        else
            error('Jacobian function should accept two arguments or Jacobian vector product functions should accept three arguments.')
        end
    end


    %% For any TLM integrator, we check that Y_TLM is provided and should be of
    % dimension NVAR*NTLM
    if strcmp(OPTIONS.Implementation, 'TLM')
        if ~isempty(OPTIONS.Y_TLM)
            Y_TLM = OPTIONS.Y_TLM;
            
            if size(Y_TLM, 1) ~= OPTIONS.NVAR
                error('Dimension of rows for Y_TLM should be NVAR')
            else
                OPTIONS.NTLM = size(Y_TLM, 2);
            end
        else
            error('Y_TLM should be provided in TLM integrator.')
        end
    end


    %%  For any Adjoint  integrator, we check and set
    % OPTIONS.NVAR (number of state variables),
    % OPTIONS.NP (number of parameters used in the adjoint algorithm)
    % OPTIONS.NADJ (number of cost function in the adjoint algorithm)
    % DRDP should return NADJ*NP
    % DRDY should return NADJ*NVAR
    % Jacp should return NVAR*NP
    % Lambda should return NADJ*NVAR (d(terminal costfn)/d(state))
    % Mu should return NADJ*NP (d(terminal costfn)/d(parameters))
    if strcmp(OPTIONS.Implementation, 'ADJ')
        if ~isempty(OPTIONS.Lambda)
            if ~isa(OPTIONS.Lambda,'function_handle')
                error('Lambda function should be a function handle like: OPTIONS.Lambda=(t,Y)function_name(t,Y,p1,p2,..).')
            end
            if nargin(OPTIONS.Lambda) ~= 2
                error('Lambda function should have exactly two arguments like: OPTIONS.Lambda=(t,Y)function_name(t,Y,p1,p2,..).')
            end
            
            try
                Lambda = OPTIONS.Lambda(t0,Y0);
            catch ME
                error(strcat('Calling Lambda function with parameters (t0, Y0) returned the following error:', ME.message));
            end
            
            if size(Lambda, 1) ~= OPTIONS.NVAR
                error('Lambda function should return a NVAR*NADJ matrix')
            else
                OPTIONS.NADJ = size(Lambda, 2);
            end
        else
            error('Lambda function should be provided like: OPTIONS.Lambda=(t,Y)function_name(t,Y,p1,p2,..).')
        end
        
        if ~isempty(OPTIONS.Mu)
            if ~isa(OPTIONS.Mu,'function_handle')
                error('Mu function should be a function handle like: OPTIONS.Mu=(t,Y)function_name(t,Y,p1,p2,..).')
            end
            if nargin(OPTIONS.Mu) ~= 2
                warning('Mu function should have exactly two arguments like: OPTIONS.Mu=(t,Y)function_name(t,Y,p1,p2,..).')
            end
            
            try
                Mu = OPTIONS.Mu(t0, Y0);
            catch ME
                error(strcat('Calling Mu function with parameters (t0, Y0) returned the following error:', ME.message));
            end
            
            if size(Mu, 1) ~= OPTIONS.NADJ
                error('Mu function should return a NADJ*NP matrix')
            else
                OPTIONS.NP = size(Mu, 2);
            end
        end
        
        if ~isempty(OPTIONS.Jacp)
            if ~isa(OPTIONS.Jacp,'function_handle')
                error('Jacp function should be a function handle like: OPTIONS.Jacp=(t,Y)function_name(t,Y,p1,p2,..).')
            end

            if nargin(OPTIONS.Jacp) ~= 2
                error('Jacp function should have exactly two arguments like: OPTIONS.Jacp=(t,Y)function_name(t,Y,p1,p2,..).')
            end
            
            try
                Jacp = OPTIONS.Jacp(t0, Y0);
            catch ME
                error(strcat('Calling Jacp function with parameters (t0, Y0) returned the following error:', ME.message));
            end
            
            if size(Jacp, 1) ~= OPTIONS.NVAR && size(Jacp, 2) ~= OPTIONS.NP
                error('Jacp function  should return a NVAR*NP matrix')
            end
        end
        
        if ~isempty(OPTIONS.DRDP)
            if ~isa(OPTIONS.DRDP,'function_handle')
                error('DRDP function should be a function handle like: OPTIONS.DRDP=(t,Y)function_name(t,Y,p1,p2,..).')
            end
            if nargin(OPTIONS.DRDP) ~= 2
                warning('DRDP function should have exactly two arguments like: OPTIONS.DRDP=(t,Y)function_name(t,Y,p1,p2,..).')
            end
            
            try
                DRDP = OPTIONS.DRDP(t0, Y0);
            catch ME
                error(strcat('Calling DRDP function with parameters (t0, Y0) returned the following error:', ME.message));
            end
            
            if size(DRDP, 1) ~= OPTIONS.NADJ && size(DRDP, 2) ~= OPTIONS.NP
                error('DRDP function needs to return a  NADJ*NP matrix')
            end
        end
        
        if ~isempty(OPTIONS.DRDY)
            if ~isa(OPTIONS.DRDY ,'function_handle')
                error('DRDY function should be a function handle like: OPTIONS.DRDY=(t,Y)function_name(t,Y,p1,p2,..).')
            end

            if nargin(OPTIONS.DRDY) ~= 2
                warning('DRDY should have exactly two arguments like: OPTIONS.DRDY=(t,Y)function_name(t,Y,p1,p2,..).')
            end
            
            try
                DRDY = OPTIONS.DRDY(t0, Y0);
            catch ME
                error(strcat('Calling DRDY function with parameters (t0, Y0) returned the following error:', ME.message));
            end
            
            if size(DRDY, 1) ~= OPTIONS.NADJ && size(DRDY, 2) ~= OPTIONS.NVAR
                error('DRDY function shoud return a NADJ*NVAR matrix')
            end
        end
    end
end
