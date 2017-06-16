function [ OPTIONS ] = Input_Dimension(Y0,OdeFunction,OPTIONS)
%% For any integrator

if  isempty(Y0)
    error('The Initial vector for the ODE system should not be empty')
end

if size(Y0,2)~=1
    error('Initial vector for the ODE system should be a column vector')
else
    OPTIONS.NVAR = size(Y0, 1);
end

if  isempty(OdeFunction)
    error('The ODE system should be provided')
end

% Check the number of arguments the function accepts
if nargin(OdeFunction) ~= 2
    warning('ODE function should have exactly two arguments.')
end

t0=0;    
Yp0=OdeFunction(t0,Y0);
if size(Yp0,2)~=1
  error('The ODE system should return a column vector')
end
if size(Yp0,1)~=OPTIONS.NVAR
  error('Dimension of the OdeFunction needs to be NVAR')   
end

%% Integrator that requires the Jacobian of the ODE system

if  ~isempty(OPTIONS.Jacobian)
    if nargin(OPTIONS.Jacobian) ~= 2
    warning('The Jacobian function should have exactly two arguments.')
    end
    
    Jac=OPTIONS.Jacobian(t0,Y0);
    if size(Jac,2)~=OPTIONS.NVAR
     error('Dimension of the Jacobian needs to be NVAR*NVAR')   
    end
    if size(Jac,1)~=OPTIONS.NVAR
     error('Dimension of the Jacobian needs to be NVAR*NVAR')   
    end
end

%% TLM integrator
if  ~isempty(OPTIONS.Y_TLM)
    Y_TLM=OPTIONS.Y_TLM;
    if size(Y_TLM,1)~=OPTIONS.NVAR
     error('Dimension of rows for Y_TLM should be NVAR')   
    else
    OPTIONS.NTLM = size(Y_TLM, 2);
    end
end

%% Adjoint Integrator 
if  ~isempty(OPTIONS.Jacp)
    if nargin(OPTIONS.Jacp) ~= 2
    warning('Jacp should have exactly two arguments.')
    end
    Jacp=OPTIONS.Jacp(t0,Y0);
    if size(Jacp,1)~=OPTIONS.NVAR
     error('Dimension of the Jacobian needs to be NVAR*NVAR')   
    else
    OPTIONS.NP = size(Jacp, 2);
    end
end

if  ~isempty(OPTIONS.DRDP)
    if nargin(OPTIONS.DRDP) ~= 2
    warning('DRDP should have exactly two arguments.')
    end
    DRDP=OPTIONS.DRDP(t0,Y0);
    if size(DRDP,2)~=OPTIONS.NP
     error('Dimension of DRDP needs to be NADJ*NP')   
    end
    OPTIONS.NADJ = size(DRDP, 1);
end

if  ~isempty(OPTIONS.DRDY)
    if nargin(OPTIONS.DRDY) ~= 2
    warning('DRDY should have exactly two arguments.')
    end
    DRDY=OPTIONS.DRDY(t0,Y0);
    if size(DRDY,2)~=OPTIONS.NVAR
     error('Dimension of DRDY needs to be NADJ*NVAR')   
    end
    if size(DRDY,1)~=OPTIONS.NADJ
     error('Dimension of DRDY needs to be NADJ*NVAR')   
    end
end 
 
end


