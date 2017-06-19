function [ OPTIONS ] = Input_Dimension(Y0,OdeFunction,OPTIONS)

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

%% For any integrator, we check that Y0 is a NVAR*1 vector

if  isempty(Y0)
    error('The Initial vector Y0 for the ODE system should be provided')
end

if size(Y0,2)~=1 &&  size(Y0,2)~=1
    error('The Initial Y0 for the ODE system should be a vector and cannot be a matrix')
end

if strcmp(OPTIONS.Implementation,'FWD') 
    if size(Y0,2)~=1
        Y0=transpose(Y0);
    end
else    
    if size(Y0,2)~=1
        error('Initial vector for the ODE system should be a column vector')
    end
end
OPTIONS.NVAR = size(Y0, 1);

%% For any integrator, we check that the ODE function returns a NVAR*1 vector
if  isempty(OdeFunction)
    error('The ODE system should be provided')
end

% Check the number of arguments the function accepts
if nargin(OdeFunction) ~= 2
    warning('ODE function should have exactly two arguments.')
end

t0=0;    
Yp0=OdeFunction(t0,Y0);
if size(Yp0,2)~=1 && size(Yp0,1)~=OPTIONS.NVAR
  error('The ODE system should return a column vector NVAR*1')
end
 

%% For any integrator that requires the Jacobian of the ODE system, we check
% it returns a NVAR*NVAR matrix

if  ~isempty(OPTIONS.Jacobian)
    if nargin(OPTIONS.Jacobian) ~= 2
    warning('The Jacobian function should have exactly two arguments.')
    end
    
    Jac=OPTIONS.Jacobian(t0,Y0);
    if (size(Jac,2)~=OPTIONS.NVAR  || size(Jac,1)~=OPTIONS.NVAR)
     error('Dimension of the Jacobian should be be NVAR*NVAR')   
    end
end

%% For any TLM integrator, we check that Y_TLM is provided and should be of 
% dimension NVAR*NTLM

if strcmp(OPTIONS.Implementation,'TLM')
    if  ~isempty(OPTIONS.Y_TLM)
        Y_TLM=OPTIONS.Y_TLM;
        if size(Y_TLM,1)~=OPTIONS.NVAR
         error('Dimension of rows for Y_TLM should be NVAR')   
        else
        OPTIONS.NTLM = size(Y_TLM, 2);
        end
    else
        error('Y_TLM should be provided') 
    end
else
    OPTIONS.Y_TLM=[];
end
    

%%  For any Adjoint  integrator, we check and set
% OPTIONS.NVAR (number of state variables),
% OPTIONS.NP (number of parameters used in the adjoint algorithm)
% OPTIONS.NADJ (number of cost function in the adjoint algorithm)
% DRDP should return NADJ*NP
% DRDY should return NADJ*NVAR
% Jacp should return NVAR*NP

if strcmp(OPTIONS.Implementation,'ADJ')
    
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
    else
        error('Jacp should be provided') 
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
    else
        error('DRDP should be provided') 
    end
    

    if  ~isempty(OPTIONS.DRDY)
        if nargin(OPTIONS.DRDY) ~= 2
        warning('DRDY should have exactly two arguments.')
        end
        DRDY=OPTIONS.DRDY(t0,Y0);
        if (size(DRDY,2)~=OPTIONS.NVAR || size(DRDY,1)~=OPTIONS.NADJ)
        error('Dimension of DRDY needs to be NADJ*NVAR')    
        end
      else
        error('DRDY should be provided') 
    end 
else
    OPTIONS.NP =[];
    OPTIONS.NADJ =[];
    OPTIONS.DRDY=[];
    OPTIONS.DRDP=[];
    OPTIONS.Jacp=[];
end
    
end


