function [w, stats] = phipmPaul(t, A, u, tol, symm, m)
% PHIPM - Evaluates a linear combinaton of the phi functions
%         evaluated at tA acting on vectors from u, that is
%
%         w = phi_0(tA) u(:, 1) + t phi_1(tA) u(:, 2) +
%             t^2 phi_2(tA) u(:, 3) + ...
%
%         The evaluation expresses eveything in terms of the highest
%         order phi function and evaluates the action of this on a
%         vector using a Krylov technique and then computes w using
%         the recurrence relation.
%
%         The size of the Krylov subspace is changed dynamically
%         during the integration. The Krylov subspace is computed
%         using Arnoldi if A is non-symmetric and Lancozs if A is
%         symmetric.
%
% PARAMETERS:
%   t    - constant value represent time.
%   A    - the matrix argument of the phi functions.
%   u    - the matrix with columns representing the vectors to be
%          multiplied by the phi functions.
%   tol  - the convergence tolarance required.
%   symm - true if the matrix A is symmetric.
%   m    - an estimate of the appropriate Krylov size.
%
% RETURNS:
%   w        - the linear combination of the phi functions
%              evaluated at tA acting on the vectors from u.
%   stats(1) - number of substeps
%   stats(2) - number of rejected steps
%   stats(3) - number of Krylov steps
%   stats(4) - number of matrix exponentials

persistent V int

%V = [];
%int = 0;

% n is the size of the original problem
% p is the number of phi functions
[n, p] = size(u);
Aisahandle = isa(A, 'function_handle');
if ~Aisahandle
    nnze = nnz(A);
else
    nnze = 10 * n; % wild guess
end;

% Add extra column of zeros if p=1
if p == 1
    p = 2;
    u = [u, zeros(size(u))];
end

% Check inputs
if nargin < 6
    m = 1;
    if nargin < 5
        if max(A - A') < 100 * eps
            symm = true;
        else
            symm = false;
        end
        if nargin < 4
            tol = 1.0e-7;
            if nargin < 3
                error('phipm:NotEnoughInputs',...
                      'Not enough input arguments.');
            end
        end
    end
end

% Krylov parameters
mmax = 500;
m_new = m;

% Preallocate matrices
% V = zeros(n, mmax + 1);
% int = zeros(n, p);
if isempty(V) && isempty(int)
    V = zeros(n, mmax + 1);
    int = zeros(n, p);
end

if numel(int) ~= n * p
    int = zeros(n, p);
end

if numel(V) ~= n * (mmax + 1)
    V = zeros(n, mmax + 1);
    int = zeros(n, p);
end

% Initializing the variables
step = 0;
krystep = 0;
ireject = 0;
reject = 0;
exps = 0;
happy = 0;
%sgn = sign(t);
sgn = sign(t(1));
t_now = 0;
t_out = abs(t(end));
j = 0;
sizet = size(t);
numSteps = sizet(2);

% Compute and initial starting approximation for the timestep
tau = t_out;

% Setting the safety factors and tolarance requirements
gamma = 0.8;
delta = 1.2;

% Used for toeplitz trick
cidx = (0:p-1)';
ridx = p:-1:1;
idx = cidx(:,ones(p,1)) + ridx(ones(p,1),:);

% Initial condition
w = zeros(length(u(:,1)), numSteps);
w(:,1) = u(:, 1);
oldm = NaN; oldtau = NaN; omega = NaN;
orderold = true; kestold = true;

l = 1;
while t_now < t_out

    % Compute necessary starting information
    if j == 0

        % Initialize the matrices V and H
        H = zeros(mmax + p, mmax + p);
        x = [zeros(1, p - 1), cumprod([1, t_now ./ (1:p - 1)])];
        up = u * x(idx);  % Code inspired by toeplitz.m
        up2 = u * x(idx);

        % Compute the update factors
        int(:, 1) = w(:,l);
        for i = 1:p - 1

            % Determine if matrix free
            if ~Aisahandle
                int(:, i + 1) = A * int(:, i) + up(:, i + 1);
            else
                int(:, i + 1) = A(int(:, i)) + up(:, i + 1);
            end;

        end

        % Normalize initial vector
        beta = norm(int(:, end));
        if beta == 0

            % Multiplying with a zero vector, hence result is zero
            % Finish all in one step
            reject = reject + ireject;
            step = step + 1;
            tau = t_out - t_now;
            w(:,l) = w(:,l) + int(:, 2:p - 1) * cumprod(tau * 1 ./ (1: p - 2)');
            break;

        end;

        % The first Krylov basis vector
        V(:, 1) = int(:, end) ./ beta;

    end;

    % Check if matrix is symmetric
    if symm

        % Symmetric use the Lanczos process
        while j < m

            % Determine if matrix free
            j = j + 1;
            if ~Aisahandle
                vv = A * V(:, j);
            else
                vv = A(V(:, j));
            end;
            H(j, j) = V(:, j)' * vv;
            if j == 1,
                vv = vv - H(j, j) * V(:, j);
            else
                vv = vv - H(j - 1, j) * V(:, j - 1) - H(j, j) * V(:, j);
            end
            krystep = krystep + 1;
            s = norm(vv);

            % Happy breakdown
            if s < tol
                happy = 1;
                tau = t_out - t_now;
                break;
            end;
            H(j + 1, j) = s;
            H(j, j + 1) = s;
            V(:, j + 1) = vv ./ s;

        end;
        
        % Keep a record of H
        H2 = H;
        H(m, m + 1) = 0;

        % Matrix is not symmetric
    else

        % Not symmetric use the Arnoldi process
        while j < m

            % Determine if matrix free
            j = j + 1;
            if ~Aisahandle
                vv = A * V(:, j);
            else
                vv = A(V(:, j));
            end;
            for i = 1:j
                H(i, j) = V(:, i)' * vv;
                vv = vv - H(i, j) * V(:, i);
            end
            krystep = krystep + 1;
            s = norm(vv);

            % Happy breakdown
            if s < tol 
                happy = 1;
                tau = t_out - t_now;
                break;
            end;
            H(j + 1, j) = s;
            V(:, j + 1) = vv ./ s;

        end

        % Keep a record of H
        H2 = H;

    end

    % We use the vector e1 in the computations
    H(1, j + 1) = 1;

    % Construct the augmented matrix
    for i = 1:p - 1
        H(j + i, j + i + 1) = 1;
    end
    h = H(j + 1, j);
    H(j + 1, j) = 0;

    % Compute the exponential of the augmented matrix
    [F, hnorm] = expmhigham(sgn * tau * H(1:j + p, 1:j + p));
    exps = exps + 1;

    % Local truncation error estimation
    %err = abs(beta * h * F(j, j + p));
    err = norm(beta * h * F(j, j + p - 1) * V(:, j+1));

    % Error per unit step
    oldomega = omega;
    omega = t_out * err / (tau * tol);

    % Estimate order
    if m == oldm && tau ~= oldtau && ireject >= 1
        order = max(1, log(omega/oldomega) / log(tau/oldtau));
        orderold = false;
    elseif orderold || ireject == 0
        orderold = true;
        order = j/4;
    else
        orderold = true;
    end;
    % Estimate k
    if m ~= oldm && tau == oldtau && ireject >= 1
        kest = max(1.1, (omega/oldomega) ^ (1/(oldm-m)));
        kestold = false;
    elseif kestold || ireject == 0
        kestold = true;
        kest = 2;
    else
        kestold = true;
    end;

    % This if statement is the main difference between fixed and
    % variable m
    oldtau = tau; oldm = m;
    if happy == 1

        % Happy breakdown; wrap up
        omega = 0;
        tau_new = tau;
        m_new = m;

    elseif j == mmax && omega > delta

        % Krylov subspace to small and stepsize to large
        tau_new = tau * (omega / gamma) ^ (-1 / order);

    else

        % Determine optimal tau and m
        tau_opt = tau * (omega / gamma) ^ (-1 / order);
        m_opt = max(1, ceil(j + log(omega / gamma) / log(kest)));
        nom = 5 + max(log(hnorm), 0) / log(2); % number of mult's in expm

        if symm

            % Cost of Lanczos; a factor of 2 has been ignored
            cost1 = ((j + p) * nnze + 3 * (j + p) * n ...
                     + nom * (j + p - 1)^3) * ceil((t_out-t_now) / tau_opt);
            cost2 = ((m_opt + p) * nnze + 3*(m_opt + p) * n ...
                     + nom * (m_opt + p - 1)^3) * ceil((t_out-t_now) / tau);
        else

            % Cost of Arnoldi
            cost1 = ((j + p) * nnze + (j^2 + 3 * p + 2) * n ...
                     + nom * (j + p - 1)^3) * ceil((t_out-t_now) / tau_opt);
            cost2 = ((m_opt + p) * nnze + (m_opt^2 + 3 * p + 2) * n ...
                     + nom * (m_opt + p - 1)^3) * ceil((t_out-t_now) / tau);
        end;

        % Determine whether to vary tau or m
        if cost1 < cost2
            tau_new = tau_opt;
            m_new = m;
        else
            m_new = m_opt;
            tau_new = tau;
        end;

    end;

    % Check error against target
    if omega <= delta

        blownTs = 0;
        nextT = t_now + tau;
        for k = l:numSteps
            if t(k) < nextT
                blownTs = blownTs + 1;
            end
        end

        % Copy current w to w we continue with.
        w(:,l + blownTs) = w(:,l);

        for k = 0:blownTs - 1
            tauPhantom = t(l+k) - t_now;
            F2 = expmhigham(sgn * tauPhantom * H(1:j + p, 1:j + p));
            up2 = w(:,l+blownTs) + int(:, 2:p - 1) * cumprod(tauPhantom * 1 ./ (1: p - 2)');
            F2(j + 1, j + p - 1) = h * F2(j, j + p);
            w(:,l+k) = beta * V(:, 1:j + 1) * F2(1:j + 1, j + p - 1) + up2;
        end

        % Advance l.
        l = l + blownTs;

        % Yep, got the required tolerance; update
        reject = reject + ireject;
        step = step + 1;
        up = w(:,l) + int(:, 2:p - 1) * cumprod(tau * 1 ./ (1: p - 2)');

        % Using the corrected quantity
        F(j + 1, j + p - 1) = h * F(j, j + p);
        w(:,l) = beta * V(:, 1:j + 1) * F(1:j + 1, j + p - 1) + up;

        % Update t
        %fprintf('Step = %d, basisVectors = %d\n', step, j);
        t_now = t_now + tau;
        j = 0;
        ireject = 0;

    else

        % Nope, try again
        H = H2;
        ireject = ireject + 1;

    end;

    % Safety factors for tau
    tau = min(t_out - t_now, max(tau/5, min(5*tau, tau_new)));

    % Safety factors for m
    m = max(1, min(mmax, max(floor(3/4*m), min(m_new, ceil(4/3*m)))));

end

phiHan = find(max(abs(u)));
if isempty(phiHan)
    phiHan = length(u);
else
    phiHan = phiHan-1;
end

for l = 1:numSteps
    w(:,l) = w(:,l)*(1/t(l).^phiHan);
end
stats = [step, reject, krystep, exps];
