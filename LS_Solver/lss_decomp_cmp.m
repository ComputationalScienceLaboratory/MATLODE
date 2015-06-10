%%
%  Authored by Tony D'Augustine, Adrian Sandu, and Hong Zhang.
%  Computational Science Laboratory, Virginia Tech.
%  ©2015 Virginia Tech Intellectual Properties, Inc.
%
function [ ISING, e2 ] = lss_decomp_cmp( NVAR, Alpha, Beta, fjac )
% Based off of lapack_decomp_cmp( alpha, beta, ising ) [LS_Solver.F90]
    
    cmplx = Alpha + Beta*1i;
    e2 = -fjac + eye(NVAR,NVAR)*cmplx;
    
%     condNum = cond(e2);
%     if ( condNum >= 10e4 )
%         ISING = true; % ISING -> 1
%     else
%         ISING = false; % ISING -> 0
%     end

    condNum = rcond(e2);
    if ( condNum <= eps )
        ISING = true; % ising -> 1
    else
        ISING = false; % ising -> 0
    end
    

return;

