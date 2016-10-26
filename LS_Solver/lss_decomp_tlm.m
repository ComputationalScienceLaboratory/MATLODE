%%
%  Authored by Tony D'Augustine, Adrian Sandu, and Hong Zhang.
%  Computational Science Laboratory, Virginia Tech.
%  ©2015 Virginia Tech Intellectual Properties, Inc.
%
function [ ISING e_tlm ] = lss_decomp_tlm( NVAR, hgamma, ISING, fjac1 )
% Based off of lapack_decomp_tlm( hgamma, ising )

    e_tlm = zeros(NVAR,NVAR);

    for j=1:NVAR
        for i=1:NVAR
            e_tlm(i,j) = -fjac1(i,j);
        end
        e_tlm(j,j) = e_tlm(j,j) + hgamma;
    end
    
    % Use cond(e_tlm) to determine if matrix ( hgamma-fjac1 ) is singular
    condNum = cond(e_tlm);
    if ( condNum >= 10e5 )
        ISING = true; % ising -> 1
    else
        ISING = false; % ising -> 0
    end

return;

