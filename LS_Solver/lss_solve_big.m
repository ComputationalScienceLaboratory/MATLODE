%%
%  Authored by Tony D'Augustine, Adrian Sandu, and Hong Zhang.
%  Computational Science Laboratory, Virginia Tech.
%  ©2015 Virginia Tech Intellectual Properties, Inc.
%
function [ Xout ] = lss_solve_big( Transp, Zbig, ax_big, ISING )
%LSS_SOLVE_BIG Summary of this function goes here
%   Detailed explanation goes here

    Xout = ax_big\Zbig;
    Xout = Xout';


end

