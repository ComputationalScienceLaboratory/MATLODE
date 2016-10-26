%% rRMS
%
%%
%  Authored by Tony D'Augustine, Adrian Sandu, and Hong Zhang.
%  Computational Science Laboratory, Virginia Tech.
%  ©2015 Virginia Tech Intellectual Properties, Inc.
%
function [ err ] = rRMS( y, y_ref )

    n = size(y,1);
    err = 0;
    threshold = 1e-10; %<----- temporary
    
    for i=1:n
        err = err + ((y(i)-y_ref(i))/max( max(abs(y(i)),abs(y_ref(i))),threshold))^2; 
    end
    err = sqrt(err/n);


return;

