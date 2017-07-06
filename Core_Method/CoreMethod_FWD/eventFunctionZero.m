function [theta, component] = eventFunctionZero(T, Y, Z, H, EventFunction, components)
    
    temptheta = ones(size(components));
    for i = 1:length(components)
        myFun = @(theta)eventComponent(theta, T, Y, Z, H, components(i), EventFunction);
        temptheta(i) = fzero(myFun, [0 1]);
    end
    [theta, j] = min(temptheta);
    component = components(j);

return

function y = interpolatedY(theta, Y, Z, H)

%     b = zeros(size(rkDOb,1), 1);
% 
%     for i = 1:size(rkDOb,1)
%        b(i) = polyval(rkDOb(i,:), theta);
%     end
    y = Y + theta*Z(:,5);

return

function f = eventComponent(theta, T, Y, Z, H, component, EventFunction)

    fun = EventFunction(T+theta*H, interpolatedY(theta, Y, Z, H));
    f = fun(component);

return






