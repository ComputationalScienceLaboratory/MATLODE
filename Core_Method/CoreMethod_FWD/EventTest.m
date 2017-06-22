function [withinTOL, signChange, theta] = ...
    EventTest(directions, eventOldVal, eventNewVal, eventTOL, eventFunction, ...
                T, Y, Z, H)

N = length(eventOldVal);
tolComponents = zeros(size(eventOldVal));
componentsOfInterest = [];
signChange = false;
eventOldSign = sign(eventOldVal);
eventNewSign = sign(eventNewVal);

eventSignChanges = eventOldSign.*eventNewSign;

for i = 1:N
    
    if (abs(eventNewVal(i)) < eventTOL && eventOldSign(i) == -directions(i))
        tolComponents(i) = true;
        
    end
    
    if (eventSignChanges(i) < 0 && eventOldSign(i) == -directions(i))
        componentsOfInterest(end+1) = i;
        signChange = true;
    end
    
end

[theta, component] = eventFunctionZero(T, Y, Z, H, ...
                    eventFunction, componentsOfInterest); 
                

withinTOL = tolComponents(component);
if( withinTOL )
    keyboard
end


return