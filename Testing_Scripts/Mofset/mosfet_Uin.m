function Uin = mosfet_Uin(t)

% define a piecewise function for Uin - the input voltage
% Replicates figure 6
if (t < 5)
    Uin = 0;
elseif (t >= 5 && t < 10)
    Uin = t - 5;
elseif (t >= 10 && t < 15)
    Uin = 5;
elseif (t >= 15 && t < 17)
    Uin = -2.5 * t + 42.5;
else
    Uin = 0;
end

end