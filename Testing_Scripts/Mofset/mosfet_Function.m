function Udot = mosfet_Function(t, U, Uop, U0, ~, fTilde)

Udot = [Uop - U(1) - fTilde(mosfet_Uin(t), U(1), U0);
    Uop - U(2:end) - fTilde(U(1:end-1), U(2:end), U0)];
return;