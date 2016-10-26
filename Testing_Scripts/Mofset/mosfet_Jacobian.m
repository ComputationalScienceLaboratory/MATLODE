function jac = mosfet_Jacobian(t, U, ~, U0, UT, ~)

mainDiag = [2 * min(U(1) + UT - mosfet_Uin(t), 0) - 1;
    2 * min(U(2:end) + UT - U(1:end-1), 0) - 1];

Utemp = UT - U(1:end-1); % UT - U_(k-1)
lowerDiag = [abs(U(2:end) + Utemp) - abs(U0 + Utemp) - U(2:end) + U0; 0];

jac = spdiags([lowerDiag mainDiag], [-1 0], length(U), length(U));
end