function [ jacV ] = arenstorfOrbit_JacobianVector( ~, y, v )

    % Model Parameters
    mu = 0.012277471;
    mu_hat = 1 - mu;

    % Analytical Jacobian
    jac(1,1) = 0;
    jac(1,2) = 0;
    jac(1,3) = 1;
    jac(1,4) = 0;

    jac(2,1) = 0;
    jac(2,2) = 0;
    jac(2,3) = 0;
    jac(2,4) = 1;

    jac(3,1) = (3*mu*(2*mu_hat - 2*y(1))*(mu_hat - y(1)))/(2*((mu_hat - y(1))^2 + y(2)^2)^(5/2)) - mu/((mu_hat - y(1))^2 + y(2)^2)^(3/2) - mu_hat/((mu + y(1))^2 + y(2)^2)^(3/2) + (3*mu_hat*(2*mu + 2*y(1))*(mu + y(1)))/(2*((mu + y(1))^2 + y(2)^2)^(5/2)) + 1;
    jac(3,2) = (3*mu_hat*y(2)*(mu + y(1)))/((mu + y(1))^2 + y(2)^2)^(5/2) - (3*mu*y(2)*(mu_hat - y(1)))/((mu_hat - y(1))^2 + y(2)^2)^(5/2);
    jac(3,3) = 0;
    jac(3,4) = 2;

    jac(4,1) = (3*mu_hat*y(2)*(2*mu + 2*y(1)))/(2*((mu + y(1))^2 + y(2)^2)^(5/2)) - (3*mu*y(2)*(2*mu_hat - 2*y(1)))/(2*((mu_hat - y(1))^2 + y(2)^2)^(5/2));
    jac(4,2) = (3*mu*y(2)^2)/((mu_hat - y(1))^2 + y(2)^2)^(5/2) - mu/((mu_hat - y(1))^2 + y(2)^2)^(3/2) - mu_hat/((mu + y(1))^2 + y(2)^2)^(3/2) + (3*mu_hat*y(2)^2)/((mu + y(1))^2 + y(2)^2)^(5/2) + 1;
    jac(4,3) = -2;
    jac(4,4) = 0;

    % Jacobian vector product
    jacV = jac*v;


end

