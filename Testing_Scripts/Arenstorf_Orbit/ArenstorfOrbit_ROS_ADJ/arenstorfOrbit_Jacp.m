function fpjac= arenstorfOrbit_Jacp( ~,y )

mu = 0.012277471;
mu_hat = 1 - mu;

fpjac(1,1) = 0;
fpjac(1,2) = 0;

fpjac(2,1) = 0;
fpjac(2,2) = 0;

fpjac(3,1) = (mu_hat - y(1))/((mu_hat - y(1))^2 + y(2)^2)^(3/2) - mu_hat/((mu + y(1))^2 + y(2)^2)^(3/2) + (3*mu_hat*(2*mu + 2*y(1))*(mu + y(1)))/(2*((mu + y(1))^2 + y(2)^2)^(5/2));
fpjac(3,2) = mu/((mu_hat - y(1))^2 + y(2)^2)^(3/2) - (mu + y(1))/((mu + y(1))^2 + y(2)^2)^(3/2) - (3*mu*(2*mu_hat - 2*y(1))*(mu_hat - y(1)))/(2*((mu_hat - y(1))^2 + y(2)^2)^(5/2));

fpjac(4,1) = (3*mu_hat*y(2)*(2*mu + 2*y(1)))/(2*((mu + y(1))^2 + y(2)^2)^(5/2)) - y(2)/((mu_hat - y(1))^2 + y(2)^2)^(3/2);
fpjac(4,2) = (3*mu*y(2)*(2*mu_hat - 2*y(1)))/(2*((mu_hat - y(1))^2 + y(2)^2)^(5/2)) - y(2)/((mu + y(1))^2 + y(2)^2)^(3/2);

return;

