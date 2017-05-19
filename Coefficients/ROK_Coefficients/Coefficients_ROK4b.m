function [ rokMethod, rokELO, rokS, rokName ] = Coefficients_ROK4b( ROK4B )

global alpha gamma c b bhat ROK_E gam

    % Rosenbrock method number
    rokMethod = ROK4B;

    % Number of stages
    rokS = 6;
    
    % rosELO  = estimator of local order - the minimum between the
    %           main and the embedded scheme orders plus 1
    rokELO = 4.0d0;

    % Method name
    rokName = 'ROK4b';

    alpha = zeros(6,6);
    gamma = zeros(6,6);

    alpha(2,1) = 1;
    alpha(3,1) = 15919/30000;
    alpha(3,2) = -(919/30000);
    alpha(4,1) = 161/180;
    alpha(4,2) = 1/18;
    alpha(4,3) = 1/20;
    alpha(5,1) = 443/600;
    alpha(5,2) = -73/600;
    alpha(5,3) = 1/3;
    alpha(5,4) = 1/20;
    alpha(6,1) = -(6772995349881/69875766435800);
    alpha(6,2) = -73/600;
    alpha(6,3) = 5479567934713/5240682482685;
    alpha(6,4) = 1208940757253/6987576643580;
    alpha(6,5) = 0;

    gamma(2,1) = -(4195163/183800);
    gamma(3,1) = -(79658501/1148750);
    gamma(3,2) = -(919/30000);
    gamma(4,1) = 148771649/367600;
    gamma(4,2) = 1/18;
    gamma(4,3) = 1/20;
    gamma(5,1) = -343/600;
    gamma(5,2) = -73/600;
    gamma(5,3) = 1/3;
    gamma(5,4) = 1/20;
    gamma(6,1) = 55256869267543/209627299307400;
    gamma(6,2) = -73/600;
    gamma(6,3) = -(661926537641/1746894160895);
    gamma(6,4) = -(102036618579/1397515328716);
    gamma(6,5) = 0;

    c(1) = 0;
    c(2) = sum(alpha(2,:));
    c(3) = sum(alpha(3,:));
    c(4) = sum(alpha(4,:));
    c(5) = sum(alpha(5,:));
    c(6) = sum(alpha(6,:));

    gam = 31/100;
    b(1) = 1/6;
    b(2) = -73/300;
    b(3) = 2/3;
    b(4) = 1/10;
    b(5) = 0;
    b(6) = 31/100;
    bhat(1) = 1/6;
    bhat(2) = -73/300;
    bhat(3) = 2/3;
    bhat(4) = 1/10;
    bhat(5) = 31/100;
    bhat(6) = 0;

    ROK_E = bhat - b;

return;

