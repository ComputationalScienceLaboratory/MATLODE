function [ expkMethod, expkELO, expkS, expkName ] = Coefficients_EXPK4( EXPK4 )

global alpha gamma c b bhat e gam

    % Rosenbrock method number
    expkMethod = EXPK4;

    % Number of stages
    expkS = 4;
    
    % rosELO  = estimator of local order - the minimum between the
    %           main and the embedded scheme orders plus 1
    expkELO = 4.0d0;

    % Method name
    expkName = 'EXPK4';

    alpha = zeros(4,4); gamma = zeros(4,4);
    alpha(2,1) = 1;
    alpha(3,1) = 41/80;
    alpha(3,2) = 1/80;
    alpha(4,1) = 1/4;
    alpha(4,2) = 1/12;
    alpha(4,3) = 1/6;
    gamma(2,1) = 7/8;
    gamma(3,1) = 1/16;
    gamma(3,2) = 0;
    gamma(4,1) = -1/32;
    gamma(4,2) = 1/24;
    gamma(4,3) = -5/12;
    b = zeros(4,1);
    b(1) = 1/6;
    b(2) = 1/6;
    b(3) = 0;
    b(4) = 2/3;
    gam = 0.572816062482134855408001384976768340931514124329888624090;
    bhat = zeros(4,1);
    bhat(1) = 8/3;
    bhat(2) = 1;
    bhat(3) = -8/3;
    bhat(4) = 0;
    c(1) = 0;
    c(2) = sum(alpha(2,:));
    c(3) = sum(alpha(3,:));
    c(4) = sum(alpha(4,:));


    e = bhat - b;

return;

