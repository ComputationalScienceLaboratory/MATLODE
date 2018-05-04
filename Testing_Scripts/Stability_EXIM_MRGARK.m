%% Stability Testing for Explicit-Explicit Multirate

function Stability_EXIM_MRGARK(test)

%Using grayscott as test problem

%Small problem size
% model = csl.odetestproblems.grayscott.presets.Random(10);

%Default Size
%model = csl.odetestproblems.grayscott.presets.Random;

model = csl.odetestproblems.brusselator.presets.Canonical;
% %Large Size
% model = csl.odetestproblems.grayscott.presets.Random(120);
% 
% 
% %Use Brusselator
% model = csl.odetestproblems.brusselator.presets.Canonical;

%model.TimeSpan = [0 0.25];
%model.TimeSpan = [0 100];

[~, yRef] = ode15s(model.F, model.TimeSpan, model.Y0, odeset('AbsTol', 100*eps, 'RelTol', 100*eps, 'Jacobian', model.Jacobian));

disp('done ref');

MethodOrder = 2;

NLF = model.FNonLinear;
LF = model.FLinear;
TIMESPAN = model.TimeSpan;
Y0 = model.Y0;
YREFEND = yRef(end,:);
JACOBIAN = model.JacobianNonLinear;



if nargin == 0
    test = 'all';
end

switch test
    case 'err_diff'
        err_diff(LF, NLF, TIMESPAN, Y0, YREFEND, MethodOrder, JACOBIAN);
%         fprintf('Last time val of reference sol was %d\n', t(end));
    case 'stability'
        test_stability(LF, NLF, TIMESPAN, Y0, YREFEND, MethodOrder, JACOBIAN);
    case 'convergence'
        test_convergence(LF, NLF, TIMESPAN, Y0, YREFEND, MethodOrder, JACOBIAN);
    otherwise
        err_diff(LF, NLF, TIMESPAN, Y0, YREFEND, MethodOrder, JACOBIAN);
        test_convergence(LF, NLF, TIMESPAN, Y0, YREFEND, MethodOrder, JACOBIAN);
        %test_stability(LF, NLF, TIMESPAN, Y0, YREFEND, MethodOrder, JACOBIAN);
end

end

%% Test Explicit Explicit MRGARK Error Difference
function err_diff(LF, NLF, TIMESPAN, Y0, YREFEND, MethodOrder, JACOBIAN)

    fprintf('Running MRGARK with method order %d and default settings\n', MethodOrder);
    
    [~, yOrder] = MATLODE_EXIM_MRGARK_FWD_Integrator(LF, NLF, TIMESPAN, Y0, MATLODE_OPTIONS('Method', MethodOrder, 'Jacobian', JACOBIAN, 'LU', 1));
%     fprintf('Time End for MRGARK was %d\n', time(end));
    Err_Norm = norm(YREFEND - yOrder(end, :))/ norm(YREFEND);
    fprintf('Norm of y(end:, ) of ode15s and EE_MRGARK order %d with default tolerences = %d\n', MethodOrder,  Err_Norm);
end

%% Check for Convergence
function test_convergence(LF, NLF, TIMESPAN, Y0, YREFEND, MethodOrder, JACOBIAN)
    disp('Running Convergence Test');
    tdiff = diff(TIMESPAN);
    HS = tdiff./(2.^(10:15));
    Err_Norm = zeros(1, numel(HS));
    disp(numel(HS));
    for Hi = 1:numel(HS)
        H = HS(Hi);
        
        [~, yStep] = MATLODE_EXIM_MRGARK_FWD_Integrator(LF, NLF, TIMESPAN, Y0, MATLODE_OPTIONS('Method', MethodOrder, ...
            'InitialM', 2, 'RadiusM', 0, 'Hmax', H, 'Hmin', H, 'HStart', H, 'Jacobian', JACOBIAN, 'NewtonMaxIt', inf, 'NewtonTol', sqrt(eps)));
        Err_Norm(Hi) = norm(YREFEND - yStep(end, :))/ norm(YREFEND);
        fprintf('Error with fixed step size of %d is %d\n', H, Err_Norm(Hi)) 
    end
    
    a = polyfit(log(HS), log(Err_Norm), 1);
    fprintf('Convergence with method order %d is %d\n', MethodOrder, a(1));
    clf;
    hold all;
    plot(HS, Err_Norm);
    plot(HS, exp(polyval(a, log(HS))));
    
    ax = gca;
    ax.XScale = 'log';
    ax.YScale = 'log';
end

%% Test for Stability
function test_stability(LF, NLF, TIMESPAN, Y0, YREFEND, MethodOrder, JACOBIAN)
    tdiff = diff(TIMESPAN);
    Hinitial = tdiff / (2 ^ 11);
    fprintf('Running Stability Test with initial step-size of %d\n', Hinitial);
    fac = 1.5;
    while fac > 1.01
        
        [~, yStep] = MATLODE_EXIM_MRGARK_FWD_Integrator(LF, NLF, TIMESPAN, Y0, MATLODE_OPTIONS('Method', MethodOrder, ...
            'InitialM', 2, 'RadiusM', 0, 'Hmax', Hinitial, 'Hmin', Hinitial, 'Hstart', Hinitial, 'Jacobian', JACOBIAN));
        Err_Norm = norm(YREFEND - yStep(end, :))/ norm(YREFEND);
        %     disp(fac);
        if Err_Norm >= 0.99 || isnan(Err_Norm)
            Hinitial = Hinitial / fac;
            fac = 1 + (fac - 1) * 0.75;
        else
            Hinitial = Hinitial * fac;
        end   
    end
    fprintf('Method of order %d broke at step-size of %d\n', MethodOrder, Hinitial);
end
%% END OF TESTS