%% Stability Testing for Explicit-Explicit Multirate

function Stability_EXEX_MRGARK(test)
%Using grayscott as test problem
%model = csl.odetestproblems.grayscott.presets.Random(10, 1);
model = csl.odetestproblems.brusselator.presets.Canonical;

[~, yRef] = ode15s(model.F, model.TimeSpan, model.Y0, odeset('AbsTol', 100*eps, 'RelTol', 100*eps));

MethodOrder = 2;

LF = model.FLinear;
NLF = model.FNonLinear;
TIMESPAN = model.TimeSpan;
Y0 = model.Y0;
YREFEND = yRef(end,:);

if nargin == 0
    test = 'all';
end

switch test
    case 'err_diff'
        err_diff(LF, NLF, TIMESPAN, Y0, YREFEND, MethodOrder);
    case 'stability'
        test_stability(LF, NLF, TIMESPAN, Y0, YREFEND, MethodOrder);
    case 'convergence'
        test_convergence(LF, NLF, TIMESPAN, Y0, YREFEND, MethodOrder);
    otherwise
        err_diff(LF, NLF, TIMESPAN, Y0, YREFEND, MethodOrder);
        test_convergence(LF, NLF, TIMESPAN, Y0, YREFEND, MethodOrder);
        %test_stability(LF, NLF, TIMESPAN, Y0, YREFEND, MethodOrder);
end

end

%% Test Explicit Explicit MRGARK Error Difference
function err_diff(LF, NLF, TIMESPAN, Y0, YREFEND, MethodOrder)

    fprintf('Running MRGARK with method order %d and default settings\n', MethodOrder);
    [~, yOrder] = MATLODE_EE_MRGARK_FWD_Integrator(LF, NLF, TIMESPAN, Y0, MATLODE_OPTIONS('Method', MethodOrder));
    Err_Norm = norm(YREFEND - yOrder(end, :))/ norm(YREFEND);
    fprintf('Norm of y(end:, ) of ode15s and EE_MRGARK order %d with default tolerences = %d\n', MethodOrder,  Err_Norm);
end

%% Check for Convergence
function test_convergence(LF, NLF, TIMESPAN, Y0, YREFEND, MethodOrder)
    disp('Running Convergence Test');
    tdiff = diff(TIMESPAN);
    HS = tdiff./(2.^(10:15));
    Err_Norm = zeros(1, numel(HS));
    
    for Hi = 1:numel(HS)
        H = HS(Hi);
        
        [~, yStep] = MATLODE_EE_MRGARK_FWD_Integrator(LF, NLF, TIMESPAN, Y0, MATLODE_OPTIONS('Method', MethodOrder, ...
            'InitialM', 1, 'RadiusM', 0, 'Hmax', H, 'Hmin', H));
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
function test_stability(LF, NLF, TIMESPAN, Y0, YREFEND, MethodOrder)

    searchRangeLow = 1e-6;
    searchRangeHigh = diff(TIMESPAN);
    Hinitial = (searchRangeLow + searchRangeHigh) / 2;
    fprintf('Running Stability Test with initial step-size of %d\n', Hinitial);
    tol = 1e-12;
    numItr = ceil(log2((searchRangeHigh - searchRangeLow)/tol));

    for i = 1:numItr
        
        [~, yStep] = MATLODE_EE_MRGARK_FWD_Integrator(LF, NLF, TIMESPAN, Y0, MATLODE_OPTIONS('Method', MethodOrder, ...
            'InitialM', 2, 'RadiusM', 0, 'Hmax', Hinitial, 'Hmin', Hinitial, 'Hstart', Hinitial));
        Err_Norm = norm(YREFEND - yStep(end, :))/ norm(YREFEND);
        
        if  Err_Norm >= 1 || isnan(Err_Norm)
            searchRangeHigh = Hinitial;
            Hinitial = (searchRangeLow + searchRangeHigh) / 2;
        else
            searchRangeLow = Hinitial;
            Hinitial = (searchRangeLow + searchRangeHigh) / 2;
        end
        
    end
    fprintf('Method of order %d broke at step-size of %d\n', MethodOrder, Hinitial);
end
%% END OF TESTS