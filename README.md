
# MATLODE V2.0.0beta1

##Example Code
*Using ODE Test Problems

```matlab
problem = otp.hires.presets.Canonical;
integrator = matlode.rk.erk.DormanPrince;

options.ErrNorm = matlode.errnorm.HistNorm.errEstimate(1e-6, 1e-6);
options.StepSizeController = matlode.stepsizecontroller.Gustafson;

sol = integrator.integrate(problem.Rhs.F, problem.TimeSpan, problem.Y0, options);
```

[insert other useful stuff]