
# MATLODE V2.0.0beta1

##Example Code
*Using ODE Test Problems

```matlab
problem = otp.hires.presets.Canonical;
integrator = matlode.rk.erk.ClassicRK4;
options = integrator.matlodeSets('StepsN', 100000);
sol = integrator.integrate(problem.Rhs.F, problem.TimeSpan, problem.Y0);
```
OR
```matlab
problem = otp.hires.presets.Canonical;
integrator = matlode.rk.erk.ClassicRK4;
sol = integrator.integrate(problem.Rhs.F, problem.TimeSpan, problem.Y0, 'StepsN', 100000);
```
OR
```matlab
problem = otp.hires.presets.Canonical;
integrator = matlode.rk.erk.ClassicRK4;
options = integrator.matlodeSets('StepsN', 100000);
sol = integrator.integrate(problem.Rhs.F, problem.TimeSpan, problem.Y0, options, 'AbsTol', 1e-8);
```

[insert other useful stuff]