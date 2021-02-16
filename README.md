# Bifurcation Diagram
This is a MATLAB function designed to numerically compute and draw  bifurcation diagrams. Functionality is currently limited to 1st order 1D functions, but will expand as time goes on.

Last updated 2/16/2021.

## Bifurcations

 Consider a function of the form below, where `x` is the state variable , `r` is the parameter.
$$
\frac{dx}{dt} = f(x;r)
$$
Then the bifurcation diagram is the solution to:
$$
f(x_e,r) = 0
$$
The stability of the fixed points thus found is given by
$$
\frac{df}{dx}(x,r)\large |_{x_e} = \begin{cases} < 0 & \text{stable} \\ = 0 & \text{half-stable} \\ > 0 & \text{unstable} \end{cases}
$$

## Use Guide

The main function to use is `bifurcation1D()`. In order to use it, you must first define the dynamical system as a MATLAB function of the form:

```matlab
[dx, fprime] = dynSys(x,r)
```

Where the output `dx` is the value of `f(x,r)` and `fprime` is `df/dx`. In the current implementation, the function must be **vectorized** for values of `x`, though not for `r`. 

Once the dynamic system is defined, you may call one of the following commands to generate the bifurcation diagram. `xLim` and `rLim` represent limits on `x` and `r` respectively, such that `-xLim < x < xLim` and `-rLim < r < rLim`. Default values are used if these are not specified.

```matlab
bifurcation1D(@dynSys)

bifurcation1D(@dynSys,xLim)

bifurcation1D(@dynSys,xLim,rLim)
```

------

Developed in conjunction with taking ENM512 - Nonlinear Dynamics and Chaos at the University of Pennsylvania in Spring 2021