# SincFun.jl

[![Build Status](https://travis-ci.org/MikaelSlevinsky/Sincfun.jl.svg?branch=master)](https://travis-ci.org/MikaelSlevinsky/Sincfun.jl)

This package is a work in progress with the goal of implementing fast algorithms relating to cardinal Sinc expansions of functions and solutions to equations such as ODEs and integral equations, among others.

The basic usage of the package is like so:

```julia
using SincFun
```

We use the sincfun constructor on a function, potentially with the domain and precision specified as well. Without specification, the domain is [-1,1] with double precision.

```julia
f(x) = exp(x);
sf = sincfun(f);
x = linspace(-1.0,1.0,101);
println(norm(f(x)-sf[x])," ",sum(sf)," ",length(sf))
```

The constructor also allows other domains and precisions and can also incorporate singularities.

```julia
f(x) = exp(x);g(x) = f(x)./(1.-x.^2).^(4//5).*log(1.+x);

sf = sincfun(f,Finite(-1.0,1.0,-0.8,-0.8,0.0,1.0));
x = linspace(-0.999,0.999,101);
println(norm(g(x)-sf[x])," ",sum(sf)," ",length(sf))

sf = sincfun(f,Finite(big(-1.0),big(1.0),BigFloat("-0.8"),BigFloat("-0.8"),big(0.0),big(1.0)));
x = linspace(BigFloat("-0.999"),BigFloat("0.999"),101);
println(norm(g(x)-sf[x])," ",sum(sf)," ",length(sf))
```
```julia
f(x) = 1./(x.^2.+1);
sf = sincfun(f,Infinite1(BigFloat));
x = linspace(big(-10.0),big(10.0),101);
println(norm(f(x)-sf[x])," ",sum(sf)," ",length(sf))
```
