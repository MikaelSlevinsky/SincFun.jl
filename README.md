# SincFun.jl

[![Build Status](https://travis-ci.org/MikaelSlevinsky/Sincfun.jl.svg?branch=master)](https://travis-ci.org/MikaelSlevinsky/Sincfun.jl)

This package is a work in progress with the goal of implementing fast algorithms relating to cardinal Sinc expansions of functions and solutions to equations such as ODEs and integral equations, among others.

The basic usage of the package is like so:

```julia
using SincFun
```

We use the sincfun constructor on a function, potentially with the domain and precision specified as well. Without specification, the domain is [-1,1] with double precision.

```julia
f(x) = sin(x);
sf = sincfun(f)
sfc = cumsum(sf)
g(x) = cos(-1)-cos(x);
sg = sincfun(g)
norm(sfc-sg)
```

```julia
f(x) = x*exp(x)*cospi(5x)/(1+x^2);
sf = sincfun(f)
roots(sf)
norm(sf[roots(sf)])
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

There is also preliminary support for the Hilbert transform:

```julia
sf = sincfun(x->x/(x^2+1)^2,Infinite1());
sg = sincfun(x->(1-x^2)/2/(x^2+1)^2,Infinite1());
sh=hilbert(sf);
x = linspace(-5.0,5.0,1001);
println(norm(sg[x]-sh[x]))
```
