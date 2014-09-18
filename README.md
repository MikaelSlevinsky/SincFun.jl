# Sincfun

[![Build Status](https://travis-ci.org/MikaelSlevinsky/Sincfun.jl.svg?branch=master)](https://travis-ci.org/MikaelSlevinsky/Sincfun.jl)

This package is a work in progress with the goal of implementing fast algorithms relating to cardinal Sinc expansions of functions and solutions to equations such as ODEs and integral equations, among others.

The basic usage of the package is like so:


	using Sincfun


We use the sincfun constructor on a function, potentially with the domain and precision specified as well. Without specification, the domain is [-1,1] with double precision.


	f(x) = exp(x)./(x.^2.+1)
	sf = sincfun(f)
	x = linspace(-1.0,1.0,101)
	norm(f(x)-sf[x])


