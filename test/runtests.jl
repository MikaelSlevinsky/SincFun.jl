using SincFun
using Base.Test

x = linspace(-1,1,101)

sf = sincfun(exp)

@test norm(sf[x]-exp(x)) < 100eps(eltype(sf))

println("Time for 101 samples should be <0.002 seconds.")
@time sf[x]

@test_approx_eq sum(sf) 2sinh(1)

for p in [2,4,6]
    @test_approx_eq norm(sf,p) (2/p*sinh(p))^(1/p)
end

sfs = sincfun(x->sinpi(5x))
sfc = sincfun(x->cospi(5x))

@test norm(sfs[x]-5π*cumsum(sfc)[x]+sinpi(5)) < 200eps(eltype(sfs))

@test norm(roots(sfc) - linspace(-0.9,0.9,10)) < 200eps(eltype(sfc))

println("Time for 10 roots should be <0.80 seconds.")
@time roots(sfc)


h = ConformalMap(1.0,0.1.^(1:4))

@test norm(h[inv(h,x)]-x) < 10eps()

println("Time for 101 samples and inverses of a ConformalMap should be < 0.0002 seconds.")
@time h[inv(h,x)]
