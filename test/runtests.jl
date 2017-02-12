using SincFun
using Base.Test

x = linspace(-1,1,101)

sf = sincfun(exp)

@test norm(sf[x]-exp.(x)) < 100eps(eltype(sf))

println("Time for 101 samples should be <0.002 seconds.")
@time sf[x]

@test sum(sf) ≈ 2sinh(1)

for p in [2,4,6]
    norm(sf,p) ≈ (2/p*sinh(p))^(1/p)
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

println("Testing sinc differentiation and indefinite integration matrices.")

k=collect(-5.0:5.0)
@test sinc.(-1,k) == [-0.020107164191308535,0.025030330116344923,-0.03309323761827199,0.04858833320985967,-0.08948987223608373,0.5,1.0894898722360837,0.9514116667901403,1.033093237618272,0.9749696698836551,1.0201071641913084]
@test sinc.(0,k) == [0.0,0.0,0.0,0.0,0.0,1.0,-0.0,-0.0,-0.0,-0.0,-0.0]
@test sinc.(1,k) == [0.2,-0.25,0.3333333333333333,-0.5,1.0,0.0,-1.0,0.5,-0.3333333333333333,0.25,-0.2]
@test sinc.(2,k) == [0.08,-0.125,0.2222222222222222,-0.5,2.0,-3.289868133696453,2.0,-0.5,0.2222222222222222,-0.125,0.08]

k=collect(-500.0:500.0)
@test all(Bool[norm(sinc.(i,k)-convert(Vector{Float64},sinc.(i,convert(Vector{BigFloat},k))),Inf) < gamma(i+1)eps() for i=2:20])
