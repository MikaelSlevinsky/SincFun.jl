module SincFun

include("helper.jl")
include("Domains.jl")

export sincfun, hilbert, domain

#
# A sincfun represents a function by
#
# f(x) = ω(ϕ^{-1}(x))/ϕ'(ϕ^{-1}(x))
#       * Σ_{j=-n}^{+n} (-1)^j f(ϕ(jh))*ϕ'(jh)/(ϕ^{-1}(x)-jh)
#       / Σ_{j=-n}^{+n} (-1)^j ω(jh)/(ϕ^{-1}(x)-jh)
#
# where ϕ is a double exponential variable transformation
# and ω is a double exponential weight.
#

type sincfun{D,T}
    n::Integer
    h::T
    fϕv::Vector{T}
    ϕpv::Vector{T}
    ωv::Vector{T}
    ωscale::T
    ωβ::T
    jh::Vector{T}
    domain::D
end

function sincfun{T<:Number}(f::Function,domain::Domain{T})
    n=2^2
    intold,intnew = Inf,-one(T)
    while n < 2^14 && abs(intold-intnew) > abs(intnew)*eps(T)^(1/3)
        n *= 2
        jh = h(T,n)*[-n:n]
        sinhv,coshv = sinh(jh)*π/2,cosh(jh)*π/2
        fϕv = T[f(ψ(domain,z)) for z in sinhv]
        ϕpv = ψp(domain,sinhv).*coshv
        test = abs(fϕv.*ϕpv)
        cutoff = !(!isinf(test).*!isnan(test))
        fϕv[cutoff],ϕpv[cutoff] = zeros(T,sum(cutoff)),zeros(T,sum(cutoff))
        intold,intnew = intnew,h(T,n)*sum(fϕv.^2.*ϕpv)
    end
    jh = h(T,n)*([-n:n-1]+one(T)/2)
    sinhv,coshv = sinh(jh)*π/2,cosh(jh)*π/2
    fϕv = interlace([fϕv,T[f(ψ(domain,z)) for z in sinhv]])
    ϕpv = interlace([ϕpv,ψp(domain,sinhv).*coshv])
    test = abs(fϕv.*ϕpv)
    cutoff = !(!isinf(test).*!isnan(test))
    fϕv[cutoff],ϕpv[cutoff] = zeros(T,sum(cutoff)),zeros(T,sum(cutoff))
    ωscale = maximum(test)
    tru = div(findfirst(reverse(interlace2(h(T,n)/2*fϕv.^2.*ϕpv)) .> ωscale*eps(T)^3),2)
    fϕv,ϕpv,jh = fϕv[tru+1:4n-tru+1],ϕpv[tru+1:4n-tru+1],h(T,n)/2*[-2n+tru:2n-tru]
    ωβ = -2log(eps(T))/(cosh(jh[end])-one(T))
    ωv = ωscale*ω(ωβ,jh)
    return sincfun{typeof(domain),T}(length(fϕv),h(T,n)/2,fϕv,ϕpv,ωv,ωscale,ωβ,jh,domain)
end
sincfun(f::Function) = sincfun(f,Finite())

Base.length(sf::sincfun) = sf.n
Base.eltype{D,T}(sf::sincfun{D,T}) = T
domain(sf::sincfun) = sf.domain


## Evaluation by the barycentric formula

function Base.getindex{D<:Domain,S<:Number,T<:Number}(sf::sincfun{D,S},x::T)
    xc = convert(promote_type(S,T),x)
    z = ψinv(sf.domain,xc)
    barycentric(sf,xc)*singularities(sf.domain,z)
end
Base.getindex{D<:Domain,S<:Number,T<:Number}(sf::sincfun{D,S},x::Vector{T}) = promote_type(S,T)[sf[xk] for xk in x]
Base.getindex{D<:Domain,S<:Number,T<:Number,N}(sf::sincfun{D,S},x::Array{T,N}) = reshape(sf[vec(x)],size(x))


function barycentric{D<:Domain,T<:Number}(sf::sincfun{D,T},x::T)
    t = asinh(2ψinv(sf.domain,x)/π)
    if t < sf.jh[1]
        return envelope(sf,sf.jh[1])*sf.fϕv[1]*sf.ϕpv[1]/sf.ωv[1]
    elseif t > sf.jh[end]
        return envelope(sf,sf.jh[end])*sf.fϕv[end]*sf.ϕpv[end]/sf.ωv[end]
    else
        idx = findfirst(sf.jh,t)
        if idx == 0
            common = t-sf.jh
            num,den = sf.fϕv.*sf.ϕpv./common,sf.ωv./common
            valN,valD = sum(num[1:2:end])-sum(num[2:2:end]),sum(den[1:2:end])-sum(den[2:2:end])
            return envelope(sf,t)*valN/valD
        else
            return envelope(sf,t)*sf.fϕv[idx]*sf.ϕpv[idx]/sf.ωv[idx]
        end
    end
end
envelope{D<:Domain,T<:Number}(sf::sincfun{D,T},t::T) = sf.ωscale*ω(sf.ωβ,t)/ψp(sf.domain,sinh(t)*π/2)/(cosh(t)*π/2)



# Algebra

for op in (:+,:-,:*,:.*)
    @eval begin
        function $op{D<:Domain,T<:Number}(c::Number,sf::sincfun{D,T})
            sf1 = deepcopy(sf)
            sf1.fϕv = $op(convert(T,c),sf.fϕv)
            return sf1
        end
        $op{D<:Domain,T<:Number}(sf::sincfun{D,T},c::Number) = $op(c,sf)

        function $op{D<:Domain,T<:Number}(sf1::sincfun{D,T},sf2::sincfun{D,T})
            #TODO: sf.domain = $op(sf1.domain,sf2.domain)
            sincfun(x->$op(sf1[x],sf2[x]),sf1.domain)
        end
    end
end

# real, imag, conj, abs2

for op in (:(Base.real),:(Base.imag),:(Base.conj),:(Base.abs2))
    @eval begin
        $op{D<:Domain,T<:Real}(sf::sincfun{D,T}) = deepcopy(sf)
    end
end

for op in (:(Base.real),:(Base.imag))
    @eval begin
        function $op{D<:Domain,T<:Real}(sf::sincfun{D,Complex{T}})
            sincfun{typeof(sf.domain),T}(sf.n,real(sf.h),$op(sf.fϕv),real(sf.ϕpv),real(sf.ωv),real(sf.ωscale),real(sf.ωβ),real(sf.jh),sf.domain)
        end
    end
end

function Base.conj{D<:Domain,T<:Real}(sf::sincfun{D,Complex{T}})
    sf1 = deepcopy(sf)
    sf1.fϕv = conj(sf1.fϕv)
    return sf1
end

function Base.abs2{D<:Domain,T<:Real}(sf::sincfun{D,Complex{T}})
    rsf,isf = real(sf),imag(sf)
    return rsf*rsf + isf*isf
end

# sum, dot, cumsum, norm, and diff.

function Base.sum{D<:Domain,T<:Number}(sf::sincfun{D,T})
    sinhv = sinh(sf.jh)*π/2
    singv = singularities(sf.domain,sinhv)
    sf.h*sum(sf.fϕv.*sf.ϕpv.*singv)
end
Base.dot{D<:Domain,T<:Number}(sf1::sincfun{D,T},sf2::sincfun{D,T}) = sum(conj(sf1)*sf2)

function Base.cumsum{D<:Domain,T<:Number}(sf::sincfun{D,T})
    SM = sinc(-1,one(T)*[-sf.n+1:sf.n-1])
    sf1 = deepcopy(sf)
    temp = sf.h*sf.fϕv.*sf.ϕpv
    [sf1.fϕv[i] = sum(SM[i+sf.n-1:-1:i].*temp) for i=1:sf.n]
    return sf1
end

function Base.norm{D<:Domain,T<:Number}(sf::sincfun{D,T},p::Int)
    if isodd(p) error("Only even integer norms supported.") end
    sinhv = sinh(sf.jh)*π/2
    singv = singularities(sf.domain,sinhv)
    abs(sf.h*sum((sf.fϕv.*singv).^p.*sf.ϕpv))^(one(T)/p)
end
Base.norm{D<:Domain,T<:Number}(sf::sincfun{D,T}) = norm(sf,2)

function hilbert{D<:Domain,T<:Number}(sf::sincfun{D,T})
    sinhv = sinh(sf.jh)*π/2
    SM = (sinc(0,one(T)/2/sf.h*(sf.jh.-sf.jh')).*(sf.jh.-sf.jh')).^2./(ψ(sf.domain,sinhv).-ψ(sf.domain,sinhv)')
    [SM[i,i] = zero(T) for i=1:sf.n]
    sf1 = deepcopy(sf)
    singv = singularities(sf.domain,sinhv)
    temp = -sf.fϕv.*sf.ϕpv.*singv*π/2/sf.h
    sf1.fϕv = SM*temp
    return sf1
end

#=
function Base.diff{D<:Domain,T<:Number}(sf::sincfun{D,T})
    SM = sinc(1,one(T)*[-sf.n+1:sf.n-1])
    #SM = sinc(1,one(T)*[-(sf.n-1)/2:(sf.n-1)/2].-[-(sf.n-1)/2:(sf.n-1)/2]')
    sf1 = deepcopy(sf)
    temp = sf.fϕv./(sf.h*sf.ϕpv)
    [sf1.fϕv[i] = sum(SM[i+sf.n-1:-1:i].*temp) for i=1:sf.n]
    sf1.fϕv .*= sf.ϕpv
    #sf1.fϕv = SM*temp
    return sf1
end
=#
#=
function Base.diff{D<:Domain,T<:Number}(sf::sincfun{D,T})
    SM = sinc(1,one(T)*[-sf.n:sf.n])
    sinhv,coshv = convert(T,π)/2*sinh(sf.jh),convert(T,π)/2*cosh(sf.jh)
    temp = ψp(sf.domain,sinhv).*coshv
    phippoverphip = sinhv./coshv - 2tanh(sinhv).*coshv
    sf1 = deepcopy(sf)
    [sf1.fϕpv[i] = sum(SM[i:i+sf.n-1].*(sf.fϕpv/sf.h)) for i=1:sf.n]
    sf1.fϕpv = (sf1.fϕpv-sf.fϕpv.*phippoverphip)#./temp
    #sf1.fϕv = sf1.fϕv./(sf.ϕpv*sf.h)
    #sf1.ϕpv = ones(T,length(sf.ϕpv))
    return sf1
end
=#

include("roots.jl")
include("fftBigFloat.jl")
include("KrylovMethods.jl")
include("steig.jl")
include("sinc.jl")
include("ConformalMap.jl")
include("SincMatrix.jl")
include("BandMatrix.jl")

end #module
