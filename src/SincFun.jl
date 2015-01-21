module SincFun

import Base

Base.eps{T<:Number}(::Type{T}) = eps(real(T))
Base.real{T<:Real}(::Union(Type{T},Type{Complex{T}})) = T

include("Domains.jl")

export sincfun

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
        jh = h(n,T)*[-n:n]
        sinhv,coshv = sinh(jh)*π/2,cosh(jh)*π/2
        fϕv = T[f(domain.ψ(z)) for z in sinhv]
        ϕpv = domain.ψp(sinhv).*coshv
        test = abs(fϕv.*ϕpv)
        cutoff = !(!isinf(test).*!isnan(test))
        fϕv[cutoff],ϕpv[cutoff] = zeros(T,sum(cutoff)),zeros(T,sum(cutoff))
        intold,intnew = intnew,h(n,T)*sum(fϕv.^2.*ϕpv)
    end
    jh = h(n,T)*([-n:n-1]+one(T)/2)
    sinhv,coshv = sinh(jh)*π/2,cosh(jh)*π/2
    fϕv = interlace([fϕv,T[f(domain.ψ(z)) for z in sinhv]])
    ϕpv = interlace([ϕpv,domain.ψp(sinhv).*coshv])
    test = abs(fϕv.*ϕpv)
    cutoff = !(!isinf(test).*!isnan(test))
    fϕv[cutoff],ϕpv[cutoff] = zeros(T,sum(cutoff)),zeros(T,sum(cutoff))
    ωscale = maximum(test)
    tru = div(findfirst(reverse(interlace2(h(n,T)/2*fϕv.^2.*ϕpv)) .> ωscale*eps(T)^3),2)
    fϕv,ϕpv,jh = fϕv[tru+1:4n-tru+1],ϕpv[tru+1:4n-tru+1],h(n,T)/2*[-2n+tru:2n-tru]
    ωβ = -2log(eps(T))/(cosh(jh[end])-one(T))
    ωv = ωscale*ω(ωβ,jh)
    return sincfun{typeof(domain),T}(length(fϕv),h(n,T)/2,fϕv,ϕpv,ωv,ωscale,ωβ,jh,domain)
end
sincfun(f::Function) = sincfun(f,Finite())

## Evaluation by the barycentric formula

function barycentric{D<:Domain,T<:Number}(sf::sincfun{D,T},x::T)
    t = asinh(2sf.domain.ψinv(x)/π)
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

function Base.getindex{D<:Domain,T<:Number,T1<:Number}(sf::sincfun{D,T},x::T1)
    xc = convert(promote_type(T,T1),x)
    z = sf.domain.ψinv(xc)
    barycentric(sf,xc)*singularities(sf.domain,z)
end
Base.getindex{D<:Domain,T<:Number,T1<:Number}(sf::sincfun{D,T},x::Vector{T1}) = T[sf[xk] for xk in x]

# Helper routines

ω{T<:Number}(β::T,t::T) = exp(-β*(cosh(t)-one(T)))
@vectorize_2arg Number ω

h{T<:Number}(n::Integer,::Type{T}) = log(convert(T,π)*n)/n

envelope{D<:Domain,T<:Number}(sf::sincfun{D,T},t::T) = sf.ωscale*ω(sf.ωβ,t)/sf.domain.ψp(sinh(t)*π/2)/(cosh(t)*π/2)

function interlace{T<:Number}(u::Vector{T})
    n = length(u)
    v = zeros(T,n)
    if isodd(n)
        m = div(n,2) # n = 2m+1
        v[1] = u[1]
        for i = 1:m
          v[2i],v[2i+1] = u[i+m+1],u[i+1]
        end
    else
        m = div(n,2) # n = 2m
        for i=1:m
            v[2i-1],v[2i] = u[i+m],u[i]
        end
    end
    v
end

function interlace2{T<:Number}(u::Vector{T})
    n = length(u)
    v = zeros(T,n)
    if isodd(n)
        m = div(n,2) # n = 2m+1
        v[1] = u[m+1]
        for i = 1:m
          v[2i],v[2i+1] = u[m+i+1],u[m-i+1]
        end
    else
        m = div(n,2) # n = 2m
        for i=1:m
            v[2i-1],v[2i] = u[m-i+1],u[m+i]
        end
    end
    v
end

# Algebra

for op in (:+,:-,:*,:.*)
    @eval begin
        function $op{D<:Domain,T<:Number}(sf::sincfun{D,T},c::Number)
            sf1 = deepcopy(sf)
            sf1.fϕv = $op(sf.fϕv,convert(T,c))
            return sf1
        end
        function $op{D<:Domain,T<:Number}(c::Number,sf::sincfun{D,T})
            sf1 = deepcopy(sf)
            sf1.fϕv = $op(convert(T,c),sf.fϕv)
            return sf1
        end
        function $op{D<:Domain,T<:Number}(sf1::sincfun{D,T},sf2::sincfun{D,T})
            #TODO: sf.domain = $op(sf1.domain,sf2.domain)
            sincfun(x->$op(sf1[x],sf2[x]),sf1.domain)
        end
    end
end

# real, imag, conj

for op in (:(Base.real),:(Base.imag),:(Base.conj))
    @eval begin
        function $op{D<:Domain,T<:Number}(sf::sincfun{D,T})
            sf1 = deepcopy(sf)
            sf1.fϕv = $op(sf.fϕv)
            return sf1
        end
    end
end

# sum, cumsum, norm, dot, and diff.

function Base.sum{D<:Domain,T<:Number}(sf::sincfun{D,T})
    sinhv = sinh(sf.jh)*π/2
    singv = singularities(sf.domain,sinhv)
    sf.h*sum(sf.fϕv.*sf.ϕpv.*singv)
end

function Base.cumsum{D<:Domain,T<:Number}(sf::sincfun{D,T})
    SM = Sinc(-1,one(T)*[-sf.n+1:sf.n-1])
    #SM = Sinc(-1,one(T)*[-(sf.n-1)/2:(sf.n-1)/2].-[-(sf.n-1)/2:(sf.n-1)/2]')
    sf1 = deepcopy(sf)
    temp = sf.h*sf.fϕv.*sf.ϕpv
    [sf1.fϕv[i] = sum(SM[i+sf.n-1:-1:i].*temp) for i=1:sf.n]
    #sf1.fϕv = SM*temp
    return sf1
end

function Base.dot{D<:Domain,T<:Number}(sf1::sincfun{D,T},sf2::sincfun{D,T})
    sum(conj(sf1)*sf2)
end

function Base.norm{D<:Domain,T<:Number}(sf::sincfun{D,T})
    sinhv = sinh(sf.jh)*π/2
    singv = singularities(sf.domain,sinhv)
    sqrt(abs(sf.h*sum((sf.fϕv.*singv).^2.*sf.ϕpv)))
end

#=
function Base.diff{D<:Domain,T<:Number}(sf::sincfun{D,T})
    SM = Sinc(1,one(T)*[-sf.n+1:sf.n-1])
    #SM = Sinc(1,one(T)*[-(sf.n-1)/2:(sf.n-1)/2].-[-(sf.n-1)/2:(sf.n-1)/2]')
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
    SM = Sinc(1,one(T)*[-sf.n:sf.n])
    sinhv,coshv = convert(T,π)/2*sinh(sf.jh),convert(T,π)/2*cosh(sf.jh)
    temp = sf.domain.ψp(sinhv).*coshv
    phippoverphip = sinhv./coshv - 2tanh(sinhv).*coshv
    sf1 = deepcopy(sf)
    [sf1.fϕpv[i] = sum(SM[i:i+sf.n-1].*(sf.fϕpv/sf.h)) for i=1:sf.n]
    sf1.fϕpv = (sf1.fϕpv-sf.fϕpv.*phippoverphip)#./temp
    #sf1.fϕv = sf1.fϕv./(sf.ϕpv*sf.h)
    #sf1.ϕpv = ones(T,length(sf.ϕpv))
    return sf1
end
=#

Base.length(sf::sincfun) = sf.n

include("roots.jl")
include("fftBigFloat.jl")
include("KrylovMethods.jl")
include("steig.jl")
include("Sinc.jl")
include("SincMatrix.jl")
include("BandMatrix.jl")

end #module
