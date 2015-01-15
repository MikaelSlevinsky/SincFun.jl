module SincFun

import Base

include("Domains.jl")

export sincfun

type sincfun{D,T}
    n::Integer
    h::T
    fϕv::Vector{T}
    ϕpv::Vector{T}
    ωv::Vector{T}
    ωscale::T
    j::Vector{Int64}
    domain::D
end

ω{T<:Number}(t::T) = exp(-convert(T,π)/2*cosh(t))
@vectorize_1arg Number ω

function sincfun{T<:Number}(f::Function,domain::Domain{T})
    Tπ = convert(T,π)
    h(n) = log(Tπ*Tπ/2*n/(Tπ/2))/n

    n=2^3
    jh = h(n)*[-n:n]
    sinhv,coshv = Tπ/2*sinh(jh),Tπ/2*cosh(jh)
    fϕv = T[f(domain.ψ(z)) for z in sinhv]
    ϕpv = domain.ψp(sinhv).*coshv
    ωv = ω(jh)
    test = abs(fϕv.*ϕpv)
    cutoff = !isinf(test).*!isnan(test)
    fϕv = fϕv[cutoff]
    ϕpv = ϕpv[cutoff]
    ωscale = maxabs(fϕv.*ϕpv)
    ωv = ωscale*ωv[cutoff]
    intold,intnew = Inf,h(n)*dot(fϕv.^2,ϕpv)
    while n < 2^14 && abs(intold-intnew) > -3log(eps(T))*abs(intnew)*eps(T)
        n *= 2
        jh = h(n)*[-n:n]
        sinhv,coshv = Tπ/2*sinh(jh),Tπ/2*cosh(jh)
        fϕv = T[f(domain.ψ(z)) for z in sinhv]
        ϕpv = domain.ψp(sinhv).*coshv
        ωv = ω(jh)
        test = abs(fϕv.*ϕpv)
        cutoff = !isinf(test).*!isnan(test)
        fϕv = fϕv[cutoff]
        ϕpv = ϕpv[cutoff]
        ωscale = maxabs(fϕv.*ϕpv)
        ωv = ωscale*ωv[cutoff]
        intold,intnew = intnew,h(n)*dot(fϕv.^2,ϕpv)
    end
    j=[-n:n]
    n=length(cutoff[cutoff.==true])
    j = j[cutoff]
    return sincfun{typeof(domain),T}(n,h((n-1)/2),fϕv,ϕpv,ωv,ωscale,j,domain)
end
sincfun(f::Function) = sincfun(f,Finite())

## Evaluation by the barycentric formula

function barycentric{D<:Domain,T<:Number}(sf::sincfun{D,T},x::T)
    Tπ = convert(T,π)
    t = asinh(2/Tπ*sf.domain.ψinv(x))

    if t < sf.j[1]*sf.h
        return sf.fϕv[1]
    elseif t > sf.j[end]*sf.h
        return sf.fϕv[end]
    else
        idx = findfirst(sf.j*sf.h,t)
        if idx == 0
            valN = zero(T)
            valD = zero(T)
            for j=1:2:sf.n
                common = t-sf.j[j]*sf.h
                valN += sf.fϕv[j]*sf.ϕpv[j]/common
                valD += sf.ωv[j]/common
            end
            for j=2:2:sf.n
                common = t-sf.j[j]*sf.h
                valN -= sf.fϕv[j]*sf.ϕpv[j]/common
                valD -= sf.ωv[j]/common
            end
            return sf.ωscale*ω(t)/sf.domain.ψp(Tπ/2*sinh(t))/(Tπ/2)/cosh(t)*valN/valD
        else
            return sf.ωscale*ω(t)/sf.domain.ψp(Tπ/2*sinh(t))/(Tπ/2)/cosh(t)*sf.fϕv[idx]*sf.ϕpv[idx]/sf.ωv[idx]
        end
    end
end

function Base.getindex{D<:Domain,T<:Number,T1<:Number}(sf::sincfun{D,T},x::T1)
    xc = convert(promote_type(T,T1),x)
    z = sf.domain.ψinv(xc)
    barycentric(sf,xc)*singularities(sf.domain,z)
end
Base.getindex{D<:Domain,T<:Number,T1<:Number}(sf::sincfun{D,T},x::Vector{T1}) = T[sf[x[i]] for i=1:length(x)]

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

# sum, norm, and dot

function Base.sum{D<:Domain,T<:Number}(sf::sincfun{D,T})
    sinhv = convert(T,π)/2*sinh(sf.j*sf.h)
    singv = T[singularities(sf.domain,sinhv[j]) for j=1:length(sinhv)]
    sf.h*sum(sf.fϕv.*sf.ϕpv.*singv)
end

function Base.norm{D<:Domain,T<:Number}(sf::sincfun{D,T})
    sinhv = convert(T,π)/2*sinh(sf.j*sf.h)
    singv = T[singularities(sf.domain,sinhv[j]) for j=1:length(sinhv)]
    sqrt(abs(sf.h*sum(conj(sf.fϕv.*singv).*sf.fϕv.*singv.*sf.ϕpv)))
end

function Base.dot{D<:Domain,T<:Number}(sf1::sincfun{D,T},sf2::sincfun{D,T})
    sum(conj(sf1)*sf2)
end

Base.length(sf::sincfun) = sf.n

include("roots.jl")
include("fftBigFloat.jl")
include("KrylovMethods.jl")
include("steig.jl")
include("Sinc.jl")
include("SincMatrix.jl")
include("BandMatrix.jl")

end #module
