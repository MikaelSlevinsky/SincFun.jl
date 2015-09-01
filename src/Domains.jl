export Domain, Finite, Infinite1, Infinite2, SemiInfinite1, SemiInfinite2
export ψ, ψinv, ψp

abstract Domain{T}

# For functions on a Finite domain with algebraic and logarithmic endpoint singularities, use the numbers α, β, γ, and δ to compute the singularities accurately.

type Finite{T} <: Domain{T}
    ab::Vector{T} # endpoints of the interval
    algebraic::Vector{T} # exponents of the algebraic singularities of the endpoints
    logarithmic::Vector{T} # exponents of the logarithmic singularities of the endpoints
end
Finite{T}(a::T,b::T,α::T,β::T,γ::T,δ::T) = Finite{T}([a,b],[α,β],[γ,δ]) # Algebraic and logarithmic endpoint singularities
Finite{T}(a::T,b::T,α::T,β::T) = Finite(a,b,α,β,zero(T),zero(T)) # Only algebraic endpoint singularities
Finite{T}(a::T,b::T) = Finite(a,b,zero(T),zero(T)) # No endpoint singularities
Finite{T}(::Type{T}) = Finite(-one(T),one(T)) # Canonical domain
Finite() = Finite(Float64) # Canonical precision

for dom in (:(Infinite1),:(Infinite2),:(SemiInfinite1),:(SemiInfinite2))
    @eval begin
        type $dom{T} <: Domain{T}
        end
        $dom{T}(::Type{T}) = $dom{T}()
        $dom() = $dom(Float64)
    end
end

# Every domain must have three functions: ψ, ψinv, and ψp.

ψ(d::Finite,t) = (d.ab[1]+d.ab[2])/2+(d.ab[2]-d.ab[1])/2*tanh(t)
ψinv(d::Finite,t) = atanh((2t-d.ab[1]-d.ab[2])/(d.ab[2]-d.ab[1]))
ψp(d::Finite,t) = (d.ab[2]-d.ab[1])/2*sech(t).^2

ψ(d::Infinite1,t) = sinh(t)
ψinv(d::Infinite1,t) = asinh(t)
ψp(d::Infinite1,t) = cosh(t)

ψ(d::Infinite2,t) = t
ψinv(d::Infinite2,t) = t
ψp(d::Infinite2,t) = 0*t+1

ψ(d::SemiInfinite1,t) = log1p(exp(t))
ψinv(d::SemiInfinite1,t) = log(expm1(t))
ψp(d::SemiInfinite1,t) = 1./(1+exp(-t))

ψ(d::SemiInfinite2,t) = exp(t)
ψinv(d::SemiInfinite2,t) = log(t)
ψp(d::SemiInfinite2,t) = exp(t)

# singularities is a function that should have an override for every domain that has endpoint singularities.

singularities{S,T<:Number}(domain::Domain{S},z::T) = one(promote_type(S,T))
singularities{S,T<:Number}(domain::Domain{S},z::Array{T}) = map(z->singularities(domain,z),z)

function singularities{S,T<:Number}(domain::Finite{S},z::T)
    α,β = domain.algebraic
    γ,δ = domain.logarithmic
    #(2/(exp(-2z)+1))^α*(2/(exp(2z)+1))^β*log(2/(exp(-2z)+1))^γ*log(2/(exp(2z)+1))^δ
    exp2z = exp(2z)
    temp1,temp2 = 2/(1/exp2z+1),2/(exp2z+1)
    (temp1)^α*(temp2)^β*log(temp1)^γ*log(temp2)^δ
end
