# Every domain must have three functions: ψ, ψinv, and ψp.

export Domain,Finite,Infinite1,Infinite2,SemiInfinite1,SemiInfinite2

abstract Domain{T}

# For functions on a Finite domain with algebraic and logarithmic endpoint singularities, use the numbers α, β, γ, and δ to compute the singularities in a stable way.

type Finite{T} <: Domain{T}
    ψ::Function
    ψinv::Function
    ψp::Function
    ab::Vector{T} # endpoints of the interval
    algebraic::Vector{T} # exponents of the algebraic singularities of the endpoints
    logarithmic::Vector{T} # exponents of the logarithmic singularities of the endpoints
end
Finite{T<:Number}(a::T,b::T,α::T,β::T,γ::T,δ::T) = Finite{T}(t->(a+b)/2+(b-a)/2*tanh(t),t->atanh((2t-a-b)/(b-a)),t->(b-a)/2*sech(t).^2,[a,b],[α,β],[γ,δ]) # Algebraic and logarithmic endpoint singularities
Finite{T<:Number}(a::T,b::T,α::T,β::T) = Finite(a,b,α,β,zero(T),zero(T)) # Only algebraic endpoint singularities
Finite{T<:Number}(a::T,b::T) = Finite(a,b,zero(T),zero(T)) # No endpoint singularities
Finite{T<:Number}(::Type{T}) = Finite(-one(T),one(T)) # Canonical domain
Finite() = Finite(Float64) # Canonical precision

for dom in (:(Infinite1),:(Infinite2),:(SemiInfinite1),:(SemiInfinite2))
    @eval begin
        type $dom{T} <: Domain{T}
            ψ::Function
            ψinv::Function
            ψp::Function
        end
        $dom() = $dom(Float64)
    end
end

Infinite1{T}(::Type{T}) = Infinite1{T}(sinh,asinh,cosh)
Infinite2{T}(::Type{T}) = Infinite2{T}(identity,identity,t->0t+1)
SemiInfinite1{T}(::Type{T}) = SemiInfinite1{T}(t->log1p(exp(t)),t->log(expm1(t)),t->1./(1+exp(-t)))
SemiInfinite2{T}(::Type{T}) = SemiInfinite2{T}(exp,log,exp)

# singularities is a function that should have an override for every domain that has endpoint singularities.

singularities{S<:Number,T<:Number}(domain::Domain{S},z::T) = one(promote_type(S,T))
singularities{S<:Number,T<:Number}(domain::Domain{S},z::Vector{T}) = promote_type(S,T)[singularities(domain,zk) for zk in z]
singularities{S<:Number,T<:Number,N}(domain::Domain{S},z::Array{T,N}) = reshape(singularities(domain,vec(z)),size(z))

function singularities{S<:Number,T<:Number}(domain::Finite{S},z::T)
    α,β = domain.algebraic
    γ,δ = domain.logarithmic
    #(2/(exp(-2z)+1))^α*(2/(exp(2z)+1))^β*log(2/(exp(-2z)+1))^γ*log(2/(exp(2z)+1))^δ
    exp2z = exp(2z)
    temp1,temp2 = 2/(1/exp2z+1),2/(exp2z+1)
    (temp1)^α*(temp2)^β*log(temp1)^γ*log(temp2)^δ
end
