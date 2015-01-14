# Every domain must have three functions: ψ, ψinv, and ψp.

export Domain,Finite,Infinite1,Infinite2,SemiInfinite1,SemiInfinite2

abstract Domain{T}

#For functions on a Finite domain with algebraic and logarithmic endpoint singularities, use the numbers α, β, γ, and δ to compute the singularities in a stable way.

type Finite{T} <: Domain{T}
    ψ::Function
    ψinv::Function
    ψp::Function
    ab::Vector{T} # endpoints of the interval
    algebraic::Vector{T} # exponents of the algebraic singularities of the endpoints
    logarithmic::Vector{T} # exponents of the logarithmic singularities of the endpoints
end
Finite{T<:Number}(a::T,b::T,α::T,β::T,γ::T,δ::T) = Finite{T}(t->(a+b)/2+(b-a)/2*tanh(t),t->atanh((2t-a-b)/(b-a)),t->(b-a)/2*sech(t).^2,[a,b],[α,β],[γ,δ]) #Algebraic and logarithmic endpoint singularities
Finite{T<:Number}(a::T,b::T,α::T,β::T) = Finite(a,b,α,β,zero(T),zero(T)) #Only algebraic endpoint singularities
Finite{T<:Number}(a::T,b::T) = Finite(a,b,zero(T),zero(T),zero(T),zero(T)) #No endpoint singularities
Finite() = Finite(-1.0,1.0) #Canonical Domain

type Infinite1{T} <: Domain{T}
    ψ::Function
    ψinv::Function
    ψp::Function
end
function Infinite1{T}(::Type{T})
    return Infinite1{T}(sinh,asinh,cosh)
end
Infinite1() = Infinite1(Float64)

type Infinite2{T} <: Domain{T}
    ψ::Function
    ψinv::Function
    ψp::Function
end
function Infinite2{T}(::Type{T})
    return Infinite2{T}(identity,identity,t->0t+1)
end
Infinite2() = Infinite2(Float64)

type SemiInfinite1{T} <: Domain{T}
    ψ::Function
    ψinv::Function
    ψp::Function
end
function SemiInfinite1{T}(::Type{T})
    return SemiInfinite1{T}(t->log(exp(t)+1),t->log(exp(t)-1),t->1./(1+exp(-t)))
end
SemiInfinite1() = SemiInfinite1(Float64)

type SemiInfinite2{T} <: Domain{T}
    ψ::Function
    ψinv::Function
    ψp::Function
end
function SemiInfinite2{T}(::Type{T})
    return SemiInfinite2{T}(exp,log,exp)
end
SemiInfinite2() = SemiInfinite2(Float64)
