export ConformalMap, differentiate

#
# A ConformalMap stores the data for h(t) = u₀sinh(t) + u₁ + u₂t + u₃ t^2 + ⋯ + uₙt^(n-1).
#
# u0 stores the sinh coefficient
# u stores the polynomial coefficients
# xpre stores the x-coordinates of the pre-images of the singularities
# O is the order of derivative
#

type ConformalMap{T}
    u0::T
    u::Vector{T}
    xpre::Vector{T}
    O::Int
end

ConformalMap{T}(u0::T,u::Vector{T},xpre::Vector{T}) = ConformalMap(u0,u,xpre,0)
ConformalMap{T}(u0::T,u::Vector{T}) = ConformalMap(u0,u,T[])

function differentiate{T}(h::ConformalMap{T})
    n = length(h.u)-1
    if n > 0
        u = Array(T,n)
        for j=1:n
            u[j] = j*h.u[j+1]
        end
        return ConformalMap(h.u0,u,h.xpre,h.O+1)
    else
        return ConformalMap(h.u0,T[],h.xpre,h.O+1)
    end
end
differentiate{T}(h::ConformalMap{T},n::Int) = n > 1 ? differentiate(differentiate(h),n-1) : n == 1 ? differentiate(h) : n == 0 ? h : throw("Cannot differentiate with order "*string(n)*".")

function evaluate{T}(h::ConformalMap{T},t::Number)
    n = length(h.u)
    ret = n > 0 ? convert(promote_type(T,typeof(t)),h.u[n]) : zero(promote_type(T,typeof(t)))
    for j = n-1:-1:1
        ret = h.u[j] + t*ret
    end
    ret += iseven(h.O) ? h.u0*sinh(t) : h.u0*cosh(t)
end

Base.getindex{T}(h::ConformalMap{T},t::Number) = evaluate(h,t)
Base.getindex{S,T<:Number}(h::ConformalMap{S},t::Array{T}) = map(t->evaluate(h,t),t)