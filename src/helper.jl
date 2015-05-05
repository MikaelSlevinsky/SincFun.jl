# Helper routines

real(x...) = Base.real(x...)
real{T<:Real}(::Type{T}) = T
real{T<:Real}(::Type{Complex{T}}) = T

eps{T<:Real}(::Type{T}) = Base.eps(T)
eps{T<:Real}(::Type{Complex{T}}) = eps(real(T))

#
# This function provides a convenient way to query or specify the BigFloat precision.
#
digits(n::Integer) = set_bigfloat_precision(int(ceil(n*log2(10))))
digits() = int(floor(get_bigfloat_precision()*log10(2)))


h{T<:Number}(::Type{T},n::Int) = log(convert(T,π)*n)/n
h(n::Int) = h(Float64,n)

ω{T<:Number}(β::T,t::T) = exp(-β*(cosh(t)-one(T)))
@vectorize_2arg Number ω


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
