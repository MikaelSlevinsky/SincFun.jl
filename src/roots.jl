# This file calculates the roots of a sincfun in barycentric form
# based on the work of P.W. Lawrence, SIAM J. MATRIX ANAL. APPL., 34:1277--1300, 2013.

## TODO: This is simply the implementation via Julia's eig() routine. It could be faster
## using the specific algorithm in the paper. Also, once a differentiation command is
## implemented, it will be possible to find roots in BigFloat precision using Newton iteration.

export roots

function complexroots{D<:Domain,T<:Float64}(sf::sincfun{D,T})
    wv,xv,fv = (-one(T)).^(sf.j).*sf.ωv, sf.domain.ψ(convert(T,π)/2*sinh(sf.j*sf.h)), sf.fϕv.*sf.ϕpv
    cutoff = abs(wv) .≥ 10eps(T)
    wv,xv,fv = wv[cutoff],xv[cutoff],fv[cutoff]
    Dm = diagm(xv)
    A = [0.0 -fv';
         wv   Dm]
    B = diagm([zero(T),ones(T,length(xv))])
    rts,V = eig(A,B)
    rts = rts[abs(rts).<Inf]
end

function roots{D<:Domain,T<:Float64}(sf::sincfun{D,T})
    rts = complexroots(sf)
    rts = sort(real(rts[imag(rts) .== zero(T)]))
end

function roots{D<:Finite,T<:Float64}(sf::sincfun{D,T})
    rts = complexroots(sf)
    rts = real(rts[imag(rts) .== zero(T)])
    a,b = sf.domain.ab
    rts = sort(rts[a + 10eps(T) .< abs(rts) .< b - 10eps(T)])
end

# This is a first attempt to use Lawrence's O(n^2) algorithm to transform the arrowhead matrix
# to upper Hessenberg.

#=
function complexroots{D<:Domain,T<:Number}(sf::sincfun{D,T})
    wv,xv,fv = (-one(T)).^(sf.j).*sf.ωv, sf.domain.ψ(convert(T,π)/2*sinh(sf.j*sf.h)), sf.fϕv.*sf.ϕpv
    cutoff = abs(wv) .≥ 10eps(T)
    xv,fv,wv = xv[cutoff],fv[cutoff],wv[cutoff]

    t,d,c = toupperhessenberg(wv,xv,fv)

    println("Done upperhessenberging.")

    A = diagm([zero(T),d])+diagm(t,1)+diagm(t,-1)
    [A[1,i] += c[i] for i=1:length(c)]

    A = convert(Matrix{Float64},A)

    B = eye(T,size(A)...)
    B[1,1] = 0.0
    B = convert(Matrix{Float64},B)
    rts,V = eig(A,B)
    rts = rts[abs(rts).<Inf]
end

function toupperhessenberg{T<:Number}(w::Vector{T},x::Vector{T},f::Vector{T})
    @assert length(w) == length(x) == length(f)
    n = length(w)
    t,d,c = zeros(T,n),zeros(T,n),zeros(T,n)
    q = zeros(T,n,n)
    q[:,0+1] = w
    t[0+1] = norm(q[:,0+1])
    q[:,0+1] /= t[0+1]
    d[0+1] = dot(q[:,0+1],x.*q[:,0+1])
    q[:,1+1] = (x-d[0+1]).*q[:,0+1]
    t[1+1] = norm(q[:,1+1])
    q[:,1+1] /= t[1+1]
    for i=1:n-2
        d[i+1] = dot(q[:,i+1],x.*q[:,i+1])
        q[:,i+1+1] = (x-d[i+1]).*q[:,i+1]-t[i+1]*q[:,i-1+1]
        t[i+1+1] = norm(q[:,i+1+1])
        q[:,i+1+1] /= t[i+1+1]
    end
    d[n] = dot(q[:,n],x.*q[:,n])
    c = [zero(T), vec(-f'*q)]
    c[1+1] -= t[0+1]
    return t,d,c
end
=#
