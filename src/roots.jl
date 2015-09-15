# This file calculates the roots of a sincfun in barycentric form
# based on the work of P.W. Lawrence, SIAM J. MATRIX ANAL. APPL., 34:1277--1300, 2013.

## TODO: This is simply the implementation via Julia's eig() routine. It could be faster
## using the specific algorithm in the paper. Also, once a differentiation command is
## implemented, it will be possible to find roots in BigFloat precision using Newton iteration.

export roots

function complexroots{D<:Domain,T<:Float64}(sf::sincfun{D,T})
    wv,xv,fv = (-one(T)).^(collect(0:sf.n-1)).*sf.ωv, ψ(sf.domain,convert(T,π)/2*sinh(sf.jh)), sf.fϕv.*sf.ϕpv
    cutoff = abs(wv) .≥ 10eps(T)
    wv,xv,fv = wv[cutoff],xv[cutoff],fv[cutoff]
    n = length(wv)
    wv,fv = wv/norm(wv),fv/norm(fv)
    A = sparse([zero(T) -fv'; wv diagm(xv)])
    s = [one(T);T[abs(fv[i]) == zero(T) ? one(T) : sqrt(abs(wv[i])/abs(fv[i])) for i=1:n]]
    S = sparse(collect(1:n+1),collect(1:n+1),s)
    B = diagm([zero(T);ones(T,n)])
    rts,V = eig(full(S\A*S),B)
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
    rts = sort(rts[a + 10sf.ωscale*eps(T) .< abs(rts) .< b - 10sf.ωscale*eps(T)])
end

# This is a first attempt to use Lawrence's O(n^2) algorithm to transform the arrowhead matrix
# to upper Hessenberg.

#=
function complexroots{D<:Domain,T<:Float64}(sf::sincfun{D,T})
    wv,xv,fv = (-one(T)).^([0:sf.n-1]).*sf.ωv, sf.jh, sf.fϕpv
    cutoff = abs(wv) .≥ 10eps(T)
    wv,xv,fv = wv[cutoff],xv[cutoff],fv[cutoff]
    #wvint,xvint,fvint = interlace2(wv),interlace2(xv),interlace2(fv)
    t,d,c = toupperhessenberg(wv,xv,fv)
    A = diagm([zero(T),d])+diagm(t,1)+diagm(t,-1)
    [A[1,i] += c[i] for i=1:length(c)]
    B = diagm([zero(T),ones(T,length(xv))])
    rts,V = eig(A,B)
    rts = rts[abs(rts).<Inf]
    rts = ψ(sf.domain,convert(T,π)/2*sinh(rts))
end

function toupperhessenberg{T<:Number}(w::Vector{T},x::Vector{T},f::Vector{T})
    @assert length(w) == length(x) == length(f)
    n = length(w)
    w,f = w/norm(w),f/norm(f)
    t,d = zeros(T,n),zeros(T,n)
    q = zeros(T,n,n)
    q[:,1] = w
    t[1] = norm(q[:,1])
    q[:,1] /= t[1]
    d[1] = dot(q[:,1],x.*q[:,1])
    q[:,2] = (x-d[1]).*q[:,1]
    t[2] = norm(q[:,2])
    q[:,2] /= t[2]
    for i=2:n-1
        d[i] = dot(q[:,i],x.*q[:,i])
        q[:,i+1] = (x-d[i]).*q[:,i]-t[i]*q[:,i-1]
        t[i+1] = norm(q[:,i+1])
        q[:,i+1] /= t[i+1]
    end
    d[n] = dot(q[:,n],x.*q[:,n])
    c = [zero(T), vec(-f'*q)]
    c[2] -= t[1]
    return t,d,c
end
=#
# Base.LinAlg.givensAlgorithm(f::T,g::T) returns cs, sn, r
