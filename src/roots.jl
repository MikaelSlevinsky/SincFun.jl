# This file calculates the roots of a sincfun in barycentric form
# based on the work of P.W. Lawrence, SIAM J. MATRIX ANAL. APPL., 34:1277--1300, 2013.

## TODO: This is simply the implementation via Julia's eig() routine. It could be faster
## using the specific algorithm in the paper. Also, once a differentiation command is
## implemented, it will be possible to find roots in BigFloat precision using Newton iteration.

export roots

function complexroots{D<:Domain,T<:Float64}(sf::sincfun{D,T})
    dv = sf.domain.ψ(convert(T,π)/2*sinh(sf.j*sf.h))
    fv,wv = sf.fϕv.*sf.ϕpv,(-one(T)).^(sf.j).*sf.ωv
    cutoff = abs(fv) .≥ 10eps(T)
    dv,fv,wv = dv[cutoff],fv[cutoff],wv[cutoff]
    Dm = diagm(dv)
    A = [0.0 -fv';
         wv   Dm]
    B = eye(T,size(A)...)
    B[1,1] = 0.0
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
