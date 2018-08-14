export sinint, cosint

function sinint{T<:BigFloat}(z::T;n::Integer=2^6)
    Tπ=convert(T,π);dDE=Tπ/2;ga=one(T);b2=one(T);h=log(Tπ*dDE*ga*n/b2)/ga/n
    phi(t::T) = log(exp(Tπ/2*sinh(t))+one(T))
    phip(t::T) = Tπ/2*cosh(t)/(one(T)+exp(-Tπ/2*sinh(t)))
    phiinv(x::T) = asinh(2/Tπ*log(expm1(x)))
    val1,val2=zero(T),zero(T)
    for j=-n:n
        jh=j*h
        phijh = phi(jh)
        temp = exp(-phijh)*phip(jh)/(phijh^2+z^2)*(one(T)-(-1)^(one(T)*j)*exp(-phiinv(abs(z))))
        val1 += temp*z
        val2 += temp*z*phijh
    end
    retval = Tπ/2 - h*(val1*cos(z)+val2*sinc(z/Tπ))
    return z == zero(T) ? zero(T) : z >= zero(T) ? retval : retval-Tπ
end

function Base.sinc{T<:Number}(n::Integer,x::T)
#This program computes the nth sinc differentiation matrices.
    val = zero(T)
    if n ≥ 0
        if x == 0
            val = isodd(n) ? zero(T) : (-1)^(n/2)*convert(T,π)^n/(n+1)
        else
            sp,cp,πk,xk,den = sinpi(x),cospi(x),one(T)/π,one(T),1
            for k=0:n
                m = mod(k,4)
                if m == 0
                    val -= sp*πk*xk/den
                elseif m == 1
                    val += cp*πk*xk/den
                elseif m == 2
                    val += sp*πk*xk/den
                else
                    val -= cp*πk*xk/den
                end
                πk,xk,den = πk*π,xk*x,den*(k+1)
            end
            val *= iseven(n) ? -den/(n+1)/xk : den/(n+1)/xk
        end
    elseif n == -1
        val = 0.5+sinint(π*x)/π
    end
    return val
end
Base.sinc{T<:Number}(n::T,x::T) = sinc(rount(Int,n),x)
