module Sincfun

import Base

export sincfun,Domain,Finite,Infinite,SemiInfinite1,SemiInfinite2

type Domain
	psi::Function
	psiinv::Function
	psip::Function
end

#For functions on a Finite domain with algebraic and logarithmic endpoint singularities, use the numbers α, β, γ, and δ to compute the singularities in a stable way.
Finite{T<:Number}(α::T,β::T,γ::T,δ::T) = Domain(t->tanh(t),t->atanh(t),t->sech(t).^2.*(2./(exp(2t).+1)).^α.*(2./(exp(-2t).+1)).^β.*log(2./(exp(2t).+1)).^γ.*log(2./(exp(-2t).+1)).^δ)
Infinite = Domain(t->sinh(t),t->asinh(t),t->cosh(t))
SemiInfinite1 = Domain(t->log(exp(t)+1),t->log(exp(t)-1),t->1./(1+exp(-t)))
SemiInfinite2 = Domain(t->exp(t),t->log(t),t->exp(t))

type sincfun{T<:Number,D<:Domain,F<:Function}
	n::Integer
	h::T
	fphiv::Vector{T}
	phipv::Vector{T}
	ωv::Vector{T}
	j::Vector{Int64}
	domain::D
	ω::F
	function sincfun{T<:Number,D<:Domain,F<:Function}(f::F,domain::D,TT::Type{T})
		Tone,Tpi = one(TT),convert(TT,pi)
		ω(t) = exp(-Tpi/2*cosh(t))
		n=2^8
		h = log(Tpi*Tpi/2*Tone*n/(Tpi/2))/Tone/n
		jh = h*[-n:n]
		sinhv,coshv = Tpi/2*sinh(jh),Tpi/2*cosh(jh)
		fphiv = f(domain.psi(sinhv))
		phipv = domain.psip(sinhv).*coshv
		ωv = ω(jh)
		intold,intnew = Inf,h*dot(fphiv.^2,phipv)
		while n < 2^16 && abs(intold-intnew) > -log(eps(TT))*intnew*eps(TT)
			n *= 2
			h = log(Tpi*Tpi/2*Tone*n/(Tpi/2))/Tone/n
			jh = h*[-n:n]
			sinhv,coshv = Tpi/2*sinh(jh),Tpi/2*cosh(jh)
			fphiv = f(domain.psi(sinhv))
			phipv = domain.psip(sinhv).*coshv
			ωv = ω(jh)
			intold,intnew = intnew,h*dot(fphiv.^2,phipv)
			println(intold," ",intnew," ",n)
		end
		j = [-n:n]
		test=abs(fphiv.*phipv)
		cutoff = !isinf(test).*!isnan(test)
		n=length(cutoff[cutoff.==true])
		println(length(cutoff)," ",length(fphiv))
		fphiv = fphiv[cutoff]
		phipv = phipv[cutoff]
		ωv = ωv[cutoff]
		j = j[cutoff]
		println(intnew," ",h*dot(fphiv.^2,phipv)," ",n)
		new(n,h,fphiv,phipv,ωv,j,domain,ω)
	end
end
sincfun{T<:Number,D<:Domain,F<:Function}(f::F,domain::D,TT::Type{T}) = sincfun{T,D,F}(f,domain,TT)
sincfun{D<:Domain,F<:Function}(f::F,domain::D) = sincfun{Float64,D,F}(f,domain,Float64)
sincfun{F<:Function}(f::F) = sincfun{Float64,Domain,F}(f,Finite(0,0,0,0),Float64)
sincfun{T<:Number,F<:Function}(f::F,TT::Type{T}) = sincfun{T,Domain,F}(f,Finite(0,0,0,0),TT)

function barycentric{T<:Number,D<:Domain,F<:Function}(sf::sincfun{T,D,F},x::T)
	Tpi = convert(T,pi)
	t = asinh(2/Tpi*sf.domain.psiinv(x))

	if t < sf.j[1]*sf.h
		return sf.fphiv[1]
	elseif t > sf.j[end]*sf.h
		return sf.fphiv[end]
	else
		idx = findfirst(sf.j*sf.h,t)
		if idx == 0
			valN = zero(T)
			valD = zero(T)
			
			for j=1:length(sf.j)
				common = (t-sf.j[j]*sf.h)*(-one(T))^(sf.j[j])
				valN += sf.fphiv[j]*sf.phipv[j]/common
				valD += sf.ωv[j]/common
			end
			
			#=
			for j=1:2:2sf.n+1
				common = t-sf.j[j]*sf.h
				valN += sf.fphiv[j]*sf.phipv[j]/common
				valD += sf.ωv[j]/common
			end
			for j=2:2:2sf.n
				common = t-sf.j[j]*sf.h
				valN -= sf.fphiv[j]*sf.phipv[j]/common
				valD -= sf.ωv[j]/common
			end
			=#
			return sf.ω(t)/sf.domain.psip(Tpi/2*sinh(t))/(Tpi/2)/cosh(t)*valN/valD
		else
			return sf.ω(t)/sf.domain.psip(Tpi/2*sinh(t))/(Tpi/2)/cosh(t)*sf.fphiv[idx]*sf.phipv[idx]/sf.ωv[idx]
		end
	end
end
Base.getindex{T<:Number,D<:Domain,F<:Function}(sf::sincfun{T,D,F},x::T) = barycentric(sf,x)
function Base.getindex{T<:Number,D<:Domain,F<:Function}(sf::sincfun{T,D,F},x::Vector{T})
	n = length(x)
	sfx = zeros(T,n)
	for i=1:n
		sfx[i] = barycentric(sf,x[i])
	end
	return sfx
end


include("fftBigFloat.jl")
include("KrylovMethods.jl")
include("steig.jl")
include("Sinc.jl")
include("SincMatrix.jl")
include("BandMatrix.jl")

end #module