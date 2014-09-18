module Sincfun

import Base

export sincfun,Domain,Finite,Infinite,SemiInfinite1,SemiInfinite2

type Domain
	psi::Function
	psiinv::Function
	psip::Function
end

Finite = Domain(t->tanh(t),t->atanh(t),t->sech(t).^2)
Infinite = Domain(t->sinh(t),t->asinh(t),t->cosh(t))
SemiInfinite1 = Domain(t->log(exp(t)+1),t->log(exp(t)-1),t->1./(1+exp(-t)))
SemiInfinite2 = Domain(t->exp(t),t->log(t),t->exp(t))

#For functions on a Finite domain with algebraic and logarithmic endpoint singularities, use the numbers α, β, γ, and δ to compute the singularities in a stable way.
type sincfun{F<:Function,D<:Domain,T<:Number}
	n::Integer
	h::T
	fphiv::Vector{T}
	phipv::Vector{T}
	ωv::Vector{T}
	j::Vector{Int64}
	domain::D
	α::T
	β::T
	γ::T
	δ::T
	ω::F
	function sincfun{F<:Function,D<:Domain,T<:Number}(f::F,domain::D,::Type{T},α::T,β::T,γ::T,δ::T)
		Tone,Tpi = one(T),convert(T,pi)
		ω(t) = exp(-Tpi/2*cosh(t))
		n=2^3
		h = log(Tpi*Tpi/2*Tone*n/(Tpi/2))/Tone/n
		jh = h*[-n:n]
		sinhv,coshv = Tpi/2*sinh(jh),Tpi/2*cosh(jh)
		fphiv = f(domain.psi(sinhv))
		phipv = domain.psip(sinhv).*coshv.*(2./(exp(2sinhv).+1)).^α.*(2./(exp(-2sinhv).+1)).^β.*log(2./(exp(2sinhv).+1)).^γ.*log(2./(exp(-2sinhv).+1)).^δ
		ωv = ω(jh)
		test = abs(fphiv.*phipv)
		cutoff = !isinf(test).*!isnan(test)
		fphiv = fphiv[cutoff]
		phipv = phipv[cutoff]
		ωv = ωv[cutoff]
		intold,intnew = Inf,h*dot(fphiv,phipv)
		while n < 2^14 && abs(intold-intnew) > -log(eps(T)^3)*intnew*eps(T)
			n *= 2
			h = log(Tpi*Tpi/2*Tone*n/(Tpi/2))/Tone/n
			jh = h*[-n:n]
			sinhv,coshv = Tpi/2*sinh(jh),Tpi/2*cosh(jh)
			fphiv = f(domain.psi(sinhv))
			phipv = domain.psip(sinhv).*coshv.*(2./(exp(2sinhv).+1)).^α.*(2./(exp(-2sinhv).+1)).^β.*log(2./(exp(2sinhv).+1)).^γ.*log(2./(exp(-2sinhv).+1)).^δ
			ωv = ω(jh)
			test = abs(fphiv.*phipv)
			cutoff = !isinf(test).*!isnan(test)
			fphiv = fphiv[cutoff]
			phipv = phipv[cutoff]
			ωv = ωv[cutoff]
			intold,intnew = intnew,h*dot(fphiv,phipv)
			println(intold," ",intnew," ",n)
		end
		j=[-n:n]
		n=length(cutoff[cutoff.==true])
		println(length(cutoff)," ",length(fphiv))
		j = j[cutoff]
		println(intnew," ",h*dot(fphiv,phipv)," ",n)
		new(n,h,fphiv,phipv,ωv,j,domain,α,β,γ,δ,ω)
	end
end
#General Domain constructors
sincfun{F<:Function,D<:Domain,T<:Number}(f::F,domain::D,::Type{T}) = sincfun{F,D,T}(f,domain,T,zero(T),zero(T),zero(T),zero(T))
sincfun{F<:Function,D<:Domain}(f::F,domain::D) = sincfun{F,D,Float64}(f,domain,Float64,0.0,0.0,0.0,0.0)

#Finite Domain constructors
sincfun{F<:Function,T<:Number}(f::F,::Type{T},α::T,β::T,γ::T,δ::T) = sincfun{F,Domain,T}(f,Finite,T,α,β,γ,δ)
sincfun{F<:Function,T<:Number}(f::F,::Type{T}) = sincfun{F,Domain,T}(f,Finite,T,zero(T),zero(T),zero(T),zero(T))
sincfun{F<:Function}(f::F,α::Float64,β::Float64,γ::Float64,δ::Float64) = sincfun{F,Domain,Float64}(f,Finite,Float64,α,β,γ,δ)
sincfun{F<:Function}(f::F) = sincfun{F,Domain,Float64}(f,Finite,Float64,0.0,0.0,0.0,0.0)

function barycentric{F<:Function,D<:Domain,T<:Number}(sf::sincfun{F,D,T},x::T)
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
			#=
			for j=1:length(sf.j)
				common = (t-sf.j[j]*sf.h)*(-one(T))^(sf.j[j])
				valN += sf.fphiv[j]*sf.phipv[j]/common
				valD += sf.ωv[j]/common
			end
			=#
			for j=1:2:sf.n
				common = t-sf.j[j]*sf.h
				valN += sf.fphiv[j]*sf.phipv[j]/common
				valD += sf.ωv[j]/common
			end
			for j=2:2:sf.n
				common = t-sf.j[j]*sf.h
				valN -= sf.fphiv[j]*sf.phipv[j]/common
				valD -= sf.ωv[j]/common
			end
			return sf.ω(t)/sf.domain.psip(Tpi/2*sinh(t))/(Tpi/2)/cosh(t)*valN/valD
		else
			return sf.ω(t)/sf.domain.psip(Tpi/2*sinh(t))/(Tpi/2)/cosh(t)*sf.fphiv[idx]*sf.phipv[idx]/sf.ωv[idx]
		end
	end
end
Base.getindex{F<:Function,D<:Domain,T<:Number}(sf::sincfun{F,D,T},x::T) = barycentric(sf,x)
function Base.getindex{F<:Function,D<:Domain,T<:Number}(sf::sincfun{F,D,T},x::Vector{T})
	n = length(x)
	sfx = zeros(T,n)
	for i=1:n
		sfx[i] = barycentric(sf,x[i])
	end
	return sfx
end
Base.length(sf::sincfun) = sf.n
Base.sum(sf::sincfun) = sf.h*dot(sf.fphiv,sf.phipv)

include("fftBigFloat.jl")
include("KrylovMethods.jl")
include("steig.jl")
include("Sinc.jl")
include("SincMatrix.jl")
include("BandMatrix.jl")

end #module