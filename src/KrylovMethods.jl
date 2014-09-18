export cg,bicgstab,Lanczos

function cg{T1<:Number,T2<:Number}(A::AbstractMatrix{T1},b::Vector{T2})
	n=length(b)
	n1,n2=size(A)
	n == n1 == n2 || throw(DimensionMismatch(""))
	T=promote_type(T1,T2)
	x=zeros(T,n)
	r=b-A*x
	p=r
	rsold=dot(r,r)
	TR=typeof(real(rsold))
	for i=1:n
		Ap=A*p
		alpha=rsold/dot(p,Ap)
		x=x+alpha*p
		r=r-alpha*Ap
		rsnew=dot(r,r)
		if real(sqrt(rsnew))<10eps(TR)
			return x
		end
		p=r+rsnew/rsold*p
		rsold=rsnew
	end
	return x
end

function cg{T1<:Number,T2<:Number,T3<:Number}(A::AbstractMatrix{T1},b::Vector{T2},P::AbstractMatrix{T3})
	n=length(b)
	n1,n2=size(A)
	n3,n4=size(P)
	n == n1 == n2 == n3 == n4 || throw(DimensionMismatch(""))
	T=promote_type(T1,T2,T3)
	x=P\b
	r=b-A*x
	z=P\r
	p=z
	rsold=dot(r,z)
	TR=typeof(real(rsold))
	for i=1:n
		Ap=A*p
		alpha=rsold/dot(p,Ap)
		x=x+alpha*p
		r=r-alpha*Ap
		z=P\r
		rsnew=dot(r,z)
		if real(sqrt(rsnew))<10eps(TR)
			return x
		end
		p=z+rsnew/rsold*p
		rsold=rsnew
	end
	return x
end

function bicgstab{T1<:Number,T2<:Number}(A::AbstractMatrix{T1},b::Vector{T2})
	n=length(b)
	n1,n2=size(A)
	n == n1 == n2 || throw(DimensionMismatch(""))
	T=promote_type(T1,T2)
	x=zeros(T,n)
	r=b-A*x
	rhat = r
	ρold = one(T)
	α =one(T)
	ωold = one(T)
	v=zeros(T,n)
	p=zeros(T,n)
	TR = typeof(real(α))
	for i=1:n
		ρnew = dot(rhat,r)
		β = (ρnew/ρold)*(α/ωold)
		p = r + β*(p-ωold*v)
		v = A*p
		α = ρnew/dot(rhat,v)
		s = r - α*v
		t = A*s
		ωnew = dot(t,s)/dot(t,t)
		upd = α*p + ωnew*s
		x = x + upd
		if norm(upd)<10eps(TR)
			return x
		end
		r = s - ωnew*t
		ρold = ρnew
		ωold = ωnew
	end
	return x
end

function bicgstab{T1<:Number,T2<:Number,T3<:Number}(A::AbstractMatrix{T1},b::Vector{T2},P::AbstractMatrix{T3})
	n=length(b)
	n1,n2=size(A)
	n3,n4=size(P)
	n == n1 == n2 == n3 == n4 || throw(DimensionMismatch(""))
	T=promote_type(T1,T2,T3)
	x=P\b
	r=b-A*x
	rhat = r
	ρold = one(T)
	α =one(T)
	ωold = one(T)
	v=zeros(T,n)
	p=zeros(T,n)
	TR = typeof(real(α))
	for i=1:n
		ρnew = dot(rhat,r)
		β = (ρnew/ρold)*(α/ωold)
		p = r + β*(p-ωold*v)
		y = P\p
		v = A*y
		α = ρnew/dot(rhat,v)
		s = r - α*v
		z = P\s
		t = A*z
		zt = P\t
		ωnew = dot(zt,z)/dot(zt,zt)
		upd = α*y + ωnew*z
		x = x + upd
		if norm(upd)<10eps(TR)
			return x
		end
		r = s - ωnew*t
		ρold = ρnew
		ωold = ωnew
	end
	return x
end

function Lanczos{T<:Number}(A::AbstractMatrix{T})
	n,n1 = size(A)
	n == n1 || throw(DimensionMismatch(""))
	z = Array(T,n)
	vold = zero(T)
	vnew = [one(T),zeros(T,n-1)]
	beta = zeros(T,n)
	alpha = zeros(T,n)
	for j = 1:n-1
		w = A*vnew
		alpha[j] = dot(w,vnew)
		w = w - alpha[j]*vnew - beta[j]*vold
		beta[j+1] = norm(w)
		vold = vnew
		vnew = w/beta[j+1]
	end
	w = A*vnew
	alpha[n] = dot(w,vnew)
	@time steig!(alpha,[beta[2:end],zero(T)],z,10000)
	return sort!(alpha)
 end