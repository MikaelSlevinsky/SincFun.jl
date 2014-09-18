function Base.fft!{T<:BigFloat}(data::Vector{T})
	nn=int(length(data)/2)
	Tpi=convert(T,pi)
	n=nn << 1
	j=1
	for i=1:2:n-1
		if j>i
			data[j],data[i]=data[i],data[j]
			data[j+1],data[i+1]=data[i+1],data[j+1]
		end
		m=nn
		while m >= 2 && j > m
			j -= m
			m >>= 1
		end
		j += m
	end
	mmax=2
	while n > mmax
		istep=mmax << 1
		theta=-2Tpi/mmax
		wtemp=sin(theta/2)
		wpr = -2wtemp*wtemp
		wpi=sin(theta)
		wr=one(T)
		wi=zero(T)
		for m=1:2:mmax-1
			for i=m:istep:n
				j=i+mmax
				tempr=wr*data[j]-wi*data[j+1]
				tempi=wr*data[j+1]+wi*data[j]
				data[j]=data[i]-tempr
				data[j+1]=data[i+1]-tempi
				data[i] += tempr
				data[i+1] += tempi
			end
			wr=(wtemp=wr)*wpr-wi*wpi+wr
			wi=wi*wpr+wtemp*wpi+wi
		end
        mmax=istep
	end
	return data
end

function Base.fft{T<:BigFloat}(x::Vector{Complex{T}})
	n=length(x)
	y = Array(T,2n)
	for i = 1:n
		y[2i-1] = x[i].re
		y[2i] = x[i].im
	end
	Base.fft!(y)
	return complex(y[1:2:2n-1],y[2:2:2n])
end
Base.fft{T<:BigFloat}(x::Vector{T}) = Base.fft(complex(x))

function Base.ifft{T<:BigFloat}(x::Vector{Complex{T}})
	n=length(x)
	y = Array(T,2n)
	for i = 1:n
		y[2i-1] = x[i].re
		y[2i] = -x[i].im
	end
	Base.fft!(y)
	x = complex(y[1:2:2n-1],-y[2:2:2n])/n
end