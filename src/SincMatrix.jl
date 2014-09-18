export SincMatrix

type SincMatrix{T} <: AbstractMatrix{T}
	D::Matrix{T}
	TM::Matrix{T}
	TFFT::Matrix{Complex{T}}
	function SincMatrix{T}(D::Matrix{T},TM::Matrix{T})
		TMn,TMm = size(TM)
		Dn,Dm = size(D)
		TMn == 2Dn-1 && TMm == Dm || throw(DimensionMismatch(""))
		n,m=int((Dn-1)/2),Dm
		pad = zeros(T,4n-1)
		TFFT = zeros(Complex{T},8n,m)
		for i =1:m
			TFFT[:,i] = fft([TM[:,i],pad])
		end
		new(D,TM,TFFT)
	end
end
SincMatrix{T}(D::Matrix{T},TM::Matrix{T}) = SincMatrix{T}(D,TM)

function Base.ndims(M::SincMatrix)
	return 2
end
function Base.size(A::SincMatrix)
	An,Am = size(A.D)
	return An,An
end
function Base.size(A::SincMatrix,d::Integer)
	An,Am = size(A.D)
	return d == 1 ? An : An
end
function Base.getindex{T}(A::SincMatrix{T},k::Int64,j::Int64)
	n,m=size(A.D)
	temp = zero(T)
	for i=1:m
		temp += A.D[k,i]*A.TM[n+k-j,i]
	end
	return temp
end

function (*){T1<:Number,T2<:Number}(A::SincMatrix{T1},x::Vector{T2})
	nx = int((length(x)-1)/2)
	n,m = size(A.D)
	n = int((n-1)/2)
	n == nx || throw(DimensionMismatch(""))
	T = promote_type(T1,T2)
	xfft = fft([x,zeros(T,6n-1)])
	temp = zeros(T,2n+1)
	for i=1:m
		temp += A.D[:,i].*real(ifft(A.TFFT[:,i].*xfft))[2n+1:4n+1]
	end
	return temp
end