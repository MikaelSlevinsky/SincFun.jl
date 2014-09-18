export BandMatrix

type BandMatrix{T} <: AbstractMatrix{T}
	B::Matrix{T}
	m::Integer
	n::Integer
	function BandMatrix{T}(A::AbstractMatrix{T},m::Integer)
		n,An1 = size(A)
		B = zeros(T,n,2m+1)
		for k=1:m
			for j=1:k+m
				B[k,j-k+m+1] = A[k,j]
			end
		end
		for k=m+1:n-m
			for j=k-m:k+m
				B[k,j-k+m+1] = A[k,j]
			end
		end
		for k=n-m+1:n
			for j=k-m:n
				B[k,j-k+m+1] = A[k,j]
			end
		end
		new(B,m,n)
	end
end
BandMatrix{T}(A::AbstractMatrix{T},m::Integer) = BandMatrix{T}(A,m)


function Base.ndims(B::BandMatrix)
	return 2
end
function Base.size(B::BandMatrix)
	return B.n,B.n
end
function Base.size(B::BandMatrix,d::Integer)
	return B.n
end

function Base.getindex{T}(B::BandMatrix{T},k::Int64,j::Int64)
	abs(k-j) > B.m ? zero(T) : B.B[k,j-k+B.m+1]
end

function Base.setindex!{T}(B::BandMatrix{T},val::T,k::Int64,j::Int64)
	if abs(k-j) <= B.m
		B.B[k,j-k+B.m+1] = val
	end
end

function (\){T1<:Number,T2<:Number}(B::BandMatrix{T1},b::Vector{T2})
	A = deepcopy(B)
	x = deepcopy(b)
	bandelim!(A,x)
	return x
end

function bandelim!{T1<:Number,T2<:Number}(A::BandMatrix{T1},b::Vector{T2})
	m,n = A.m,A.n
	n == length(b) || throw(DimensionMismatch(""))
	T=promote_type(T1,T2)
	b = convert(Vector{T},b)
	for i=2:m
		for j=0:m-1
			coef = A[i+j,i-1]/A[i-1,i-1]
			for k=i-m-1:m
				A[i+j,i+k] -= coef*A[i-1,i+k]
			end
			b[i+j] -= coef*b[i-1]
		end
	end
	for i=m+1:n-m
		for j=0:m-1
			coef = A[i+j,i-1]/A[i-1,i-1]
			for k=-m:m
				A[i+j,i+k] -= coef*A[i-1,i+k]
			end
			b[i+j] -= coef*b[i-1]
		end
	end
	for i=n-m+1:n
		for j=0:n-i
			coef = A[i+j,i-1]/A[i-1,i-1]
			for k=-m:n-i
				A[i+j,i+k] -= coef*A[i-1,i+k]
			end
			b[i+j] -= coef*b[i-1]
		end
	end
	
	for i=n:-1:m+1
		b[i] /=A[i,i]
		A[i,i] = one(T)
		for j=1:m
			b[i-j] -= A[i-j,i]*b[i]
			A[i-j,i] = zero(T)
		end
	end
	for i=m:-1:1
		b[i] /=A[i,i]
		A[i,i] = one(T)
		for j=1:i-1
			b[i-j] -= A[i-j,i]*b[i]
			A[i-j,i] = zero(T)
		end
	end
end