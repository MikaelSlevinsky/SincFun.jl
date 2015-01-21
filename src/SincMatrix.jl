export SincMatrix

#
# A SincMatrix is a sum of products of Diagonal and Toeplitz matrices.
# S = D1*T1 + D2*T2 + ... + Dm*Tm
# Diagonal and Toeplitz matrices can be stored as vectors. Therefore,
# m diagonal or Toeplitz matrices can be stored as a single matrix.
#

type SincMatrix{T} <: AbstractMatrix{T}
    D::Matrix{T}
    TM::Matrix{T}
    TFFT::Matrix{Complex{T}}
end

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
    SincMatrix(D,TM,TFFT)
end

SincMatrix{T}(D::Matrix{T},TM::Vector{T}) = SincMatrix(D,reshape(TM,length(TM),1))
SincMatrix{T}(D::Vector{T},TM::Matrix{T}) = SincMatrix(reshape(D,length(D),1),TM)
SincMatrix{T}(D::Vector{T},TM::Vector{T}) = SincMatrix(reshape(D,length(D),1),reshape(TM,length(TM),1))


Base.ndims(A::SincMatrix) = 2
Base.size(A::SincMatrix) = size(A.D)[1],size(A.D)[1]
Base.size(A::SincMatrix,d::Integer) = d == 1 || d == 2 ? size(A.D)[1] : throw(DimensionMismatch(""))

function Base.getindex{T}(A::SincMatrix{T},k::Int64,j::Int64)
    n,m=size(A.D)
    sum(A.D[k,1:m]*A.TM[n+k-j,1:m])
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
