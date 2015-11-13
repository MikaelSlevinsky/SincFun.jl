function steig!{T<:AbstractFloat}(d::Array{T,1}, e::Array{T,1}, z::Array{T,1}, maxits::Integer)
    #
    # Finds the eigenvalues and first components of the normalised
    # eigenvectors of a symmetric tridiagonal matrix by the implicit
    # QL method.
    #
    # d[i]   On entry, holds the ith diagonal entry of the matrix.
    #        On exit, holds the ith eigenvalue.
    #
    # e[i]   On entry, holds the [i+1,i] entry of the matrix for
    #        i = 1, 2, ..., n-1.  (The value of e[n] is not used.)
    #        On exit, e is overwritten.
    #
    # z[i]   On exit, holds the first component of the ith normalised
    #        eigenvector associated with d[i].
    #
    # maxits The maximum number of QL iterations.
    #
    # Martin and Wilkinson, Numer. Math. 12: 377-383 (1968).
    # Dubrulle, Numer. Math. 15: 450 (1970).
    # Handbook for Automatic Computation, Vol ii, Linear Algebra,
    #        pp. 241-248, 1971.
    #
    # This is a modified version of the Eispack routine imtql2.
    #
    n = length(z)
    z[1] = 1
    z[2:n] = 0
    e[n] = 0

    if n == 1 # Nothing to do for a 1x1 matrix.
        return
    end
    for l = 1:n
        for j = 1:maxits
            # Look for small off-diagonal elements.
            m = n
            for i = l:n-1
                if abs(e[i]) <= eps(T) * ( abs(d[i]) + abs(d[i+1]) )
                    m = i
                    break
                end
            end
            p = d[l]
            if m == l
                continue
            end
            if j == maxits
                msg = @sprintf("No convergence after %d iterations", j)
                msg *= " (try increasing maxits)"
                error(msg)
            end
            # Form shift
            g = ( d[l+1] - p ) / ( 2 * e[l] )
            r = hypot(g, one(T))
            g = d[m] - p + e[l] / ( g + copysign(r, g) )
            s = one(T)
            c = one(T)
            p = zero(T)
            for i = m-1:-1:l
                f = s * e[i]
                b = c * e[i]
                if abs(f) <  abs(g)
                    s = f / g
                    r = hypot(s, one(T))
                    e[i+1] = g * r
                    c = one(T) / r
                    s *= c
                else
                    c = g / f
                    r = hypot(c, one(T))
                    e[i+1] = f * r
                    s = one(T) / r
                    c *= s
                end
                g = d[i+1] - p
                r = ( d[i] - g ) * s + 2 * c * b
                p = s * r
                d[i+1] = g + p
                g = c * r - b
                # Form first component of vector.
                f = z[i+1]
                z[i+1] = s * z[i] + c * f
                z[i]   = c * z[i] - s * f
            end # loop over i
            d[l] -= p
            e[l] = g
            e[m] = zero(T)
        end # loop over j
    end # loop over l
end
