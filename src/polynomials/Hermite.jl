export Hermite

"""
    Hermite{<:Number}(coeffs::AbstractVector, var=:x)

Hermite polynomials https://en.wikipedia.org/wiki/Hermite_polynomials (the physicsts' kind)


Construct a polynomial from its coefficients
"""
struct Hermite{T <: Number} <: AbstractPolynomial{T}
    coeffs::Vector{T}
    var::Symbol
    function Hermite{T}(coeffs::AbstractVector{T}, var::SymbolLike=:x) where {T <: Number}
        length(coeffs) == 0 && return new{T}(zeros(T, 1), var)
        last_nz = findlast(!iszero, coeffs)
        last = max(1, last_nz === nothing ? 0 : last_nz)
        return new{T}(coeffs[1:last], var)
    end
end

@register Hermite

function Hermite{T}(coeffs::OffsetArray{T,1, Array{T,1}}, var::SymbolLike=:x) where {T <: Number}
    hs = zeros(T, 1 + lastindex(coeffs))
    hs[1 .+ (firstindex(coeffs):lastindex(coeffs))] = coeffs.parent
    Hermite(hs, var)
end

"""
    (::Hermite)(x)

Evaluate the Hermite polynomial at 'x'.

"""
function (he::Hermite{T})(x::S) where {T,S}
    R = promote_type(T, S)
    length(he) == 0 && return zero(R)
    length(he) == 1 && return R(he[0])
    h0 = he[end - 1]
    h1 = he[end]
    @inbounds for i in lastindex(he) - 2:-1:0
        h0, h1 = he[i] - 2(i + 1)*h1, h0 + 2x * h1
    end
    return R(h0 + 2x * h1)
end


# println(Hermite(-2.0))
