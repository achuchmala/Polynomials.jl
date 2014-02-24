# Poly type manipulations

module Polynomial
#todo: sparse polynomials?

export Poly, polyval, polyint, polyder, poly, roots, degree
export polydir #Deprecated

import Base: length, endof, getindex, setindex!, copy, zero, one, convert
import Base: string, show, *, /, //, -, +, ==, divrem, rem, eltype

eps{T}(::Type{T}) = convert(T,0)
eps{F<:FloatingPoint}(x::Type{F}) = Base.eps(F)

immutable Poly{T<:Number}
    a::Vector{T}
    nzfirst::Int #for effiencicy, track the first non-zero index
    var::Char
    function Poly(a::Vector{T}, var::Char)
        nzfirst = 0 #find and chop leading zeros
        for i = 1:length(a)
            if abs(a[i]) > 2*abs(eps(T)) #abs ensures valid comparison for complex
                break
            end
            nzfirst = i
        end
        new(a, nzfirst, var)
    end
end

#Poly{T<:Number}(a::Vector{T}) = Poly{T}(a, 'x')
Poly{T<:Number}(a::Vector{T}, var::Char='x') = Poly{T}(a, var)

convert{T}(::Type{Poly{T}}, p::Poly) = Poly(convert(Vector{T}, p.a), p.var)
promote_rule{T, S}(::Type{Poly{T}}, ::Type{Poly{S}}) = Poly{promote_type(T, S)}
eltype{T}(::Poly{T}) = T

length(p::Poly) = length(p.a)-p.nzfirst
endof(p::Poly) = length(p)
deg(p::Poly) = length(p) - 1

getindex(p::Poly, i) = p.a[i+p.nzfirst]
setindex!(p::Poly, v, i) = (p.a[i+p.nzfirst] = v)

copy(p::Poly) = Poly(copy(p.a[1+p.nzfirst:end]), p.var)

zero{T}(p::Poly{T}) = Poly([zero(T)], p.var)
zero{T}(::Type{Poly{T}}) = Poly([zero(T)])
one{T}(p::Poly{T}) = Poly([one(T)], p.var)
one{T}(::Type{Poly{T}}) = Poly([one(T)])

function _getcoefstr(value::Real)
    #Helper func, nicely format coefficient number for poly
    coefstr = @sprintf("%.4f", abs(value))
    #If the coefstr ends in 4 zeros, trim them (high precision not needed)
    if coefstr[end-3:] == "0000"
        coefstr = coefstr[1:end-5]
    end
    return coefstr
end

function _getcoefstr(value::Complex)
    #Helper func, nicely format coefficient number for poly
    realstr = @sprintf("%.4f", abs(value.re))
    imagstr = @sprintf("%.4f", abs(value.im))
    sgn = sign(value.im) < 0 ? "-" : "+"
    #If the coefstr ends in 4 zeros, trim them (high precision not needed)
    if realstr[end-3:] == "0000"
        realstr = realstr[1:end-5]
    end
    if imagstr[end-3:] == "0000"
        imagstr = imagstr[1:end-5]
    end
    #Create the coefstr for the current coefficient
    if imagstr == "0"
        coefstr = realstr
    else
        coefstr = "[$realstr $sgn $(imagstr)im]"
    end
    return coefstr
end

function string(p::Poly)
    #Convert polynomial into a string

    #Initalize the outstring to "0"
    outstr = "0"

    #Compute the number of coefficients
    N = length(p)

    for k=1:(N)
        #Get current coefficient in string form
        coefstr = _getcoefstr(p[k])
        #Power of current coefficient
        power = N-k
        if power == 0
            if coefstr != "0"
                newstr = coefstr
            else
                if k == 1
                    newstr = "0"
                else
                    newstr = ""
                end
            end
        elseif power == 1
            if coefstr == "0"
                newstr = ""
            elseif coefstr == "1"
                newstr = p.var
            else
                newstr = "$coefstr $(p.var)"
            end
        else
            if coefstr == "0"
                newstr = ""
            elseif coefstr == "1"
                newstr = "$(p.var)^$power"
            else
                newstr = "$coefstr $(p.var)^$power"
            end
        end

        if k > 1
            if newstr != ""
                if real(p[k]) < 0
                    outstr = "$outstr - $newstr"
                else
                    outstr = "$outstr + $newstr"
                end
            end
        elseif k == 1 && newstr != "" && real(p[k]) < 0
            outstr = "-$newstr"
        else
            outstr = newstr
        end
    end
    return outstr
end

function show(io::IO, p::Poly)
    n = length(p)
    print(io,"Poly($(string(p)))")
end

*(c::Number, p::Poly) = Poly(c * p.a[1+p.nzfirst:end], p.var)
*(p::Poly, c::Number) = Poly(c * p.a[1+p.nzfirst:end], p.var)
/(p::Poly, c::Number) = Poly(p.a[1+p.nzfirst:end] / c, p.var)
-(p::Poly) = Poly(-p.a[1+p.nzfirst:end], p.var)

-(p::Poly, c::Number) = +(p, -c)
+(c::Number, p::Poly) = +(p, c)
function +(p::Poly, c::Number)
    if length(p) < 1
        return Poly([c,], p.var)
    else
        p2 = copy(p)
        p2.a[end] += c;
        return p2;
    end
end
function -(c::Number, p::Poly)
    if length(p) < 1
        return Poly([c,], p.var)
    else
        p2 = -p;
        p2.a[end] += c;
        return p2;
    end
end

function +{T,S}(p1::Poly{T}, p2::Poly{S})
    if p1.var != p2.var
        error("Polynomials must have same variable")
    end
    R = promote_type(T,S)
    n = length(p1)
    m = length(p2)
    if n > m
        a = Array(R, n)
        for i = 1:m
            a[n-m+i] = p1[n-m+i] + p2[i]
        end
        for i = 1:n-m
            a[i] = p1[i]
        end
    else
        a = Array(R, m)
        for i = 1:n
            a[m-n+i] = p1[i] + p2[m-n+i]
        end
        for i = 1:m-n
            a[i] = p2[i]
        end
    end
    Poly(a, p1.var)
end

function -{T,S}(p1::Poly{T}, p2::Poly{S})
    if p1.var != p2.var
        error("Polynomials must have same variable")
    end
    R = promote_type(T,S)
    n = length(p1)
    m = length(p2)
    if n > m
        a = Array(R, n)
        for i = 1:m
            a[n-m+i] = p1[n-m+i] - p2[i]
        end
        for i = 1:n-m
            a[i] = p1[i]
        end
    else
        a = Array(R, m)
        for i = 1:n
            a[m-n+i] = p1[i] - p2[m-n+i]
        end
        for i = 1:m-n
            a[i] = -p2[i]
        end
    end
    Poly(a, p1.var)
end

function *{T,S}(p1::Poly{T}, p2::Poly{S})
    if p1.var != p2.var
        error("Polynomials must have same variable")
    end
    R = promote_type(T,S)
    n = length(p1)
    m = length(p2)
    if n == 0 || m == 0
        return Poly(R[], p1.var)
    end
    a = zeros(R, n+m-1)
    for i = 1:length(p1)
        for j = 1:length(p2)
            a[i+j-1] += p1[i] * p2[j]
        end
    end
    Poly(a, p1.var)
end

function divrem{T, S}(num::Poly{T}, den::Poly{S})
    if num.var != den.var
        error("Polynomials must have same variable")
    end
    m = length(den)
    if m == 0
        throw(DivideError())
    end
    R = typeof(one(T)/one(S))
    n = length(num)
    deg = n-m+1
    if deg <= 0
        return zero(Poly{R}), convert(Poly{R}, num)
    end
    d = zeros(R, n)
    q = zeros(R, deg)
    r = zeros(R, n)
    r[:] = num.a[1+num.nzfirst:end]
    for i = 1:deg
        quot = r[i] / den[1]
        q[i] = quot
        if i > 1
            d[i-1] = 0
            r[i-1] = 0
        end
        for j = 1:m
            k = i+j-1
            elem = den[j]*quot
            d[k] = elem
            r[k] -= elem
        end
    end
    return Poly(q, num.var), Poly(r, num.var)
end
/(num::Poly, den::Poly) = divrem(num, den)[1]
rem(num::Poly, den::Poly) = divrem(num, den)[2]

function ==(p1::Poly, p2::Poly)
    if length(p1) != length(p2)
        return false
    elseif p1.var != p2.var
        return false
    else
        return p1.a[1+p1.nzfirst:end] == p2.a[1+p2.nzfirst:end]
    end
end

function polyval{T}(p::Poly{T}, x::Number)
    R = promote_type(T, typeof(x))
    lenp = length(p)
    if lenp == 0
        return zero(R)
    else
        y = convert(R, p[1])
        for i = 2:lenp
            y = p[i] + x.*y
        end
        return y
    end
end

polyval(p::Poly, v::AbstractVector) = map(x->polyval(p, x), v)

function polyint{T}(p::Poly{T}, k::Number=0)
    n = length(p)
    R = typeof(one(T)/1)
    a2 = Array(R, n+1)
    for i = 1:n
        a2[i] = p[i] / (n-i+1)
    end
    a2[end] = k
    Poly(a2, p.var)
end

@deprecate polydir polyder

function polyder{T}(p::Poly{T})
    n = length(p)
    if n > 0
        a2 = Array(T, n-1)
        for i = 1:n-1
            a2[i] = p[i] * (n-i)
        end
    else
        a2 = zeros(T, 0)
    end
    Poly(a2, p.var)
end

# create a Poly object from its roots
function poly{T}(r::AbstractVector{T}, var='x')
    n = length(r)
    c = zeros(T, n+1)
    c[1] = 1
    for j = 1:n
        c[2:j+1] = c[2:j+1]-r[j]*c[1:j]
    end
    return Poly(c, var)
end
poly(A::Matrix, var='x') = poly(eig(A)[1], var)

roots{T}(p::Poly{Rational{T}}) = roots(convert(Poly{promote_type(T, Float64)}, p))

# compute the roots of a polynomial
function roots{T}(p::Poly{T})
    R = promote_type(T, Float64)
    num_zeros = 0
    if length(p) == 0
        return zeros(R, 0)
    end
    while abs(p[end-num_zeros]) <= 2*eps(T)
        if num_zeros == length(p)-1
            return zeros(R, 0)
        end
        num_zeros += 1
    end
    n = length(p)-num_zeros-1
    if n < 1
        return zeros(R, length(p)-1)
    end
    companion = zeros(R, n, n)
    a0 = p[end-num_zeros]
    for i = 1:n-1
        companion[1,i] = -p[end-num_zeros-i] / a0
        companion[i+1,i] = 1;
    end
    companion[1,end] = -p[1] / a0
    D,V = eig(companion)
    r = zeros(eltype(D),length(p)-1)
    r[1:n] = 1./D
    return r
end

end # module Poly
