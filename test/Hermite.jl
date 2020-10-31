@testset "Construction" for coeff in [
    Int64[1, 1, 1, 1],
    Float32[1.5, 3.0, -2.0],
    ComplexF64[1.5 + 2.0im, -3.0 + 2.5im],
    [1 // 7, -2 // 3, 2 // 2]
]
    p = Hermite(coeff)
    @test p.coeffs == coeff
    @test coeffs(p) == coeff
    @test degree(p) == length(coeff) - 1
    @test p.var == :x
    @test size(p) == size(coeff)
    @test size(P, 1) == size(coeff, 1)
    @test typeof(p).parameters[1] == eltype(coeff)
    @test eltype(p) == eltype(coeff)
end

@testset "Values" begin
    h1 = Hermite([1.5, 3.0, -1.0])  # 3.5 + 6x -4x^2
    x = -2:1:2
    expected = [-24.5, -6.5, 3.5,  5.5, -0.5]
    println(h1(x))
    @test h1(x) == expected
end
