# include("mpc.jl")

using Test

function test_xref(debug=false)
  @testset "test_xref" begin
    times = LinRange(0, 10, 101)
    dt = 0.1
    n = 13
    m = 4
    Nmpc = 20
    x0 = [
      0.0,  10.0,  10.0,  
      1.0,  0.0,  0.0,  0.0,
      0.0,  0.0,  0.0,
      0.0,  0.0,  0.0
    ]
    xeq = [
      0.0,  0.0,  10.0,  
      1.0,  0.0,  0.0,  0.0,
      0.0,  0.0,  0.0,
      0.0,  0.0,  0.0
    ]

    Xref = [copy(x0) for i=1:length(times)]

    mpc1 = OSQPController(n, m, Nmpc, length(Xref))
    mpc1.times .= times
    mpc1.Xref .= Xref

    update_xref!(mpc1, x0, times[1], dt, xeq, times[end])
    @test mpc1.Xref[1][2] ≈ 10.0 atol=1e-6
    @test mpc1.Xref[2][2] ≈ 9.9 atol=1e-6
    @test mpc1.Xref[end][2] ≈ 0.0 atol=1e-6
    @test mpc1.Xref[end-1][2] ≈ 0.1 atol=1e-6

    xeq = [
      0.0,  0.0,  15.0,  
      1.0,  0.0,  0.0,  0.0,
      0.0,  0.0,  0.0,
      0.0,  0.0,  0.0
    ]
    k = get_k(mpc1, 5.0)
    end_k = get_k(mpc1, 7.5)
    update_xref!(mpc1, Xref[k], times[k], dt, xeq, times[end_k])
    @test mpc1.Xref[1][2] ≈ 10.0 atol=1e-6
    @test mpc1.Xref[2][2] ≈ 9.9 atol=1e-6
    @test mpc1.Xref[end][2] ≈ 0.0 atol=1e-6
    @test mpc1.Xref[end-1][2] ≈ 0.0 atol=1e-6
    @test mpc1.Xref[1][3] ≈ 10.0 atol=1e-6
    @test mpc1.Xref[2][3] ≈ 10.0 atol=1e-6
    @test mpc1.Xref[end_k][3] ≈ 15.0 atol=1e-6
    @test mpc1.Xref[end_k-1][3] ≈ 14.8 atol=1e-6
    @test mpc1.Xref[end][3] ≈ 15.0 atol=1e-6
    @test mpc1.Xref[end-1][3] ≈ 15.0 atol=1e-6

    if debug
      for i=1:length(times)
        println(mpc1.Xref[i])
      end
    end
  end
end