@testset "stability" begin
    rules = BindingRules(
        [1 3 2 3;
         2 1 4 1;
         2 2 3 2;
         3 1 4 1],
        UnitTriangle)   # 16 structures, largest has 5 particles
    asys = AssemblySystem(rules, EntropyModel(TreeLike()))
    percut(rxn) = inv(length(cut(rxn)))
    net = ReactionNetwork(asys; fwdkernel=percut)

    ns, nstr = nspecies(asys), nstructures(asys)
    ξ = [fill(-3.0, ns); fill(2.0, nbonds(asys))]
    N = @view compositionmatrix(asys)[:, 1:ns]

    # the basis defaults to DensityBasis, positionally or omitted
    @test stabilitymatrix(net, ξ) == stabilitymatrix(net, ξ, DensityBasis())
    @test stabilityjacobian(net, ξ) == stabilityjacobian(net, ξ, DensityBasis())

    for basis in (DensityBasis(), SymmetricBasis())
        symmetric = basis isa SymmetricBasis

        S = stabilitymatrix(net, ξ, basis)
        @test size(S) == (nstr, nstr)
        @test all(isfinite, S)
        symmetric && @test S ≈ S'

        # the in-place method must agree with the allocating one, and `scale` is a plain prefactor
        @test stabilitymatrix!(fill(NaN, nstr, nstr), net, ξ, basis) ≈ S
        @test stabilitymatrix(net, ξ, basis; scale=3) ≈ 3S

        # Each conserved particle count is a left null vector of S. The symmetric basis is inv(D)*S*D,
        # so there the null vectors pick up a factor of D rather than losing one.
        D = symmetric ? sqrt.(densities(asys, ξ)) : ones(nstr)
        Ncons = N .* D
        @test maximum(abs, Ncons' * S) < 1e-10 * maximum(abs, Ncons) * maximum(abs, S)

        # The analytic ξ-derivative must agree with differentiating the matrix itself
        J = stabilityjacobian(net, ξ, basis)
        @test size(J) == (nstr, nstr, length(ξ))
        @test stabilityjacobian!(fill(NaN, nstr, nstr, length(ξ)), net, ξ, basis) ≈ J
        @test stabilityjacobian(net, ξ, basis; scale=3) ≈ 3J
        fd = ForwardDiff.jacobian(x -> vec(stabilitymatrix(net, x, basis)), ξ)
        @test J ≈ reshape(fd, nstr, nstr, length(ξ))

        # nothing is cached, so a later `rate!` is picked up rather than going stale
        rate!(net; fwdkernel=Returns(2.0))
        @test !isapprox(stabilitymatrix(net, ξ, basis), S)
        rate!(net; fwdkernel=percut)
        @test stabilitymatrix(net, ξ, basis) ≈ S
    end

    # the two bases differ by a similarity transform, so they share a spectrum
    D = Diagonal(sqrt.(densities(asys, ξ)))
    Sρ = stabilitymatrix(net, ξ, DensityBasis())
    Ssym = stabilitymatrix(net, ξ, SymmetricBasis())
    @test Ssym ≈ D \ Sρ * D
    @test sort(real(eigvals(Sρ))) ≈ sort(eigvals(Symmetric(Ssym)))

    # deactivated reactions drop out of the sum
    quiet = ReactionNetwork(asys; fwdkernel=Returns(0.0), bwdkernel=Returns(0.0))
    @test !any(quiet.active)
    @test iszero(stabilitymatrix(quiet, ξ))
    @test iszero(stabilityjacobian(quiet, ξ))

    @test_throws ArgumentError stabilitymatrix(net, [0.0, 0.0])

    # the slowest mode is the least negative eigenvalue that is not one of the `ns` zero modes
    λs = sort(real(eigvals(Sρ)))
    @test count(<(1e-10), abs.(λs)) == ns    # one zero mode per conserved particle count
    τ = correlationtime(net, ξ)
    @test τ ≈ -inv(λs[nstr-ns])
    @test τ > 0

    # the two bases share a spectrum, so either matrix gives the same time
    @test correlationtime(Sρ; nconserved=ns) ≈ τ
    @test correlationtime(Symmetric(Ssym); nconserved=ns) ≈ τ

    # S is linear in the rates, so scaling them scales the rate of relaxation
    @test correlationtime(net, ξ; scale=2) ≈ τ / 2

    # a network that cannot react has nothing but conserved quantities
    @test correlationtime(quiet, ξ) == Inf
    @test correlationtime(zeros(4, 4); nconserved=1) == Inf

    # undercounting the zero modes selects one of them.
    @test -1e-10 < λs[nstr-ns+1] < 1e-10
    @test correlationtime(net, ξ; nconserved=ns - 1) == Inf
    @test correlationtime(Sρ; nconserved=ns - 1) == Inf
    @test isfinite(correlationtime(net, ξ; nconserved=ns - 1, rtol=1e-20))

    @test_throws ArgumentError correlationtime(Sρ; nconserved=nstr)

    # Nothing here is defined off equilibrium, and equilibrium is only a steady state under detailed
    # balance, which is also what makes the symmetric basis symmetric.
    @test isdetailedbalanced(net) && issymmetric(Ssym)
    driven = ReactionNetwork(asys; fwdkernel=Returns(1.0), bwdkernel=Returns(4.0))
    @test !isdetailedbalanced(driven)
    for basis in (DensityBasis(), SymmetricBasis())
        @test_throws ArgumentError stabilitymatrix(driven, ξ, basis)
        @test_throws ArgumentError stabilitymatrix!(zeros(nstr, nstr), driven, ξ, basis)
        @test_throws ArgumentError stabilityjacobian(driven, ξ, basis)
        @test_throws ArgumentError stabilityjacobian!(zeros(nstr, nstr, length(ξ)), driven, ξ, basis)
    end
    @test_throws ArgumentError correlationtime(driven, ξ)

    # rating it back into balance makes it usable again, since nothing about the check is cached
    rate!(driven; fwdkernel=percut)
    @test isdetailedbalanced(driven)
    @test stabilitymatrix(driven, ξ) ≈ Sρ
end;
