scales = Dict(
    "K_0"=> (tDigest.K_0(), 1e-15),
    "K_1"=> (tDigest.K_1(), 2e-5),
    "K_2"=> (tDigest.K_2(), 1e-12),
    "K_3"=> (tDigest.K_3(), 1e-11)
)

for (name, s) in scales
    if name == "K_0"
        continue
    end
    scale, _ = s
    @info "Invariant test" name

    # demonstrate invariant as we add data
    m = tDigest.MergingDigest(50, scale)
    tDigest.checkWeights(m)
    tDigest.fit!(m, rand(1))
    tDigest.checkWeights(m)
    tDigest.fit!(m, rand(10))
    tDigest.checkWeights(m)
    tDigest.fit!(m, rand(100))
    tDigest.checkWeights(m)
    tDigest.fit!(m, rand(100))
    tDigest.checkWeights(m)
    tDigest.fit!(m, rand(1000))
    tDigest.checkWeights(m)
    tDigest.fit!(m, rand(1000_000))
    tDigest.checkWeights(m)

    # demonstrate invariant for merges
    m1 = tDigest.MergingDigest(50, scale)
    m2 = tDigest.MergingDigest(50, scale)
    tDigest.checkWeights(m)

    tDigest.fit!(m1, rand(1))
    tDigest.fit!(m2, rand(1_000))
    tDigest.merge!(m1, m2)
    tDigest.checkWeights(m1)

    m1 = tDigest.MergingDigest(50, scale)
    m2 = tDigest.MergingDigest(50, scale)
    tDigest.fit!(m1, rand(1_000))
    tDigest.fit!(m2, rand(1_000))
    tDigest.merge!(m1, m2)
    tDigest.checkWeights(m1)

    m1 = tDigest.MergingDigest(50, scale)
    m2 = tDigest.MergingDigest(50, scale)
    tDigest.fit!(m1, rand(1_000))
    tDigest.merge!(m1, m2)
    tDigest.checkWeights(m1)


    m1 = tDigest.MergingDigest(50, scale)
    m2 = tDigest.MergingDigest(50, scale)
    tDigest.fit!(m1, rand(1_000))
    tDigest.fit!(m2, rand(1_000))
    tDigest.merge!(m1, m2)
    tDigest.checkWeights(m1)

    m1 = tDigest.MergingDigest(50, scale)
    m2 = tDigest.MergingDigest(50, scale)
    tDigest.fit!(m1, rand(1_000_000))
    tDigest.fit!(m2, rand(1_000_000))
    tDigest.merge!(m1, m2)
    tDigest.checkWeights(m1)
end
