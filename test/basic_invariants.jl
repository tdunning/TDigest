scales = Dict(
    "K_0"=> (TDigest.K_0(), 1e-15),
    "K_1"=> (TDigest.K_1(), 2e-5),
    "K_2"=> (TDigest.K_2(), 1e-12),
    "K_3"=> (TDigest.K_3(), 1e-11)
)

for (name, s) in scales
    if name == "K_0"
        continue
    end
    scale, _ = s
    @info "Invariant test" name

    # demonstrate invariant as we add data
    m = TDigest.MergingDigest(50, scale)
    TDigest.checkWeights(m)
    TDigest.fit!(m, rand(1))
    TDigest.checkWeights(m)
    TDigest.fit!(m, rand(10))
    TDigest.checkWeights(m)
    TDigest.fit!(m, rand(100))
    TDigest.checkWeights(m)
    TDigest.fit!(m, rand(100))
    TDigest.checkWeights(m)
    TDigest.fit!(m, rand(1000))
    TDigest.checkWeights(m)
    TDigest.fit!(m, rand(1000_000))
    TDigest.checkWeights(m)

    # demonstrate invariant for merges
    m1 = TDigest.MergingDigest(50, scale)
    m2 = TDigest.MergingDigest(50, scale)
    TDigest.checkWeights(m)

    TDigest.fit!(m1, rand(1))
    TDigest.fit!(m2, rand(1_000))
    TDigest.merge!(m1, m2)
    TDigest.checkWeights(m1)

    m1 = TDigest.MergingDigest(50, scale)
    m2 = TDigest.MergingDigest(50, scale)
    TDigest.fit!(m1, rand(1_000))
    TDigest.fit!(m2, rand(1_000))
    TDigest.merge!(m1, m2)
    TDigest.checkWeights(m1)

    m1 = TDigest.MergingDigest(50, scale)
    m2 = TDigest.MergingDigest(50, scale)
    TDigest.fit!(m1, rand(1_000))
    TDigest.merge!(m1, m2)
    TDigest.checkWeights(m1)


    m1 = TDigest.MergingDigest(50, scale)
    m2 = TDigest.MergingDigest(50, scale)
    TDigest.fit!(m1, rand(1_000))
    TDigest.fit!(m2, rand(1_000))
    TDigest.merge!(m1, m2)
    TDigest.checkWeights(m1)

    m1 = TDigest.MergingDigest(50, scale)
    m2 = TDigest.MergingDigest(50, scale)
    TDigest.fit!(m1, rand(1_000_000))
    TDigest.fit!(m2, rand(1_000_000))
    TDigest.merge!(m1, m2)
    TDigest.checkWeights(m1)
end
