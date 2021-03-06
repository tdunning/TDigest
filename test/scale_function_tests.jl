using CSV, DataFrames

# check the bounding evaluator used in the scale functions
let domain = randn(100),
    range = Vector(),
    f = x -> push!(range, x)

    TDigest.limit.(f, domain, -1, 1)
    @test minimum(domain) < -1 && minimum(range) ≥ -1
    @test maximum(domain) > 1 && maximum(range) ≤ 1
end

# test that the scale functions have internal consistency and provide
# a reliable round trip
scales = Dict(
    "K_0"=> (TDigest.K_0(), 1e-15),
    "K_1"=> (TDigest.K_1(), 2e-5),
    "K_2"=> (TDigest.K_2(), 1e-12),
    "K_3"=> (TDigest.K_3(), 1e-11)
)

for (name, s) in scales
    scale, tolerance = s
    @info "scale function round trip" name
    for δ in [10, 30, 100, 300, 1000]
        for n in [10, 1000, 1_000_000, 1000_000_000]
            norm = TDigest.normalizer(scale, δ, n)
            let f1 = q -> TDigest.q_scale(scale, TDigest.k_scale(scale, q, δ, n), δ, n),
                f2 = q -> TDigest.q_scale(scale, TDigest.k_scale(scale, q, norm), norm),
                step1 = q -> TDigest.max_step(scale, q, δ, n),
                step2 = q -> TDigest.max_step(scale, q, norm),
                q = [[0, 1, 1e-10, 1-1e-10]; rand(10000)]

                # the two forms of forward scale are compatible
                k1 = TDigest.k_scale.(Ref(scale), q, norm)
                k2 = TDigest.k_scale.(Ref(scale), q, δ, n)
                @test maximum(abs.(k2 .- k1)) ≈ 0 atol = tolerance

                # the inverse functions are compatible
                q1 = TDigest.q_scale.(Ref(scale), k1, δ, n)
                q2 = TDigest.q_scale.(Ref(scale), k1, norm)
                @test maximum(abs.(q1 .- q2)) ≈ 0 atol = 1e-12

                # and the round trip is clean
                @test maximum(abs.(f1.(q) .- q)) ≈ 0 atol = tolerance
                @test maximum(abs.(f2.(q) .- q)) ≈ 0 atol = tolerance

                # check the max step as well 
                @test maximum(abs.(step1.(q) .- step2.(q))) ≈ 0 atol = tolerance
                @test all(f1.(min.(1, q + step1.(q))) .- f1.(q) .≤ 1)
                @test all(f1.(min.(1, q + step2.(q))) .- f1.(q) .≤ 1)

                @test all(f1.(q) .- f1.(max.(0, q - step1.(q))) .≤ 1)
                @test all(f1.(q) .- f1.(max.(0, q - step2.(q))) .≤ 1)
            end
        end
    end
end

# check basic shape of scale function
for (name, s) in scales
    scale, tolerance = s
    @info "scale function shape" name

    q = range(0, 1, 10_001)
    k = TDigest.k_scale.(Ref(scale), q, 100, 10000)

    @test issorted(k)
    d2k = diff(diff(k))
end
    
@info "Check scale function against data reference"
ref = CSV.read("ref.csv", DataFrame, comment="#")
for row in eachrow(ref)
    q = row[:q]
    for (name, s) in scales
        scale, _ = s
        @test TDigest.k_scale(scale, q, 100, 10_000) ≈ row[name] atol = 1e-5
    end
end
