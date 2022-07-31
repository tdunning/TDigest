using TDigest
using Markdown

begin
    @info "small digest I"
    data = [15.0, 20.0, 32.0, 60.0]
    digest = TDigest.MergingDigest(200)
    TDigest.fit!(digest, data)

    @test TDigest.checkWeights(digest)
    @test 20 ≈ TDigest.quantile(digest, 0.4) atol=1e-10
    @test 20 ≈ TDigest.quantile(digest, 0.25) atol=1e-10
    @test 15 ≈ TDigest.quantile(digest, 0.25 - 1e-10) atol=1e-10
    @test 20 ≈ TDigest.quantile(digest, 0.5 - 1e-10) atol=1e-10
    @test 32 ≈ TDigest.quantile(digest, 0.5) atol=1e-10
end

begin
    @info "small digest II"
    data = [
                245, 246, 247.249, 240, 243, 248, 250, 241, 244, 245, 245, 247, 243, 242, 241,
                50100, 51246, 52247, 52249, 51240, 53243, 59248, 59250, 57241, 56244, 55245,
                56245, 575247, 58243, 51242, 54241
    ]

    digest = TDigest.MergingDigest(50)
    TDigest.fit!(digest, data)
    @test TDigest.checkWeights(digest)
    @test quantile_ref(data, 0.5) == TDigest.quantile(digest, 0.5)
end


md"""
See https://github.com/tdunning/t-digest/issues/114

The problem in that issue seems to have been due to adding samples with non-unit weight.
This resulted in a violation of the t-digest invariant.

The question in the issue about the origin of the shuffle still applies.

testQuantile
"""
begin
    @info "many repeated values and stable sorting"
    δ = 100
    samples = [1.0, 2.0, 2.0, 2.0, 2.0, 2.0, 3.0, 3.0, 3.0, 3.0, 3.0, 4.0, 5.0, 6.0, 7.0]
    for i in 1:2
        hist1 = TDigest.MergingDigest(δ)
        let data = []
            for j in 1:100
                append!(data, samples)
                TDigest.fit!(hist1, samples)
            end
            @test TDigest.checkWeights(hist1)
            
            hist2 = TDigest.MergingDigest(δ)
            if hist1.sketch[begin].count > 1 || hist1.sketch[end].count > 1
                @error "starting with first or last centroid too large" hist1.sketch hist1.mergeCount
            end
            TDigest.compress(hist1)
            @test TDigest.checkWeights(hist1)
            if hist1.sketch[begin].count > 1 || hist1.sketch[end].count > 1
                @error "starting with first or last centroid too large" hist1.sketch hist1.mergeCount
            end
            TDigest.merge!(hist2, hist1)
            @test TDigest.checkWeights(hist2)
            
            sort!(data)
            
            TDigest.compress(hist2)
            x1 = TDigest.quantile(hist1, 0.5)
            x2 = TDigest.quantile(hist2, 0.5)
            @test quantile_ref(data, 0.5) ≈ x1 atol=0.2
            @test x1 ≈ x2 atol=0.01
        end
    end
end

md"""
Brute force test that cdf and quantile give reference behavior in digest made up of all singletons.

singletonQuantiles()
"""
begin
    @info "small digest III" 
    data = collect(0:19)
    digest = TDigest.MergingDigest(100)
    TDigest.fit!(digest, data)
    @test TDigest.checkWeights(digest)

    for x in LinRange(TDigest.minimum(digest) - 0.1, TDigest.maximum(digest) + 0.1, 20000)
        @test cdf_ref(data, x) == TDigest.cdf(digest, x)
    end

    for q in LinRange(0, 1, 1001)
        @test quantile_ref(data, q) == TDigest.quantile(digest, q)
    end
end


md"""
Verifies behavior involving interpolation (or lack of same, really) between singleton centroids.

singleSingleRange()
"""
begin
    @info "small digest IV"
    digest = TDigest.MergingDigest(100)
    TDigest.fit!(digest, 1);
    TDigest.fit!(digest, 2);
    TDigest.fit!(digest, 3);
    @test TDigest.checkWeights(digest)

    # verify the cdf is a step between singletons
    @test 0.5 / 3.0 == TDigest.cdf(digest, 1)
    @test 1 / 3.0 == TDigest.cdf(digest, 1 + 1e-10)
    @test 1 / 3.0 == TDigest.cdf(digest, 2 - 1e-10)
    @test 1.5 / 3.0 == TDigest.cdf(digest, 2)
    @test 2 / 3.0 == TDigest.cdf(digest, 2 + 1e-10)
    @test 2 / 3.0 == TDigest.cdf(digest, 3 - 1e-10)
    @test 2.5 / 3.0 == TDigest.cdf(digest, 3)
    @test 1.0 == TDigest.cdf(digest, 3 + 1e-10)
end

md"""
Tests cases where min or max is not the same as the extreme centroid which has weight>1. In these cases min and
max give us a little information we wouldn't otherwise have.

singletonAtEnd()
"""
begin
    @info "singleton at end"
    digest = TDigest.MergingDigest(100)
    TDigest.fit!(digest, [1, 2, 3])
    @test TDigest.checkWeights(digest)

    @test 1 == TDigest.minimum(digest)
    @test 3 == TDigest.maximum(digest)
    @test 3 == length(digest.sketch)
    @test 0 == TDigest.cdf(digest,0)
    @test 0 == TDigest.cdf(digest, 1 - 1e-9)
    @test 0.5 / 3 ≈ TDigest.cdf(digest, 1) atol=1e-10
    @test 1.0 / 3 ≈ TDigest.cdf(digest, 1 + 1e-10) atol=1e-10
    @test 2.0 / 3 ≈ TDigest.cdf(digest, 3 - 1e-9)
    @test 2.5 / 3 ≈ TDigest.cdf(digest, 3)
    @test 1.0 ≈ TDigest.cdf(digest, 3 + 1e-9)

    TDigest.fit!(digest, 1)
    @test TDigest.checkWeights(digest)
    @test cdf_ref([1, 1, 2, 3], 1) ≈ TDigest.cdf(digest, 1)

    # normally min == mean[1] because weight[1] == 1
    # we can force this not to be true for testing
    digest = TDigest.MergingDigest(10, TDigest.K_0())
    data = [1.0, 1, 2, 3]
    TDigest.fit!(digest, [1, 1, 2, 3])
    @test TDigest.checkWeights(digest)

    for i in 1:100
        append!(data, [1.0, 2, 3])
        TDigest.fit!(digest, [1, 2, 3])
    end
    @test TDigest.checkWeights(digest)
    
    # This sample would normally be added to the first cluster that already exists
    # but there is special code in place to prevent that centroid from ever
    # having weight of more than one
    # As such, near q=0, cdf and quantiles
    # should reflect this single sample as a singleton
    TDigest.fit!(digest, 0)
    push!(data, 0)
    sort!(data)
    
    @test length(digest.sketch) > 0

    @test 0.0 == TDigest.minimum(digest)
    first_centroid = digest.sketch[1]
    @test 1 == first_centroid.count

    @test 0 ≈ TDigest.cdf(digest, 0 - 1e-9)
    @test cdf_ref(data, 0) ≈ TDigest.cdf(digest, 0) atol=1e-10
    @test cdf_ref(data, 1e-9) ≈ TDigest.cdf(digest, 1e-9) atol=1e-10

    @test 0 ≈ TDigest.quantile(digest, 0) atol=0
    let q = 0.5 / length(digest.sketch)
        @test quantile_ref(data, q) ≈ TDigest.quantile(digest, q) atol=0
    end
    let q = 1.0 / length(digest.sketch)
        @test quantile_ref(data, q - 1e-10) ≈ TDigest.quantile(digest, q - 1e-10) atol=0
        @test quantile_ref(data, q) ≈ TDigest.quantile(digest, q) atol=0
    end
    @test first_centroid.mean == 0.0

    TDigest.fit!(digest, 4)
    @test TDigest.checkWeights(digest)

    push!(data, 4)
    sort!(data)


    let x = TDigest.maximum(digest)
        last_centroid = digest.sketch[end]
        @test 1.0 == last_centroid.count
        @test 4 == last_centroid.mean
        @test cdf_ref(data, x + 1e-9) ≈ TDigest.cdf(digest, x + 1e-9) atol=0
        @test cdf_ref(data, x) ≈ TDigest.cdf(digest, x)
        @test cdf_ref(data, x - 1e-9) ≈ TDigest.cdf(digest, x - 1e-9)
    end
                      
    let x = TDigest.minimum(digest)
        @test cdf_ref(data, x) ≈ TDigest.cdf(digest, x) atol=0
    end

    @test 4.0 == TDigest.quantile(digest, 1)
    let ϵ = 0.5 / length(digest.sketch)
        @test quantile_ref(data, 1 - ϵ - 1e-10) == TDigest.quantile(digest, 1 - ϵ - 1e-10) 
        @test quantile_ref(data, 1 - ϵ) == TDigest.quantile(digest, 1 - ϵ) 
        @test quantile_ref(data, 1 - ϵ + 1e-10) == TDigest.quantile(digest, 1 - ϵ + 1e-10) 

        @test quantile_ref(data, 1 - 2ϵ - 1e-10) == TDigest.quantile(digest, 1 - 2ϵ - 1e-10) 
        @test quantile_ref(data, 1 - 2ϵ) == TDigest.quantile(digest, 1 - 2ϵ) 
        @test quantile_ref(data, 1 - 2ϵ + 1e-10) == TDigest.quantile(digest, 1 - 2ϵ + 1e-10) 
    end
end

md"""
The example from issue https://github.com/tdunning/t-digest/issues/167
"""
begin
    @info "many repeated values"
    d = TDigest.MergingDigest(100)
    data = []
    for i in 1:2
        TDigest.fit!(d, 9000)
        append!(data, 9000)
    end
    @test TDigest.checkWeights(d)
    for i in 1:11
        TDigest.fit!(d, 3000)
        append!(data, 3000)
    end
    @test TDigest.checkWeights(d)
    for i in 1:26
        TDigest.fit!(d, 1000)
        append!(data, 1000)
    end
    @test TDigest.checkWeights(d)
    @test quantile_ref(data, 0.9) == TDigest.quantile(d, 0.9)
    @test quantile_ref(data, 0.95) == TDigest.quantile(d, 0.95)
end
