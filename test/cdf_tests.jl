
@info "CDF tests"

cdf_ref(data::Vector{<:Number}, x::Number) = (sum(data .< x) + sum(data .== x) / 2) / length(data)

m = tDigest.MergingDigest(100)
data = [1.0,2.0,3.0,5.0]
tDigest.fit!(m, data)
@test length(m.sketch) == length(data)

@test tDigest.cdf(m, 0) == 0
@test tDigest.cdf(m, 10) == 1

for v in data
    for x in [prevfloat(v), v, nextfloat(v)]
        @test tDigest.cdf(m, x) == cdf_ref(data, x)
    end
end
