
@info "CDF tests"

cdf_ref(data::Vector{<:Number}, x::Number) = (sum(data .< x) + sum(data .== x) / 2) / length(data)
function quantile_ref(data::Vector{<:Number}, q::Number)
    sort!(data)
    if q ≤ 0
        return data[1]
    elseif q ≥ 1
        return data[end]
    end
    return data[UInt32(ceil(q * length(data)))]
end

m = tDigest.MergingDigest(100)
data = [1.0,2.0,3.0,5.0]
tDigest.fit!(m, data)
@test length(m.sketch) == length(data)

@test tDigest.cdf(m, 0) == 0
@test tDigest.cdf(m, 10) == 1

@test tDigest.quantile(m, 0.0) == minimum(data)
@test tDigest.quantile(m, 1.0) == maximum(data)

for v in data
    steps = [prevfloat(v), v, nextfloat(v)]
    for x in steps
        @test tDigest.cdf(m, x) == cdf_ref(data, x)
    end
    qx = tDigest.cdf.(Ref(m), steps)
    @test tDigest.quantile(m, qx[1]) == v
    @test tDigest.quantile(m, qx[2]) == v
    @test tDigest.quantile(m, prevfloat(qx[3])) == v
end

function discrepancy(m, data, q)
    c1 = tDigest.cdf.(Ref(m), q)
    c2 = cdf_ref.(Ref(data), q)
    c1 - c2
end
