
@info "CDF tests"

m = TDigest.MergingDigest(100)
data = [1.0,2.0,3.0,5.0]
TDigest.fit!(m, data)
@test length(m.sketch) == length(data)

@test TDigest.cdf(m, 0) == 0
@test TDigest.cdf(m, 10) == 1

@test TDigest.quantile(m, 0.0) == minimum(data)
@test TDigest.quantile(m, 1.0) == maximum(data)

for v in data
    steps = [prevfloat(v), v, nextfloat(v)]
    for x in steps
        @test TDigest.cdf(m, x) == cdf_ref(data, x)
    end
    qx = TDigest.cdf.(Ref(m), steps)
    @test TDigest.quantile(m, qx[1]) == v
    @test TDigest.quantile(m, qx[2]) == v
    @test TDigest.quantile(m, prevfloat(qx[3])) == v
end

function discrepancy(m, data, q)
    c1 = TDigest.cdf.(Ref(m), q)
    c2 = cdf_ref.(Ref(data), q)
    c1 - c2
end
