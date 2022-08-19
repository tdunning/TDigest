# Licensed to Ted Dunning under one or more
# contributor license agreements.  See the NOTICE file distributed with
# this work for additional information regarding copyright ownership.
# The ASF licenses this file to You under the Apache License, Version 2.0
# (the "License"); you may not use this file except in compliance with
# the License.  You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

using Test, Distributions
using TDigest

include("./SimpleCompression.jl")


# testEmpty
begin
    @info "Static histogram structure tests"
    x = TDigest.LogHistogram(1.0, 100.0)
    bins = TDigest.extrema(x)
    @test 1.0 ≈ bins[1]
    @test 100.0 ≈ bins[2]
    TDigest.lowerBound(x, length(x.counts)) ≥ 100.0
    @test length(x.counts) == 50
end

#
#     testFitToLog 
# 
# The point of this test is to make sure that the floating point representation
# can be used as a quick approximation of log_2
# 
begin
    @info "Test internal log approximation"
    # 4 bits, worst case is mid octave
    lowerBound = 1 / 16.0 * sqrt(2)
    open("log-fit.csv", "w") do out
        scale = 2^52
        x = 0.001
        print(out, "x,y1,y2\n")        
        while x < 10
            xz = reinterpret(UInt64, x)
            # the magic 0x3ff is the offset for the floating point exponent
            v1 = xz / scale - 0x3ff
            v2 = log2(x)
            write(out, "$x,$v1,$v2\n")
            @test v2 - v1 > 0
            @test v2 - v1 < lowerBound
            x *= 1.02;
        end
    end
end


begin
    @info "testCompression"
    n = 1_000_000
    h = TDigest.LogHistogram(1e-3, 10.0)

    for _ in 1:n
        TDigest.fit!(h, rand())
    end

    compressed = compress!(Vector{UInt64}(), h.counts)
    @test length(compressed) < 45

    uncompressed = uncompress!(Vector{UInt64}(), compressed);
    @test length(uncompressed) ≥ length(h.counts)
    @test all(h.counts .== uncompressed[1:length(h.counts)])
end

begin
    @info "approxLog2"
    let x1 = 1e-6
        for i in 1:1000
            @test log(x1) / log(2) ≈ TDigest.approxLog2(x1) atol=0.01
            x1 *= 1.0 + π / 100.0
        end
        @test x1 > 1e6
    end
end

begin
    @info "round trip"
    for x = 0.001:1e-3:100
        log = TDigest.approxLog2(x);
        roundTrip = TDigest.pow2(log);
        @test x ≈ roundTrip atol=1e-13
    end
end

begin
    h = TDigest.LogHistogram(0.001, 100.0)
    @test length(h.counts) == 122
    @test TDigest.lowerBound(h, 91) ≈ 5.332317239079148
    @test TDigest.lowerBound(h, 100) ≈ 12.501048545139563
    @test TDigest.bucket(h, 5.0) == 90

    e = extrema((k -> TDigest.lowerBound(h, k+1) / TDigest.lowerBound(h, k)).(1:121))
    @test e[1] > 1.08
    @test e[2] < 1.12
end

begin
    @info "counts"
    N = Normal(4,1)
    h = TDigest.LogHistogram(0.001, 100.0)
    TDigest.fit!.(Ref(h), 4 .+ randn(10000))

    expected(n) = (cdf(Normal(4,1), TDigest.lowerBound(h, n+1)) - cdf(Normal(4,1), TDigest.lowerBound(h, n)))*10000
    std_error(n) = let p = (cdf(Normal(4,1), TDigest.lowerBound(h, n+1)) - cdf(Normal(4,1), TDigest.lowerBound(h, n)))
        return sqrt(p * (1-p) * 10000)
    end
    actual = abs.(h.counts - expected.(1:length(h.counts)))

    @test sum(actual .> 2 * std_error.(1:length(h.counts))) < 10
    @test sum(actual .> 3 * std_error.(1:length(h.counts))) < 3
end
    



