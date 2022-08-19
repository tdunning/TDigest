#
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


abstract type Histogram end

"""
Non-linear histogram that uses floating point representation plus a
quadratic correction to bin width to achieve tighter fit to the ideal
log2 sizing. Samples are assumed to be positive only.
"""
struct LogHistogram{K} <: Histogram
    min
    max
    logFactor
    logOffset
    counts::Vector{K}

    LogHistogram(min, max) = LogHistogram{Int64}(min, max, 0.1)

    function LogHistogram{K}(min,  max,  epsilonFactor) where K <: Integer
        logFactor, logOffset = _setupFactors(min, max, epsilonFactor)
        new{K}(min, max, logFactor, logOffset, _setupBins(K, logFactor, logOffset, min, max))
    end
end

"""
Non-linear histogram that uses floating point representation plus a
quadratic correction to bin width to achieve tighter fit to the ideal
log2 sizing. This is just like a `LogHistogram`(@Ref) except that signed
values are handled by keeping two arrays, one for negative values and
one for positive values.
"""
struct SignedLogHistogram{K} <: Histogram
    min
    max
    logFactor
    logOffset
    positive::Vector{K}
    negative::Vector{K}

    SignedLogHistogram(min, max) = LogHistogram{Int64}(min, max, 0.1, false)

    function SignedLogHistogram{K}( min,  max,  epsilonFactor) where K <: Integer
        logFactor, logOffset = _setupFactors(min, max, epsilonFactor)
        new{K}(min, max, logFactor, logOffset, _setupBins(K, logFactor, logOffset, min, max), _setupBins(K, logFactor, logOffset, min, max))
    end
end

getBounds(this::Histogram) = [lowerBound(this, i) for i in 1:length(this.counts)]
extrema(this::Histogram) = (this.min, this.max)

fit!(this::LogHistogram, v) = this.counts[bucket(this, v)] += 1

function fit!(this::SignedLogHistogram, v)
    if v ≥ 0
        this.positive[bucket(this, v)] += 1
    else
        this.negative[bucket(this, -v)] += 1
    end
end

function cdf(this::LogHistogram, x::Number)
    N = sum(this.counts)
    buckets = length(this.counts)
    bounds = lowerBound.(Ref(x), 1:buckets+1)
    marks = sum(bounds[1:end] .< x)

    (sum(marks) - marks[1]/2 - marks[end]/2) / N
end


function quantile(this::LogHistogram, q::Number)
    s = cumsum(this.counts) / sum(this.counts) .< q
    n = findfirst(s)
    (lowerBound(this, n) + lowerBound(this, n+1)) / 2
end

"""
Merges a bunch of histograms into this one. By enforcing exact conformity, this can be done with vector addition.
"""
function merge!(this::SignedLogHistogram, others::AbstractVector{SignedLogHistogram})
    for other in others
        if conformal(this, other)
            this.positive += other.positive
            this.negative += other.negative
        else
            throw(Argument("Can only merge histograms with identical bounds and precision"))
        end
    end
end


function merge!(this::LogHistogram, others::AbstractVector{LogHistogram}) 
    for other in others
        if conformal(this, other)
            this.counts .= this.counts + other.counts
        else
            throw(ArgumentException("Can only merge histograms with identical bounds and precision"))
        end
    end
end

function conformal(this::K, that::K) where K <: Histogram
    typeof(this) == typeof(that) &&
        extrema(this) == extrema(that) &&
        length(this.counts) == length(that.counts) &&
        this.allowSignedValues == that.allowSignedValues
end


################
# Internal functions

bucketIndexa(this::Histogram, x::Number) = bucketIndex(this.logFactor, this.logOffset, x)
bucketIndex(logFactor, logOffset, x) = Int32(ceil(approxLog2(x) * logFactor - logOffset))

lowerBound(this::Histogram, k) = pow2((k - 1 + this.logOffset) / this.logFactor)

function bucket(this::Histogram, x)
    if x ≤ this.min
        return 1
    elseif x ≥ this.max
        return length(this.counts)
    else
        return bucketIndex(this.logFactor, this.logOffset, x)
    end
end

"""
Internal function to set up scale factors for a log histogram

Returns `logFactor` and `logOffset` for direct inclusion in a constructor call.
"""
function _setupFactors(min, max, epsilonFactor)
    logFactor = log(2.0) / log1p(epsilonFactor)
    logOffset = approxLog2(min) * logFactor

    if max ≤ 2 * min
        throw(ArgumentException("Illegal/nonsensical min, max ($min, $max)"))
    end
    if min ≤ 0 || max ≤ 0
        throw(ArgumentException("Min and max must be positive"))
    end
    if epsilonFactor < 1e-6 || epsilonFactor > 0.5
        throw(ArgumentException(
            "Unreasonable number of bins per decade $epsilonFactor. Expected value in range [1e-6,0.5]"))
    end
    (logFactor, logOffset)
end

"""
Internal function to set up logarithmic bins
"""
function _setupBins(K, logFactor, logOffset, min, max) 
    binCount = bucketIndex(logFactor, logOffset, max) + 1;
    if binCount > 10000
        throw(ArgumentError("Excessive number of bins $binCount resulting from min,max = $min, $max"))
    end
    return zeros(K, binCount)
end

# Internal numerical functions

"""
Approximates log_2(value) by abusing floating point hardware. The floating point exponent
is used to get the integer part of the log. The mantissa is then adjusted with a second order
polynomial to get a better approximation. The error is bounded to be less than ±0.01 and is
zero at every power of two (which also implies the approximation is continuous).

@param value The argument of the log
@return log_2(value) (within an error of about ± 0.01)
"""
function approxLog2(value)
    valueBits = reinterpret(Int64, Float64(value))
    exponent = (Int64(valueBits & 0x7ff0_0000_0000_0000) >>> 52) - 1024
    m = reinterpret(Float64, (valueBits & 0x800fffffffffffff) | 0x3ff0000000000000)
    return m * (2 - (1.0 / 3) * m) + exponent - (2.0 / 3.0)
end

"""
Computes an approximate value of 2^x. This is done as an exact inverse of `approxLog2`(@Ref) so
that bin boundaries can be computed exactly.
"""
function pow2(x) 
    exponent = floor(x) - 1
    x = x - exponent
    m = 3 - sqrt(7 - 3 * x)
    return 2^(exponent + 1) * m
end
