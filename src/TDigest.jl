module TDigest

using Markdown
using LinearAlgebra



struct Centroid{T, K}
    mean::T
    count::K
end
defalg(v::AbstractArray{<:Centroid}) = MergeSort
Base.isless(a::Centroid{T, K}, b::Centroid{T, K}) where {T, K} = a.mean < b.mean

include("scale.jl")

md"""
    MergingDigest(compression[, scaleFunction])

Maintains a t-digest consisting of "centroids", each of which is
characterized by the mean of all points that have contributed to
that centroid as well as the total number of points that have
contributed. Points are added to centroids in such a way that centroids
nearly the minimum and maximum values seen are kept small while
centroids near the median value are allowed to be large.

This is done by collecting new points in a buffer as centroids with
unit weight and then occasionally sorting and merging old and new
centroids while maintaining the size limits. 

The overall t-digest is bounded in size, but can be used to approximate the
distribution of all of the data seen so far with the approximation being
particularly accurate near the tails.

Adding samples can be very fast because the cost of sorting and merging is
amortized over many insertions.

The virtues of this kind of t-digest implementation include:

* No allocation is required after initialization
* The data structure automatically compresses existing centroids when possible
* No dynamic allocation is necessary

The current implementation uses a single buffer to keep both
previously merged and retained centroids together with new samples.

New values are added to a digest using [`TDigest.fit!`](@ref). 

"""
mutable struct MergingDigest{T<:Number, K<:Number}
    function MergingDigest{T1,K1}(public::Real, private::Real,
                                  scale::ScaleFunction, maxSize::Int,
                                  logData::Bool) where {T1 <: Number, K1 <: Number} 
        new{T1, K1}(public, private, scale,
                    maxSize, 
                    logData, Vector{Vector{T1}}(), 0,
                    Vector{Centroid{T1, K1}}(),
                    0, true, 0, false, true)
    end
    
    publicCompression::Float64
    privateCompression::Float64
    scale::ScaleFunction

    # max number of centroids before we merge
    maxSize::Int

    # these hold a record of all samples we have seen
    logData::Bool
    log::Vector{Vector{T}}
    
    # sum_i weight[i]
    totalWeight::K

    # current centroids, guaranteed to meet the t-digest invariant
    # also, sorted in order of mean
    sketch::Vector{Centroid{T, K}}

    mergeCount::Int32

    # if true, alternate upward and downward merge passes
    useAlternatingSort

    # this records how many centroids survived the last merge
    watermark::Int

    # what order is the watermarked part of the digest?
    isReversed

    # if true, use higher working value of compression during
    # construction, then reduce on presentation
    useTwoLevelCompression
end

function MergingDigest(compression) 
    MergingDigest(compression, K_3())
end

function MergingDigest(compression, scale::ScaleFunction) 
    MergingDigest{Float64, Float64}(compression, 5*compression, true, scale)
end

function MergingDigest{T,K}(compression) where {T, K} 
    MergingDigest{T, K}(compression, K_3())
end

function MergingDigest{T,K}(compression, scale::ScaleFunction) where {T, K} 
    MergingDigest{T, K}(compression, 5*compression, true, scale)
end

function MergingDigest{T,K}(compression, maxPending, useTwoLevelCompression, scale) where {T,K}
    if compression < 10
        compression = 10
    end
    maxSize = 2*compression + max(50, maxPending)

    publicCompression = compression
    if useTwoLevelCompression
        compression = sqrt(maxSize / (2*compression)) * compression
    end
    MergingDigest{T, K}(publicCompression, compression, scale, maxSize, false)
end

weight(c::Centroid) = c.count

function max_step(digest::MergingDigest, q::Number, private=true)
    compression = private ? digest.privateCompression : digest.publicCompression
    max_step(digest.scale, q, compression, length(digest))
end    

function fit!(digest::MergingDigest{T, K}, vals::AbstractVector{<:Number}) where {T, K}
    if length(digest.sketch) > 0 && (digest.sketch[begin].count > 1 || digest.sketch[end].count > 1)
        @error "First or last centroid is too large" digest.sketch digest.mergeCount
        end
    if length(vals) > 10_000
        for first = 1:10_000:length(vals)
            last = min(length(vals), first + 9_999)
            fit!(digest, vals[first:last])
        end
    else
        if any(isnan.(vals))
            throw(ArgumentError("Cannot add NaN to t-digest"))
        end
        
        digest.totalWeight += length(vals)
        if length(vals) + length(digest.sketch) > digest.maxSize
            tmp = copy(digest.sketch)
            append!(tmp, [Centroid{T, K}(x, 1) for x in vals])
            if digest.logData
                log = copy(digest.log)
                append!(log, [[x] for x in vals])
            else
                log = digest.log
            end
            
            # merging on a temporary copy of the data avoids excess expansion of the sketch itself 
            mergeNewValues!(digest, tmp, log, digest.watermark, false, digest.privateCompression)
            
            copy!(digest.sketch, copy(tmp))
            copy!(digest.log, log)
        else
            append!(digest.sketch, [Centroid{T, K}(x, 1) for x in vals])
        end
    end
    if length(digest.sketch) > 0 && (digest.sketch[begin].count > 1 || digest.sketch[end].count > 1)
        @error "First or last centroid is too large" digest.sketch digest.mergeCount
    end
end

function fit!(digest::MergingDigest{T, K}, x) where {T, K}
    if isnan(x)
        throw(ArgumentError("Cannot add NaN to t-digest"))
    end
    push!(digest.sketch, Centroid{T, K}(x, 1))
    digest.totalWeight += 1
    if digest.logData
        push!(digest.log, [x])
    end
    if length(digest.sketch) > digest.maxSize
        mergeNewValues!(digest, false, digest.privateCompression)
    end
    return
end

function merge!(digest::MergingDigest, other::MergingDigest)
    if digest.logData && !other.logData
        throw(ArgumentError("Can't merge a digest that hasn't logged samples to one that has"))
    end

    digest.totalWeight += other.totalWeight
    if length(digest.sketch) + length(other.sketch) > digest.maxSize
        tmp = copy(digest.sketch)
        append!(tmp, other.sketch)

        log = copy(digest.log)
        append!(log, other.log)

        mergeNewValues!(digest, tmp, log, true, digest.privateCompression)

        copy!(digest.sketch, tmp)
        if digest.logData
            copy!(digest.log, log)
        end
    else
        append!(digest.sketch, other.sketch)
        if digest.logData
            append!(digest.log, other.log)
        end
    end
end

function mergeNewValues!(digest::MergingDigest) 
    mergeNewValues!(digest, digest.sketch, digest.log, digest.watermark, true, digest.privateCompression)
end

function mergeNewValues!(digest::MergingDigest, force::Bool, compression) 
    mergeNewValues!(digest, digest.sketch, digest.log, digest.watermark, force, compression)
end

md"""
Merge the clusters of a t-digest as much as possible. This has no effect 
if not enough data has been added since the last merge unless the merge
is forced.

The order of the clusters after this operation is either perfectly ascending 
or perfectly descending except if the merge is forced, then the result is
always in ascending order.
"""
function mergeNewValues!(digest::MergingDigest, sketch::Vector, log::Vector, watermark, force::Bool, compression) 
    if length(sketch) ≤ 1
        # nothing to do ... one element can't even be out of order
        return
    end
    if force || length(sketch) > digest.maxSize
        # note that we run the merge in reverse every other
        # merge to avoid left-to-right bias in merging
        reverseOrder = !force && digest.useAlternatingSort && digest.mergeCount % 2 == 1

        if digest.isReversed != reverseOrder
            # reversing the watermarked part of the sketch keeps the order of established centroids stable
            # this is very important with repeated values because otherwise the digest invariant can go
            # out of whack when centroids with identical means get reordered
            reverse!(@view(sketch[1:digest.watermark]))
        end
        if sketch[begin].count > 1 || sketch[end].count > 1
            @error "first or last centroid too large" sketch digest.isReversed reverseOrder digest.mergeCount
            throw(AssertionError("Boundary centroid not singleton"))
        end
        if digest.watermark > 0
            if sketch[digest.watermark].count > 1
                @show digest.sketch
                @error "watermarked centroid too large" sketch reverseOrder digest.mergeCount digest.watermark
                  throw(AssertionError("Boundary centroid not singleton"))
            end
        end
     
        if digest.logData
            order = sortperm(sketch, rev=reverseOrder)
            permute!(log, order)
            permute!(sketch, order)
        else
            sort!(sketch, rev=reverseOrder)
        end
        digest.isReversed = reverseOrder
        
        if !reverseOrder && length(sketch) < compression
            # we just sorted it so the sketch is a trivial t-digest
            # this return can only happen when we forced the sort,
            # but it is still worth checking
            digest.watermark = length(sketch)
            return
        end

        digest.mergeCount += 1

        # now pass through the data merging clusters where possible
        # but start partway in because we won't merge the first cluster
        total = digest.totalWeight
        norm = normalizer(digest.scale, compression, total)

        length(sketch) ≥ 2 || throw(AssertionError("Impossibly enough, sketch is too short")) 

        # w_so_far is total count up to and including `sketch[to]`
        # k0 is the scale up to but not including `sketch[to]`
        w_so_far = sketch[1].count
        k0 = k_scale(digest.scale, w_so_far / total, norm)
        w_so_far += sketch[2].count
        to = 2
        from = 3
        # the limiting weight is computed by finding a diff of 1 in scale
        limit = total * q_scale(digest.scale, k0 + 1, norm)
        while from <= length(sketch)
            from > to || throw(AssertionError("Impossible: from ≤ to"))
            
            dw = sketch[from].count
            if w_so_far + dw > limit || from == length(sketch)
                # can't merge due to size or due to forcing singleton at end
                to += 1
                k0 = k_scale(digest.scale, w_so_far / total, norm)
                limit = total * q_scale(digest.scale, k0 + 1, norm)

                if to < from
                    sketch[to] = sketch[from]
                    if digest.logData
                        log[to] = log[from]
                    end
                end
                from += 1
            else
                if to == 1
                    throw(AssertionError("merging first cluster!"))
                end
                # but we can merge other clusters if small enough
                sketch[to] = merge(sketch[to], sketch[from])
                if digest.logData
                    append!(log[to], log[from])
                end
                from += 1
            end
            w_so_far += dw
        end
        resize!(sketch, to)
        digest.watermark = length(sketch)
        
        length(sketch) < compression || begin
            @error "Merging was ineffective" compression sketch
            throw(AssertionError("Merging was ineffective"))
        end

        if sketch[begin].count > 1 || sketch[end].count > 1
            @error "First or last centroid is too large" sketch reverseOrder digest.mergeCount
        end
    end
end


md"""
Returns the minimum sample added to a digest
"""
function minimum(digest::MergingDigest)
    mergeNewValues!(digest)
    digest.sketch[1].mean
end


md"""
Returns the maximum sample added to a digest
"""
function maximum(digest::MergingDigest)
    mergeNewValues!(digest)
    digest.sketch[end].mean
end

md"""
Private function the combines two clusters
"""
function merge(a::Centroid{T, K}, b::Centroid{T, K}) where {T, K}
    total = a.count + b.count
    if total == 0
        Centroid{T, K}(0, 0)
    else
        μ = a.count / total
        Centroid{T, K}(μ * a.mean + (1-μ) * b.mean, total)
    end
end

md"""
Scans a digest to verify that the digest invariant is satisfied without 
compressing or sorting the data in the digest.
"""
function checkWeights(digest::MergingDigest)
    if length(digest) == 0
        return true
    end

    sizeof(digest.sketch) / sizeof(digest.sketch[1]) ≤ digest.maxSize ||
        throw(AssertionError("Digest sketch is oversized"))

    sum(map(c->c.count, digest.sketch)) ≈ digest.totalWeight ||
        throw(AssertionError("Digest has lost track of size"))

    if digest.logData
        sum(map(length, digest.log)) ≈ digest.totalWeight ||
            throw(AssertionError("Digest has lost track of samples"))
    else
        length(digest.log) == 0 ||
            throw(AssertionError("Digest has shouldn't have logged samples"))
    end
            
    if digest.isReversed
        tmp = sort(reverse(digest.sketch))
    else
        tmp = sort(digest.sketch)
    end
    
    (tmp[1].count == 1 && tmp[end].count == 1) ||
        @error "debug" tmp[1] tmp[end]

    scale = digest.scale
    norm = normalizer(digest.scale, digest.publicCompression, digest.totalWeight)
    q1 = 0
    k1 = k_scale(scale, q1 / digest.totalWeight, norm)
    
    i = 0
    for c in tmp
        i += 1
        q2 = q1 + c.count / digest.totalWeight
        k2 = k_scale(scale, q2, norm)

        c.count == 1 || k2 - k1 ≤ 1 || begin
            @show digest
            @error "Weight is too large" i c q1 q2 k1 k2
            throw(AssertionError("Weight is too large"))
        end

        q1 = q2
        k1 = k2
    end
    return true
end

md"""
Merges any pending inputs and compresses the data down to the public setting.
Note that this typically loses a bit of precision and thus isn't a thing to
be doing all the time. It is best done only when we want to persist the digest.
"""
compress(digest::MergingDigest) = mergeNewValues!(digest, true, digest.publicCompression)

Base.length(digest::MergingDigest) = length(digest.sketch)
         
md"""
Approximate the value of the empirical CDF at a value of `x`.
"""
function cdf(digest::MergingDigest, x)
    if isnan(x) || isinf(x)
        throw(ArgumentError("Invalid value: $x"))
    end
    mergeNewValues!(digest, true, digest.privateCompression)

    if length(digest) == 0
        # no data to examine
        return NaN
    elseif length(digest) == 1
        # exactly one centroid, should have max==min
        v = digest.sketch[1].center
        if x < v
            return 0.0
        elseif x > v
            return 1.0
        else
            return 0.5
        end
    else
        n = length(digest)
        min, max = first(digest.sketch).mean, last(digest.sketch).mean
        total = digest.totalWeight
        
        if x < min
            return 0.0
        elseif x == min
            w = 0.0
            for c in digest.sketch
                if c.mean > x
                    break
                end
                w += c.count
            end
            return w / 2.0 / total
        end

        if x > max
            return 1.0
        elseif x == max
            return 1 - 0.5 / total
        end
        min < x < max || throw(AssertionError("Can't happen"))

        # we know that there are at least two centroids and mean[0] < x < mean[n-1]
        # that means that there are either one or more consecutive centroids all at x
        # or there are consecutive centroids, c0 < x < c1
        weightSoFar = 0
        n = length(digest)
        for i in 1:n
            c1 = digest.sketch[i]
            c2 = digest.sketch[i + 1]
            if x == c1.mean
                dw = 0
                while i < n && c1.mean == x
                    dw += c1.count
                    i += 1
                    c1 = digest.sketch[i]
                end
                return (weightSoFar + dw/2) / total
            elseif c1.mean ≤ x < c2.mean
                # landed between centroids ... check for floating point madness
                if c2.mean - c1.mean > 0
                    # note how we handle singleton centroids here
                    # the point is that for singleton centroids, we know that their entire
                    # weight is exactly at the centroid and thus shouldn't be involved in
                    # interpolation
                    leftExcludedW = 0
                    rightExcludedW = 0
                    if c1.count == 1
                        if c2.count == 1
                            # two singletons means no interpolation
                            # left singleton is in, right is out
                            return (weightSoFar + 1) / total
                        else
                            leftExcludedW = 0.5
                        end
                    elseif c2.count == 1
                        rightExcludedW = 0.5
                    end
                    dw = (c1.count + c2.count) / 2

                    # can't have double singleton (handled that earlier)
                    dw > 1 || throw(AssertionError("Can't have double singleton"))
                    (leftExcludedW + rightExcludedW) <= 0.5 || 
                        throw(AssertionError("Can't have double singleton"))

                    # adjust endpoints for any singleton
                    left = c1.mean
                    right = c2.mean

                    dwNoSingleton = dw - leftExcludedW - rightExcludedW

                    # adjustments have only limited effect on endpoints
                    dwNoSingleton > dw / 2 || throw(AssertionError("Can't have excess effect of singletons"))

                    base = weightSoFar + c1.count / 2 + leftExcludedW
                    return (base + dwNoSingleton * (x - left) / (right - left)) / total
                else
                    # this is simply caution against floating point madness
                    # it is conceivable that the centroids will be different
                    # but too near to allow safe interpolation
                    dw = (weight[it] + weight[it + 1]) / 2
                    return (weightSoFar + dw) / totalWeight
                end
            else
                weightSoFar += c1.count
            end
        end
        throw(AssertionError("Can't happen ... loop fell through"))
    end
end

"""
Estimates the quantile that corresponds to a given value of q.

This is trickier than you might think because we have to deal with
loss of information where centroids are due to more than one sample
but we want to be as accurate as possible where a centroid has only
a single sample. We also want to handle centroids with two samples
well.
"""
function quantile(digest::MergingDigest, q::Number)
    if q < 0 || q > 1
            throw(ArgumentError("q should be in [0,1], got $q"))
    end
    mergeNewValues!(digest, true, digest.privateCompression)

    if length(digest.sketch) == 0
        # no centroids means no data, no way to get a quantile
        return NaN
    elseif length(digest.sketch) == 1 
        # with one data point, all quantiles lead to Rome
        return digest.sketch[1].mean
    end

    # we know that there are at least two centroids now
    n = length(digest.sketch)

    # if all of the values were stored in a sorted array, index would be the offset we are interested in
    totalWeight = digest.totalWeight
    index = q * totalWeight

    (digest.sketch[1].count == 1 && digest.sketch[end].count == 1) ||
        throw(AssertionError("Ill-formed digest has too large weight in first or last centroid"))

    # beyond the boundaries, we return min or max
    # the first centroid must have unit weight so this will make it moot
    if index < 1
        return digest.sketch[1].mean
    end

    if index > totalWeight - 1
        return digest.sketch[end].mean
    end


    1 ≤ index ≤ totalWeight - 1 ||
        throw(AssertionError("Impossible value of index"))

    # in between extremes we interpolate between centroids
    weightSoFar = digest.sketch[1].count / 2
    for i in 1:(length(digest.sketch)-1)
        dw = (digest.sketch[i].count + digest.sketch[i + 1].count) / 2;
        if weightSoFar + dw > index
            # centroids i and i+1 bracket our current point
            # this has to happen before loop completes because index ≤ totalWeight - 1

            # check for unit weight
            leftUnit = 0.0
            if digest.sketch[i].count == 1
                if index - weightSoFar < 0.5
                    # within the singleton's sphere
                    return digest.sketch[i].mean
                else
                    leftUnit = 0.5
                end
            end

            rightUnit = 0.0
            if digest.sketch[i + 1].count == 1
                if weightSoFar + dw - index <= 0.5
                    # no interpolation needed near singleton
                    return digest.sketch[i + 1].mean
                end
                rightUnit = 0.5
            end
            z1 = index - weightSoFar - leftUnit
            z2 = weightSoFar + dw - index - rightUnit
            return weightedAverageSorted(digest.sketch[i].mean, z2, digest.sketch[i + 1].mean, z1)
        end
        weightSoFar += dw
    end
    # if this loop completes then weightSoFar ≈ totalWeight - 1/2
    # but at the start, we knew that index ≤ totalWeight - 1
    @error "Fell out of the loop" index totalWeight
    throw(AssertionError("Fell out the bottom of the loop (which should be impossible)"))
end

"""
Same as {@link #weightedAverageSorted(double, double, double, double)} but flips
the order of the variables if `x2` is greater than
`x1`.
"""
function weightedAverage(x1, w1, x2, w2)
    if x1 ≤ x2
        return weightedAverageSorted(x1, w1, x2, w2)
    else
        return weightedAverageSorted(x2, w2, x1, w1)
    end
end

"""
Compute the weighted average between `x1` with a weight of
`w1` and `x2` with a weight of `w2`.
This expects `x1` to be less than or equal to `x2`
and is guaranteed to return a number in `[x1, x2]`. An
explicit check is required since this isn't guaranteed with floating-point
numbers.
"""
function weightedAverageSorted(x1, w1, x2, w2)
    x1 ≤ x2 || throw(AssertionError("Out of order values"))
    x = (x1 * w1 + x2 * w2) / (w1 + w2)
    return max(x1, min(x, x2))
end

            
md"""
How many bytes will it take to hold this t-digest in verbose form?
"""
function byteSize(digest::MergingDigest) 
    compress(digest)
    # format code, compression(float), buffer-size(int), temp-size(int), #centroids-1(int),
    # then two doubles per centroid
    return length(digest.sketch) * 16 + 32
end

md"""
How many bytes will it take to hold this t-digest in default form?
"""
function smallByteSize(digest::MergingDigest) 
    compress(digest)
    # format code(int), compression(float), buffer-size(short), temp-size(short), #centroids-1(short),
    # then two floats per centroid
    return length(digest.sketch) * 8 + 30;
end


md"""
Reads a MergingDigest
"""
read(io::IO, digest::MergingDigest)  = read(io, MergingDigest{Float64, Float64})

md"""
Reads a MergingDigest
"""
function read(io::IO, digest::MergingDigest{T, K}) where {T, K}
    encoding = read(io, UInt32)
    if encoding == Encoding.VERBOSE_ENCODING.code
        min = read(io, Float64)
        max = read(io, Float64)
        δ = read(io, Float64)
        n = read(io, UInt32)
        r = MergingDigest{T, K}(δ)
        total = 0.0
        if n > 0
            push!(r.sketch, Centroid{T, K}(min, 1))
            for i in 1:(n-2)
                w = read(io, Float64)
                m = read(io, Float64)
                push!(r.sketch, Centroid{T, K}(m, w))
                total += r.sketch[end].count
            end
            if n > 1
                push!(r.sketch, Centroid{T, K}(max, 1))
            end
        end
        r.totalWeight = total
        return r
    elseif encoding == Encoding.SMALL_ENCODING.code
        min = read(io, Float64)
        max = read(io, Float64)
        δ = read(io, Float32)
        n = read(io, UInt16)
        tmpSize = read(io, UInt16)
        r = MergingDigest{T, K}(compression, tmpSize - δ, true, K_3)
        total = 0.0
        if n > 0
            push!(r.sketch, Centroid{T, K}(min, 1))
            for i in 1:(n-2)
                w = read(io, Float32)
                m = read(io, Float32)
                push!(r.sketch, Centroid{T, K}(m, w))
                total += r.sketch[end].count
            end
            if n > 1
                push!(r.sketch, Centroid{T, K}(max, 1))
            end
        end
        r.totalWeight = total
        return r
    else
        throw(AssertionError("Invalid format for serialized t-digest"))
    end
end

function write(io::IO, digest::MergingDigest{T, K}) where {T, K}
    if n > 1
        digest.sketch[1].count == 1 ||
            throw(AssertionError("First centroid weight is > 1"))
        digest.sketch[max].count == 1 ||
            throw(AssertionError("Last centroid weight is > 1"))
        min = Float32(digest.sketch[1].mean)
        max = Float32(digest.sketch[end].mean)
    else
        min, max = NaN, NaN
    end
    compress(digest)
    write(io, Encoding.SMALL_ENCODING.code)    # 4
    write(io, min)
    write(io, max)
    write(io, Float32(digest.publicCompression))
    write(io, Int16(ceil(digest.publicCompression)))
    write(io, Int16(ceil(digest.maxSize - digest.publicCompression)))
    write(io, Int16(length(digest.sketch)))

    for i = 2:(length(digest.sketch)-1)
        write(io, Float32(digest.sketch[i].count))
        write(io, Float32(digest.sketch[i].mean))
    end
end

md"""
Writes a t-digest using 64-bit floats
"""
function writeVerbose(io::IO, digest::MergingDigest{T, K}) where {T, K}
    compress(digest)
    write(io, Encoding.VERBOSE_ENCODING.code)
    if n > 1
        digest.sketch[1].count == 1 ||
            throw(AssertionError("First centroid weight is > 1"))
        digest.sketch[max].count == 1 ||
            throw(AssertionError("Last centroid weight is > 1"))
        min = Float64(digest.sketch[1].mean)
        max = Float64(digest.sketch[end].mean)
    else
        min, max = NaN, NaN
    end

    write(io, min)
    write(io, max)
    write(io, Float64(digest.publicCompression))
    write(io, Int32(lastUsedCell))
    for i = 2:(length(digest.sketch)-1)
        write(io, Float64(digest.sketch[i].count))
        write(io, Float64(digest.sketch[i].mean))
    end
end

end # module
