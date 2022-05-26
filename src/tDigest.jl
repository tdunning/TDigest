module tDigest

using Markdown
using LinearAlgebra



struct Centroid{T, K}
    mean::T
    count::K
end
defalg(v::AbstractArray{<:Centroid}) = QuickSort
Base.isless(a::Centroid{T, K}, b::Centroid{T, K}) where {T, K} = a.mean < b.mean

include("scale.jl")

md""" 
Maintains a t-digest by collecting new points in a buffer that
is then sorted occasionally and merged into a sorted array that
contains previously computed centroids.

This can be very fast because the cost of sorting and merging is
amortized over several insertion. If we keep N centroids total and
have the input array is k long, then the amortized cost is something
like

$N/k + log k$

These costs even out when $N/k = log k$.  Balancing costs is often a
good place to start in optimizing an algorithm.  For different values
of compression factor, the following table shows estimated asymptotic
values of N and suggested values of k:

| Compression | N | k |
| 50 | 78 | 25 |
| 100 | 157 | 42 |
| 200 | 314 | 73 |
Sizing considerations for t-digest

The virtues of this kind of t-digest implementation include:

* No allocation is required after initialization
* The data structure automatically compresses existing centroids when possible
* No Java object overhead is incurred for centroids since data is kept in primitive arrays

The current implementation takes the liberty of using ping-pong
buffers for implementing the merge resulting in a substantial memory
penalty, but the complexity of an in place merge was not considered as
worthwhile since even with the overhead, the memory cost is less than
40 bytes per centroid which is much less than half what the
AVLTreeDigest uses and no dynamic allocation is required at all.
"""
struct MergingDigest{T, K} 
    function MergingDigest{T1,K1}(public::Real, private::Real,
                                  scale::ScaleFunction, maxSize::Int,
                                  logData::Bool) where {T1, K1} 
        new{T1, K1}(public, private, scale,
                    maxSize, 
                    logData, Vector{Vector{T1}}(), [0],
                    Vector{Centroid{T1, K1}}(),
                    [0], true, true)
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
    totalWeight::Vector{K}

    # current centroids, guaranteed to meet the t-digest invariant
    # also, sorted in order of mean
    sketch::Vector{Centroid{T, K}}

    mergeCount::Vector{Int32}

    # if true, alternate upward and downward merge passes
    useAlternatingSort

    # if true, use higher working value of compression during
    # construction, then reduce on presentation
    useTwoLevelCompression
end

function MergingDigest(compression) 
    MergingDigest(compression, K_3())
end

function MergingDigest(compression, scale::ScaleFunction) 
    MergingDigest{Float64, Int64}(compression, 5*compression, true, scale)
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

function max_step(digest::MergingDigest, q::Number, private=true)
    compression = private ? digest.privateCompression : digest.publicCompression
    max_step(digest.scale, q, compression, length(digest))
end    

function fit!(digest::MergingDigest{T, K}, vals::AbstractVector{<:T}) where {T, K}
    if any(isnan.(vals))
        throw(ArgumentError("Cannot add NaN to t-digest"))
    end
    append!(digest.sketch, [Centroid{T, K}(x, 1) for x in vals])
    digest.totalWeight[1] += length(vals)
    if digest.logData
        append!(digest.log, [[x] for x in vals])
    end
    if length(digest) > digest.maxSize
        mergeNewValues!(digest, false, digest.privateCompression)
    end
end

function fit!(digest::MergingDigest{T, K}, x) where {T, K}
    if isnan(x)
        throw(ArgumentError("Cannot add NaN to t-digest"))
    end
    push!(digest.sketch, Centroid{T, K}(x, 1))
    digest.totalWeight[1] += 1
    if digest.logData
        push!(digest.log, [x])
    end
    if length(digest) > digest.maxSize
        mergeNewValues!(digest, false, digest.privateCompression)
    end
end

function merge!(digest::MergingDigest, other::MergingDigest)
    append!(digest.sketch, other.sketch)
    if digest.logData
        if !other.logData
            throw(ArgumentError("Can't merge a digest that hasn't logged samples to one that has"))
        end
        append!(digest.log, other.log)
    end
end

md"""
Merge the clusters of a t-digest as much as possible. This has no effect 
if not enough data has been added since the last merge unless the merge
is forced.

The order of the clusters after this operation is either perfectly ascending 
or perfectly descending except if the merge is forced, then the result is
always in ascending order.
"""
function mergeNewValues!(digest::MergingDigest, force::Bool, compression) 
    if length(digest.sketch) < 2
        # nothing to do
        return
    end
    if force || length(digest.sketch) > digest.maxSize
        # note that we run the merge in reverse every other
        # merge to avoid left-to-right bias in merging
        reverse = !force && digest.useAlternatingSort && digest.mergeCount[1] % 2 == 1
        if digest.logData
            order = sortperm(digest.sketch, rev=reverse)
            digest.log = digest.log[order]
            digest.sketch = digest.sketch[order]
        else
            sort!(digest.sketch, rev=reverse)
        end
        digest.mergeCount[1] += 1

        # now pass through the data merging clusters where possible
        # but start partway in because we won't merge the first cluster
        total = digest.totalWeight[1]
        norm = normalizer(digest.scale, compression, total)

        length(digest.sketch) ≥ 2 || throw(AssertionError("Impossible")) 

        # w_so_far is total count up to and including `sketch[to]`
        # k0 is the scale up to but not including `sketch[to]`
        w_so_far = digest.sketch[1].count
        k0 = k_scale(digest.scale, w_so_far / total, norm)
        w_so_far += digest.sketch[2].count
        to = 2
        from = 3
        # the limiting weight is computed by finding a diff of 1 in scale 
        limit = total * q_scale(digest.scale, k0 + 1, norm)
        while from <= length(digest.sketch)
            from > to || throw(AssertionError("from ≤ to"))
            
            dw = digest.sketch[from].count
            if w_so_far + dw > limit || from == length(digest.sketch)
                # can't merge due to size or due to forcing singleton at end
                to += 1
                k0 = k_scale(digest.scale, w_so_far / total, norm)
                limit = total * q_scale(digest.scale, k0 + 1, norm)

                if to < from
                    digest.sketch[to] = digest.sketch[from]
                    if digest.logData
                        digest.log[to] = digest.sketch[from]
                    end
                end
                from += 1
            else
                # but we can merge other clusters if small enough
                digest.sketch[to] = merge(digest.sketch[to], digest.sketch[from])
                if digest.logData
                    append!(digest.sketch[to], digest.sketch[from])
                end
                from += 1
            end
            w_so_far += dw
        end
        resize!(digest.sketch, to)
        length(digest.sketch) < compression ||
            @error "Merging was ineffective" digest.sketch
    end
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
        return
    end

    tmp = sort(digest.sketch)
    
    (tmp[1].count == 1 && tmp[end].count == 1) ||
        @error "debug" tmp[1] tmp[end]

    scale = digest.scale
    norm = normalizer(digest.scale, digest.privateCompression, digest.totalWeight[1])
    q1 = 0
    k1 = k_scale(scale, q1, norm)
    
    
    i = 0
    for c in tmp
        i += 1
        q2 = q1 + c.count
        k2 = k_scale(scale, q2, norm)

        c.count == 1 || k2 - k1 ≤ 1 || begin
            @show tmp[max(1, i-3):min(end, i+3)]
            @warn "Weight is too large" i c q1 q2 k1 k2
            throw(AssertionError("Weight is too large $i"))
        end

        q1 = q2
        k1 = k2
    end
end

md"""
Merges any pending inputs and compresses the data down to the public setting.
Note that this typically loses a bit of precision and thus isn't a thing to
be doing all the time. It is best done only when we want to show results to
the outside world.
"""
compress(digest::MergingDigest) = mergeNewValues(digest, true, digest.publicCompression)

Base.length(digest::MergingDigest) = length(digest.sketch)
         
md"""
Approximate the value of the empirical CDF at a value of `x`.
"""
function cdf(digest::MergingDigest, x)
    if isnan(x) || isinf(x)
        throw(ArgumentError("Invalid value: $x"))
    end
    mergeNewValues(digest, true, digest.privateCompression)

    if length(digest) == 0
        # no data to examine
        return NaN
    elseif length(digest) == 1
        # exactly one centroid, should have max==min
        v = digest.sketch[1].center
        if x < v
            return 0
        elseif x > v
            return 1
        else
            return 0.5
        end
    else
        n = length(digest)
        min, max = first(digest.sketch), last(digest.sketch)
        total = digest.totalWeight[1]
        
        if x < min
            return 0
        elseif x == min
            return 0.5 / total
        end

        if x > max
            return 1
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
            if x == c.mean
                dw = 0
                while i < n && digest.sketch[i].mean == x
                    dw += digest.sketch[i].count
                    i += 1
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
                            return (weightSoFar + 1) / totalWeight
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
                    return (base + dwNoSingleton * (x - left) / (right - left)) / totalWeight
                else
                    # this is simply caution against floating point madness
                    # it is conceivable that the centroids will be different
                    # but too near to allow safe interpolation
                    dw = (weight[it] + weight[it + 1]) / 2
                    return (weightSoFar + dw) / totalWeight
                end
            else
                weightSoFar += weight[it];
            end
        end
        throw(AssertionError("Can't happen ... loop fell through"))
    end
end
"""
    @Override
    public double quantile(double q) {
        if (q < 0 || q > 1) {
            throw new IllegalArgumentException("q should be in [0,1], got " + q);
        }
        mergeNewValues();

        if (lastUsedCell == 0) {
            # no centroids means no data, no way to get a quantile
            return Double.NaN;
        } else if (lastUsedCell == 1) {
            # with one data point, all quantiles lead to Rome
            return mean[0];
        }

        # we know that there are at least two centroids now
        int n = lastUsedCell;

        # if values were stored in a sorted array, index would be the offset we are interested in
        final double index = q * totalWeight;

        # beyond the boundaries, we return min or max
        # usually, the first centroid will have unit weight so this will make it moot
        if (index < 1) {
            return min;
        }

        # if the left centroid has more than one sample, we still know
        # that one sample occurred at min so we can do some interpolation
        if (weight[0] > 1 && index < weight[0] / 2) {
            # there is a single sample at min so we interpolate with less weight
            return min + (index - 1) / (weight[0] / 2 - 1) * (mean[0] - min);
        }

        # usually the last centroid will have unit weight so this test will make it moot
        if (index > totalWeight - 1) {
            return max;
        }

        # if the right-most centroid has more than one sample, we still know
        # that one sample occurred at max so we can do some interpolation
        if (weight[n - 1] > 1 && totalWeight - index <= weight[n - 1] / 2) {
            return max - (totalWeight - index - 1) / (weight[n - 1] / 2 - 1) * (max - mean[n - 1]);
        }

        # in between extremes we interpolate between centroids
        double weightSoFar = weight[0] / 2;
        for (int i = 0; i < n - 1; i++) {
            double dw = (weight[i] + weight[i + 1]) / 2;
            if (weightSoFar + dw > index) {
                # centroids i and i+1 bracket our current point

                # check for unit weight
                double leftUnit = 0;
                if (weight[i] == 1) {
                    if (index - weightSoFar < 0.5) {
                        # within the singleton's sphere
                        return mean[i];
                    } else {
                        leftUnit = 0.5;
                    }
                }
                double rightUnit = 0;
                if (weight[i + 1] == 1) {
                    if (weightSoFar + dw - index <= 0.5) {
                        # no interpolation needed near singleton
                        return mean[i + 1];
                    }
                    rightUnit = 0.5;
                }
                double z1 = index - weightSoFar - leftUnit;
                double z2 = weightSoFar + dw - index - rightUnit;
                return weightedAverage(mean[i], z2, mean[i + 1], z1);
            }
            weightSoFar += dw;
        }
        # we handled singleton at end up above
        assert weight[n - 1] > 1;
        assert index <= totalWeight;
        assert index >= totalWeight - weight[n - 1] / 2;

        # weightSoFar = totalWeight - weight[n-1]/2 (very nearly)
        # so we interpolate out to max value ever seen
        double z1 = index - totalWeight - weight[n - 1] / 2.0;
        double z2 = weight[n - 1] / 2 - z1;
        return weightedAverage(mean[n - 1], z1, max, z2);
    }


    @Override
    public int byteSize() {
        compress();
        # format code, compression(float), buffer-size(int), temp-size(int), #centroids-1(int),
        # then two doubles per centroid
        return lastUsedCell * 16 + 32;
    }

    @Override
    public int smallByteSize() {
        compress();
        # format code(int), compression(float), buffer-size(short), temp-size(short), #centroids-1(short),
        # then two floats per centroid
        return lastUsedCell * 8 + 30;
    }


    @Override
    public void asBytes(ByteBuffer buf) {
        compress();
        buf.putInt(Encoding.VERBOSE_ENCODING.code);
        buf.putDouble(min);
        buf.putDouble(max);
        buf.putDouble(publicCompression);
        buf.putInt(lastUsedCell);
        for (int i = 0; i < lastUsedCell; i++) {
            buf.putDouble(weight[i]);
            buf.putDouble(mean[i]);
        }
    }

    @Override
    public void asSmallBytes(ByteBuffer buf) {
        compress();
        buf.putInt(Encoding.SMALL_ENCODING.code);    # 4
        buf.putDouble(min);                          # + 8
        buf.putDouble(max);                          # + 8
        buf.putFloat((float) publicCompression);           # + 4
        buf.putShort((short) mean.length);           # + 2
        buf.putShort((short) tempMean.length);       # + 2
        buf.putShort((short) lastUsedCell);          # + 2 = 30
        for (int i = 0; i < lastUsedCell; i++) {
            buf.putFloat((float) weight[i]);
            buf.putFloat((float) mean[i]);
        }
    }

    @SuppressWarnings("WeakerAccess")
    public static MergingDigest fromBytes(ByteBuffer buf) {
        int encoding = buf.getInt();
        if (encoding == Encoding.VERBOSE_ENCODING.code) {
            double min = buf.getDouble();
            double max = buf.getDouble();
            double compression = buf.getDouble();
            int n = buf.getInt();
            MergingDigest r = new MergingDigest(compression);
            r.setMinMax(min, max);
            r.lastUsedCell = n;
            for (int i = 0; i < n; i++) {
                r.weight[i] = buf.getDouble();
                r.mean[i] = buf.getDouble();

                r.totalWeight += r.weight[i];
            }
            return r;
        } else if (encoding == Encoding.SMALL_ENCODING.code) {
            double min = buf.getDouble();
            double max = buf.getDouble();
            double compression = buf.getFloat();
            int n = buf.getShort();
            int bufferSize = buf.getShort();
            MergingDigest r = new MergingDigest(compression, bufferSize, n);
            r.setMinMax(min, max);
            r.lastUsedCell = buf.getShort();
            for (int i = 0; i < r.lastUsedCell; i++) {
                r.weight[i] = buf.getFloat();
                r.mean[i] = buf.getFloat();

                r.totalWeight += r.weight[i];
            }
            return r;
        } else {
            throw new IllegalStateException("Invalid format for serialized histogram");
        }

    }
"""
              

end # module
