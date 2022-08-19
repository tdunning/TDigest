
cdf_ref(data::Vector{<:Number}, x::Number) = (sum(data .< x) + sum(data .== x) / 2) / length(data)

function quantile_ref(data::Vector, q::Number)
    sort!(data)
    if q ≤ 0
        return data[1]
    elseif q ≥ 1
        return data[end]
    end
    return data[UInt32(1 + floor(q * length(data)))]
end

"""
Verify that weights are nondecreasing until we pass the median
"""
function check_pattern(weights, median_weight)
    w_so_far = 0
    sign = 1
    max_weight = -1
    for w in w_so_far
        @test sign * (w-max_weight) ≥ 0
        if w > max_weight
            max_weight = w
        end
        if w_so_far ≤ median_weight ≤ w_so_far + w
            @test w == max_weight
            sign = -1
        end
        w_so_far += w
    end
end

