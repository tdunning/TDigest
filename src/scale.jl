struct ScaleFunction
    qcn_to_k::Function
    qn_to_k::Function
    kcn_to_q::Function
    kn_to_q::Function
    max_qcn::Function
    max_qn::Function
    normalizer::Function
end

k_scale(scale::ScaleFunction, q::Number, norm::Number) = scale.qn_to_k(q, norm)
k_scale(scale::ScaleFunction, q::Number, δ::Number, n::Number) = scale.qcn_to_k(q, δ, n)
q_scale(scale::ScaleFunction, k::Number, norm::Number) = scale.kn_to_q(k, norm)
q_scale(scale::ScaleFunction, k::Number, δ::Number, n::Number) = scale.kcn_to_q(k, δ, n)
max_step(scale::ScaleFunction, q::Number, norm::Number) = scale.max_qn(q, norm)
max_step(scale::ScaleFunction, q::Number, δ::Number, n::Number) = scale.max_qcn(q, δ, n)
normalizer(scale::ScaleFunction, δ, n) = scale.normalizer(δ, n)

limit(f, x, x0, x1) = f(max(x0, min(x1, x)))

md"Generates uniform cluster sizes. Used for comparison only."
K_0 = ScaleFunction(
         (q, compression, n) -> compression * q / 2,
         (q, normalizer) -> normalizer * q,
         (k, compression, n) -> 2 * k / compression,
         (k, normalizer) -> k / normalizer,
         (q, compression, n) -> 2 / compression,
         (q, normalizer) -> 1 / normalizer,
         (compression, n) -> compression / 2)

md"""
Generates cluster sizes proportional to sqrt(q*(1-q)). This gives 
constant relative accuracy if accuracy is proportional to squared 
cluster size. It is expected that K_2 and K_3 will give better 
practical results.
"""
K_1 = ScaleFunction(
    (q, compression, n) -> limit(
        qx -> compression * asin(2 * qx - 1) / (2 * pi),
        q, 1e-15, 1e15),
    
    (q, normalizer) -> limit(
        qx -> normalizer * asin(2 * qx - 1),
        q, 1e-15, min(1 - 1e-15, q)),
    
    (k, compression, n) -> limit(
        kx -> (sin(kx * (2 * pi / compression)) + 1) / 2,
        k, -compression/4, min(compression/4, k)),
    
    (k, normalizer) -> limit(
        k -> (sin(k/normalizer) + 1) / 2,
        k, -π/2 * normalizer, π/2 * normalizer),
    
    (q, compression, n) -> limit(
        q -> 2 * sin(pi / compression) * sqrt(q * (1 - q)),
        q, 0, 1),
    
    (q, normalizer) -> limit(
        q -> 2 * sin(0.5 / normalizer) * sqrt(q * (1 - q)),
        q, 0, 1),
    
    (compression, n) -> compression / (2 * pi)
)

md"""
    Generates cluster sizes proportional to sqrt(q*(1-q)) but avoids computation 
    of asin in the critical path by using an approximate version.
    """
K_1_FAST = ScaleFunction(
    (q, compression, n) -> limit(
        qx -> compression * fast_asin(2 * qx - 1) / (2 * pi),
        q, 1e-15, 1e15),
    
    (q, normalizer) -> limit(
        qx -> normalizer * fast_asin(2 * qx - 1),
        q, 1e-15, min(1 - 1e-15, q)),
    
    (k, compression, n) -> limit(
        kx -> (sin(kx * (2 * pi / compression)) + 1) / 2,
        k, -compression/4, min(compression/4, k)),
    
    (k, normalizer) -> limit(
        k -> (sin(k/normalizer) + 1) / 2,
        k, -π/2 * normalizer, π/2 * normalizer),
    
    (q, compression, n) -> limit(
        q -> 2 * sin(pi / compression) * sqrt(q * (1 - q)),
        q, 0, 1),
    
    (q, normalizer) -> limit(
        q -> 2 * sin(0.5 / normalizer) * sqrt(q * (1 - q)),
        q, 0, 1),
    
    (compression, n) -> compression / (2 * pi)
)

md"""
    Generates cluster sizes proportional to q*(1-q). This makes tail error 
    bounds tighter than for K_1. The use of a normalizing function results 
    in a strictly bounded number of clusters no matter how many samples.
    """
K_2 = let Z = (n, compression) -> 4 * log(n / compression) + 24
    ScaleFunction(
        (q, compression, n) -> begin
            if n <= 1
                if q <= 0
                    return -10
                elseif q >= 1
                    return 10
                else
                    return 0
                end
            end
            
            return limit(
                q -> compression * log(q / (1 - q)) / Z(compression, n),
                q, 1e-15, 1 - 1e-15)
        end,
        
        (q, normalizer) -> limit(
            q -> log(q / (1 - q)) * normalizer,
            q, 1e-15, 1 - 1e-15),
        
        (k, compression, n) -> let w = exp(k * Z(compression, n) / compression)
            w / (1 + w)
        end,
        
        (k, normalizer) -> let w = exp(k / normalizer)
            w / (1 + w)
        end,
        
        (q, compression, n) -> Z(compression, n) * q * (1 - q) / compression,
        
        (q, normalizer) -> q * (1 - q) / normalizer,
        
        (compression, n) -> compression / Z(compression, n)
    )
end

md"""
    Generates cluster sizes proportional to min(q, 1-q). This makes tail error 
    bounds tighter than for K_1 or K_2. The use of a normalizing function 
    results in a strictly bounded number of clusters no matter how many samples.
    """
K_3 = let Z = (compression, n) -> 4 * log(n / compression) + 21
    ScaleFunction(
        (q, compression, n) -> limit(
            q -> if q <= 0.5
                return compression * log(2 * q) / Z(compression, n)
            else
                return -compression * log(2 * (1-q)) / Z(compression, n)
            end,
            q, 1e-15, 1 - 1e-15),
        
        (q, normalizer) -> limit(
            q -> if q <= 0.5
                return log(2 * q) * normalizer
            else
                return -log(2 * (1-q)) * normalizer
            end,
            q, 1e-15, 1 - 1e-15),
        
        (k, compression, n) ->
            if k <= 0
                return exp(k * Z(compression, n) / compression) / 2
            else
                return 1 - exp(-k * Z(compression, n) / compression) / 2
            end,
        
        (k, normalizer) ->
            if k <= 0
                return exp(k / normalizer) / 2
            else
                return 1 - exp(-k / normalizer) / 2
            end,
        (q, compression, n) -> Z(compression, n) * min(q, 1 - q) / compression,
        
        (q, normalizer) -> min(q, 1 - q) / normalizer,
        
        (compression, n) -> compression / Z(compression, n)
    )
end
