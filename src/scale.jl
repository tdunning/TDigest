abstract type ScaleFunction end
md"""
This function defines a subtype of ScaleFunction
and then sets up the type specialized methods for `k_scale`, `q_scale`, 
`max_step` and `normalizer`.
"""
function defineScaleFunction(name; Z = 1,
                             qcn_to_k::Expr, qn_to_k::Expr,
                             kcn_to_q::Expr, kn_to_q::Expr,
                             max_qcn::Expr, max_qn::Expr,
                             normalizer::Expr)
    eval(quote
             struct $name <: ScaleFunction end 
             k_scale(scale::$name, q::Number, norm::Number) = $qn_to_k
             k_scale(scale::$name, q::Number, compression::Number, n::Number) =
                 let Z = $Z
                     $qcn_to_k
                 end
             q_scale(scale::$name, k::Number, norm::Number) = $kn_to_q
             q_scale(scale::$name, k::Number, compression::Number, n::Number) =
                 let Z = $Z
                     $kcn_to_q
                 end
             max_step(scale::$name, q::Number, norm::Number) = $max_qn
             max_step(scale::$name, q::Number, compression::Number, n::Number) =
                 let Z = $Z
                     $max_qcn
                 end
             normalizer(scale::$name, compression::Number, n::Number) =
                 let Z = $Z
                     $normalizer
                 end
         end)
end

limit(f, x, x0, x1) = f(max(x0, min(x1, x)))

md"Generates uniform cluster sizes. Used for comparison only."
defineScaleFunction(:K_0,
                    qcn_to_k = :(compression * q / 2),
                    qn_to_k = :(norm * q),
                    kcn_to_q = :(2 * k / compression),
                    kn_to_q = :(k / norm),
                    max_qcn = :(2 / compression),
                    max_qn = :(1 / norm),
                    normalizer = :(compression / 2))

md"""
Generates cluster sizes proportional to sqrt(q*(1-q)). This gives 
constant relative accuracy if accuracy is proportional to squared 
cluster size. It is expected that K_2 and K_3 will give better 
practical results.
"""
defineScaleFunction(:K_1,
                    qcn_to_k = :(limit(
                        qx -> compression * asin(2 * qx - 1) / (2 * pi),
                        q, 1e-15, 1e15)),
    
                    qn_to_k = :(limit(
                        qx -> norm * asin(2 * qx - 1),
                        q, 1e-15, min(1 - 1e-15, q))),
                    
                    kcn_to_q = :(limit(
                        kx -> (sin(kx * (2 * pi / compression)) + 1) / 2,
                        k, -compression/4, min(compression/4, k))),
                    
                    kn_to_q = :(limit(
                        k -> (sin(k/norm) + 1) / 2,
                        k, -π/2 * norm, π/2 * norm)),
                    
                    max_qcn = :(limit(
                        q -> 2 * sin(pi / compression) * sqrt(q * (1 - q)),
                        q, 0, 1)),
                    
                    max_qn = :(limit(
                        q -> 2 * sin(0.5 / norm) * sqrt(q * (1 - q)),
                        q, 0, 1)),
                    
                    normalizer = :(compression / (2 * pi))
                    )

md"""
Generates cluster sizes proportional to q*(1-q). This makes tail error 
bounds tighter than for K_1. The use of a normalizing function results 
in a strictly bounded number of clusters no matter how many samples.
"""
defineScaleFunction(:K_2, Z = :(4 * log(n / compression) + 24),
                    qcn_to_k = quote
                        if n <= 1
                            if q <= 0
                                return -10
                            elseif q >= 1
                                return 10
                            else
                                return 0
                            end
                        else
                            return limit(
                                q -> compression * log(q / (1 - q)) / Z,
                                q, 1e-15, 1 - 1e-15)
                        end
                    end,
                        
                    qn_to_k = quote
                        limit(
                            q -> log(q / (1 - q)) * norm,
                            q, 1e-15, 1 - 1e-15)
                    end,
                    kcn_to_q = quote
                        let w = exp(k * Z / compression)
                            w / (1 + w)
                        end
                    end,
                    kn_to_q = quote
                        let w = exp(k / norm)
                            w / (1 + w)
                        end
                    end,
                    max_qcn = :(Z * q * (1 - q) / compression),
                    max_qn = :(q * (1 - q) / norm),
                    normalizer = :(compression / Z)
                    )

md"""
Generates cluster sizes proportional to min(q, 1-q). This makes tail error 
bounds tighter than for K_1 or K_2. The use of a normalizing function 
results in a strictly bounded number of clusters no matter how many samples.
"""
defineScaleFunction(:K_3, Z = :(4 * log(n / compression) + 21),
                    qcn_to_k = quote
                        limit(
                            q -> if q <= 0.5
                                return compression * log(2 * q) / Z
                            else
                                return -compression * log(2 * (1-q)) / Z
                            end,
                            q, 1e-15, 1 - 1e-15)
                    end,
                    
                    qn_to_k = quote
                        limit(
                            q -> if q <= 0.5
                                return log(2 * q) * norm
                            else
                                return -log(2 * (1-q)) * norm
                            end,
                            q, 1e-15, 1 - 1e-15)
                    end,
                    kcn_to_q = quote
                        if k <= 0
                            return exp(k * Z / compression) / 2
                        else
                            return 1 - exp(-k * Z / compression) / 2
                        end
                    end,
                    kn_to_q = quote
                        if k <= 0
                            return exp(k / norm) / 2
                        else
                            return 1 - exp(-k / norm) / 2
                        end
                    end,
                    max_qcn = :(Z * min(q, 1 - q) / compression),
                    max_qn = :(min(q, 1 - q) / norm),
                    normalizer = :(compression / Z)
                    )

