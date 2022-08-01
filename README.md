# TDigest

This is a native Julia implementation of the t-digest that is at near parity with the Java
reference implementation. The code is considerably simpler than
the Java version, but the functionality should be essentially identical. The Julia version
has not had any optimization for speed, but appears to already be roughly on par with the Java
version (adding individual samples takes roughly 50-100ns).

The current plan is to prepare this package for inclusion in the `OnlineStats` package.

# How to Use It

First, construct a digest:
```julia
digest = TDigest.MergingDigest(200)
```
The single parameter in this example is the compression factor which bounds the size of the serialized form of the digest. If you only want to estimate extreme quantiles (99-th percentile or more) then a value of 100 is a good choice. If you want to estimate medians as well, 200 is a good value. 

Next, add data one sample at a time
```julia
x = randn()
TDigest.fit!(digest, x)
```
Or in a great mass, all at once
```julia
bunch_of_x = randn(1_000_000)
TDigest.fit!(digest, bunch_of_x)
```

Finally, find how how extreme a value is relative to the observed distribution of your data
```julia
julia> TDigest.cdf(digest, 3)
0.9986316818404886
```
Or find out a threshold such that 99% of the observed data is less than the threshold
```julia
julia> TDigest.quantile(digest, 0.99)
2.3244557159639254
julia> TDigest.quantile.(Ref(digest), [0, 0.5, 0.95, 0.99, 0.999])
5-element Vector{Float64}:
 -4.7915381780277615
  0.00048690715119342943
  1.641346672479079
  2.3244557159639254
  3.0941552936402497
```

You can see the internals of the t-digest by examining the `sketch` component of the digest. Here we plot the cumulative distribution of the data summarized by our digest.
```julia
```



# Things left to do

This package looks nice right now, but it is by no means done. In case you have a hankering to hack, please help with the following:

1. Help me find some remaining issues with maintenance of the t-digest invariants
2. Add appropriate exports
3. Comment on the API design
4. Get some more documentation and expamples written especially regarding anomaly detection
5. Add CI/CD hooks using Github actions
