# Benchmarks for DiscretePersistentHomologyTransform


##
- ran PHT(s1,200) on BoundaryCSV/1.csv
  params: BenchmarkTools.Parameters
    seconds: Float64 5.0
    samples: Int64 10000
    evals: Int64 1
    overhead: Float64 0.0
    gctrial: Bool true
    gcsample: Bool false
    time_tolerance: Float64 0.05
    memory_tolerance: Float64 0.01
  times: Array{Float64}((4,)) [1.485798792e9, 1.496669502e9, 1.515156549e9, 1.516813454e9]
  gctimes: Array{Float64}((4,)) [2.48715413e8, 2.52564369e8, 2.54381421e8, 2.62019311e8]
  memory: Int64 2871305648
  allocs: Int64 2208281

BenchmarkTools.Trial: 
  memory estimate:  2.67 GiB
  allocs estimate:  2208276
  --------------
  minimum time:     1.389 s (16.32% GC)
  median time:      1.633 s (18.52% GC)
  mean time:        1.590 s (18.26% GC)
  maximum time:     1.704 s (18.36% GC)
  --------------
  samples:          4
  evals/sample:     1

- ran PHT(s1,50) on BoundaryCSV/1.csv 
BenchmarkTools.Trial
  params: BenchmarkTools.Parameters
    seconds: Float64 5.0
    samples: Int64 10000
    evals: Int64 1
    overhead: Float64 0.0
    gctrial: Bool true
    gcsample: Bool false
    time_tolerance: Float64 0.05
    memory_tolerance: Float64 0.01
  times: Array{Float64}((13,)) [3.81759842e8, 3.82844893e8, 3.89405846e8, 3.89499383e8, 3.91372933e8, 3.91654021e8, 3.94706806e8, 3.94808695e8, 3.96448971e8, 3.9715521e8, 3.98389384e8, 4.01807528e8, 4.0570488e8]
  gctimes: Array{Float64}((13,)) [7.4618649e7, 6.7723043e7, 7.2615237e7, 6.9714693e7, 7.4674826e7, 7.4510973e7, 7.1063262e7, 7.6989363e7, 6.8776731e7, 7.3407929e7, 7.2627434e7, 7.4844979e7, 7.072253e7]
  memory: Int64 718323168
  allocs: Int64 551831
BenchmarkTools.Trial: 
  memory estimate:  685.05 MiB
  allocs estimate:  551831
  --------------
  minimum time:     381.760 ms (19.55% GC)
  median time:      394.707 ms (18.40% GC)
  mean time:        393.504 ms (18.42% GC)
  maximum time:     405.705 ms (17.43% GC)
  --------------
  samples:          13
  evals/sample:     1
