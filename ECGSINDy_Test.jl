using DataDrivenDiffEq
using ModelingToolkit
using OrdinaryDiffEq
using LinearAlgebra
using Plots
using CSV
using DataFrames
using SignalDecomposition
using DynamicalSystems

#=
I want to add respiration data as a control input and see what that gets 
in terms of a model. Idea here is that breathing adjusts heart rate
=#

println("creating library")

#Create the library
@parameters t
@variables x(t);
u = [x];
polys = [];
for i ∈ 0:5;
    push!(polys, u[1]^i);
end

basis = Basis(polys, u);

println("Setting up the problem")
#Create the model 
sampler = DataSampler(Split(ratio = 0.8),Batcher(n=10, shuffle = true, repeated =true))
λs = exp10.(-10:0.1:-1)
opt = STLSQ(λs)
#problem = ContinuousDataDrivenProblem(ecgDataMatrixDownSampled1', timeDownSampled1, GaussianKernel())
#problem = ContinuousDataDrivenProblem(ecgDataMatrixDownSampled1', timeDownSampled1, TriangularKernel())
#problem = ContinuousDataDrivenProblem(ecgDataMatrixDownSampled1', timeDownSampled1, EpanechnikovKernel())
#problem = ContinuousDataDrivenProblem(ecgDataMatrixDownSampled1', timeDownSampled1, SilvermanKernel())
problem = ContinuousDataDrivenProblem(smoothedEcgData', timeDownSampled1, EpanechnikovKernel())

println("Solving the problem")
#solve the problem
res = solve(problem, basis, opt, sampler = sampler)

#Display model results
println(res)
println(result(res))
println(res.parameters)
