using DataDrivenDiffEq
using ModelingToolkit
using OrdinaryDiffEq
using LinearAlgebra
using Plots
using CSV
using DataFrames
using SignalDecomposition
using DynamicalSystems
using Distributions
using Printf

#=
I want to add respiration data as a control input and see what that gets 
in terms of a model. Idea here is that breathing adjusts heart rate
=#

println("Setting up the problem")

#Create the model 
problemDMD = ContinuousDataDrivenProblem(smoothedEcgData', timesFiltered)

sampler = DataSampler(Split(ratio = 0.7),Batcher(n=10, shuffle = false, repeated =false))
λs = exp10.(-10:0.1:-1)
opt = STLSQ(λs)

resDMD = solve(problemDMD, DMDSVD(), opt, sampler = sampler)

dump(resDMD)