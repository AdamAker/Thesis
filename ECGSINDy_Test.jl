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

println("creating library")


#Create the library
@parameters t
@variables x(t);
u = [x];
polys = [];
for i ∈ 0:1;
    for j∈0:5;
        push!(polys, u[1].^i*exp(-u[1].^(2*j)));
    end
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
problem = ContinuousDataDrivenProblem(smoothedEcgData', timesFiltered)

#=
Plot the problem and it's derivative
=#

println("Plot the problem and its Derivative")

plot(timesFiltered, problem.X', label="Signal", title="ECG (mV) vs time(s)")
xlabel!("Time (s)")
ylabel!("ECG (mV)")
#plot!(timesFiltered, problem.DX', label="Derivative")

println("Solving the problem using Ensemble SINDy")
#solve the problem
maxiter = 100
modelParameters = []
performanceL₂norm = Float64[]
performanceAIC = Float64[]
performanceR² = Float64[]
for i∈1:maxiter
    local res = solve(problem, basis, opt, sampler = sampler)
    resultsDict = metrics(res)
    push!(modelParameters, res.parameters)
    push!(performanceL₂norm, resultsDict[:L₂][1])
    push!(performanceAIC, resultsDict[:AIC][1])
    push!(performanceR², resultsDict[:R²][1])
end

#=
Calculate Statistics of model fit (i.e. average and standard deviation)
and fit it to a normal distribution
=#

#These are the performace metrics fit to the statistics
fitL₂norm = fit(Normal, performanceL₂norm)
μL₂norm = fitL₂norm.μ[1]
strμL₂norm = @sprintf("%.2e",μL₂norm)
σL₂norm = fitL₂norm.σ[1]
strσL₂norm = @sprintf("%.2e",σL₂norm)
metricplt1 = scatter(performanceL₂norm, (1/(fitL₂norm.σ*sqrt(2*π))).*exp.(-(1/(2*fitL₂norm.σ)^2)*(performanceL₂norm.-fitL₂norm.μ).^2))
plot!(x->pdf(fitL₂norm,x))

fitAIC = fit(Normal, performanceAIC)
μAIC = fitAIC.μ[1]
strμAIC = @sprintf("%.2e",μAIC)
σAIC = fitAIC.σ[1]
strσAIC = @sprintf("%.2e",σAIC)
metricplt2 = scatter(performanceAIC, (1/(fitAIC.σ*sqrt(2*π))).*exp.(-(1/(2*fitAIC.σ)^2)*(performanceAIC.-fitAIC.μ).^2))
plot!(x->pdf(fitAIC,x))

fitR² = fit(Normal, performanceR²)
μR² = fitR².μ[1]
strμR² = @sprintf("%.2e",μR²)
σR² = fitR².σ[1]
strσR² = @sprintf("%.2e",σR²)
metricplt3 = scatter(performanceR², (1/(fitR².σ*sqrt(2*π))).*exp.(-(1/(2*fitR².σ)^2)*(performanceR².-fitR².μ).^2))
plot!(x->pdf(fitR²,x))

plot(
    metricplt1, metricplt2, metricplt3, 
    layout = (3,1), label = ["L₂ Error" "fit" "AIC" "fit" "R²" "fit"], 
    title = ["L₂ Error: μ = "*strμL₂norm*", σ = "*strσL₂norm "AIC: μ = "*strμAIC*", σ = "*strσAIC "R²: μ = "*strμR²*", σ = "*strσR²]
    )

#Model parameters fit to the statistics

#Display model results
#=
println(res)
println(result(res))
println(res.parameters)
=#