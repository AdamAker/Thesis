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
@parameters t;
@variables x[1:2];
x = collect(x);
h = Num[ exp.(-x[1].^(2)) ;exp.(-x[2].^(2));polynomial_basis(x,5)];


basis = Basis(h,x);

println("Setting up the problem")
#Create the model 
sampler = DataSampler(Split(ratio = 0.8),Batcher(n=5, shuffle = true, repeated =true))
λs = exp10.(-10:0.1:-1)
opt = STLSQ(λs)
#problem = ContinuousDataDrivenProblem(ecgDataMatrixDownSampled1', timeDownSampled1, GaussianKernel())
#problem = ContinuousDataDrivenProblem(ecgDataMatrixDownSampled1', timeDownSampled1, TriangularKernel())
#problem = ContinuousDataDrivenProblem(ecgDataMatrixDownSampled1', timeDownSampled1, EpanechnikovKernel())
#problem = ContinuousDataDrivenProblem(ecgDataMatrixDownSampled1', timeDownSampled1, SilvermanKernel())
problem = ContinuousDataDrivenProblem(smoothedEcgData', timesFiltered)
problem2D = ContinuousDataDrivenProblem([problem.X;problem.DX], timesFiltered)

phasePlot = plot(problem.X',problem.DX', label = "Attractor",title = "Phase Diagram")
xlabel!("Time Series")
ylabel!("Derivative")
savefig(phasePlot, savedir*"phasePlot.png")

#=
Plot the problem and it's derivative
=#

println("Plot the problem, its Derivative, and 2nd Derivative")

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
    local res = solve(problem2D, basis, opt, sampler = sampler)
    
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

performancePlots = plot(
    metricplt1, metricplt2, metricplt3, 
    layout = (3,1), label = ["L₂ Error" "fit" "AIC" "fit" "R²" "fit"], 
    title = ["L₂ Error: μ = "*strμL₂norm*", σ = "*strσL₂norm "AIC: μ = "*strμAIC*", σ = "*strσAIC "R²: μ = "*strμR²*", σ = "*strσR²]
    )

savefig(performancePlots, savedir*"performancePlots.png")

#Model parameters fit to the statistics

#=First seach modelParameters to find the model with the most parameters, then 
we'll make a matrix of size "maxiter x maxparamsize" so that each column of
the matrix will be a a series of values for a parameter.=#

maxparamsize = 0
for i∈1:maxiter
    paramsize = length(modelParameters[i]) 
    if paramsize > maxparamsize
        global maxparamsize = paramsize
    end
end

modelParamMatrix = zeros(maxiter,maxparamsize)

for i∈1:maxparamsize
    for j∈1:maxiter
        try
            modelParamMatrix[j,i] = modelParameters[j][i]
        catch
            modelParamMatrix[j,i] = 0
        end
    end
end

paramPlots = []
for i∈1:maxparamsize
    fitparam = fit(Normal, modelParamMatrix[:,i])
    μparam = fitparam.μ[1]
    strμparam = @sprintf("%.2e",μparam)
    σparam = fitparam.σ[1]
    strσparam = @sprintf("%.2e",σparam)
    pplt = scatter(modelParamMatrix[:,i], (1/(fitparam.σ*sqrt(2*π))).*exp.(-(1/(2*fitparam.σ)^2)*(modelParamMatrix[:,i].-fitparam.μ).^2))
    plot!(x->pdf(fitparam,x))
    title!("Parameter "*string(i)*": μ = "*strμparam*" σ = "*strσparam)
    ylabel!("Normalized PDF")
    xlabel!("Parameter "*string(i)*" value")
    push!(paramPlots, pplt)
    savefig(paramPlots[i], savedir*"parametersPlot"*string(i)*".png")
end



#Display model results
#=
println(res)
println(result(res))
println(res.parameters)
=#