using DataDrivenDiffEq
using ModelingToolkit
using OrdinaryDiffEq
using LinearAlgebra
using Plots
using CSV
using DataFrames
using SignalDecomposition

#Load the Data into a dataframe
ecgDataFrame = CSV.read("/home/adam/Desktop/ThesisData/2020_06_04_T02_U00T_EEG01_EEGAccelTimetable.csv", DataFrame);

variables = names(ecgDataFrame);

findall(occursin.("ECG",variables));
startind = 51001;
stopind = 51500;
#format the data into a vector
ecgData = ecgDataFrame.ECG[startind:stopind];

DecgData = ecgData;

#approximate derivative by choppin off last data point of "DecgData" and first of "ecgData"
popfirst!(ecgData);
pop!(DecgData);

#Create the library
@variables x;
u = [x];
polys = [];
for i ∈ 0:3;
    push!(polys, u[1]^i);
end

basis = Basis(polys, u);

#Now setup the optimiser
maxiter = 100;
c_error = 1e-3;

#Create the model 
ecgDataMatrix = ecgData[:,:]
ecgDataMatrix = ecgDataMatrix'
DecgDataMatrix = DecgData[:,:]
DecgDataMatrix = DecgDataMatrix'
problem = DirectDataDrivenProblem(ecgDataMatrix,DecgDataMatrix, name = :Test)

res = solve(problem,basis,STLSQ())
println(res)
println(result(res))
println(res.parameters)