using DataDrivenDiffEq
using ModelingToolkit
using OrdinaryDiffEq
using LinearAlgebra
using Plots
using CSV
using DataFrames
using SignalDecomposition

println("running")

#Load the Data into a dataframe
ecgDataFrame = CSV.read("/home/adam/Desktop/ThesisData/2020_06_04_T02_U00T_EEG01_EEGAccelTimetable.csv", DataFrame);

println("loading data into dataframe")

#Isolate the ECG dataset from the dataframe
#findall(occursin.("ECG",variables))
#variables = names(ecgDataFrame);

#select the data range, these are indicies not times
startind = 50001;
stopind = 55000;

println("formating data")

#format the data into a a "1 x length(ecgData)" Matrix type
ecgDataMatrix = ecgDataFrame.ECG[startind:stopind]
ecgDataMatrix = ecgDataMatrix'

println("Scaling data")

#Scale the data to be in mV
ecgDataMatrix = ecgDataMatrix.*.001

println("Creating the vector containing sample times")

#create a matrix of times
NumberofSamples = length(ecgDataMatrix); #Total number of samples in the ecgData
Samplerate = 500; #samples taken per second
TotalSeconds = NumberofSamples/Samplerate; #Length in seconds of time series

#setup the time matrix, same dimensions as ecgDataMatrix
time = Float64[0.0];
for i∈2:length(ecgDataMatrix)
    push!(time,time[i-1]+1/Samplerate)
end

smoothedData = QuadraticSpline(vec(ecgDataMatrix),vec(time))

smoothedData.
