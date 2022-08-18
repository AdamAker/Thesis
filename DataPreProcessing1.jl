using DataDrivenDiffEq
using ModelingToolkit
using OrdinaryDiffEq
using LinearAlgebra
using Plots
using CSV
using DataFrames
using SignalDecomposition
using DynamicalSystems

savedir = "/Users/adamaker/Desktop/Research/Thesis/ThesisReports/"

println("running")

#Load the Data into a dataframe
#Linux machine source
#ecgDataFrame = CSV.read("/home/adam/Desktop/ThesisData/2020_06_04_T02_U00T_EEG01_EEGAccelTimetable.csv", DataFrame);
#MacOS machine source
ecgDataFrame = CSV.read("/Users/adamaker/Desktop/Research/Thesis/ThesisData/2020_06_04_T02_U00T_EEG01_EEGAccelTimetable.csv", DataFrame)

println("loading data into dataframe")

#Isolate the ECG dataset from the dataframe
#findall(occursin.("ECG",variables))
#variables = names(ecgDataFrame);

#select the data range, these are indicies not times
startind = 50451;
stopind = 50951;
#startind = 50451;
#stopind = 52951;
#startind = 1;
#stopind = length(ecgDataFrame.ECG)
println("formating data")

#format the data into a a "1 x length(ecgData)" Matrix type
ecgDataMatrix = ecgDataFrame.ECG[startind:stopind]
ecgDataMatrix = ecgDataMatrix'
respDataMatrix = ecgDataFrame.Resp[startind:stopind]
respDataMatrix = respDataMatrix'

println("Scaling data")

#Scale the data to be in mV
ecgDataMatrix = ecgDataMatrix.*.001
respDataMatrix = respDataMatrix.*.001


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

#=
Applying a filter to the data which calculates the averge of a window
of windowSize    
=#

println("Filtering Data")


ecgDataMatrixFiltered = Float64[]
timesFiltered = Float64[]
i=1
ε = .25
windowSize = 50
while i+windowSize <= length(ecgDataMatrix)
    windowAvg = 0
    for j∈i:i+windowSize
        windowAvg = windowAvg+ecgDataMatrix[j];
    end
    windowAvg = windowAvg/windowSize;
    if abs(ecgDataMatrix[i]-windowAvg)<ε
        push!(ecgDataMatrixFiltered, windowAvg)
        push!(timesFiltered, time[i])
        global i += 1;
    else
        push!(ecgDataMatrixFiltered, ecgDataMatrix[i])
        push!(timesFiltered, time[i])
        global i += 1;
    end
end

#Smoothing the Downsampled data using cubic splines
splinedData = CubicSpline(vec(ecgDataMatrixFiltered),vec(timesFiltered))

smoothedEcgData = splinedData[1:floor(Int,length(splinedData)/2),1]

#=
println("Downsampling data")


#Down sample to data to try to remove noise
ecgDataMatrixDownSampled1 = Float64[]
timeDownSampled1 = Float64[]
i=0
while 8i+1<length(ecgDataMatrix)
    push!(ecgDataMatrixDownSampled1, ecgDataMatrix[8i+1])
    push!(timeDownSampled1, time[8i+1])
    global i += 1;
end


#Smoothing the Downsampled data using cubic splines
splinedData = CubicSpline(vec(ecgDataMatrixDownSampled1),vec(timeDownSampled1))

smoothedEcgData = splinedData[1:floor(Int,length(splinedData)/2),1]
println("Averaging data")
=#
#=
#Simple pointwise average of the data to try to remove noise
ecgDataMatrixSimpleAvg = Float64[]
timeSimpleAvg = Float64[]
i=1
while i<length(ecgDataMatrix)
    push!(ecgDataMatrixSimpleAvg, (ecgDataMatrix[i]+ecgDataMatrix[i+1])/2)
    push!(timeSimpleAvg, (time[i]+time[i+1])/2)
    global i += 1;
end
=#

println("plotting data")

#plot the data
dataPlot = plot(time, ecgDataMatrix', label="Original", title="ECG (mV) vs time(s)")
xlabel!("Time (s)")
ylabel!("ECG (mV)")
plot!(timesFiltered, ecgDataMatrixFiltered, label="Filtered")
savefig(dataPlot, savedir*"dataPlot.png")