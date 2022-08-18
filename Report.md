---
author: "Adam Aker"
title: "Test Report"
date: "28th July 2022"
---



# Section 1: Data Pre-processing

```julia
f = (x) -> x^2
```

```julia; echo=false
using DataDrivenDiffEq
using ModelingToolkit
using OrdinaryDiffEq
using LinearAlgebra
using Plots
using CSV
using DataFrames
using SignalDecomposition
using DynamicalSystems

startind = 50451;
stopind = 52951;
Samplerate = 500;
ε = .25;
windowSize = 50;

ecgDataFrame = CSV.read("/Users/adamaker/Desktop/Research/Thesis/ThesisData/2020_06_04_T02_U00T_EEG01_EEGAccelTimetable.csv", DataFrame)

ecgDataMatrix = ecgDataFrame.ECG[startind:stopind];
ecgDataMatrix = ecgDataMatrix';
respDataMatrix = ecgDataFrame.Resp[startind:stopind];
respDataMatrix = respDataMatrix';
ecgDataMatrix = ecgDataMatrix.*.001;
respDataMatrix = respDataMatrix.*.001;

NumberofSamples = length(ecgDataMatrix);
TotalSeconds = NumberofSamples/Samplerate;

time = Float64[0.0];
for i∈2:length(ecgDataMatrix)
    push!(time,time[i-1]+1/Samplerate)
end

ecgDataMatrixFiltered = Float64[];
timesFiltered = Float64[];
i=1;
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

splinedData = CubicSpline(vec(ecgDataMatrixFiltered),vec(timesFiltered));

smoothedEcgData = splinedData[1:floor(Int,length(splinedData)/2),1];

DataPlot = plot(time, ecgDataMatrix', label="Original", title="ECG (mV) vs time(s)")
xlabel!("Time (s)")
ylabel!("ECG (mV)")
plot!(timesFiltered, ecgDataMatrixFiltered, label="Filtered")
```

Plot of the Data

```{julia}
DataPlot
```

