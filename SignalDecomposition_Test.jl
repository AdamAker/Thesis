using SignalDecomposition
using Plots
using CSV
using DataFrames

#Load the Data into a dataframe
ecgDataFrame = CSV.read("/home/adam/Desktop/ThesisData/2020_06_04_T02_U00T_EEG01_EEGAccelTimetable.csv", DataFrame);

#Isolate the ECG dataset from the dataframe
variables = names(ecgDataFrame);

#select the data range, these are indicies not times
startind = 51001;
stopind = 51500;

#format the data into a a "1 x length(ecgData)" Matrix type
ecgData = ecgDataFrame.ECG[startind:stopind];
ecgData = ecgData[:,:];
ecgDataMatrix = ecgData';

#decompose the data using SignalDecomposition.jl
k=1 #past delay
ℓ=1 #forward delay
w=5 #Theiler window
τ=10 #delay time
ε=0.01 #radius of neighborhood in embedded space

x,r = decompose(ecgDataMatrix ,ExtremelySimpleNL(k,ℓ,w,τ,ε))
