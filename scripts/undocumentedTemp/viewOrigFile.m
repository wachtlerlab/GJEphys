% clear workspace
clear;
% add path to CED code
if isempty(getenv('CEDS64ML'))
    setenv('CEDS64ML', '/home/ajay/Documents/MATLAB/matson/CEDS64ML');
end
cedpath = getenv('CEDS64ML');
addpath(cedpath);
% load ceds64int.dll
CEDS64LoadLib( cedpath );
% Open a file
fileName = uigetfile('*.smr');
fhand1 = CEDS64Open( fileName );
if (fhand1 <= 0); unloadlibrary ceds64int; return; end
timebase = CEDS64TimeBase( fhand1 );
sampRate = 1 / timebase;
% Read signals 
voltageSignal = CED64ReadWaveF(fhand1, 1);
vibrationSignal = CED64ReadWaveF(fhand1, 3);
% Read calibrations
voltageCalibration = CED64ReadWaveF(fhand1, 2)
vibrationCalibration = CED64ReadWaveF(fhand1, 4)
% start stop
startStop = CED64ReadWaveF(fhand1, 5)

figure(); hold on;

for ind = 1:ceil(CED64MaxTime(fhand1))
    clf()
    t = (ind - 1) + (0: 7 :sampRate) * timebase;
    subplot(311)
    voltSig = voltageSignal((ind - 1) * sampRate: 7: ind * sampRate);
    plot(t, voltSig, 'b-')
    subplot(312)
    vibSig = vibrationSignal((ind - 1) * sampRate: 7: ind * sampRate);
    plot(t, vibSig, 'g-')
    ch = input('Continue(y/n):', 's');
    if ch == 'n'
        break
    end
end