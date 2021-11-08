% Hardcode the data from the Chylinski paper
% Coding it here because the spreadsheet is too funky to read in
% Is there a better way to do this?
% Save as a .mat file to be able to read it into the parameter fitting
% script easy

% Data from figure 2b
% Cre negative
% Do I need any of this for fitting
% xEmpty = [2 2 2 3 3 3 4 4 4 5 5 5 6 6 6 10 10 10 13 13 13];
% yEmpty = [0.64 0.47 0.56 0.33 0.33 0.43 0.28 0.14 0.38 0.24 0.23 0.39 0.61 0.46 0.59 1.23 1.16 0.72 0.49 1.33 1.26];
% 
% xUnmodified = [2 2 2 3 3 3 4 4 4 5 5 5 6 6 6 10 10 10 13 13 13];
% yUnmodified = [];
% 
% xSwitchOff = [2 2 2 3 3 3 4 4 4 5 5 5 6 6 6 10 10 10 13 13 13];
% ySwitchOff = [];
% 
% xSingleLoxP = [2 2 2 3 3 3 4 4 4 5 5 5 6 6 6 10 10 10 13 13 13];
% ySingleLoxP = [];

% Cre positive
% Which of these do I need for the fitting?
xEmpty = [2 2 2 3 3 3 4 4 4 5 5 5 6 6 6 10 10 10 13 13 13];
yEmpty = [0.12 0.12 0.069 0.13 0.14 0.1 0.078 0.068 0.059 0.082 0.074 0.081 0.21 0.13 0.2 0.37 0.22 0.23 0.19 0.1 0.16];

xUnmodified = [2 2 2 3 3 3 4 4 4 5 5 5 6 6 6 10 10 10 13 13 13];
yUnmodified = [0.5 0.51 0.47 22.6 22.3 22.7 71.8 71.2 71.2 88.4 88.6 88.7 92.5 92.2 92.7 96.5 96.5 97.3 94.6 96.2 97.6];

xSwitchOff = [2 2 2 3 3 3 4 4 4 5 5 5 6 6 6 10 10 10 13 13 13];
ySwitchOff = [0.07 0.12 0.13 0.4 0.42 0.42 0.86 1.06 1.07 1.15 1.53 1.4 1.53 1.45 1.21 1.46 1.48 1.56 1.31 1.59 1.67];

xSingleLoxP = [2 2 2 3 3 3 4 4 4 5 5 5 6 6 6 10 10 10 13 13 13];
ySingleLoxP = [0.18 0.22 0.2 5.46 5.79 5.88 45.4 45.4 44.9 70.9 70.6 71 81.7 81.6 81.4 92.1 92.6 93 92.8 93.2 93.3];

% Check
plot(xEmpty, yEmpty, 'ko')
