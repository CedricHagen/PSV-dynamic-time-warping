%% DTWPSV - Dynamic Time Warping of Paleomagnetic Secular Variation Data
% Example set-up
clear all

% Select routines to run
timewarp = 1; % Warp PSV candidate record to target record (1 = run; 0 = skip)
significancetest = 1; % Test significance of dynamic time warping (1 = run; 0 = skip)

% Name your target curve (name_T) and candidate curve (name_C)
% Note: should be same name as folder in PSV_data directory
name_T = '2269'; % Target Curve Name
name_C = '2322'; % Candidate Curve Name

% Import your data (note: 3 columns--(1)depth, (2)dec, (3)inc)
psv_T = xlsread('PSV_Data/2269/NNA_2269_4dtw.xlsx');
psv_C = xlsread('PSV_Data/2322/NNA_2322_4dtw.xlsx');
% psv_C = xlsread('PSV_Data/U1305/NNA_U1305_4dtw.xlsx');

% match vector lengths
% nd = [min(psv_C(:, 1)):(max(psv_C(:, 1))-min(psv_C(:, 1)))/length(psv_T(:, 1)):max(psv_C(:, 1))]';
% [x y z] = xyz(psv_C(:, 2), psv_C(:, 3), psv_C(:, 2)*0+1);
% tnd = interp1(psv_C(:, 1), [x y z], nd);
% [tdec tinc ~] = decinc(tnd(:, 1), tnd(:, 2), tnd(:, 3));
% psv_Co = psv_C;
% psv_C = [nd tdec tinc];

% Select g and edge values for dynamic timewarping and significance test parameters
% Note significance test is time consuming and may be best to run one
% g-value at a time
g = [1]; % Suggested [0.98:0.01:1.01]
edge = [.1]; % Suggested [0.01:0.01:0.15]
nitt = 2; % Number of Sythetic Timeseries Suggested 10000
speed = 1; % 1 = Drops unnessesary calculations on Target curve to increase speed
gensyn = 1; % 1 = generate synthetics (only needs to be done once)

%--- Do Not Edit Below ----
% Run routines
if timewarp == 1
    DTWPSV_main(psv_T, psv_C, name_T, name_C, g, edge);
end

if significancetest == 1
    DTWPSV_Stat(psv_T, psv_C, name_T, name_C, g, edge, nitt, speed, gensyn);
end