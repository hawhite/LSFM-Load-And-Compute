%%  Author and License Header
%   Author: Hamilton White, Postdoctoral Researcher, Brigham and Women's
%   Hospital, Harvard Medical School, and Boston University
%   Copyright (c) 2024 - present: Hamilton White
%
%   License: CC BY-NC-ND
%   (https://creativecommons.org/licenses/by-nc-nd/4.0/)
%   
%   Distributed via: 
%   Permanent DOI: 
%
%   You need:
%       -A function to get the variables out of the structures
%       -To automate the computation of each element of the structure per file
%       -To normalize fonts, size, shape of figures, elements in each panel to
%       the same size/shape/color
%       -A function to differentiate noisy neural data. Here we use a
%       custom implementation ("diffall" command) of TVRegDiff from Rick Chartrand: https://github.com/JeffreyEarly/GLNumericalModelingKit/tree/master/Matlab
%           --> Loops over each row of neural data and computed the energy
%           regularized differential
%       -A function to normalize the data (by row) to the lowest 1% of
%       datapoints: herein "NormalizeFlnrm"
%
%   Version History:
%       - v1.0 [20260201]: Original publically released version
%       -



clc;clear;close all;
Comps = [];

fileList = dir("*.mat");
for i = 1:size(fileList,1)
    tmp = load(fileList(i).name); % Load data
    Comps(i).("fn") = fileList(1).name; % Capture filename
    Comps(i).("input") = tmp.activity.GCaMP; % Add raw data to Comps
    Comps(i).("wholetrace_input") = NormalizeFlnrm(Comps(i).input); % Normalize data
    %%% Insert code to reject stimulation flashes (will saturate analyzed data), only if they exist. Do not just use a
    %%% static threshold --> not sufficient overall
    Comps(i).("input_diff") = diffall(Comps(i).input); % Diff of raw input
    Comps(i).("wholetrace_diff") = diffall(Comps(i).wholetrace_input); % Diff of Normalized data, raw data for "diffs" in my own code
    
    % Correlations
    Comps(i).("corrs") = corr(Comps(i).input_diff','rows','complete'); % Compute correlations, reject nan, size = [neurons x neurons]
    tril_mat = tril(Comps(i).corrs,-1); tril_mat(tril_mat==0)=nan;
    Comps(i).("corrs_tril") = tril_mat(:);
end

%% Plots of first set (row) of data in Comps set
figure(1);clf;
imagesc(Comps(1).wholetrace_input);

figure(2);clf;
edges = linspace(-1,1,20);
histogram(Comps(1).input_diff,'Normalization','probability','BinEdges',edges,DisplayStyle='stairs');
