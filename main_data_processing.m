%% Assignment 1 | 0. Preliminary step: processing of the data
% artifact removal and computation of useful metrics

clear
close all
clc

%% Add useful paths to Matlab path

% Add biosig and eeglab toolboxes
%
%

% Add functions' path
addpath(genpath('./functions'))

%% Data loading 

% Load Laplacian matrix
load laplacian16.mat

% Parameters
wlength = 0.5;         % [s] length of the external window 
pshift  = 0.25;        % [s] shift of the internal window 
wshift  = 0.0625;      % [s] shift of the external window 
mlength = 1;           % length of the window
freq_range = 4:2:48;   % selected frequencies
winconv = 'backward';  % window type

%% Extracting subjects' folders

folders = dir(fullfile('./dataset/a*')); % taking the folders of each subject
num_subj = size(folders,1); % number of subjects

%% Processing of the gdf files for each subject

for j=1:num_subj % for each subject

    % Extract subject folder
    folder_name = folders(j).name;

    process_all_gdf(wlength, wshift, pshift, mlength, lap, freq_range, winconv, folder_name);

end