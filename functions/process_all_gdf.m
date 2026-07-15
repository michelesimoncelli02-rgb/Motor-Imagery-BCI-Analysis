function process_all_gdf(wlength, wshift, pshift, mlength, lap, freq_range, winconv, folder)
%
% process_all_gdf.m
% This function processes all .gdf files in the folder.
%
% Input:
% · wlength: [s] length of the external window
% · wshift: [s] shift of the internal window 
% · pshift: [s] shift of the external window 
% · mlength: length of the window
% · lap: Laplacian filter
% · freq_range: range of frequencies to select
% · winconv: window type
% . folder: name of the folder containing gdf files
%
% Output:
% · .mat file
%
% Processing steps for each .gdf file:
%   - loading of every .gdf found in the folder
%   - artifact removal
%   - Laplacian filter application
%   - PSD computation with the function proc_spectrogram()
%   - meaningful frequencies selection
%   - recomputation of event's information
%   - saving
%

% Importing the required biosig package
%
%

% Add dataset's path
addpath(genpath('./dataset'))

% Extracting the .gdf files in the folder
files = dir(fullfile(['./dataset/' folder],'*.gdf')); 

for k = 1:numel(files)

    fname = files(k).name;
    fprintf('Processing %s...\n', fname);

    % Load GDF
    [data, h] = sload(fname);     % data: samples × channels

    %% Preprocessing - artifact removal

    threshold = 300; % [microV]
    
    if contains(fname, 'offline') % we can apply it only to offline files
        idx_fix = find(h.EVENT.TYP == 786);
        idx_cue = find(h.EVENT.TYP == 771 | h.EVENT.TYP == 773);
        idx_feed = find(h.EVENT.TYP == 781);
        
        for t = 1:length(idx_fix)
            t_start = h.EVENT.POS(idx_fix(t));
            t_end   = h.EVENT.POS(idx_feed(t)) + h.EVENT.DUR(idx_feed(t)) - 1;
            
            if max(abs(data(t_start:t_end, :)), [], 'all') > threshold
                % set to zero all event type of that trial
                h.EVENT.TYP(idx_fix(t)) = 0;
                h.EVENT.TYP(idx_feed(t)) = 0;
                h.EVENT.TYP(idx_cue(t)) = 0;
                
                fprintf('Trial %d completely removed \n', t);
            end
        end
    end
    
    %% PSD computation and event recalculation
    % Apply Laplacian
    data = data(:,1:end-1) * lap;

    % Compute PSD
    [PSD, f] = proc_spectrogram(data, wlength, wshift, pshift, h.SampleRate, mlength); % [windows x frequencies x channels]

    % Select useful frequencies
    [~, idx] = arrayfun(@(ff) min(abs(f - ff)), freq_range); % nearest match, because the exact values may not exist
    PSD_sel = PSD(:, idx, :);
    f_sel = f(idx);

    % Convert event positions to PSD windows
    wPOS = proc_pos2win(h.EVENT.POS, wshift*h.SampleRate, winconv, wlength*h.SampleRate);
    
    % Convert event DUR to window durations
    dur_sec = h.EVENT.DUR / h.SampleRate;   % samples -> seconds
    wDUR = ceil(dur_sec / wshift);       % seconds -> window units

    %% Save
    
    % Parent directory of the current working directory
    parentDir = fullfile(pwd);

    parentFolder = 'processed_dataset';     % common parent
    subjectFolder = folder;                 % per-subject folder

    subDir = fullfile(parentDir, parentFolder, subjectFolder);
    
    if ~exist(subDir, 'dir')
        mkdir(subDir);
    end
    
    outName = [fname(1:end-4)];

    % Build full file path (add .mat extension if needed)
    filePath = fullfile(subDir, [outName '.mat']);
    
    % Save variables
    save(filePath, 'PSD', 'PSD_sel', 'f', 'f_sel', 'wPOS', 'wDUR', 'data', 'h');

end

disp('Done.');
