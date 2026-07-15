%% Assignment 1 | 1. Grand average analysis 
% on the whole population and on representative subjects

clear
close all
clc

%% Add utils path

% add biosig and eeglab toolboxes
%
%

% Add functions' path
addpath(genpath('./functions'))

% Add dataset's path
addpath(genpath('./processed_dataset'))

% Load channels' matrix and denomination
load('chanlocs16.mat');
ch = {"Fz", "FC3", "FC1", "FCz", "FC2", "FC4", "C3", ...
    "C1", "Cz", "C2", "C4", "CP3", "CP1", "CPz", "CP2", "CP4"};

%% Extract all the subjects' subfolders
% Extracting the .gdf files in the folder
folders = dir(fullfile('./processed_dataset','a*')); % taking subjects' folders
num_subj = size(folders, 1); % number of subjects

%% ERD/ERS computation for each subject

population = struct;

for j=1:num_subj % for each subject

    % Extract subject folder
    folder_name = folders(j).name;

    % Extract all gdf files in the folder
    dataDir = fullfile(pwd, 'processed_dataset', folder_name);
    files = dir(fullfile([dataDir, '\*offline*.mat']));

    % concatenating all GDF files of the j-th subject
    s=[]; EVENT=struct; EVENT.TYP=[]; EVENT.DUR=[]; EVENT.POS=[]; PSD=[]; EVENT.wDUR=[]; EVENT.wPOS=[];

    for i=1:size(files,1) % for each file of the j-th subject

        % load the processed file
        file = load(fullfile(files(i).folder, files(i).name));

        % Load useful data
        s_i = file.data;             % [samples × channels]
        h_i = file.h;   
        
        EVENT.TYP = [EVENT.TYP; h_i.EVENT.TYP];
        EVENT.DUR = [EVENT.DUR; h_i.EVENT.DUR];
        EVENT.POS = [EVENT.POS; size(s,1)+h_i.EVENT.POS];
        EVENT.wDUR = [EVENT.wDUR; file.wDUR];
        EVENT.wPOS = [EVENT.wPOS; size(PSD,1)+file.wPOS];

        s = [s; s_i];
        PSD = [PSD; file.PSD_sel];    % [window x frequencies x channels]
        
    end
    
    % sampling rate
    Sr = h_i.SampleRate;

    % extracting the sizes
    [N, M] = size(s);  % [# samples, # channels]

    %% Filtering 
    
    % Filtering but without applying the log-transform
    % since the data will be normalized with the ERD/ERS
    
    % Butterworth filters
    [b_mi, a_mi] = butter(5, 2*[8 12]/Sr);
    [b_beta, a_beta] = butter(5, 2*[18 22]/Sr);
    % checking their stability
    % fvtool(b_mi, a_mi, 'Fs', Sr)
    % fvtool(b_beta, a_beta, 'Fs', Sr)
    % applying the filters
    s_mi = filtfilt(b_mi, a_mi, s);
    s_beta = filtfilt(b_beta, a_beta, s);
    
    % squaring
    s_mi2 = (s_mi).^2;
    s_beta2 = (s_beta).^2;
    
    % applying a moving average
    s_mi2_avg = movavg(s_mi2, Sr);
    s_beta2_avg = movavg(s_beta2, Sr);
    
    % final signals after the pre-processing steps
    s_filt_mi = s_mi2_avg;
    s_filt_beta = s_beta2_avg;
    
    %% Trial extraction
    
    % 771 = both feet
    % 773 = both hands
    % extracting the trial labels
    Ck = EVENT.TYP(or(EVENT.TYP==771, EVENT.TYP==773));
    
    % trials' starting points
    strt = EVENT.POS(EVENT.TYP==786); % starting from the fixation cross event
    w_strt = EVENT.wPOS(EVENT.TYP==786); % starting from the fixation cross event
    nT = length(strt);  % # trials
    % trials' end points
    endt = EVENT.POS(EVENT.TYP==781) + EVENT.DUR(EVENT.TYP==781) - 1;
    w_endt = EVENT.wPOS(EVENT.TYP==781) + EVENT.wDUR(EVENT.TYP==781) - 1;
    
    % selecting the length of the trials as the minumum one
    l_tr = min(endt-strt);
    w_l_tr = min(w_endt-w_strt);
    
    % initializing the matrices
    [~, nF, nC] = size(PSD);
    Activity = zeros(w_l_tr, nF, nC, nT); % [windows (in a trials) x frequencies x channels x trials]
    trials_mi = zeros(l_tr, M, nT);
    trials_beta = zeros(l_tr, M, nT);
    
    for i = 1:nT
        % signal indexes
        ind = (strt(i):strt(i)+l_tr-1);
        w_ind = (w_strt(i):w_strt(i)+w_l_tr-1);
        % filling the matrices
        Activity(:,:,:,i) = PSD(w_ind,:,:);
        trials_mi(:,:,i) = s_filt_mi(ind,:);
        trials_beta(:,:,i) = s_filt_beta(ind,:);
    end
    
    % length of the fixation period (computed as l_tr)
    l_f = min(EVENT.DUR(EVENT.TYP==786)) - 1;
    w_l_f = min(EVENT.wDUR(EVENT.TYP==786)) - 1;
    
    % length of the continuous feedback period (computed as l_tr)
    l_cf = min(EVENT.DUR(EVENT.TYP==781)) - 1;

    %% ERD/ERS computation

    % we consider the fixation period as the reference period
    % and the continuous feedback period as activity period

    % whole band
    Reference = repmat(mean(Activity(1:w_l_f,:,:,:)), [size(Activity,1), 1,1,1]);
    ERD = log(Activity./Reference); % [windows x frequencies x channels x trials]

    % mi band
    % We compute the mean of the fixation data for every column (channel)
    % and replicate it to match the size of trials_mi for ERD/ERS computation.
    Reference_mi = repmat(mean(trials_mi(1:l_f,:,:)), [size(trials_mi, 1) 1 1]);
    ERD_mi = 100 * (trials_mi - Reference_mi)./ Reference_mi;
        
    % beta band
    % We compute the mean of the fixation data for every column (channel)
    % and replicate it to match the size of trials_beta for ERD/ERS computation.
    Reference_beta = repmat(mean(trials_beta(1:l_f,:,:)), [size(trials_beta, 1) 1 1]);
    ERD_beta = 100 * (trials_beta - Reference_beta)./ Reference_beta;

    % saving subject's matrices in population's structure
    population(j).ERD_771 = ERD(:,:,:,Ck==771);
    population(j).ERD_773 = ERD(:,:,:,Ck==773);
    population(j).ERD_mi = ERD_mi;
    population(j).ERD_beta = ERD_beta;

    %% ERD Ref/Act - mi band
    
    % mi band - useful for spatial visualization
    population(j).ERD_Ref_771_mi = mean(mean(ERD_mi(1:l_f, :, Ck == 771), 3), 1);
    population(j).ERD_Act_771_mi = mean(mean(ERD_mi(end-l_cf:end, :, Ck == 771), 3), 1);
    population(j).ERD_Ref_773_mi = mean(mean(ERD_mi(1:l_f, :, Ck == 773), 3), 1);
    population(j).ERD_Act_773_mi = mean(mean(ERD_mi(end-l_cf:end, :, Ck == 773), 3), 1);
    

    %% ERD Ref/Act - beta band
    
    % beta band - useful for spatial visualization
    population(j).ERD_Ref_771_beta = mean(mean(ERD_beta(1:l_f, :, Ck == 771), 3), 1);
    population(j).ERD_Act_771_beta = mean(mean(ERD_beta(end-l_cf:end, :, Ck == 771), 3), 1);
    population(j).ERD_Ref_773_beta = mean(mean(ERD_beta(1:l_f, :, Ck == 773), 3), 1);
    population(j).ERD_Act_773_beta = mean(mean(ERD_beta(end-l_cf:end, :, Ck == 773), 3), 1);

    %% averaging the ERD/ERS across trials for both classes
    
    % used for temporal visualization
    % average values for both classes
    population(j).avg_ERD_mi_771 = mean(population(j).ERD_mi(:,:,Ck==771), 3);
    population(j).avg_ERD_beta_771 = mean(population(j).ERD_beta(:,:,Ck==771), 3);
    population(j).avg_ERD_mi_773 = mean(population(j).ERD_mi(:,:,Ck==773), 3);
    population(j).avg_ERD_beta_773 = mean(population(j).ERD_beta(:,:,Ck==773), 3);
    % standard error values for both classes
    population(j).se_ERD_mi_771 = std(population(j).ERD_mi(:,:,Ck==771), 0, 3)/sqrt(nT/2);
    population(j).se_ERD_beta_771 = std(population(j).ERD_beta(:,:,Ck==771), 0, 3)/sqrt(nT/2);
    population(j).se_ERD_mi_773 = std(population(j).ERD_mi(:,:,Ck==773), 0, 3)/sqrt(nT/2);
    population(j).se_ERD_beta_773 = std(population(j).ERD_beta(:,:,Ck==773), 0, 3)/sqrt(nT/2);

    % time axis
    population(j).time = 0:1/Sr:(l_tr-1)/Sr;  % [0:end] sec

end

%% Spatial visualization - topoplot (for each subject)

for j=1:length(population)
    % mi band 
    
    figure
    sgtitle(['Subject ' num2str(j)], 'Color', 'k')
    
    subplot(221)
    t = title('ERD/ERS during activity | \mu band | both feet');
    t.Color = 'k';
    h = [];
    topoplot(squeeze(population(j).ERD_Act_771_mi), chanlocs16);
    h = [h, gca];
    set(h, 'clim', [-50 100]);
    cb = colorbar;
    cb.Color = 'k';
    
    subplot(222)
    t = title('ERD/ERS during activity | \mu band | both hands');
    t.Color = 'k';
    h = [];
    topoplot(squeeze(population(j).ERD_Act_773_mi), chanlocs16);
    h = [h, gca];
    set(h, 'clim', [-50 100]);
    cb = colorbar;
    cb.Color = 'k';


    % beta band
    
    subplot(223)
    t = title('ERD/ERS during activity | \beta band | both feet');
    t.Color = 'k';
    h = [];
    topoplot(squeeze(population(j).ERD_Act_771_beta), chanlocs16);
    h = [h, gca];
    set(h, 'clim', [-50 100]);
    cb = colorbar;
    cb.Color = 'k';
    
    subplot(224)
    t = title('ERD/ERS during activity | \beta band | both hands');
    t.Color = 'k';
    h = [];
    topoplot(squeeze(population(j).ERD_Act_773_beta), chanlocs16);
    h = [h, gca];
    set(h, 'clim', [-50 100]);
    cb = colorbar;
    cb.Color = 'k';

end

%% Temporal visualization (for each subject)

% selecting a meaningful channel
channel = 9; 

for j=1:length(population) 

    time = population(j).time;
    
    % plot
    figure
    sgtitle(['Subject ' num2str(j) ' | Channel ' ch{channel}], 'Color', 'w')
    
    subplot(121)
    hold on
    plot(time, population(j).avg_ERD_mi_771(:,channel), 'g')
    plot(time, population(j).avg_ERD_mi_773(:,channel), 'r')
    plot(time, population(j).avg_ERD_mi_771(:,channel)-population(j).se_ERD_mi_771(:,channel), 'g:', 'Linewidth', 0.2)
    plot(time, population(j).avg_ERD_mi_771(:,channel)+population(j).se_ERD_mi_771(:,channel), 'g:', 'Linewidth', 0.2)
    plot(time, population(j).avg_ERD_mi_773(:,channel)-population(j).se_ERD_mi_773(:,channel), 'r:', 'Linewidth', 0.2)
    plot(time, population(j).avg_ERD_mi_773(:,channel)+population(j).se_ERD_mi_773(:,channel), 'r:', 'Linewidth', 0.2)
    xline(3) % 3 sec
    xline(4) % 4 sec
    xlabel('Time [s]')
    ylabel('ERD/ERD [%]')
    title('ERD in \mu band | Mean +/- SE')
    legend('both feet', 'both hands')
    
    subplot(122)
    hold on
    plot(time, population(j).avg_ERD_beta_771(:,channel), 'g')
    plot(time, population(j).avg_ERD_beta_773(:,channel), 'r')
    plot(time, population(j).avg_ERD_beta_771(:,channel)-population(j).se_ERD_beta_771(:,channel), 'g:', 'Linewidth', 0.2)
    plot(time, population(j).avg_ERD_beta_771(:,channel)+population(j).se_ERD_beta_771(:,channel), 'g:', 'Linewidth', 0.2)
    plot(time, population(j).avg_ERD_beta_773(:,channel)-population(j).se_ERD_beta_773(:,channel), 'r:', 'Linewidth', 0.2)
    plot(time, population(j).avg_ERD_beta_773(:,channel)+population(j).se_ERD_beta_773(:,channel), 'r:', 'Linewidth', 0.2)
    xline(3) % 3 sec
    xline(4) % 4 sec
    xlabel('Time [s]')
    ylabel('ERD/ERD [%]')
    title('ERD in \beta band | Mean +/- SE')
    legend('both feet', 'both hands')

end

%% Grand average analysis
% grand average metrics' computation

GA_ERD_Ref_771_mi = mean(cat(1,population.ERD_Ref_771_mi),1);
GA_ERD_Act_771_mi = mean(cat(1,population.ERD_Act_771_mi),1);
GA_ERD_Ref_773_mi = mean(cat(1,population.ERD_Ref_773_mi),1);
GA_ERD_Act_773_mi = mean(cat(1,population.ERD_Act_773_mi),1);

GA_ERD_Ref_771_beta = mean(cat(1,population.ERD_Ref_771_beta),1);
GA_ERD_Act_771_beta = mean(cat(1,population.ERD_Act_771_beta),1);
GA_ERD_Ref_773_beta = mean(cat(1,population.ERD_Ref_773_beta),1);
GA_ERD_Act_773_beta = mean(cat(1,population.ERD_Act_773_beta),1);

%% Spatial visualization - topoplot (grand average) 

% mi band
figure
sgtitle('Grand Average | Topoplots | \mu band', 'Color', 'k')
subplot(221)
t = title('ERD/ERS during fixation | \mu band | both feet');
t.Color = 'k';
h = [];
topoplot(squeeze(GA_ERD_Ref_771_mi), chanlocs16);
h = [h, gca];
set(h, 'clim', [-50 100]);
cb = colorbar;
cb.Color = 'k';

subplot(222)
t = title('ERD/ERS during activity | \mu band | both feet');
t.Color = 'k';
h = [];
topoplot(squeeze(GA_ERD_Act_771_mi), chanlocs16);
h = [h, gca];
set(h, 'clim', [-50 100]);
cb = colorbar;
cb.Color = 'k';

subplot(223)
t = title('ERD/ERS during fixation | \mu band | both hands');
t.Color = 'k';
h = [];
topoplot(squeeze(GA_ERD_Ref_773_mi), chanlocs16);
h = [h, gca];
set(h, 'clim', [-50 100]);
cb = colorbar;
cb.Color = 'k';

subplot(224)
t = title('ERD/ERS during activity | \mu band | both hands');
t.Color = 'k';
h = [];
topoplot(squeeze(GA_ERD_Act_773_mi), chanlocs16);
h = [h, gca];
set(h, 'clim', [-50 100]);
cb = colorbar;
cb.Color = 'k';

% beta band
figure
sgtitle('Grand Average | Topoplots | \beta band', 'Color', 'k')
subplot(221)
t = title('ERD/ERS during fixation | \beta band | both feet');
t.Color = 'k';
h = [];
topoplot(squeeze(GA_ERD_Ref_771_beta), chanlocs16);
h = [h, gca];
set(h, 'clim', [-50 100]);
cb = colorbar;
cb.Color = 'k';

subplot(222)
t = title('ERD/ERS during activity | \beta band | both feet');
t.Color = 'k';
h = [];
topoplot(squeeze(GA_ERD_Act_771_beta), chanlocs16);
h = [h, gca];
set(h, 'clim', [-50 100]);
cb = colorbar;
cb.Color = 'k';

subplot(223)
t = title('ERD/ERS during fixation | \beta band | both hands');
t.Color = 'k';
h = [];
topoplot(squeeze(GA_ERD_Ref_773_beta), chanlocs16);
h = [h, gca];
set(h, 'clim', [-50 100]);
cb = colorbar;
cb.Color = 'k';

subplot(224)
t = title('ERD/ERS during activity | \beta band | both hands');
t.Color = 'k';
h = [];
topoplot(squeeze(GA_ERD_Act_773_beta), chanlocs16);
h = [h, gca];
set(h, 'clim', [-50 100]);
cb = colorbar;
cb.Color = 'k';


%% Grand average metrics' re-computation for selected subjects
% In this section one can try different combinations of subjects to be
% removed from the analysis, to check how the metrics evolve.
% One can try to remove outliers to achieve better outcomes, or to
% focus on outliers to explore possible reasons behind poor performances.

new_population = population;
new_population([1,2]) = []; % to be personalized

new_GA_ERD_Ref_771_mi = mean(cat(1,new_population.ERD_Ref_771_mi),1);
new_GA_ERD_Act_771_mi = mean(cat(1,new_population.ERD_Act_771_mi),1);
new_GA_ERD_Ref_773_mi = mean(cat(1,new_population.ERD_Ref_773_mi),1);
new_GA_ERD_Act_773_mi = mean(cat(1,new_population.ERD_Act_773_mi),1);

new_GA_ERD_Ref_771_beta = mean(cat(1,new_population.ERD_Ref_771_beta),1);
new_GA_ERD_Act_771_beta = mean(cat(1,new_population.ERD_Act_771_beta),1);
new_GA_ERD_Ref_773_beta = mean(cat(1,new_population.ERD_Ref_773_beta),1);
new_GA_ERD_Act_773_beta = mean(cat(1,new_population.ERD_Act_773_beta),1);

%% Spatial visualization - topoplot (grand average) 

% mi band
figure
sgtitle('Grand Average | Topoplots | \mu band', 'Color', 'k')
subplot(221)
t = title('ERD/ERS during fixation | \mu band | both feet');
t.Color = 'k';
h = [];
topoplot(squeeze(new_GA_ERD_Ref_771_mi), chanlocs16);
h = [h, gca];
set(h, 'clim', [-50 100]);
cb = colorbar;
cb.Color = 'k';

subplot(222)
t = title('ERD/ERS during activity | \mu band | both feet');
t.Color = 'k';
h = [];
topoplot(squeeze(new_GA_ERD_Act_771_mi), chanlocs16);
h = [h, gca];
set(h, 'clim', [-50 100]);
cb = colorbar;
cb.Color = 'k';

subplot(223)
t = title('ERD/ERS during fixation | \mu band | both hands');
t.Color = 'k';
h = [];
topoplot(squeeze(new_GA_ERD_Ref_773_mi), chanlocs16);
h = [h, gca];
set(h, 'clim', [-50 100]);
cb = colorbar;
cb.Color = 'k';

subplot(224)
t = title('ERD/ERS during activity | \mu band | both hands');
t.Color = 'k';
h = [];
topoplot(squeeze(new_GA_ERD_Act_773_mi), chanlocs16);
h = [h, gca];
set(h, 'clim', [-50 100]);
cb = colorbar;
cb.Color = 'k';

% beta band
figure
sgtitle('Grand Average | Topoplots | \beta band', 'Color', 'k')
subplot(221)
t = title('ERD/ERS during fixation | \beta band | both feet');
t.Color = 'k';
h = [];
topoplot(squeeze(new_GA_ERD_Ref_771_beta), chanlocs16);
h = [h, gca];
set(h, 'clim', [-50 100]);
cb = colorbar;
cb.Color = 'k';

subplot(222)
t = title('ERD/ERS during activity | \beta band | both feet');
t.Color = 'k';
h = [];
topoplot(squeeze(new_GA_ERD_Act_771_beta), chanlocs16);
h = [h, gca];
set(h, 'clim', [-50 100]);
cb = colorbar;
cb.Color = 'k';

subplot(223)
t = title('ERD/ERS during fixation | \beta band | both hands');
t.Color = 'k';
h = [];
topoplot(squeeze(new_GA_ERD_Ref_773_beta), chanlocs16);
h = [h, gca];
set(h, 'clim', [-50 100]);
cb = colorbar;
cb.Color = 'k';

subplot(224)
t = title('ERD/ERS during activity | \beta band | both hands');
t.Color = 'k';
h = [];
topoplot(squeeze(new_GA_ERD_Act_773_beta), chanlocs16);
h = [h, gca];
set(h, 'clim', [-50 100]);
cb = colorbar;
cb.Color = 'k';


%% Spectrogram analysis for data exploration

% neurophysiologically meaningful channels
chans = [7 9 11];  % [C3 Cz C4]

% defining the time vector
wshift  = 0.0625;      % [s] shift of the external window used in processing
time = (0:w_l_tr-1) * wshift;   % time axis in seconds

% visualize the ERD/ERS averaged across trials for the two MI classes
figure
sgtitle('Grand Average | Spectrograms', 'Color', 'w')

for i = 1:length(chans)
    
    % selected channel
    chan = chans(i);
    if chan == 7
        ch_name = 'Channel C3';
    elseif chan == 9
        ch_name = 'Channel Cz';
    elseif chan == 11
        ch_name = 'Channel C4';
    end

    subplot(2,3,i)
    tmp = cat(4,population.ERD_773);
    ERD_ch_773_avg = mean(tmp(:,:,chan,:),4)';
    % we have to transpose it, since it has l_tr (time) on y and nF (frequency) on x
    imagesc(time, file.f_sel, ERD_ch_773_avg)
    axis xy;
    colormap(hot(256));
    clim([-1.4 0.4]);
    colorbar;
    xlabel('Time [s]');
    ylabel('Frequency [Hz]');
    title([ch_name " | Both hand"]);
    xline(3, 'k')
    xline(4, 'k')

    subplot(2,3,i+3)
    tmp = cat(4,population.ERD_771);
    ERD_ch_771_avg = mean(tmp(:,:,chan,:),4)';
    % we have to transpose it, since it has l_tr (time) on y and nF (frequency) on x
    imagesc(time, file.f_sel, ERD_ch_771_avg)
    axis xy;
    colormap(hot(256));
    clim([-1.4 0.4]);
    colorbar;
    xlabel('Time [s]');
    ylabel('Frequency [Hz]');
    title([ch_name " | Both feet"]);
    xline(3, 'k')
    xline(4, 'k')
end





