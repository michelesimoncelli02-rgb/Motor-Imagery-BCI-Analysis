%% Assignment 1 | 2. Calibration
% features selection and extraction and model calibration on offline data

clearvars
close all
clc

%% Adding the function folder to the pathway
addpath(genpath([pwd '/functions']))

%% Channel labels
ch = {"Fz", "FC3", "FC1", "FCz", "FC2", "FC4", "C3", ...
    "C1", "Cz", "C2", "C4", "CP3", "CP1", "CPz", "CP2", "CP4"};

%% Features extraction and model training for each subject

% Extract all the subjects' subfolders
folders = dir(fullfile([pwd '/processed_dataset'],'a*')); % taking all the .gdf files
num_subj = size(folders, 1); % number of subjects

for j = 1:num_subj % for each subject

    fprintf('Working on %s...\n', folders(j).name);

    % Extract all the offline data files in subject's folder
    dataDir = fullfile(pwd, 'processed_dataset', folders(j).name);
    files   = dir(fullfile(dataDir, '*offline*.mat'));
    num_of  = numel(files); % number of offline files per subject
    subj_code = folders(j).name(1:end-length('_micontinuous'));

    % Initialize per-subject variables
    PSD  = [];
    CFbk = [];
    Ck   = [];
    Pk   = [];
    nWE  = 0;
    FS_subj_j = []; % FS per run

    figure('Name',['Subject ' subj_code ' – Runs']);

    % Loop over runs
    % for each offline run we compute the Fisher Score
    for k = 1:num_of % for each offline file

        file = load(fullfile(files(k).folder, files(k).name));

        % Concatenate PSDs on the selected frequencies along the window dimension
        PSD = cat(1, PSD, file.PSD_sel); % [windows x frequencies x channels]
        % selected frequencies
        f_sel = file.f_sel; % [Hz]

        % computing the Fisher Score and useful arrays
        [FS_k, CFbk_k, Ck_k, nWE_k, Pk_k] = psd_to_fscore(file);

        FS_subj_j(:,:,k) = FS_k;

        % Accumulate labels (window-level)
        CFbk = [CFbk; CFbk_k];
        Ck   = [Ck;   Ck_k];
        Pk   = [Pk;   Pk_k];
        nWE  =  nWE + nWE_k;

        % Visualization per run
        subplot(1, num_of, k)
        imagesc(FS_k)
        axis square
        colorbar
        yticks(1:numel(ch))
        yticklabels(ch)
        xticks(1:numel(f_sel))
        xticklabels(f_sel)
        xlabel('Frequency [Hz]')
        ylabel('Channel')
        title(['Run ' num2str(k)])

    end

    %% Subject-average Fisher Score
    FS_subj_avg = mean(FS_subj_j, 3);
    FS_subj(:,:,j) = FS_subj_avg;

    % visualization of the average FS
    figure('Name',['Subject ' subj_code ' – Average FS'])
    imagesc(FS_subj_avg)
    axis square
    colorbar
    yticks(1:numel(ch))
    yticklabels(ch)
    xticks(1:numel(f_sel))
    xticklabels(f_sel)
    xlabel('Frequency [Hz]')
    ylabel('Channel')

    %% Feature selection
    % getting the 5 largest values in FS and their indeces
    [max_feat, ind_feat] = maxk(FS_subj_avg(:), 5);
    % from these, we take only the ones above the 75% of the highest one
    ind_feat = ind_feat(max_feat>0.75*max_feat(1));
    % and we extract the respective indeces in the FS matrix
    [row_feat, col_feat] = ind2sub(size(FS_subj_avg), ind_feat);

    % selected features
    sel_feat = [col_feat, row_feat];  % [frequencies x channels]

    fprintf('Selected features:\n');
    for i = 1:size(sel_feat,1)
        fprintf('  %s @ %d Hz\n', ch{sel_feat(i,2)}, f_sel(sel_feat(i,1)));
    end
    fprintf('\n')

    %% Feature extraction
    % number of total windows of the concatenated PSDs
    nW = size(PSD,1);
    % number of selected features
    n_sel_feat = size(sel_feat,1);

    % data for the identified selected features
    F = zeros(nW, n_sel_feat);
    for i = 1:n_sel_feat
        F(:,i) = log(PSD(:, sel_feat(i,1), sel_feat(i,2)) + eps);
    end

    %% Model training
    LabelIdx = (CFbk == 1);

    F_sel  = F(LabelIdx,:);
    Ck_sel = Ck(LabelIdx);
    
    Model = fitcdiscr(F_sel, Ck_sel, 'DiscrimType','quadratic');      

    %% Model evaluation
    Gk = predict(Model, F_sel); 

    %% Metrics computation
    % Single sample accuracy on trainset
    acc = (Pk == Gk);
    ssa_overall = 100*sum(acc)/nWE;
    ssa_bh = 100 * mean(acc(Pk == 773));
    ssa_bf = 100 * mean(acc(Pk == 771));

    figure('Name', ['Subject ' subj_code ' - SSA'])
    b = bar({'overall', 'both hands', 'both feet'}, [ssa_overall, ssa_bh, ssa_bf]);
    b.Labels = b.YData; % may not work for all Matlab version, in case of error just comment this line
    grid on
    title('Single sample accuracy on train set')
    ylim([0 100])
    ylabel('accuracy [%]')

    %% Savings

    outDir = fullfile(pwd, 'processed_dataset', folders(j).name);

    % the decoder
    outName_decoder = fullfile(outDir, 'decoder.mat');
    if ~exist(outName_decoder, 'file')
        save(outName_decoder, 'Model')
    end

    % the extracted features
    outName_feature = fullfile(outDir, 'sel_feat.mat');
    if ~exist(outName_feature, 'file')
        save(outName_feature, 'sel_feat')
    end

end

%% Average FS across subjects

FS_pop = mean(FS_subj, 3);

figure('Name', 'Population Average FS')
imagesc(FS_pop)
axis square
colorbar
yticks(1:numel(ch))
yticklabels(ch)
xticks(1:numel(f_sel))
xticklabels(f_sel)
xlabel('Frequency [Hz]')
ylabel('Channel')




