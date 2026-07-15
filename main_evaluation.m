%% Assignment 1 | 3. Evaluation
% model evaluation on online data

clearvars
close all
clc

%% Adding the function folder to the pathway
addpath(genpath([pwd '/functions']))

%% Channel labels
ch = {"Fz", "FC3", "FC1", "FCz", "FC2", "FC4", "C3", ...
    "C1", "Cz", "C2", "C4", "CP3", "CP1", "CPz", "CP2", "CP4"};

%% Features computation and model evaluation for each subject

% Extract all the subjects' subfolders
folders = dir(fullfile([pwd '/processed_dataset'],'a*')); % taking all the .gdf files
num_subj = size(folders, 1); % number of subjects

for j = 1:num_subj % for each subject

    fprintf('Working on %s...\n', folders(j).name);

    % Extract all the online data files in subject's folder
    dataDir = fullfile(pwd, 'processed_dataset', folders(j).name);
    files   = dir(fullfile(dataDir, '*online*.mat'));
    num_onf  = numel(files); % number of online files per subject
    subj_code = folders(j).name(1:end-length('_micontinuous'));

    % Initialize per-subject variables
    PSD  = [];
    CFbk = [];
    Ck   = [];
    Pk   = [];
    nWE  = 0;
    FS_subj_j = []; % FS per run
    EVENT=struct; EVENT.TYP=[]; EVENT.wDUR=[]; EVENT.wPOS=[];

    % figure('Name',['Subject ' subj_code ' – Runs']);

    % Loop over runs
    % for each online run we compute the Fisher Score
    for k = 1:num_onf % for each online file

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
        EVENT.TYP =  [EVENT.TYP; file.h.EVENT.TYP];
        EVENT.wDUR = [EVENT.wDUR; file.wDUR];
        EVENT.wPOS = [EVENT.wPOS; size(PSD,1)+file.wPOS];

        % Visualization per run
        subplot(1, num_onf, k)
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

    %% Feature extraction
    % loading the features selected from the offline calibration
    load(fullfile(dataDir,'sel_feat.mat'), 'sel_feat');

    % number of total windows of the concatenated PSDs
    nW = size(PSD,1);
    % number of selected features
    n_sel_feat = size(sel_feat,1);

    % data for the identified selected features
    F = zeros(nW, n_sel_feat);
    for i = 1:n_sel_feat
        F(:,i) = log(PSD(:, sel_feat(i,1), sel_feat(i,2)) + eps);
    end

    %% Model loading and evaluation

    load(fullfile(dataDir, 'decoder.mat'));

    LabelIdx = (CFbk == 1);
    F_sel  = F(LabelIdx,:);

    [Gk,pp] = predict(Model, F_sel); % [predicted class, posterior probability]
    % pp: [prob of 771, prob of 773]

    %% Metrics computation
    % Single sample accuracy on testset
    acc = (Pk == Gk);
    ssa_overall = 100*sum(acc)/nWE;
    ssa_bh = 100 * mean(acc(Pk == 773));
    ssa_bf = 100 * mean(acc(Pk == 771));

    figure('Name', ['Subject ' subj_code ' - SSA'])
    b = bar({'overall', 'both hands', 'both feet'}, [ssa_overall, ssa_bh, ssa_bf]);
    b.Labels = b.YData; % may not work for all Matlab version, in case of error just comment this line
    grid on
    title('Single sample accuracy on test set')
    ylim([0 100])
    ylabel('accuracy [%]')

    %% Exponential accumulation framework

    % new cue points (so new trial)
    strt = EVENT.wPOS(EVENT.TYP==771 | EVENT.TYP==773); 
    nT = length(strt);  % number of (total) trials
    
    % indexes (nW scale) of events (nWE scale)
    indEv = find(LabelIdx);
    
    D = 0.5 * ones(size(pp,1), 2); % [prob of 771, prob of 773]
    % two balanced classes -> the default value is 0.5
    
    % thresholds for trial based classification
    ths = [0.2 0.8];
    
    % initialize the vector of classified trials ...
    tr_pred = zeros(nT, 1);
    % ... and its index
    ind_tr = 1;

    % initialize the vector of time needed for giving a command...
    comm_time = zeros(nT, 1);
    % ... and its indeces
    ind_comm = [];
    
    % defining the integration parameter
    alpha = 0.92;
    
    % initialize class assignation array
    classes = [];
    
    % initialize the array of indexes of where new trials start
    ind_new_tr = zeros(nT,1);
    ind_new_tr(1) = 1;
    
    for wId = 2:size(pp,1)
        % Is the first sample of a new trial?
        % |-YES: Reset the current evidence (i.e., D(wId) to [0.5 0.5])
        % |-NO:  Keep integrating the value
        if indEv(wId) == indEv(wId-1)+1 
            D(wId,:) = D(wId - 1,:) * alpha + pp(wId,:) * (1 - alpha);
        
        % decision making
            if D(wId,1) <= ths(1) 
                classes = [classes 773];
                ind_comm = [ind_comm indEv(wId)];
            elseif D(wId,1) >= ths(2)
                classes = [classes 771];
                ind_comm = [ind_comm indEv(wId)];
            end
    
        else % we have moved to the next trial, so D is resetted to 0.5
            % we store the indexes of new trials' start
            ind_new_tr(ind_tr+1) = wId;

            % we classify the trial
            if ~isempty(classes)
                tr_pred(ind_tr) = classes(1); % we consider only the first label (first reached threshold)
            end % otherwise the trial is classified as 0 (rejected)

            % we compute after how many samples a threshold is reached
            if ~isempty(ind_comm)
                comm_time(ind_tr) = ind_comm(1)-indEv(ind_new_tr(ind_tr))-1; % we consider only the first label (first reached threshold)
            end % otherwise the time is set to 0

            % we move to the next trial
            ind_tr = ind_tr + 1;
            % and we reset the classes and commands arrays
            classes = [];
            ind_comm = [];
        end
    end
    % and we classify the last trial
    if ~isempty(classes)
        tr_pred(ind_tr) = classes(1); % we consider only the first label (first reached threshold)
    end % otherwise the trial is classified as 0 (rejected)
    
    %% Visualization
    % plotting one random trial
    trial_to_plot = 29;
    ind_trial_to_plot = ind_new_tr(trial_to_plot):ind_new_tr(trial_to_plot+1)-1;
    
    figure('Name', ['Subject ' subj_code ' - Acc. Framework'])
    hold on
    plot(pp(ind_trial_to_plot,1), 'wo')
    plot(D(ind_trial_to_plot,1), 'y', 'Linewidth', 2)
    title(['Trial ' num2str(trial_to_plot) ' with accumulation framework'])
    yline(ths(1), 'r')
    yline(ths(2), 'g')
    legend('pp', 'D', '773', '771')
    axis([1 length(ind_trial_to_plot) 0 1])
    
    %disp(['Predicted class: ' num2str(tr_pred(trial_to_plot))]) % commented for the sake of output readability
    
    %% Trial accuracy with and without rejection
    % ground truth
    true_labels = EVENT.TYP(EVENT.TYP==771 | EVENT.TYP==773);
    
    % indexes of rejected trials
    ind_rej = tr_pred==0;
    
    acc_trial = 100 * sum(tr_pred==true_labels)/nT;
    acc_trial_rej = 100 * sum(tr_pred(~ind_rej)==true_labels(~ind_rej))/length(tr_pred(~ind_rej));
    
    figure('Name', ['Subject ' subj_code ' - Trial acc.'])
    b = bar({'no rejection', 'rejection'}, [acc_trial acc_trial_rej]);
    b.Labels = b.YData; % may not work for all Matlab version, in case of error just comment this line
    grid on
    ylim([0 100])
    ylabel('accuracy [%]')

    %% Average time to deliver a command
    
    % shift of the external window during PSD computation (main_data_processing.m)
    wshift  = 0.0625;      % [s] 

    % computed not considering the rejected trials
    comm_time_avg = mean(comm_time(comm_time ~= 0))*wshift; % [s]

    %% Classification report

    % single sample chance level
    ss_chance_lvl = 100*max(sum(Pk==771), sum(Pk==773))/length(Pk);

    % trial chance level (with rejection)
    tr_chance_lvl = 100*max(sum(true_labels(~ind_rej)==771), sum(true_labels(~ind_rej)==773))/length(true_labels(~ind_rej));

    % display
    fprintf('  Single sample Chance level:              %4.2f %%\n', ss_chance_lvl)
    fprintf('  Single sample Accuracy (overall):        %4.2f %%\n', ssa_overall)
    fprintf('  Single sample Accuracy (both hands):     %4.2f %%\n', ssa_bh)
    fprintf('  Single sample Accuracy (both feet):      %4.2f %%\n', ssa_bf)
    fprintf('  Trial Chance level (w/ rejection):       %4.2f %%\n', tr_chance_lvl)
    fprintf('  Trial Accuracy (w/o rejection):          %4.2f %%\n', acc_trial)
    fprintf('  Trial Accuracy (w/ rejection):           %4.2f %%\n', acc_trial_rej)
    fprintf('  Average time to deliver a command:       %4.2f s\n',  comm_time_avg)
    fprintf('\n')

end

