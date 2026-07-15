function [FS, CFbk, Ck, nWE, Pk] = psd_to_fscore(file)
%
% data_to_fscore computes the Fischer's score of a file containing PSD data
%
%               [FS, CFbk, Ck, nWE] = psd_to_fscore(file)
%
% Input:
% • file: file processed with proc_all_gdf
% 
% Output:
% • FS: Fischer's score matrix [channel x frequency]
% • CFbk: continuous feedback labels
% • Ck: trials' label
% • nWE: number of extracted windows
% • Pk: extracted trials' label
%

% PSD extraction
PSD = file.PSD_sel; % only on the selected frequencies [window x frequency x channel]
[nW, nF, nC] = size(PSD);

% events extraction
EVENT.TYP = file.h.EVENT.TYP;
EVENT.DUR = file.wDUR;
EVENT.POS = file.wPOS;

% continuous feedback labels
CFbk = zeros(nW, 1);
% continuous feedbacks' starting points
strcf = EVENT.POS(EVENT.TYP==781); % [both feet, both hands]
% continuous feedbacks' end points
endcf = EVENT.POS(EVENT.TYP==781) + EVENT.DUR(EVENT.TYP==781) - 1;
% getting rid of the rest trials
strcf(EVENT.TYP(EVENT.TYP==771 | EVENT.TYP==773 | EVENT.TYP==783)==783) = [];
endcf(EVENT.TYP(EVENT.TYP==771 | EVENT.TYP==773 | EVENT.TYP==783)==783) = [];
nT = length(strcf);  % number of useful trials
% filling CFbk
for n = 1:nT
    CFbk(strcf(n):endcf(n)) = 1;
end
% we use this vector to highlight the periods of interest: from the cue to
% the end of the continuous feedback

% 771 = both feet
% 773 = both hands
% extracting the labels of each trials
tr_lab = EVENT.TYP(EVENT.TYP==771 | EVENT.TYP==773);
Ck = zeros(nW, 1);
% trials' starting points
strt = EVENT.POS(EVENT.TYP==786); % starting from the fixation cross event
% trials' end points
endt = EVENT.POS(EVENT.TYP==781) + EVENT.DUR(EVENT.TYP==781) - 1;
% getting rid of the rest trials
strt(EVENT.TYP(EVENT.TYP==771 | EVENT.TYP==773 | EVENT.TYP==783)==783) = [];
endt(EVENT.TYP(EVENT.TYP==771 | EVENT.TYP==773 | EVENT.TYP==783)==783) = [];
% filling Ck with the respective code of each trial
for i = 1:nT
    ind = strt(i):endt(i);
    Ck(ind) = tr_lab(i);
end

% extracting the PSD only in the periods of interest
P = PSD(CFbk == 1, :, :); 
% extracting the respective trials' label
Pk = Ck(CFbk == 1);

% number of extracted windows
nWE = size(P,1);

% reshape the data in the format [windows x features]
% where feature = (channel-frequency) pair

% Initializing the PSD reshaped matrix
P_rshp = zeros(nWE, nF*nC);
% each row is a window, each column is a feature
for i = 1:nWE
    rshp = P(i,:,:);
    P_rshp(i,:) = reshape(rshp, 1, nF*nC);
end
% each row of P_rshp is fill with all the frequencies of a channel,
% repeated for every channel (e.g., [f1ch1 f2ch1 ... fnch1 f1ch2 f2ch2 ...])

% F-score computation
FS = abs(mean(P_rshp(Pk==771,:))-mean(P_rshp(Pk==773,:)))...
    ./sqrt(var(P_rshp(Pk==771,:))+var(P_rshp(Pk==773,:)));
% reshaping it
FS = reshape(FS, nF, nC)';

