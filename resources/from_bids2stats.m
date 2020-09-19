% BIDS Tools / EEGLAB / LIMO EEG 
% data analysis of Wakeman and Henson 2015 data
% Arnaud Delorme & Cyril Pernet


%% Import
% start EEGLAB
clear
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;
 
% call BIDS tool BIDS
studypath       = 'XXX\WakemanHenson_Faces\eeg';
[STUDY, ALLEEG] = pop_importbids(filepath, 'bidsevent','on','bidschanloc','on', 'studyName','Face_detection');
ALLEEG          = pop_select( ALLEEG, 'nochannel',{'EEG061','EEG062','EEG063','EEG064'});
CURRENTSTUDY    = 1; 
EEG             = ALLEEG; 
CURRENTSET      = 1:length(EEG);


%% Preprocessing

% Remove bad channels
EEG = pop_clean_rawdata( EEG,'FlatlineCriterion',5,'ChannelCriterion',0.8,...
    'LineNoiseCriterion',4,'Highpass',[0.25 0.75] ,...
    'BurstCriterion','off','WindowCriterion','off','BurstRejection','off',...
    'Distance','Euclidian','WindowCriterionTolerances','off' );
 
% Rereference using average reference
EEG = pop_reref( EEG,[],'interpchan',[]);
 
% Run ICA and flag artifactual components using IClabel
for s=1:size(EEG,2)
    EEG(s) = pop_runica(EEG(s), 'icatype','runica','concatcond','on','options',{'pca',EEG(s).nbchan-1});
    EEG(s) = pop_iclabel(EEG(s),'default');
    EEG(s) = pop_icflag(EEG(s),[NaN NaN;0.8 1;0.8 1;NaN NaN;NaN NaN;NaN NaN;NaN NaN]);
    EEG(s) = pop_subcomp(EEG(s), find(EEG(s).reject.gcompreject), 0);
end
 
% clear data using ASR - just the bad epochs
EEG = pop_clean_rawdata( EEG,'FlatlineCriterion','off','ChannelCriterion','off',...
    'LineNoiseCriterion','off','Highpass','off','BurstCriterion',20,...
    'WindowCriterion',0.25,'BurstRejection','on','Distance','Euclidian',...
    'WindowCriterionTolerances',[-Inf 7] );

% Extract data epochs (no baseline removed)
EEG    = pop_epoch( EEG,{'famous_new','famous_second_early','famous_second_late','scrambled_new','scrambled_second_early','scrambled_second_late','unfamiliar_new','unfamiliar_second_early','unfamiliar_second_late'},...
    [-0.5 1] ,'epochinfo','yes');
EEG    = eeg_checkset(EEG);
EEG    = pop_saveset(EEG, 'savemode', 'resave');
ALLEEG = EEG;
 
% update study & compute single trials
STUDY        = std_checkset(STUDY, ALLEEG);
[STUDY, EEG] = std_precomp(STUDY, EEG, {}, 'savetrials','on','interp','on','recompute','on',...
    'erp','on','erpparams', {'rmbase' [-200 0]}, 'spec','off', 'ersp','off','itc','off');
eeglab redraw

%% Statitiscal analysis
% to restart the analysis from here - simply reload the STUDY see pop_loadstudy

chanlocs = [STUDY.filepath filesep 'limo_gp_level_chanlocs.mat'];

% two-way ANOVA faces * repetition
% -------------------------------

% 1st level analysis
STUDY = std_makedesign(STUDY, ALLEEG, 1, 'name','FaceRepetition','delfiles','off','defaultdesign','off',...
    'variable1','type','values1',{'famous_new','famous_second_early','famous_second_late','scrambled_new','scrambled_second_early','scrambled_second_late','unfamiliar_new','unfamiliar_second_early','unfamiliar_second_late'},...
    'vartype1','categorical','subjselect',{'sub-002','sub-003','sub-004','sub-005','sub-006','sub-007','sub-008','sub-009','sub-010','sub-011','sub-012','sub-013','sub-014','sub-015','sub-016','sub-017','sub-018','sub-019'});
[STUDY, EEG] = pop_savestudy( STUDY, EEG, 'savemode','resave');
STUDY  = pop_limo(STUDY, ALLEEG, 'method','WLS','measure','daterp','timelim',[-50 650],'erase','on','splitreg','off','interaction','off');

% 2nd level analysis
mkdir([STUDY.filepath filesep '2-ways-ANOVA'])
cd([STUDY.filepath filesep '2-ways-ANOVA'])
limo_random_select('Repeated Measures ANOVA',chanlocs,'LIMOfiles',...
    {[STUDY.filepath filesep 'LIMO_Face_detection' filesep 'Beta_files_FaceRepetition_GLM_Channels_Time_WLS.txt']},...
    'analysis_type','Full scalp analysis','parameters',{[1 2 3],[4 5 6],[7 8 9]},...
    'factor names',{'face','repetition'},'type','Channels','nboot',1000,'tfce',0,'skip design check','yes');
limo_eeg(5) % print signitifanct results and create LIMO.data.timevect

% make topoplots
LIMO  = load('LIMO.mat'); LIMO = LIMO.LIMO;
stats = load('Rep_ANOVA_Main_effect_1_face.mat');
stats = stats.(cell2mat(fieldnames(stats)));
stats = squeeze(stats(:,:,1)); % keep F values
cc    = limo_color_images(stats);
opt   = {'electrodes','on','maplimits','maxmin','verbose','off','colormap', cc};
% cluster 1 starts at 140ms ends at 424ms, maximum F values 64.1281 at 280ms on channel EEG017
% cluster 2 starts at 440ms ends at 648ms, maximum F value 17.6071 at 616ms on channel EEG057
time = [49 59 84 120 124 168 176];
figure('Name','Topoplot Main face effect')
for t = 1:length(time)
    subplot(2,ceil(length(time)/2),t);
    topoplot(squeeze(stats(:,time(t),:)),LIMO.data.chanlocs,opt{:});
    title(sprintf('%g ms',LIMO.data.timevect(time(t))))
end

stats = load('Rep_ANOVA_Main_effect_2_repetition.mat');
stats = stats.(cell2mat(fieldnames(stats)));
stats = squeeze(stats(:,:,1)); % keep F values
cc    = limo_color_images(stats);
opt   = {'electrodes','on','maplimits','maxmin','verbose','off','colormap', cc};
% cluster starts at 232ms ends at 648ms, maximum F value 51.3596 at 612ms channel EEG045
time = [72 167 176];
figure('Name','Topoplot Main repetition effect')
for t = 1:length(time)
    subplot(1,length(time),t);
    topoplot(squeeze(stats(:,time(t),:)),LIMO.data.chanlocs,opt{:});
    title(sprintf('%g ms',LIMO.data.timevect(time(t))))
end

% compute the mean ERPs to visualize differences
Files = [STUDY.filepath filesep 'LIMO_' STUDY.filename(1:end-6) filesep ...
    'LIMO_files_FaceRepetition_GLM_Channels_Time_WLS.txt'];
parameters = [1 2 3];
savename1  = [pwd filesep 'famous_faces.mat'];
limo_central_tendency_and_ci(Files, parameters, chanlocs, 'Weighted mean', 'Mean', [],savename1) 
parameters = [4 5 6];
savename2  = [pwd filesep 'srambled_faces.mat'];
limo_central_tendency_and_ci(Files, parameters, chanlocs, 'Weighted mean', 'Mean', [],savename2) 
parameters = [7 8 9];
savename3  = [pwd filesep 'unfamiliar_faces.mat'];
limo_central_tendency_and_ci(Files, parameters, chanlocs, 'Weighted mean', 'Mean', [],savename3) 
limo_add_plots([savename1(1:end-4) '_Mean_of_Weighted mean.mat'],...
    [savename2(1:end-4) '_Mean_of_Weighted mean.mat'],[savename3(1:end-4) '_Mean_of_Weighted mean.mat'],...
    'channel',50); title('Face type at channel 50')
% get a measure of raw effect size based on famous faces peak
tmp = load([savename1(1:end-4) '_Mean_of_Weighted mean.mat']);
[~,peaktime]=min(tmp.Data.mean(50,:,2));
tmp = load([savename1(1:end-4) '_single_subjects_Weighted mean.mat']);
effect_size = tmp.Data.data;
tmp = load([savename3(1:end-4) '_single_subjects_Weighted mean.mat']);
effect_size = (effect_size + tmp.Data.data)./2; % mean of faces
tmp = load([savename2(1:end-4) '_single_subjects_Weighted mean.mat']);
effect_size = effect_size - tmp.Data.data;  % faces vs scrambled
TM = limo_trimmed_mean(squeeze(effect_size),0, 0.05);
fprintf('Faces vs. Scrambled @ peak = %g uV CI=[%g %g]\n',TM(50,peaktime,2),TM(50,peaktime,3),TM(50,peaktime,1))

parameters = [1 4 7];
savename1  = [pwd filesep 'first_time.mat'];
limo_central_tendency_and_ci(Files, parameters, chanlocs, 'Weighted mean', 'Mean', [],savename1) 
parameters = [2 5 8];
savename2  = [pwd filesep 'second_time.mat'];
limo_central_tendency_and_ci(Files, parameters, chanlocs, 'Weighted mean', 'Mean', [],savename2) 
parameters = [3 6 9];
savename3  = [pwd filesep 'third_time.mat'];
limo_central_tendency_and_ci(Files, parameters, chanlocs, 'Weighted mean', 'Mean', [],savename3) 
limo_add_plots([savename1(1:end-4) '_Mean_of_Weighted mean.mat'],...
    [savename2(1:end-4) '_Mean_of_Weighted mean.mat'],[savename3(1:end-4) '_Mean_of_Weighted mean.mat'],...
    'channel',45); title('Repetition effect at channel 45')
% get a measure of raw effect size based on 2nd repetition peak
tmp = load([savename2(1:end-4) '_Mean_of_Weighted mean.mat']);
[~,peaktime]=max(tmp.Data.mean(45,:,2));
tmp = load([savename1(1:end-4) '_single_subjects_Weighted mean.mat']);
effect_size = tmp.Data.data;
tmp = load([savename3(1:end-4) '_single_subjects_Weighted mean.mat']);
effect_size = (effect_size + tmp.Data.data)./2; % mean of 1st and 3rd presentation
tmp = load([savename2(1:end-4) '_single_subjects_Weighted mean.mat']);
effect_size = effect_size - tmp.Data.data;  % 1st/3rd vs 2nd
TM = limo_trimmed_mean(squeeze(effect_size),0, 0.05);
fprintf('2nd vs. 1st&3rd presentation @ peak = %g uV CI=[%g %g]\n',TM(45,peaktime,2),TM(45,peaktime,1),TM(45,peaktime,3))


%% one-way ANOVA on time distance between trials for each face condition
% ----------------------------------------------------------------------

% 1st level analysis

STUDY = std_makedesign(STUDY, ALLEEG, 2, 'name','Face_time','delfiles','off','defaultdesign','off',...
    'variable1','face_type','values1',{'famous','scrambled','unfamiliar'},'vartype1','categorical',...
    'variable2','time_dist','values2',[],'vartype2','continuous',...
    'subjselect',{'sub-002','sub-003','sub-004','sub-005','sub-006','sub-007','sub-008','sub-009','sub-010','sub-011','sub-012','sub-013','sub-014','sub-015','sub-016','sub-017','sub-018','sub-019'});
[STUDY, EEG] = pop_savestudy( STUDY, EEG, 'savemode','resave');
STUDY = pop_limo(STUDY, ALLEEG, 'method','WLS','measure','daterp','timelim',[-50 650],'erase','on','splitreg','on','interaction','off');

% 2nd level analysis
mkdir([STUDY.filepath filesep 'derivatives' filesep '1-way-ANOVA'])
cd([STUDY.filepath filesep 'derivatives' filesep '1-way-ANOVA'])
limo_random_select('Repeated Measures ANOVA',chanlocs,'LIMOfiles',...
    {[STUDY.filepath filesep 'LIMO_Face_detection' filesep 'Beta_files_Face_time_GLM_Channels_Time_WLS.txt']},...
    'analysis_type','Full scalp analysis','parameters',{4 5 6},...
    'factor names',{'face','repetition'},'type','Channels','nboot',1000,'tfce',0,'skip design check','yes');

% add contrast
load('LIMO.mat')
LIMO.contrast{1}.C = [1 -05 -05 0];             % set contrast, note the last column 0 for constant`
LIMO.contrast{1}.V = 'F';                       % always F test for repeated measure ANOVA`
limo_contrast([pwd filesep 'Yr.mat'], LIMO, 3); % run the contrast
limo_contrast([pwd filesep 'Yr.mat'], LIMO, 4)  % apply bootstrap
save LIMO LIMO

% make topoplots
stats = load('ess_1.mat');
stats = stats.(cell2mat(fieldnames(stats)));
stats = squeeze(stats(:,:,end-1)); % keep F values
[peakchannel,peaktime]=ind2sub(size(stats),find(stats == max(stats(:))));
cc    = limo_color_images(stats);
opt   = {'electrodes','on','maplimits','maxmin','verbose','off','colormap', cc};
% cluster 1 starts at 536ms ends at 576ms, maximum 15.1112 @ 556ms channel EEG039
% cluster 2 starts at 548ms ends at 584ms, maximum 9.58689 @ 572ms channel EEG055
time = [148 153 158 151 157 160];
timevect = LIMO.data.start:1/LIMO.data.sampling_rate*1000:LIMO.data.end;
figure('Name','Topoplot effect of time famous faces more than others')
for t = 1:length(time)
    subplot(2,ceil(length(time)/2),t);
    topoplot(squeeze(stats(:,time(t),:)),LIMO.data.chanlocs,opt{:});
    title(sprintf('%g ms',timevect(time(t))))
end

% extract regression parameters to make 2D plots at max channels
data = load('Yr.mat'); data = data.Yr([39 55],:,:,:);
[mean_reg39,ci_reg39]=data_plot(squeeze(data(1,153,:,:)),'estimator','mean');
title('Average regression coef - channel 39 @556ms')
[mean_reg55,ci_reg55]=data_plot(squeeze(data(2,157,:,:)),'estimator','mean');
title('Average regression coef - channel 55 @572ms')
