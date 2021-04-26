% BIDS Tools / EEGLAB / LIMO EEG 
% data analysis of Wakeman and Henson 2015 data
% Scripted by Cyril Pernet 
% ---------------------------------------------
%%                Spectral Analyses
% ----------------------------------------------

%% Import
% start EEGLAB
% clear
% [ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;
%  
% % call BIDS tool BIDS
% studypath        = 'XXX\WakemanHenson_Faces\eeg';
% [STUDY, ALLEEG] = pop_importbids(filepath, 'bidsevent','on','bidschanloc','on', 'studyName','Face_detection');
% ALLEEG = pop_select( ALLEEG, 'nochannel',{'EEG061','EEG062','EEG063','EEG064'});
% CURRENTSTUDY = 1; EEG = ALLEEG; CURRENTSET = [1:length(EEG)];
% 
% 
%% Preprocessing
% 
% % Remove bad channels
% EEG = pop_clean_rawdata( EEG,'FlatlineCriterion',5,'ChannelCriterion',0.8,...
%     'LineNoiseCriterion',4,'Highpass',[0.25 0.75] ,...
%     'BurstCriterion','off','WindowCriterion','off','BurstRejection','off',...
%     'Distance','Euclidian','WindowCriterionTolerances','off' );
%  
% % Rereference using average reference
% EEG = pop_reref( EEG,[],'interpchan',[]);
%  
% % Run ICA and flag artifactual components using IClabel
% for s=1:size(EEG,2)
%     EEG(s) = pop_runica(EEG(s), 'icatype','runica','concatcond','on','options',{'pca',EEG(s).nbchan-1});
%     EEG(s) = pop_iclabel(EEG(s),'default');
%     EEG(s) = pop_icflag(EEG(s),[NaN NaN;0.8 1;0.8 1;NaN NaN;NaN NaN;NaN NaN;NaN NaN]);
%     EEG(s) = pop_subcomp(EEG(s), find(EEG(s).reject.gcompreject), 0);
% end
%  
% % clear data using ASR - just the bad epochs
% EEG = pop_clean_rawdata( EEG,'FlatlineCriterion','off','ChannelCriterion','off',...
%     'LineNoiseCriterion','off','Highpass','off','BurstCriterion',20,...
%     'WindowCriterion',0.25,'BurstRejection','on','Distance','Euclidian',...
%     'WindowCriterionTolerances',[-Inf 7] );
% 
% % Extract data epochs (no baseline removed)
% EEG    = pop_epoch( EEG,{'famous_new','famous_second_early','famous_second_late','scrambled_new','scrambled_second_early','scrambled_second_late','unfamiliar_new','unfamiliar_second_early','unfamiliar_second_late'},...
%     [-0.5 1] ,'epochinfo','yes');
% EEG    = eeg_checkset(EEG);
% EEG    = pop_saveset(EEG, 'savemode', 'resave');
% ALLEEG = EEG;
% [STUDY, EEG] = pop_savestudy( STUDY, EEG, 'savemode','resave');

% to restart the analysis from here - simply reload the STUDY see pop_loadstudy
studyfullname       = 'F:\WakemanHenson_Faces\eeg\derivatives\Face_detection.study';
[root,std_name,ext] = fileparts(studyfullname); cd(root);       
EEG                 = eeglab;
[STUDY, ALLEEG]     = pop_loadstudy('filename', [std_name ext], 'filepath', root);
STUDY               = std_checkset(STUDY, ALLEEG);
[STUDY, EEG]        = std_precomp(STUDY, ALLEEG, {}, 'savetrials','on','interp','on','recompute','on',...
    'erp','off', 'spec','on', 'ersp','off','itc','off');
eeglab redraw

%% One way repeated measures ANOVA (Famous, Unfamiliar, Scrambled faces as conditions)
% https://github.com/LIMO-EEG-Toolbox/limo_meeg/wiki/2.-One-way-repeated-measures-ANOVA-(Famous,-Unfamiliar,-Scrambled-faces-as-conditions)
% ---------------------------------------------------------------
% 1st level analysis - specify the design
% We ignore the repetition levels using the variable 'face_type'
STUDY = std_makedesign(STUDY, ALLEEG, 1, 'name','ANOVA_Faces','delfiles','off','defaultdesign','off',...
    'variable1','face_type','values1',{'famous','scrambled','unfamiliar'},'vartype1','categorical',...
    'subjselect',{'sub-002','sub-003','sub-004','sub-005','sub-006','sub-007','sub-008','sub-009',...
    'sub-010','sub-011','sub-012','sub-013','sub-014','sub-015','sub-016','sub-017','sub-018','sub-019'});
[STUDY, EEG] = pop_savestudy( STUDY, EEG, 'savemode','resave');

% 1st level analysis - estimate parameters
STUDY  = pop_limo(STUDY, ALLEEG, 'method','WLS','measure','datspec','freqlim',[3 45],'erase','on','splitreg','off','interaction','off');

% 2nd level analysis - ANOVA on Beta parameters 1 2 3
chanlocs = [STUDY.filepath filesep 'limo_gp_level_chanlocs.mat'];
mkdir([STUDY.filepath filesep '1-way-ANOVA'])
cd([STUDY.filepath filesep '1-way-ANOVA'])
limo_random_select('Repeated Measures ANOVA',chanlocs,'LIMOfiles',...
    {[STUDY.filepath filesep 'LIMO_Face_detection' filesep 'Beta_files_ANOVA_Faces_GLM_Channels_Frequency_WLS.txt']},...
    'analysis_type','Full scalp analysis','parameters',{[1 2 3]},...
    'factor names',{'face'},'type','Channels','nboot',1000,'tfce',0,'skip design check','yes');

% add contrast famous+unfamiliar>scrambled
limo_contrast(fullfile(pwd,'Yr.mat'),fullfile(pwd,'LIMO.mat'), 3 ,[0.5 -1 0.5]); % compute a new contrast
limo_contrast(fullfile(pwd,'Yr.mat'),fullfile(pwd,'LIMO.mat'), 4);               % do the bootstrap of the last contrast

% figures
limo_eeg(5,pwd) % channel*freq imagesc for both effects and contrast
limo_display_results(3,'ess_1.mat',pwd,0.05,2,fullfile(pwd,'LIMO.mat'),0,'channels',49); % course plot
saveas(gcf, 'contrast_spectrum.fig'); close(gcf)

%% One way repeated measures ANOVA revised (Famous, Unfamiliar, Scrambled faces as 1st level contrasts)
% https://github.com/LIMO-EEG-Toolbox/limo_meeg/wiki/3.--One-way-repeated-measures-ANOVA-revised-(Famous,-Unfamiliar,-Scrambled-faces-as-1st-level-contrasts)
% Everything is modelled and contrasts are computed to merge repetition levels
% ---------------------------------------------------------------------------------
% 1st level analysis - specify the design
% Note we use the variable 'type' and use cells within a cell array to
% indicate grouping of conditions - this  means contrasts will be computed
% pooling those levels 
STUDY = std_makedesign(STUDY, ALLEEG, 2, 'name','FaceRepAll','delfiles','off','defaultdesign','off',...
    'variable1','type','values1',{{'famous_new','famous_second_early','famous_second_late'},...
    {'scrambled_new','scrambled_second_early','scrambled_second_late'},...
    {'unfamiliar_new','unfamiliar_second_early','unfamiliar_second_late'}},'vartype1','categorical',...
    'subjselect',{'sub-002','sub-003','sub-004','sub-005','sub-006','sub-007','sub-008','sub-009',...
    'sub-010','sub-011','sub-012','sub-013','sub-014','sub-015','sub-016','sub-017','sub-018','sub-019'});
[STUDY, EEG] = pop_savestudy( STUDY, EEG, 'savemode','resave');

% 1st level analysis - estimate parameters
[STUDY]      = pop_limo(STUDY, ALLEEG, 'method','WLS','measure','datspec','freqlim',[3 45],'erase','on','splitreg','off','interaction','off');

% 2nd level analysis - ANOVA on con_1, con_2, con_3
chanlocs   = [STUDY.filepath filesep 'limo_gp_level_chanlocs.mat'];
con1_files = fullfile(STUDY.filepath,['LIMO_' STUDY.filename(1:end-6)],'con_1_files_FaceRepAll_GLM_Channels_Frequency_WLS.txt');
con2_files = fullfile(STUDY.filepath,['LIMO_' STUDY.filename(1:end-6)],'con_2_files_FaceRepAll_GLM_Channels_Frequency_WLS.txt');
con3_files = fullfile(STUDY.filepath,['LIMO_' STUDY.filename(1:end-6)],'con_3_files_FaceRepAll_GLM_Channels_Frequency_WLS.txt');

mkdir([STUDY.filepath filesep '1-way-ANOVA-revised'])
cd([STUDY.filepath filesep '1-way-ANOVA-revised'])
limo_random_select('Repeated Measures ANOVA',chanlocs,'LIMOfiles', {con1_files,con2_files,con3_files},...
    'analysis_type','Full scalp analysis','parameters',{[1 1 1]},...
    'factor names',{'face'},'type','Channels','nboot',1000,'tfce',0,'skip design check','yes');

% add contrast famous+unfamiliar>scrambled
limo_contrast(fullfile(pwd,'Yr.mat'),fullfile(pwd,'LIMO.mat'), 3 ,[0.5 -1 0.5]); % compute a new contrast
limo_contrast(fullfile(pwd,'Yr.mat'),fullfile(pwd,'LIMO.mat'), 4);               % do the bootstrap of the last contrast

% figures
limo_eeg(5,pwd) % channel*freq imagesc for both effects and contrast
limo_display_results(3,'ess_1.mat',pwd,0.05,2,fullfile(pwd,'LIMO.mat'),0,'channels',49); % course plot
saveas(gcf, 'contrast_spectrum.fig'); close(gcf)

%% let's check effect sizes
% https://github.com/LIMO-EEG-Toolbox/limo_meeg/wiki/4.-Summary-statistics:-Effects-and-Effect-sizes

% the ANOVA plot shows the mean value differences corresponding to the contrast
% wich is only 2 differences - not that useful
limo_display_results(3,'Rep_ANOVA_Main_effect_1_face.mat',pwd,0.05,1,...
    fullfile(pwd,'LIMO.mat'),0,'channels',49,'sumstats','mean'); % course plot
saveas(gcf, 'Rep_ANOVA_Main_effect_spectrum.fig'); close(gcf)

% plot mean betas per condition ; this represents the data on which the
% ANOVA was conputed - as such this is the right measure
Yr       = load('Yr.mat'); % ANOVA data are channel*[freq/time]frames*subjects*conditions
chanlocs = [STUDY.filepath filesep 'limo_gp_level_chanlocs.mat'];
mkdir('average_betas'); cd('average_betas');
limo_central_tendency_and_ci(squeeze(Yr.Yr(:,:,:,1)), 'Mean',[],fullfile(pwd,'famous.mat'))
limo_central_tendency_and_ci(squeeze(Yr.Yr(:,:,:,2)), 'Mean',[],fullfile(pwd,'scrambled.mat'))
limo_central_tendency_and_ci(squeeze(Yr.Yr(:,:,:,3)), 'Mean',[],fullfile(pwd,'unfamiliar.mat'))
limo_add_plots('channel',49,fullfile(fileparts(pwd),'LIMO.mat'),...
    {fullfile(pwd,'famous_Mean.mat'), fullfile(pwd,'scrambled_Mean.mat'), fullfile(pwd,'unfamiliar_Mean.mat')})
title('mean betas channel 49'); saveas(gcf, 'Rep_ANOVA_Main_effect_Betas.fig'); close(gcf)
cd ..

% we can go further down, to the ERP (XB in our GLM) 
LIMOfiles = fullfile(STUDY.filepath,['LIMO_' STUDY.filename(1:end-6)],'LIMO_files_FaceRepAll_GLM_Channels_Frequency_WLS.txt');
mkdir('Spectra'); cd('Spectra');
limo_central_tendency_and_ci(LIMOfiles, [1 2 3], chanlocs, 'Weighted mean', 'Mean', [],fullfile(pwd,'famous.mat'))  
limo_central_tendency_and_ci(LIMOfiles, [4 5 6], chanlocs, 'Weighted mean', 'Mean', [],fullfile(pwd,'scrambled.mat'))  
limo_central_tendency_and_ci(LIMOfiles, [7 8 9], chanlocs, 'Weighted mean', 'Mean', [],fullfile(pwd,'unfamiliar.mat'))  
limo_add_plots('channel',49,fullfile(fileparts(pwd),'LIMO.mat'),...
    {fullfile(pwd,'famous_Mean_of_Weighted mean.mat'),fullfile(pwd,'scrambled_Mean_of_Weighted mean.mat'),...
    fullfile(pwd,'unfamiliar_Mean_of_Weighted mean.mat')}); title('ERPs channel 49')
saveas(gcf, 'Rep_ANOVA_Main_effect_ERPs.fig'); close(gcf)

% single subjects data were saved too - let's see the between subjects variance
% and how much our ERP reflects data
limo_add_plots('channel',49,fullfile(fileparts(pwd),'LIMO.mat'),'variable',[1:18],...
    {fullfile(pwd,'famous_Mean_of_Weighted mean.mat'),fullfile(pwd,'famous_single_subjects_Weighted mean.mat')})
title('Famous faces'); saveas(gcf, 'famous.fig'); close(gcf)
limo_add_plots('channel',49,fullfile(fileparts(pwd),'LIMO.mat'),'variable',[1:18],...
    {fullfile(pwd,'scrambled_Mean_of_Weighted mean.mat'),fullfile(pwd,'scrambled_single_subjects_Weighted mean.mat')})
title('Scrambled faces'); saveas(gcf, 'scrambled.fig'); close(gcf)
limo_add_plots('channel',49,fullfile(fileparts(pwd),'LIMO.mat'),'variable',[1:18],...
    {fullfile(pwd,'unfamiliar_Mean_of_Weighted mean.mat'),fullfile(pwd,'unfamiliar_single_subjects_Weighted mean.mat')})
title('Unfamiliar faces'); saveas(gcf, 'unfamiliar.fig'); close(gcf)

%% One sample t test (contrasting Full Faces vs Scrambled Faces at the subject level)
% https://github.com/LIMO-EEG-Toolbox/limo_meeg/wiki/5.-One-sample-t-test-(contrasting-Full-Faces-vs-Scrambled-Faces-at-the-subject-level)

% for each subject, we have a model with 9 conditions: famous 1st, 2nd,
% 3rd, scrambled 1st, 2nd, 3rd and unfamiliar 1st, 2nd, 3rd 
% --> we could use a contrast testing the interaction effect
% --> let's use limo_batch to do all the contrasts, creating
%     a contrast strucure to pass as argument in

cd(STUDY.filepath)
[~,~,contrast.LIMO_files] = limo_get_files([],[],[],...
    fullfile(STUDY.filepath,['LIMO_' STUDY.filename(1:end-6)],...
    'LIMO_files_FaceRepAll_GLM_Channels_Frequency_WLS.txt'));
contrast.mat = [0.5 0.5 0.5 -1 -1 -1 0.5 0.5 0.5];
limo_batch('contrast only',[],contrast);

% let's compute the one-sample t-test on this contrast 
chanlocs = [STUDY.filepath filesep 'limo_gp_level_chanlocs.mat'];
mkdir('one_sample'); cd('one_sample');
limo_random_select('one sample t-test',chanlocs,'LIMOfiles',...
    fullfile(STUDY.filepath,['LIMO_' STUDY.filename(1:end-6)],'con_4_files_FaceRepAll_GLM_Channels_Frequency_WLS.txt'),...
    'analysis_type','Full scalp analysis', 'type','Channels','nboot',101,'tfce',0);
limo_eeg(5,pwd)
limo_display_results(3,'one_sample_ttest_parameter_1.mat',pwd,0.05,2,...
    fullfile(pwd,'LIMO.mat'),0,'channels',49,'sumstats','mean'); % course plot
saveas(gcf, 'One_sample_spectrum.fig'); close(gcf)
    
%% Plot differences
% https://github.com/LIMO-EEG-Toolbox/limo_meeg/wiki/6.-Summary-statistics-of-differences

% since we computed the one sample t-test on differences (contrast) let's
% also make a plot of differences -- note that because we use a Bayesian
% confidence interval we test directly H1 ie that the difference is not 0;
% which differs from the one-sample t-test that tests the null, if the
% difference is 0 (an the significant effect tells you the you should
% reject the hypothesis that this is 0 - still does not prove it is!).

% trimmed mean of the contrast - same as one_sample_ttest_parameter_1 dim 1
Yr       = load('Yr.mat'); % this is the 1st level beta parameters 
chanlocs = [STUDY.filepath filesep 'limo_gp_level_chanlocs.mat'];
mkdir('average_betas'); cd('average_betas');
limo_central_tendency_and_ci(Yr.Yr, 'Trimmed mean',[],fullfile(pwd,'average.mat'))
limo_add_plots('channel',49,fullfile(fileparts(pwd),'LIMO.mat'),...
    {fullfile(pwd,'average_Trimmed_mean.mat')}); title('Average at channel 49')
saveas(gcf, 'Average.fig'); close(gcf); cd ..

% weighted mean for faces and scramble, followed by the difference
mkdir('ERPs'); cd('ERPs');
limo_central_tendency_and_ci(LIMOfiles, [1 2 3 7 8 9], chanlocs, 'Weighted mean', 'Mean', [],fullfile(pwd,'faces.mat'))  
limo_central_tendency_and_ci(LIMOfiles, [4 5 6], chanlocs, 'Weighted mean', 'Mean', [],fullfile(pwd,'scrambled.mat'))  
limo_plot_difference('faces_single_subjects_Weighted mean.mat','scrambled_single_subjects_Weighted mean',...
    'type','paired','name','faces_vs_scrambled','channel',49); % default 20% trimmed mean
saveas(gcf, 'Difference.fig'); close(gcf)

%% Two-ways Face * Repetition ANOVA
% https://github.com/LIMO-EEG-Toolbox/limo_meeg/wiki/7.-Two-way-ANOVA-(Faces-x-Repetition)

cd(STUDY.filepath)
chanlocs = [STUDY.filepath filesep 'limo_gp_level_chanlocs.mat'];
mkdir('Face-Repetition_ANOVA');cd('Face-Repetition_ANOVA')
LIMOPath = limo_random_select('Repeated Measures ANOVA',chanlocs,'LIMOfiles',...
    fullfile(STUDY.filepath,['LIMO_' STUDY.filename(1:end-6)],'Beta_files_FaceRepAll_GLM_Channels_Frequency_WLS.txt'),...
    'analysis_type','Full scalp analysis','parameters',{[1 2 3],[4 5 6],[7 8 9]},...
    'factor names',{'face','repetition'},'type','Channels','nboot',1000,'tfce',0,'skip design check','yes');

% add contrast famous>unfamiliar
limo_contrast(fullfile(pwd,'Yr.mat'),fullfile(pwd,'LIMO.mat'), 3 ,[1 1 1 0 0 0 -1 -1 -1]); % compute a new contrast
limo_contrast(fullfile(pwd,'Yr.mat'),fullfile(pwd,'LIMO.mat'), 4);                         % do the bootstrap - although here there is no effect anyway

% figures
limo_eeg(5,pwd)

%% paired t-test famous vs. unfamiliar controling for scrambled
% https://github.com/LIMO-EEG-Toolbox/limo_meeg/wiki/8.-Paired-t-test-(Famous-vs-Unfamiliar)

% let's say the research question is familiar vs unfamiliar and scrambled
% are just a control - doing the ANOVA is a little meaningless because you
% know include scrambled as a condition when in fact it's a control - using
% contrast we can buld that in

cd(STUDY.filepath)
[~,~,contrast.LIMO_files] = limo_get_files([],[],[],...
    fullfile(STUDY.filepath,['LIMO_' STUDY.filename(1:end-6)],...
    'LIMO_files_FaceRepAll_GLM_Channels_Frequency_WLS.txt'));
contrast.mat = [1 1 1 -1 -1 -1 0 0 0 ; 0 0 0 -1 -1 -1 1 1 1];
limo_batch('contrast only',[],contrast);

% note here con5 and con6 because in previous steps of the tutorial we have
% contrasts 1,2,3 from design and contrast 4 from the interaction effect
mkdir('Paired_ttest'); cd('Paired_ttest');
chanlocs = [STUDY.filepath filesep 'limo_gp_level_chanlocs.mat'];
files = {fullfile(STUDY.filepath,['LIMO_' STUDY.filename(1:end-6)],'con_5_files_FaceRepAll_GLM_Channels_Frequency_WLS.txt'), ...
    fullfile(STUDY.filepath,['LIMO_' STUDY.filename(1:end-6)],'con_6_files_FaceRepAll_GLM_Channels_Frequency_WLS.txt')};
limo_random_select('paired t-test',chanlocs,'LIMOfiles',files,...
    'analysis_type','Full scalp analysis', 'type','Channels','nboot',1000,'tfce',0);

% ---------------------------------------------------
%% let's see designs with groups
% ----------------------------------------------------
% update STUDY using median split + unknow age group
cd(STUDY.filepath)
[STUDY ALLEEG] = std_editset( STUDY, ALLEEG, 'commands',{{'index',2,'group','1'}, ...
    {'index',7,'group','1'},{'index',14,'group','1'},{'index',15,'group','1'}, ...
    {'index',16,'group','1'},{'index',17,'group','1'},{'index',1,'group','2'}, ...
    {'index',4,'group','2'},{'index',8,'group','2'},{'index',9,'group','2'}, ...
    {'index',10,'group','2'},{'index',11,'group','2'},{'index',13,'group','2'}, ...
    {'index',3,'group','3'},{'index',5,'group','3'},{'index',6,'group','3'}, ...
    {'index',12,'group','3'}, {'index',18,'group','3'}}, 'updatedat','off','rmclust','on');
[STUDY, EEG] = pop_savestudy( STUDY, EEG, 'savemode','resave');

%% Between subjects ANOVA (ie groups) with repeated measures
% https://github.com/LIMO-EEG-Toolbox/limo_meeg/wiki/9.-Between-subjects%E2%80%99-ANOVAs-with-repeated-factors

% 1st level is exactly the same as before, but because STUDY now has
% groups, LIMO will create additional text files, splitting by groups

% 'recompute' 1st level calling pop_limo (since LIMO keeps tracks of computation
% with psom, it actually computes nothing and simply creates the additional txt files) 

STUDY.currentdesign = 1; % indicates we go back to the 1st design
STUDY  = pop_limo(STUDY, ALLEEG, 'method','WLS','measure','datspec','freqlim',[3 45],'erase','on','splitreg','off','interaction','off');

% 2nd level
cd(STUDY.filepath)
mkdir('Gp-Conditions_ANOVA'); cd('Gp-Conditions_ANOVA');
chanlocs = [STUDY.filepath filesep 'limo_gp_level_chanlocs.mat'];
Files = cell(2,1); % groups in rows, repeated measures in columns  
Files{1} = fullfile(STUDY.filepath,['LIMO_' STUDY.filename(1:end-6)],'Beta_files_Gp1_ANOVA_Faces_GLM_Channels_Frequency_WLS.txt');  
Files{2} = fullfile(STUDY.filepath,['LIMO_' STUDY.filename(1:end-6)],'Beta_files_Gp2_ANOVA_Faces_GLM_Channels_Frequency_WLS.txt');  
limo_random_select('Repeated Measures ANOVA',chanlocs,'LIMOfiles',Files,...  
    'analysis_type','Full scalp analysis','parameters',{[1 2 3];[1 2 3]},...
    'factor names',{'face'},'type','Channels','nboot',1000,'tfce',0,'skip design check','yes');


%% Two samples t-test on famous faces
% https://github.com/LIMO-EEG-Toolbox/limo_meeg/wiki/10.-Two-sample-t-tests

cd(STUDY.filepath)
mkdir('Two-samples_ttest'); cd('Two-samples_ttest');
chanlocs = [STUDY.filepath filesep 'limo_gp_level_chanlocs.mat'];
Files{1} = fullfile(STUDY.filepath,['LIMO_' STUDY.filename(1:end-6)],'Beta_files_Gp1_ANOVA_Faces_GLM_Channels_Frequency_WLS.txt');  
Files{2} = fullfile(STUDY.filepath,['LIMO_' STUDY.filename(1:end-6)],'Beta_files_Gp2_ANOVA_Faces_GLM_Channels_Frequency_WLS.txt');  
limo_random_select('two-samples t-test',chanlocs,'LIMOfiles',Files,...  
   'analysis_type','Full scalp analysis', 'type','Channels','parameter',[1;1],'nboot',1000,'tfce',0);  

% --------------------------------------------------------
%% Let's see how to process data with contiunous variables
% --------------------------------------------------------

%% Regression over subjects
% https://github.com/LIMO-EEG-Toolbox/limo_meeg/wiki/11.-Regression-among-subjects  

% here we take the contrast 4 (faces>scrambled) and ask if there is an age effect
% regressing subjects' age onto the contrast

cd(STUDY.filepath)
mkdir('Regression'); cd('Regression');
chanlocs = [STUDY.filepath filesep 'limo_gp_level_chanlocs.mat'];
age_regressor = cell2mat(arrayfun(@(x) x.age,STUDY.datasetinfo,'UniformOutput',false))';  
save age_regressor age_regressor
limo_random_select('regression',chanlocs,'LIMOfiles',...  
    fullfile(STUDY.filepath,['LIMO_' STUDY.filename(1:end-6)],'con_4_files_FaceRepAll_GLM_Channels_Frequency_WLS.txt'),...
    'analysis type','Full scalp analysis', 'type','Channels','nboot',0,'tfce',0,'regressor',...
    fullfile(pwd,'age_regressor.mat'), 'zscore','Yes','skip design check','Yes') ;

% figures
limo_eeg(5,pwd)
limo_display_results(3,'Covariate_effect_1.mat',pwd,0.05,1,...
    fullfile(pwd,'LIMO.mat'),0,'channels',56,'plot3type','Modelled'); % course plot
saveas(gcf, 'contrast_spectrum.fig'); close(gcf)


%% Regression over trials
% https://github.com/LIMO-EEG-Toolbox/limo_meeg/wiki/12.-Regression-at-the-trial-level

% we use the time between each repetition of the same stimulus as continuous variable
% for a given subject we have 3 conditions (familiar faces, unfamiliar faces and scrambled faces)
% and one continuous variable (the distance between the repeat of a
% stimulus type) that we split per condition

% 1st level
STUDY = std_makedesign(STUDY, ALLEEG, 2, 'name','Face_time','delfiles','off','defaultdesign','off',...
    'variable1','face_type','values1',{'famous','scrambled','unfamiliar'},'vartype1','categorical',...
    'variable2','time_dist','values2',[],'vartype2','continuous',...
    'subjselect',{'sub-002','sub-003','sub-004','sub-005','sub-006','sub-007','sub-008','sub-009','sub-010','sub-011','sub-012','sub-013','sub-014','sub-015','sub-016','sub-017','sub-018','sub-019'});
[STUDY, EEG] = pop_savestudy( STUDY, EEG, 'savemode','resave');
% note 'splitreg' is 'on'
STUDY = pop_limo(STUDY, ALLEEG, 'method','WLS','measure','datspec','freqlim',[3 45],'erase','on','splitreg','on','interaction','off');

% 2nd level 
% now we can do a repeated measure ANOVA, but instead of the 3 conditions,
% it's between the 3 'time effect' regressors 

chanlocs   = [STUDY.filepath filesep 'limo_gp_level_chanlocs.mat'];
cd(STUDY.filepath)
mkdir([STUDY.filepath filesep '1-way-ANOVA-time'])
cd([STUDY.filepath filesep '1-way-ANOVA-time'])
limo_random_select('Repeated Measures ANOVA',chanlocs,'LIMOfiles',...
    fullfile(STUDY.filepath,['LIMO_' STUDY.filename(1:end-6)],'Beta_files_Face_time_GLM_Channels_Frequency_WLS.txt'), ...
    'analysis_type','Full scalp analysis','parameters',{[4 5 6]},...
    'factor names',{'face'},'type','Channels','nboot',1000,'tfce',0,'skip design check','yes');

% add contrast famous > unfamiliar
% = the effect of time bettwen stimuli that repeat is stronger for famous than others
limo_contrast(fullfile(pwd,'Yr.mat'),fullfile(pwd,'LIMO.mat'), 3 ,[1 0 -1]); % compute a new contrast
limo_contrast(fullfile(pwd,'Yr.mat'),fullfile(pwd,'LIMO.mat'), 4);            % do the bootstrap of the last contrast
limo_eeg(5,pwd) % channel*freq imagesc for both effects and contrast
