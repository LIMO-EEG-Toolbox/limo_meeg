% Data driven analysis to validate the Weighted Least Square minimization method
% developped in LIMO EEG. There is only three thing to do to run the script
% and fully reproduce the analysis:
% 1 - download the data @ https://openneuro.org/datasets/ds002718/versions/1.0.4
% 2 - install the tools (eeglab, limo eeg and rst_toolbox
%     https://github.com/CPernet/Robust_Statistical_Toolbox)
% 3 - define where are those downloaded data line 12
% 
% The script starts by preprocessing the data as decribed on Pernet et al.
% 2021 (https://www.frontiersin.org/articles/10.3389/fnins.2020.610388/full),
% after which 1st level statistical modelling is performed using 9
% conditions (3 types of faces * 3 repetition levels) and parameters are
% computed using OLS, WLS and IRLS.
% - high and low weight trials are then compared
% - subjects' data are (super)bootrapped under the null for OLS,WLS and IRLS to
%   determine the type 1 FWER
% - group level repeated measures ANOVA are computed (Hotelling)
% - distributions are compared: median of mean squares effects, median of
%   mean squares erros, Brown Forsythe variance test - distributions of F and
%   TFCE values are then exported as csv files, and analyzed in R (rogme package)
%
% Script created by Cyril R Pernet

% %% Indicate here were are the open neuro data
% studypath = 'F:\WakemanHenson_Faces\eeg';

%%  quick check
if ~isfolder(studypath)
    error('invalid study path')
end

if ~exist('eeglab.m','file')
    error('EEGLAB is not in your path')
end

if ~exist('pop_importbids.m','file')
    error('EEGLAB BIDS import tools are not your path')
end

if ~exist('limo_eeg.m','file')
    error('LIMO tools are not in your path')
end

if ~exist('rst_shifthd.m','file')
    warning on
    warning('RST toolbox required https://github.com/CPernet/Robust_Statistical_Toolbox')
end

%% get the data
clear variables
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;
[STUDY, ALLEEG]                   = pop_importbids(studypath, 'bidsevent','on','bidschanloc','on', 'eventtype', 'trial_type',,'studyName','Face_detection');
ALLEEG                            = pop_select( ALLEEG, 'nochannel',{'061','062','063','064'});
CURRENTSTUDY                      = 1; 
EEG                               = ALLEEG; 
CURRENTSET                        = 1:length(EEG);

%% Preprocessing
% this reproduces Pernet et al. 2021 Front. Neurosci., 11
% https://www.frontiersin.org/articles/10.3389/fnins.2020.610388/full

% Clean data - just the bad channels
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

%% Create study design and run the HLM with each method
STUDY  = std_checkset(STUDY, ALLEEG);
STUDY  = std_makedesign(STUDY, EEG, 1, 'name','STUDY.FaceRepetition','delfiles','off','defaultdesign','off','variable1','type','values1',{});
eeglab redraw

% Precompute ERP with baseline correction [-200 0]
[STUDY, EEG] = std_precomp(STUDY, EEG, {}, 'savetrials','on','interp','on','recompute','on',...
    'erp','on','erpparams', {'rmbase' [-200 0]}, 'spec','on',...
    'ersp','on','itc','on', 'specparams',{'specmode','fft','logtrials','off'});
eeglab redraw

% 1st level Analysis for each method (IRLS takes a long time)
STUDY = pop_limo(STUDY, ALLEEG, 'method','OLS','measure','daterp','timelim',[-50 650],'erase','on','splitreg','off','interaction','off');
STUDY = pop_limo(STUDY, ALLEEG, 'method','WLS','measure','daterp','timelim',[-50 650],'erase','on','splitreg','off','interaction','off');
STUDY = pop_limo(STUDY, ALLEEG, 'method','IRLS','measure','daterp','timelim',[-50 650],'erase','on','splitreg','off','interaction','off');
eeglab redraw

%% investigate weights 
% is there a consistent difference between good and bad trials among subjects? 
% show the strongest differences and topographies per subject using
% limo_ttest and bootstrapping this for inference - show side by side, raw,
% thresholded + topography + max diff + individual means
chanlocs = fullfile(studypath,['derivatives' filesep 'limo_gp_level_chanlocs.mat']);
WLS      = [studypath filesep 'derivatives' filesep 'LIMO_Face_detection' filesep 'LIMO_files_FaceRepAll_GLM_Channels_Time_WLS.txt'];
cd([studypath filesep 'derivatives' filesep 'ANOVA_WLS'])
limo_CheckWeight(WLS,chanlocs,'TestDifference','on','SingleSubjectAnalysis','on','PlotRank','off','CheckBias','off')
cd([studypath filesep 'derivatives' filesep 'ANOVA_WLS' filesep 'Weights_checking' filesep 'single_subjects'])
load maxchannels; maxchannels_wls = maxchannels; 
sub = dir('sub*'); figure('Name','WLS subjects''s high vs. low trials'); 
for s=1:length(sub)
   cd(sub(s).name); subplot(9,2,s);
   load Y1r; plot(mean(squeeze(Y1r(maxchannels_wls(s),:,:)),2),'r');
   hold on; grid on; box on; axis tight
   load Y2r; plot(mean(squeeze(Y2r(maxchannels_wls(s),:,:)),2),'k');
   set(gca,'Xticklabel',[]); ylabel(num2str(maxchannels_wls(s))); index = index + 2;   
   cd ..
end

% check metrics
for s=1:length(sub)-1
   cd(sub(s).name); 
   filename = dir('metrics*');
   t = readtable(filename.name);
   [h1(s),M1(s),CI1(:,s),p1(s)] = rst_2ttest(t.time_var_lt,t.time_var_ht,'estimator','trimmed mean','figure','off','newfig','no');
   [h2(s),M2(s),CI2(:,s),p2(s)] = rst_2ttest(t.power_lt,t.power_ht,'estimator','trimmed mean','figure','off','newfig','no');
   [h3(s),M3(s),CI3(:,s),p3(s)] = rst_2ttest(t.autocorr_lt,t.autocorr_ht,'estimator','trimmed mean','figure','off','newfig','no');
   indices{s}={{t.low_weight_trials},{t.high_weight_trials}}; cd ..
end

C = limo_color_images(18);
metric_names = {'average tSNR', 'power','autocorrelation'};
for metric = 1:3
    figure('Name',metric_names{metric});subplot(1,4,1); 
    t = readtable('mean_metrics_maxchannels.csv');
    for s=1:18
        if metric == 1
            plot([t.time_var_lt(s) t.time_var_ht(s)],'o','MarkerFaceColor',C(s,:));hold on;
            plot([t.time_var_lt(s) t.time_var_ht(s)],'color',C(s,:));
        elseif metric == 2
            plot([t.power_lt(s) t.power_ht(s)],'o','MarkerFaceColor',C(s,:));hold on;
            plot([t.power_lt(s) t.power_ht(s)],'color',C(s,:));
        elseif metric == 3
            plot([t.autocorr_lt(s) t.autocorr_ht(s)],'o','MarkerFaceColor',C(s,:));hold on;
            plot([t.autocorr_lt(s) t.autocorr_ht(s)],'color',C(s,:));
       end
    end
    title(sprintf('%s\nmax channels',metric_names{metric}),'FontSize',12);
    ylabel(metric_names{metric},'FontSize',10);
    xticklabels({'red','','black'}); grid on
    if metric == 1
        axis([0.5 2.5 4 24]);
        subplot(1,4,2);
        [~,~,CI,p]= rst_2ttest(t.time_var_lt,t.time_var_ht,'estimator','trimmed mean','figure','on','newfig','no');
        axis([0.7 1.3 -3 7]); title(sprintf('diff [%g %g]\np=%g',CI(1),CI(2),p),'FontSize',12); xlabel('');
        ylabel(sprintf('%s difference',metric_names{metric}),'FontSize',10);
    elseif metric == 2
        axis([0.5 2.5 0 510]);
        subplot(1,4,2);
        [~,~,CI,p]= rst_2ttest(t.time_var_lt,t.time_var_ht,'estimator','trimmed mean','figure','on','newfig','no');
        axis([0.7 1.3 -1 9]); ylabel(sprintf('%s difference',metric_names{metric}),'FontSize',10);
    elseif metric == 3
        axis([0.5 2.5 19 65]);
        subplot(1,4,2);
        [~,~,CI,p]= rst_2ttest(t.time_var_lt,t.time_var_ht,'estimator','trimmed mean','figure','on','newfig','no');
        axis([0.7 1.3 -4 10]); ylabel(sprintf('%s difference',metric_names{metric}),'FontSize',10);
    end
    title(sprintf('diff [%g %g]\np=%g',CI(1),CI(2),p),'FontSize',12); xlabel('');
    
   t = readtable('mean_metrics_allchannels.csv');
   for s=1:18
       if metric == 1
           subplot(1,4,3);
           plot([t.time_var_lt(s) t.time_var_ht(s)],'o','MarkerFaceColor',C(s,:));hold on;
           plot([t.time_var_lt(s) t.time_var_ht(s)],'color',C(s,:));
       elseif metric == 2
           subplot(1,4,3);
           plot([t.power_lt(s) t.power_ht(s)],'o','MarkerFaceColor',C(s,:));hold on;
           plot([t.power_lt(s) t.power_ht(s)],'color',C(s,:));
       elseif metric == 3
           subplot(1,4,3);
           plot([t.autocorr_lt(s) t.autocorr_ht(s)],'o','MarkerFaceColor',C(s,:));hold on;
           plot([t.autocorr_lt(s) t.autocorr_ht(s)],'color',C(s,:));
       end
   end
       title(sprintf('%s\nall channels',metric_names{metric}),'FontSize',12);
       ylabel(metric_names{metric},'FontSize',10);
       xticklabels({'red','black'}); grid on
   if metric == 1
       axis([0.5 2.5 3 9]);
       subplot(1,4,4);
       [~,~,CI,p]= rst_2ttest(t.time_var_lt,t.time_var_ht,'estimator','trimmed mean','figure','on','newfig','no');
       axis([0.7 1.3 -2 3]); ylabel(sprintf('%s difference',metric_names{metric}),'FontSize',10);
   elseif metric == 2
       axis([0.5 2.5 0 80]);
       subplot(1,4,4);
       [~,~,CI,p]= rst_2ttest(t.time_var_lt,t.time_var_ht,'estimator','trimmed mean','figure','on','newfig','no');
       axis([0.7 1.3 -1.5 3]); ylabel(sprintf('%s difference',metric_names{metric}),'FontSize',10);
   elseif metric == 3
       axis([0.5 2.5 35 85]);
       subplot(1,4,4);
       [~,~,CI,p]= rst_2ttest(t.time_var_lt,t.time_var_ht,'estimator','trimmed mean','figure','on','newfig','no');
       axis([0.7 1.3 -1.5 2.5]); ylabel(sprintf('%s difference',metric_names{metric}),'FontSize',10);
   end
   title(sprintf('diff [%g %g]\np=%g',CI(1),CI(2),p),'FontSize',12); xlabel('');
end

figure('Name','High vs Low trial spectra');
for s=1:length(sub)-1
    cd(sub(s).name); subplot(9,2,s)
    load Y1r; [spectra,freqs] = spectopo(Y1r, 176, 250,'plot','off');
    plot(spectra(maxchannels(s),:),'r'); hold on
    load Y2r; [spectra,freqs] = spectopo(Y2r, 176, 250,'plot','off');
    plot(spectra(maxchannels(s),:),'k');
    set(gca,'Xticklabel',[]); ylabel(num2str(maxchannels(s))); index = index + 2;
    cd ..
end


% ---------- redo for IRLS -------------------

IRLS = [studypath filesep 'derivatives' filesep 'LIMO_Face_detection' filesep 'LIMO_files_FaceRepAll_GLM_Channels_Time_IRLS.txt'];
cd([studypath filesep 'derivatives' filesep 'ANOVA_IRLS'])
limo_CheckWeight(IRLS,chanlocs,'TestDifference','on','SingleSubjectAnalysis','on','PlotRank','off','CheckBias','off')
cd([studypath filesep 'derivatives' filesep 'ANOVA_IRLS' filesep 'Weights_checking' filesep 'single_subjects'])
load maxchannels; maxchannels_irls = maxchannels;
common_channels = maxchannels_wls(find(maxchannels_wls' == maxchannels_irls'));
sub = dir('sub*'); figure('Name','IRLS subjects''s high vs. low trials'); index = 1;
for s=1:length(sub)
   cd(sub(s).name); subplot(9,2,s);
   load Y1r; plot(mean(squeeze(Y1r(maxchannels_irls(s),:,:)),2),'r');
   hold on; grid on; box on; axis tight
   load Y2r; plot(mean(squeeze(Y2r(maxchannels_irls(s),:,:)),2),'k');
   set(gca,'Xticklabel',[]); ylabel(num2str(maxchannels_irls(s)));
   cd ..
end

% check metrics
common_trials = NaN(length(sub)-1,2);
for s=1:length(sub)-1
   cd(sub(s).name); 
   filename = dir('metrics*');
   t = readtable(filename.name);
   [h1(s),M1(s),CI1(:,s),p1(s)] = rst_2ttest(t.time_var_lt,t.time_var_ht,'estimator','trimmed mean','figure','off','newfig','no');
   [h2(s),M2(s),CI2(:,s),p2(s)] = rst_2ttest(t.power_lt,t.power_ht,'estimator','trimmed mean','figure','off','newfig','no');
   [h3(s),M3(s),CI3(:,s),p3(s)] = rst_2ttest(t.autocorr_lt,t.autocorr_ht,'estimator','trimmed mean','figure','off','newfig','no');
   if any(maxchannels_irls(s) == common_channels)
       common(s,1) = length(intersect(t.low_weight_trials,cell2mat(indices{s}{1})))/min(length(t.low_weight_trials),length(cell2mat(indices{s}{1})));
       common(s,2) = length(intersect(t.high_weight_trials,cell2mat(indices{s}{2})))/min(length(t.low_weight_trials),length(cell2mat(indices{s}{2})));
   end
   cd ..
end


C = limo_color_images(18);
metric_names = {'average tSNR', 'power','autocorrelation'};
for metric = 1:3
    figure('Name',metric_names{metric});subplot(1,4,1); 
    t = readtable('mean_metrics_maxchannels.csv');
    for s=1:18
        if metric == 1
            plot([t.time_var_lt(s) t.time_var_ht(s)],'o','MarkerFaceColor',C(s,:));hold on;
            plot([t.time_var_lt(s) t.time_var_ht(s)],'color',C(s,:));
        elseif metric == 2
            plot([t.power_lt(s) t.power_ht(s)],'o','MarkerFaceColor',C(s,:));hold on;
            plot([t.power_lt(s) t.power_ht(s)],'color',C(s,:));
        elseif metric == 3
            plot([t.autocorr_lt(s) t.autocorr_ht(s)],'o','MarkerFaceColor',C(s,:));hold on;
            plot([t.autocorr_lt(s) t.autocorr_ht(s)],'color',C(s,:));
       end
    end
    title(sprintf('%s\nmax channels',metric_names{metric}),'FontSize',12);
    ylabel(metric_names{metric},'FontSize',10);
    xticklabels({'red','','black'}); grid on
    if metric == 1
        axis([0.5 2.5 4 24]);
        subplot(1,4,2);
        [~,~,CI,p]= rst_2ttest(t.time_var_lt,t.time_var_ht,'estimator','trimmed mean','figure','on','newfig','no');
        axis([0.7 1.3 -3 7]); title(sprintf('diff [%g %g]\np=%g',CI(1),CI(2),p),'FontSize',12); xlabel('');
        ylabel(sprintf('%s difference',metric_names{metric}),'FontSize',10);
    elseif metric == 2
        axis([0.5 2.5 0 510]);
        subplot(1,4,2);
        [~,~,CI,p]= rst_2ttest(t.time_var_lt,t.time_var_ht,'estimator','trimmed mean','figure','on','newfig','no');
        axis([0.7 1.3 -1 9]); ylabel(sprintf('%s difference',metric_names{metric}),'FontSize',10);
    elseif metric == 3
        axis([0.5 2.5 19 65]);
        subplot(1,4,2);
        [~,~,CI,p]= rst_2ttest(t.time_var_lt,t.time_var_ht,'estimator','trimmed mean','figure','on','newfig','no');
        axis([0.7 1.3 -4 10]); ylabel(sprintf('%s difference',metric_names{metric}),'FontSize',10);
    end
    title(sprintf('diff [%g %g]\np=%g',CI(1),CI(2),p),'FontSize',12); xlabel('');
    
   t = readtable('mean_metrics_allchannels.csv');
   for s=1:18
       if metric == 1
           subplot(1,4,3);
           plot([t.time_var_lt(s) t.time_var_ht(s)],'o','MarkerFaceColor',C(s,:));hold on;
           plot([t.time_var_lt(s) t.time_var_ht(s)],'color',C(s,:));
       elseif metric == 2
           subplot(1,4,3);
           plot([t.power_lt(s) t.power_ht(s)],'o','MarkerFaceColor',C(s,:));hold on;
           plot([t.power_lt(s) t.power_ht(s)],'color',C(s,:));
       elseif metric == 3
           subplot(1,4,3);
           plot([t.autocorr_lt(s) t.autocorr_ht(s)],'o','MarkerFaceColor',C(s,:));hold on;
           plot([t.autocorr_lt(s) t.autocorr_ht(s)],'color',C(s,:));
       end
   end
       title(sprintf('%s\nall channels',metric_names{metric}),'FontSize',12);
       ylabel(metric_names{metric},'FontSize',10);
       xticklabels({'red','black'}); grid on
   if metric == 1
       axis([0.5 2.5 3 9]);
       subplot(1,4,4);
       [~,~,CI,p]= rst_2ttest(t.time_var_lt,t.time_var_ht,'estimator','trimmed mean','figure','on','newfig','no');
       axis([0.7 1.3 -2 3]); ylabel(sprintf('%s difference',metric_names{metric}),'FontSize',10);
   elseif metric == 2
       axis([0.5 2.5 0 80]);
       subplot(1,4,4);
       [~,~,CI,p]= rst_2ttest(t.time_var_lt,t.time_var_ht,'estimator','trimmed mean','figure','on','newfig','no');
       axis([0.7 1.3 -1.5 3]); ylabel(sprintf('%s difference',metric_names{metric}),'FontSize',10);
   elseif metric == 3
       axis([0.5 2.5 35 85]);
       subplot(1,4,4);
       [~,~,CI,p]= rst_2ttest(t.time_var_lt,t.time_var_ht,'estimator','trimmed mean','figure','on','newfig','no');
       axis([0.7 1.3 -1.5 2.5]); ylabel(sprintf('%s difference',metric_names{metric}),'FontSize',10);
   end
   title(sprintf('diff [%g %g]\np=%g',CI(1),CI(2),p),'FontSize',12); xlabel('');
end

figure('Name','High vs Low trial spectra');
for s=1:length(sub)-1
    cd(sub(s).name); subplot(9,2,s)
    load Y1r; [spectra,freqs] = spectopo(Y1r, 176, 250,'plot','off');
    plot(spectra(maxchannels(s),:),'r'); hold on
    load Y2r; [spectra,freqs] = spectopo(Y2r, 176, 250,'plot','off');
    plot(spectra(maxchannels(s),:),'k');
    set(gca,'Xticklabel',[]); ylabel(num2str(maxchannels(s))); index = index + 2;
    cd ..
end

%% bootstrap subjects (under the null) to check type 1 error rates
% this takes a very long time -- not all of this is oresented in the
% article but was needed to check the type 1 error is well controlled in
% any cases (model/contrast - OLS/WLS/IRLS)

for s=1:size(STUDY.datasetinfo,2)
    fprintf('running analysis on subject %g\n',s)
    LIMO = load(fullfile(STUDY.datasetinfo(s).filepath,['FaceRepAll_GLM_Channels_Time_OLS' filesep 'LIMO.mat']));
    LIMO = LIMO.LIMO; LIMO.design.bootstrap = 2500; LIMO.design.status = 'to do';
    save(fullfile(LIMO.dir,'LIMO.mat'),'LIMO'); try rmdir([LIMO.dir filesep 'H0'],'s'); end
    limo_eeg(4,LIMO)
    LIMO = load(fullfile(STUDY.datasetinfo(s).filepath,['FaceRepAll_GLM_Channels_Time_WLS' filesep 'LIMO.mat']));
    LIMO = LIMO.LIMO; LIMO.design.bootstrap = 2500; LIMO.design.status = 'to do';
    save(fullfile(LIMO.dir,'LIMO.mat'),'LIMO'); try rmdir([LIMO.dir filesep 'H0'],'s'); end
    limo_eeg(4,LIMO)
    LIMO = load(fullfile(STUDY.datasetinfo(s).filepath,['FaceRepAll_GLM_Channels_Time_IRLS' filesep 'LIMO.mat']));
    LIMO = LIMO.LIMO; LIMO.design.bootstrap = 2500; LIMO.design.status = 'to do';
    save(fullfile(LIMO.dir,'LIMO.mat'),'LIMO'); try rmdir([LIMO.dir filesep 'H0'],'s'); end
    limo_eeg(4,LIMO)
end

if ~exist('limo_test_glmboot.m','file')
    warning on; warning('get limo validation tools to check type 1 error')
else
    for s=1:size(STUDY.datasetinfo,2)
        H0o{s}  = fullfile(STUDY.datasetinfo(s).filepath,['FaceRepAll_GLM_Channels_Time_OLS' filesep 'H0']);
        H0w{s}  = fullfile(STUDY.datasetinfo(s).filepath,['FaceRepAll_GLM_Channels_Time_WLS' filesep 'H0']);
        H0iw{s} = fullfile(STUDY.datasetinfo(s).filepath,['FaceRepAll_GLM_Channels_Time_IRLS' filesep 'H0']);
    end
end

% use 1500 for the null and 1000 for testing
[errOLS,~,maxOLS,~,clusterOLS]     = limo_test_glmboot(chanlocs,H0o, 'step_size',300,'Nboot',1000,'MinSamp',300);
[errWLS,~,maxWLS,~,clusterWLS]     = limo_test_glmboot(chanlocs,H0w, 'step_size',300,'Nboot',1000,'MinSamp',300);
[errIRLS,~,maxIRLS,~,clusterIRLS]  = limo_test_glmboot(chanlocs,H0iw,'step_size',300,'Nboot',1000,'MinSamp',300);
OLS     = cell2mat(errOLS); OLS = OLS(:,end); % pick 1500 bootstraps
tmp = cell2mat(maxOLS); OLS = [OLS tmp(:,end)];
tmp = cell2mat(clusterOLS); OLS = [OLS tmp(:,end)];
WLS     = cell2mat(errWLS); WLS = WLS(:,end);
tmp = cell2mat(maxWLS); WLS = [WLS tmp(:,end)];
tmp = cell2mat(clusterWLS); WLS = [WLS tmp(:,end)];
IRLS    = cell2mat(errIRLS); IRLS = IRLS(:,end);
[h,CI,p] = rst_pttest(OLS,WLS);

% additional validation for contrasts
for s=1:size(STUDY.datasetinfo,2)
    warning('running analysis on subject %g\n',s)
    cd(fullfile(STUDY.datasetinfo(s).filepath,'FaceRepAll_GLM_Channels_Time_OLS'))
    load LIMO; try LIMO = rmfield(LIMO,'contrast'); save LIMO LIMO; end
    delete('con_*'); delete(['H0' filesep 'H0_con_*'])
    % limo_contrast(Y, Betas, LIMO, contrast type, analysis type ,contrast)
    limo_contrast(fullfile(STUDY.datasetinfo(s).filepath,['FaceRepAll_GLM_Channels_Time_OLS' filesep 'Yr.mat']), ...
        fullfile(STUDY.datasetinfo(s).filepath,['FaceRepAll_GLM_Channels_Time_OLS' filesep 'Betas.mat']), ...
        fullfile(STUDY.datasetinfo(s).filepath,['FaceRepAll_GLM_Channels_Time_OLS' filesep 'LIMO.mat']), ...
        'T', 1, [1 -1 0 1 -1 0 1 -1 0 0]);
    limo_contrast(fullfile(STUDY.datasetinfo(s).filepath,['FaceRepAll_GLM_Channels_Time_OLS' filesep 'Yr.mat']), ...
        fullfile(STUDY.datasetinfo(s).filepath,['FaceRepAll_GLM_Channels_Time_OLS' filesep 'H0' filesep 'H0_Betas.mat']), ...
        fullfile(STUDY.datasetinfo(s).filepath,['FaceRepAll_GLM_Channels_Time_OLS' filesep 'LIMO.mat']), ...
        'T', 2);
    
    cd(fullfile(STUDY.datasetinfo(s).filepath,'FaceRepAll_GLM_Channels_Time_WLS'))
    load LIMO; try LIMO = rmfield(LIMO,'contrast'); save LIMO LIMO; end
    delete('con*'); delete(['H0' filesep 'H0_con*'])
    limo_contrast(fullfile(STUDY.datasetinfo(s).filepath,['FaceRepAll_GLM_Channels_Time_WLS' filesep 'Yr.mat']), ...
        fullfile(STUDY.datasetinfo(s).filepath,['FaceRepAll_GLM_Channels_Time_WLS' filesep 'Betas.mat']), ...
        fullfile(STUDY.datasetinfo(s).filepath,['FaceRepAll_GLM_Channels_Time_WLS' filesep 'LIMO.mat']), ...
        'T', 1, [1 -1 0 1 -1 0 1 -1 0 0]);
    limo_contrast(fullfile(STUDY.datasetinfo(s).filepath,['FaceRepAll_GLM_Channels_Time_WLS' filesep 'Yr.mat']), ...
        fullfile(STUDY.datasetinfo(s).filepath,['FaceRepAll_GLM_Channels_Time_WLS' filesep 'H0' filesep 'H0_Betas.mat']), ...
        fullfile(STUDY.datasetinfo(s).filepath,['FaceRepAll_GLM_Channels_Time_WLS' filesep 'LIMO.mat']), ...
        'T', 2);    
end

%% 2nd level analysis - ANOVAs
chanlocs = fullfile(studypath,['derivatives' filesep 'limo_gp_level_chanlocs.mat']);
mkdir([studypath filesep 'derivatives' filesep 'ANOVA_OLS'])
cd([studypath filesep derivatives' filesep 'ANOVA_OLS'])
limo_random_select('Repeated Measures ANOVA',chanlocs,'LIMOfiles',...
    {[studypath filesep 'derivatives' filesep 'LIMO_Face_detection' filesep 'Beta_files_FaceRepetition_GLM_Channels_Time_OLS.txt']},...
    'analysis_type','Full scalp analysis','parameters',{[1 2 3],[4 5 6],[7 8 9]},...
    'factor names',{'face','repetition'},'type','Channels','nboot',1000,'tfce',0,'skip design check','yes');
limo_eeg(5,pwd); % to make the channels*time figure

mkdir([studypath filesep 'derivatives' filesep 'ANOVA_WLS'])
cd([studypath filesep 'derivatives' filesep 'ANOVA_WLS'])
limo_random_select('Repeated Measures ANOVA',chanlocs,'LIMOfiles',...
    {[studypath filesep 'derivatives' filesep 'LIMO_Face_detection' filesep 'Beta_files_FaceRepetition_GLM_Channels_Time_WLS.txt']},...
    'analysis_type','Full scalp analysis','parameters',{[1 2 3],[4 5 6],[7 8 9]},...
    'factor names',{'face','repetition'},'type','Channels','nboot',1000,'tfce',0,'skip design check','yes');
limo_eeg(5,pwd)

mkdir([studypath filesep 'derivatives' filesep 'ANOVA_IRLS'])
cd([studypath filesep 'derivatives' filesep 'ANOVA_IRLS'])
limo_random_select('Repeated Measures ANOVA',chanlocs,'LIMOfiles',...
    {[studypath filesep 'derivatives' filesep 'LIMO_Face_detection' filesep 'Beta_files_FaceRepAll_GLM_Channels_Time_IRLS.txt']},...
    'analysis_type','Full scalp analysis','parameters',{[1 2 3],[4 5 6],[7 8 9]},...
    'factor names',{'face','repetition'},'type','Channels','nboot',1000,'tfce',1,'skip design check','yes');
limo_eeg(5,pwd)

%% compute the mean ERPs to visualize differences
for method = 1:3
    if method == 1
        Files = [STUDY.filepath filesep 'LIMO_' STUDY.filename(1:end-6) filesep ...
            'LIMO_files_FaceRepAll_GLM_Channels_Time_OLS.txt'];
        cd([studypath filesep 'derivatives' filesep 'ANOVA_OLS'])
    elseif method == 2
        Files = [STUDY.filepath filesep 'LIMO_' STUDY.filename(1:end-6) filesep ...
            'LIMO_files_FaceRepAll_GLM_Channels_Time_WLS.txt'];
        cd([studypath filesep 'derivatives' filesep 'ANOVA_WLS'])
    else
        Files = [STUDY.filepath filesep 'LIMO_' STUDY.filename(1:end-6) filesep ...
            'LIMO_files_FaceRepAll_GLM_Channels_Time_IRLS.txt'];
        cd([studypath filesep 'derivatives' filesep 'ANOVA_IRLS'])
    end
    parameters = [1 2 3];
    savename1  = [pwd filesep 'famous_faces.mat'];
    limo_central_tendency_and_ci(Files, parameters, chanlocs, 'Weighted mean', 'Mean', [],savename1)
    parameters = [4 5 6];
    savename2  = [pwd filesep 'scrambled_faces.mat'];
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
    fprintf('Faces vs. Scrambled @%gms = %g uV [%g %g]\n',LIMO.data.timevect(peaktime),TM(50,peaktime,2),TM(50,peaktime,3),TM(50,peaktime,1))
    fprintf('Faces vs. Scrambled @%gms = %g uV [%g %g]\n',LIMO.data.timevect(171),TM(50,171,2),TM(50,171,3),TM(50,171,1))
    
    limo_add_plots([savename1(1:end-4) '_Mean_of_Weighted mean.mat'],...
        [savename2(1:end-4) '_Mean_of_Weighted mean.mat'],[savename3(1:end-4) '_Mean_of_Weighted mean.mat'],...
        'channel',7); title('Face type at channel 007')
    % get a measure of raw effect size based on famous faces peak
    peakmethod = [85 82 80];
    peaktime   = peakmethod(method);
    fprintf('Faces vs. Scrambled @%gms = %g uV [%g %g]\n',LIMO.data.timevect(peaktime),TM(7,peaktime,2),TM(7,peaktime,3),TM(7,peaktime,1))
end


%% compare ANOVA results to understand differences between methods

% Observed data
n1 = 51; N1 = 8193;
n2 = 74; N2 = 8201;
sum(mask(:)>0) / numel(mask) * 100
% Pooled estimate of proportion
p0 = (n1+n2) / (N1+N2)
% Expected counts under H0 (null hypothesis) 5/100*numel(mask)
n10 = N1 * p0;
n20 = N2 * p0;
% Chi-square test, by hand
observed = [n1 N1-n1 n2 N2-n2];
expected = [n10 N1-n10 n20 N2-n20];
chi2stat = sum((observed-expected).^2 ./ expected)
p = 1 - chi2cdf(chi2stat,1)
       
       

F_values     = NaN([70*176,3,3]);
Cluster_mass = cell(3,3);
tfce_maps    = NaN([70*176,3,3]);
Tsquare      = NaN([70*176,3,3]);

% get data
load([studypath filesep 'derivatives\limo_gp_level_chanlocs.mat'])
for method = 1:3
    if method  == 1
        cd([studypath filesep 'derivatives\LIMO_Face_detection\2nd_level\Rep_2factors_ANOVA' filesep 'time_OLS'])
    elseif method  == 2
        cd([studypath filesep 'derivatives\LIMO_Face_detection\2nd_level\Rep_2factors_ANOVA' filesep 'time_WLS'])
    else
        cd([studypath filesep 'derivatives\LIMO_Face_detection\2nd_level\Rep_2factors_ANOVA' filesep 'time_IRLS'])
    end
    
    LIMO = load(fullfile(pwd,'LIMO.mat'));
    LIMO = LIMO.LIMO;
    for effect = length(LIMO.design.C):-1:1
        if effect == 1
            tmp = load(['tfce' filesep 'tfce_Rep_ANOVA_Main_effect_1_face.mat']);
            tfce_maps(:,method,effect) = tmp.tfce_score(:);
            tmp = load('Rep_ANOVA_Main_effect_1_face.mat');
            [~, mask] = limo_stat_values('Rep_ANOVA_Main_effect_1_face.mat',0.05,1,LIMO);
            nsig(1,method,effect) = sum(mask(:)>0);
            [~, mask] = limo_stat_values('Rep_ANOVA_Main_effect_1_face.mat',0.05,2,LIMO);
            nsig(2,method,effect) = sum(mask(:)>0);
            [~, mask] = limo_stat_values('Rep_ANOVA_Main_effect_1_face.mat',0.05,3,LIMO);
            nsig(3,method,effect) = sum(mask(:)>0);
        elseif effect == 2
            tmp = load(['tfce' filesep 'tfce_Rep_ANOVA_Main_effect_2_repetition.mat']);
            tfce_maps(:,method,effect) = tmp.tfce_score(:);
            tmp = load('Rep_ANOVA_Main_effect_2_repetition.mat');            
             [~, mask] = limo_stat_values('Rep_ANOVA_Main_effect_2_repetition.mat',0.05,1,LIMO);
            nsig(1,method,effect) = sum(mask(:)>0);
            [~, mask] = limo_stat_values('Rep_ANOVA_Main_effect_2_repetition.mat',0.05,2,LIMO);
            nsig(2,method,effect) = sum(mask(:)>0);
            [~, mask] = limo_stat_values('Rep_ANOVA_Main_effect_2_repetition.mat',0.05,3,LIMO);
            nsig(3,method,effect) = sum(mask(:)>0);
       else
            tmp = load(['tfce' filesep 'tfce_Rep_ANOVA_Interaction_Factors_12.mat']);
            tfce_maps(:,method,effect) = tmp.tfce_score(:);
            tmp = load('Rep_ANOVA_Interaction_Factors_12.mat');            
             [~, mask] = limo_stat_values('Rep_ANOVA_Interaction_Factors_12.mat',0.05,1,LIMO);
            nsig(1,method,effect) = sum(mask(:)>0);
            [~, mask] = limo_stat_values('Rep_ANOVA_Interaction_Factors_12.mat',0.05,2,LIMO);
            nsig(2,method,effect) = sum(mask(:)>0);
            [~, mask] = limo_stat_values('Rep_ANOVA_Interaction_Factors_12.mat',0.05,3,LIMO);
            nsig(3,method,effect) = sum(mask(:)>0);
        end
        F = squeeze(tmp.Rep_ANOVA(:,:,1));
        chan = diff(F,1).^2; time = diff(F,2).^2;
        Smoothness(method,effect) = sqrt(mean(chan(:))+mean(time(:)));
        F_values(:,method,effect) = F(:);
        df   = rank(LIMO.design.C{effect});
        dfe  = 18-df; % 18 subjects
        Tsquare(:,method,effect) = F_values(:,method,effect) ./ (dfe/(17*df)); % dfe/(n-1)*df
        % check clusters
        [posclusterslabelmat,nposclusters] = limo_findcluster(tmp.Rep_ANOVA(:,:,2)<0.05, channeighbstructmat, 2);
        cm = zeros(1,nposclusters);
        for C = 1:nposclusters
            cm(C) = sum(F(posclusterslabelmat==C)); % sum stat value in a cluster label
        end
        Cluster_mass{method,effect} = cm;
    end
end

% summary statistics
index = 9;
for effect = 3:-1:1
    for method = 3:-1:1
        medianT(index)       = median(Tsquare(:,method,effect));
        minT(index)          = min(Tsquare(:,method,effect));
        maxT(index)          = max(Tsquare(:,method,effect));
        medianF(index)       = median(F_values(:,method,effect));
        minF(index)          = min(F_values(:,method,effect));
        maxF(index)          = max(F_values(:,method,effect));
        medianCluster(index) = median(Cluster_mass{method,effect});
        minCluster(index)    = min(Cluster_mass{method,effect});
        maxCluster(index)    = max(Cluster_mass{method,effect});
        medianTFCE(index)    = median(tfce_maps(:,method,effect));
        minTFCE(index)       = min(tfce_maps(:,method,effect));
        maxTFCE(index)       = max(tfce_maps(:,method,effect));
        index                = index - 1;
    end
end
cd ..

t = array2table([medianT' minT' maxT' medianF' minF' maxF'...
     medianCluster' minCluster' maxCluster' medianTFCE' minTFCE' maxTFCE'],...
     'VariableNames',{'medianT' 'minT' 'maxT' 'medianF' 'minF' 'maxF'...
     'medianCluster' 'minCluster' 'maxCluster' 'medianTFCE' 'minTFCE' 'maxTFCE'},...
     'RowNames',{'Face OLS','Face WLS','Face IRLS','Rep. OLS','Rep. WLS','Rep. IRLS','Int. OLS','Int. WLS','Int. IRLS'})
 writetable(t,'distribution_summary.csv')
 
% look at empiral power
nsig = nsig/numel(mask)*100;
nsig = [round(squeeze(nsig(1,:,[1 2])))' round(squeeze(nsig(2,:,[1 2])))' round(squeeze(nsig(3,:,[1 2])))'];
t = array2table(nsig,'RowNames',{'Face','Repetition'},'VariableNames',...
    {'Uncor OLS','Uncor WLS','Uncor IRLS','Cluster OLS','Cluster WLS','Cluster IRLS',...
    'TFCE OLS','TFCE WLS','TFCE IRLS'})
writetable(t,'empirical_power.csv')

% compare effect sizes using T^2 
HotellingT = [Tsquare(:,:,1) Tsquare(:,:,2) Tsquare(:,:,3)];
[Tdiff,CI,p,alphav,h]= rst_multicompare( HotellingT,[2 1;2 3;5 4;5 6;8 7; 8 9],...
    'alphav',0.05,'estimator','median','newfig','yes');

figure
hp{1} = uipanel('position',[0   2/3 .5 1/3]);
hp{2} = uipanel('position',[0.5 2/3 .5 1/3]);
hp{3} = uipanel('position',[0   1/3 .5 1/3]);
hp{4} = uipanel('position',[0.5 1/3 .5 1/3]);
hp{5} = uipanel('position',[0   0   .5 1/3]);
hp{6} = uipanel('position',[0.5 0   .5 1/3]);
for index = 1:6
    if index == 1; A = HotellingT(:,2); B = HotellingT(:,1); 
    elseif index == 2;  A = HotellingT(:,2); B = HotellingT(:,3); 
    elseif index == 3;  A = HotellingT(:,5); B = HotellingT(:,4); 
    elseif index == 4;  A = HotellingT(:,5); B = HotellingT(:,6); 
    elseif index == 5;  A = HotellingT(:,8); B = HotellingT(:,7); 
    else, A = HotellingT(:,8); B = HotellingT(:,9);  
    end
        
    scatterhist(A(:),B(:),'kernel','on','Direction','out','parent',hp{index});
    hold on; plot(0:ceil(max(A(:))),0:ceil(max(A(:))),'--k','LineWidth',2);
    axis([-1 max(A(:)) -1 max(B(:))]); grid on;
    xlabel('WLS T^2'); 
    if index ==1 || index ==3 || index ==5
        ylabel('OLS T^2')
    else
        ylabel('IRLS T^2')
    end
    if index ==1 || index ==2
    title(sprintf('Face effect Median difference: %g [%g %g] p=%g',Tdiff(index),CI(1,index),CI(2,index),p(index)));
    elseif index ==3 || index ==4
    title(sprintf('Repetition effect Median difference: %g [%g %g] p=%g',Tdiff(index),CI(1,index),CI(2,index),p(index)));
    else
    title(sprintf('Interaction Median difference: %g [%g %g] p=%g',Tdiff(index),CI(1,index),CI(2,index),p(index)));
    end
end

% make cluster_mass a matrix
mc = max(max(cellfun(@(x) length(x), Cluster_mass))); % max size
cluster_mass = NaN(mc,3,3);
for method = 1:3
    for effect = 1:3
        cluster_mass(1:length(Cluster_mass{method,effect}),method,effect) = Cluster_mass{method,effect};
    end
end

writematrix([F_values(:,:,1) F_values(:,:,2) F_values(:,:,3)],'F_values_distrution.csv')
writematrix([tfce_maps(:,:,1) tfce_maps(:,:,2) tfce_maps(:,:,3)],'tfce_values_distrution.csv')
writematrix([cluster_mass(:,:,1) cluster_mass(:,:,2) cluster_mass(:,:,3)],'cluster_mass_values_distrution.csv')

