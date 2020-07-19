% start EEGLAB
clear variables
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;

% import BIDS
filepath        = 'F:\WakemanHenson_Faces\eeg';
[STUDY, ALLEEG] = pop_importbids(filepath, 'bidsevent','on','bidschanloc','on', 'studyName','Face_detection');
ALLEEG = pop_select( ALLEEG, 'nochannel',{'EEG061','EEG062','EEG063','EEG064'});
CURRENTSTUDY = 1; EEG = ALLEEG; CURRENTSET = 1:length(EEG);

% reorient if using previous version od the data
% EEG = pop_chanedit(EEG,'nosedir','+Y');

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

% Create study design
STUDY  = std_checkset(STUDY, ALLEEG);
STUDY  = std_makedesign(STUDY, EEG, 1, 'name','STUDY.FaceRepetition','delfiles','off','defaultdesign','off','variable1','type','values1',{});
eeglab redraw

% Precompute ERP and Spectrum measures
[STUDY, EEG] = std_precomp(STUDY, EEG, {}, 'savetrials','on','interp','on','recompute','on',...
    'erp','on','erpparams', {'rmbase' [-200 0]}, 'spec','on',...
    'ersp','on','itc','on', 'specparams',{'specmode','fft','logtrials','off'});
eeglab redraw

STUDY = pop_limo(STUDY, ALLEEG, 'method','WLS','measure','daterp','timelim',[-50 650],'erase','on','splitreg','off','interaction','off');
STUDY = pop_limo(STUDY, ALLEEG, 'method','WLS','measure','datspec','freqlim',[3 45],'erase','on','splitreg','off','interaction','off');
STUDY = pop_limo(STUDY, ALLEEG, 'method','WLS','measure','dattimef','timelim',[-50 650],'freqlim',[3 45],'erase','on','splitreg','off','interaction','off');
eeglab redraw

