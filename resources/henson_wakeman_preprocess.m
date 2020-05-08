% start EEGLAB
clear
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;

% import BIDS
filepath        = 'F:\WakemanHenson_Faces\eeg';
[STUDY, ALLEEG] = pop_importbids(filepath, 'bidsevent','on','bidschanloc','on', 'studyName','Face_detection');
ALLEEG = pop_select( ALLEEG, 'nochannel',{'EEG061','EEG062','EEG063','EEG064'});
CURRENTSTUDY = 1; EEG = ALLEEG; CURRENTSET = [1:length(EEG)];

% Clean data using ASR - here EEG contains all the datasets
EEG = pop_clean_rawdata( EEG,'FlatlineCriterion',5,'ChannelCriterion',0.8,'LineNoiseCriterion',4,'Highpass',[0.25 0.75] ,'BurstCriterion',20,'WindowCriterion',0.25,'BurstRejection','on','Distance','Euclidian','WindowCriterionTolerances',[-Inf 7] );

% Rereference using average reference
EEG = pop_reref( EEG,[],'interpchan',[]);

% Run ICA and flag artifactual components using IClabel
EEG = pop_runica(EEG, 'icatype','runica','concatcond','on','options',{'pca',16});
EEG = pop_iclabel( EEG,'default');
EEG = pop_icflag( EEG,[NaN NaN;0.9 1;0.9 1;NaN NaN;NaN NaN;NaN NaN;NaN NaN]);

% Extract data epochs (no baseline removed)
EEG = pop_epoch( EEG,{'famous_new','famous_second_early','famous_second_late','scrambled_new','scrambled_second_early','scrambled_second_late','unfamiliar_new','unfamiliar_second_early','unfamiliar_second_late'},[-0.5 1] ,'epochinfo','yes');
EEG = eeg_checkset(EEG);
EEG = pop_saveset(EEG, 'savemode', 'resave');
ALLEEG = EEG;
eeglab redraw

