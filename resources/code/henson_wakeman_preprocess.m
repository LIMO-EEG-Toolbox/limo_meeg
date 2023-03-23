% Basic Preprocessing Pipeline - adapted from Pernet & Delorme (2021)
% Ref: From BIDS-Formatted EEG Data to Sensor-Space Group Results: 
% A Fully Reproducible Workflow With EEGLAB and LIMO EEG.
% Front. Neurosci. 14:610388. doi: 10.3389/fnins.2020.610388
% <https://www.frontiersin.org/articles/10.3389/fnins.2020.610388/full>
%
% This function preprocesses Wakeman and Henson data to create ERPs, Spectrum, 
% and ERSP 

% start EEGLAB
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;

if ~exist('pop_importbids','file')
    plugin_askinstall('bids-matlab-tools',[],1);
end

if ~exist('pop_zapline_plus','file')
    plugin_askinstall('zapline-plus',[],1);
end

if ~exist('picard','file')
    plugin_askinstall('picard', 'picard', 1);
end

% call BIDS tool BIDS
[STUDY, ALLEEG] = pop_importbids(bids_folder,'bidsevent','on','bidschanloc','on',...
                  'studyName','Face_detection','outputdir', fullfile(bids_folder, 'derivatives'), ...
                  'eventtype', 'trial_type');
ALLEEG          = pop_select( ALLEEG, 'nochannel',{'EEG061','EEG062','EEG063','EEG064'});
CURRENTSTUDY    = 1;
EEG             = ALLEEG;
CURRENTSET      = 1:length(EEG);
root            = fullfile(bids_folder, 'derivatives');

% reorient if using previous version of the data
% for s=1:size(EEG,2)
%   EEG(s) = pop_chanedit(EEG(s),'nosedir','+Y');
% end

% Remove bad channels
EEG = pop_clean_rawdata( EEG,'FlatlineCriterion',5,'ChannelCriterion',0.8,...
    'LineNoiseCriterion',4,'Highpass',[0.25 0.75] ,...
    'BurstCriterion','off','WindowCriterion','off','BurstRejection','off',...
    'Distance','Euclidian','WindowCriterionTolerances','off' );

% remove line noise
for s=1:size(EEG,2)
    EEG(s) = pop_zapline_plus(EEG(s),'noisefreqs','line',...
        'coarseFreqDetectPowerDiff',4,'chunkLength',0,...
        'adaptiveNremove',1,'fixedNremove',1,'plotResults',0);
end
            
% Rereference using average reference
EEG = pop_reref( EEG,[]); % ,'interpchan',[]);

% Run ICA and flag artifactual components using IClabel
for s=1:size(EEG,2)
    try
        EEG(s) = pop_runica(EEG(s), 'icatype','picard','concatcond','on','options',{'pca',EEG(s).nbchan-1}); %#ok<*AGROW>
        EEG(s) = pop_iclabel(EEG(s),'default');
        EEG(s) = pop_icflag(EEG(s),[NaN NaN;0.8 1;0.8 1;NaN NaN;NaN NaN;NaN NaN;NaN NaN]);
        EEG(s) = pop_subcomp(EEG(s), find(EEG(s).reject.gcompreject), 0);
    catch icaerr
        suberror{s} = sprintf('ICA failed for subject %g, error %s',s,icaerr.message); %#ok<NASGU>
    end
end

% clean-up EEG and directories
if exist('suberror','var')
    error('something went wrong during ICA and IC labelling')
end

% clean data using ASR - just the bad epochs
EEG  = pop_clean_rawdata(EEG,'FlatlineCriterion','off','ChannelCriterion','off',...
     'LineNoiseCriterion','off','Highpass','off','BurstCriterion',20,...
     'WindowCriterion',0.3,'BurstRejection','on','Distance','Euclidian',...
     'WindowCriterionTolerances',[-Inf 7] );

% Extract data epochs (no baseline removed)
for s=size(EEG,2):-1:1
    try
        EEG(s)= pop_epoch(EEG(s),{'famous_new','famous_second_early','famous_second_late', ...
         'scrambled_new','scrambled_second_early','scrambled_second_late','unfamiliar_new', ...
         'unfamiliar_second_early','unfamiliar_second_late'},[-0.1 1] ,'epochinfo','yes');
    catch
        warning('subject %s failed to epoch',EEG(s).subject)
    end
end
EEG    = eeg_checkset(EEG);
EEG    = pop_saveset(EEG, 'savemode', 'resave');
ALLEEG = EEG;

% compute single trials
[STUDY, EEG]  = std_precomp(STUDY, EEG, 'channels', 'savetrials','on','interp','on','recompute','on',...
    'erp','on', 'spec','on', 'ersp','on','itc','on', ...
     'erspparams', {'freqs', [5 30], 'timelimits', [-50 650]});

[STUDY, EEG]  = std_precomp(STUDY, EEG, 'components', 'savetrials','on','interp','on','recompute','on',...
    'erp','on','erpparams', {'rmbase' [-200 0]}, 'spec','on', 'ersp','on','itc','on',...
    'erspparams', {'freqs', [5 30], 'timelimits', [-50 650]});

[STUDY, EEG] = pop_savestudy( STUDY, EEG, 'savemode','resave');
