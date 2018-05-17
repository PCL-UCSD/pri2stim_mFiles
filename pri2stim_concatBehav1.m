% pri2stim_concatBehav1.m

function pri2stim_concatBehav1(subj)
% concatenates behavioral files from each run into one file per subj, to
% line up w/ trialData files
%
% these files will be loaded and plotted later on
%
% updated 8/7/2016, TCS, added timing like in doPrt_pri2stim.m



root = load_root;

if nargin < 1
    subj = {'AI151','AP151','AR151','AS151','BA151','BB151','BC151','BF151'};
end


for ss = 1:length(subj)
   
    myfs = fullfile(root,subj{ss},sprintf('%s_Behav',subj{ss}),sprintf('%s_pri2stim_scanner1*.mat',subj{ss}));
    fprintf('Looking for: %s\n',myfs);
    
    
    myf = dir(myfs);
    fprintf('%i files found...\n',length(myf));
    
    
    all_fn = cell(length(myf),1);
    
    startidx = 1;
    
    for ff = 1:length(myf)
        
        fn = fullfile(root,subj{ss},sprintf('%s_Behav',subj{ss}),myf(ff).name);
        fprintf('Loading %s\n',fn);
        
        thisbehav = load(fn);
        all_fn{ff} = fn;
        
        if ff == 1
            % allocate variables
            nblank = length(myf)*size(thisbehav.p.conditions,1);
            
            all_r = nan(nblank,2); % run, subrun
            all_conds = nan(nblank,size(thisbehav.p.conditions,2));
            all_correct = nan(nblank,1);
            all_resp = nan(nblank,1);
            all_rt = nan(nblank,1);
            
            % WITHIN-TRIAL!!!!
            all_targOnset = nan(nblank,2); % targ, distractor onset time (for looking at hazard-type stuff)
            
            all_stimLocs{1} = nan(nblank,2);
            all_stimLocs{2} = nan(nblank,2);
            
            % WITHIN RUN
            all_cueOnset = nan(nblank,1);
            all_stimOnset = nan(nblank,1);
            all_stimOffset = nan(nblank,1);
            
        end
        
        thisidx = startidx:(startidx+size(thisbehav.p.conditions,1)-1);
        
        all_r(thisidx,1) = thisbehav.p.runNum;
        all_r(thisidx,2) = thisbehav.p.runSub;
        
        all_conds(thisidx,:) = thisbehav.p.conditions;
        all_correct(thisidx) = thisbehav.p.correct;
        all_resp(thisidx) = thisbehav.p.resp;
        all_rt(thisidx) = thisbehav.p.responseTime;
        
        all_targOnset(thisidx,:) = thisbehav.p.targOnsetFrame/60;
        
        all_stimLocs{1}(thisidx,:) = thisbehav.p.stimLocsDeg{1};
        all_stimLocs{2}(thisidx,:) = thisbehav.p.stimLocsDeg{2};
        
        % matched to how time is computed in doPrt scripts
        all_stimOnset(thisidx)  = thisbehav.p.stimStart-thisbehav.p.startExp;
        all_stimOffset(thisidx) = thisbehav.p.stimEnd  -thisbehav.p.startExp;
        all_cueOnset(thisidx)   = thisbehav.p.cueStart -thisbehav.p.startExp;
        
        startidx = thisidx(end)+1;
        
        clear fn thisbehav;
    end
    
    fn2s = fullfile(root,'pri2stim_behav',sprintf('%s_pri2stim_behav.mat',subj{ss}));
    fprintf('saving to %s\n',fn2s);
    
    save(fn2s,'all_r','all_conds','all_correct','all_resp','all_rt','all_targOnset','all_stimLocs','all_fn','all_stimOnset','all_stimOffset','all_cueOnset');
    
    
    % TODO: same for IEM_hexMap!!!!
end


return