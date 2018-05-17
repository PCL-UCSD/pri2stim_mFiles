function [bEst, scans, stimLocs,targPresent] = hexMap_extractSignal_avg(mr,whichTRs)
% extract the signal in each voxel on each trial, using trial-by-trial GLM
% also return the conditions (attn task, side, contrast level, target?
% extracts signal from mr, at times of interest after trial start
% (returns one set of vox responses per trial)
% no need to compute design matrix here
%
% adapted from wmDrop_extractSignal_raw.m by TCS 2/6/2015
% - now takes in a set of TRs ([3 4]) and returns the average of those TRs
%   on each trial within each voxel


% bEst, etc is stacked from tpt 0:nTRs-1, use tpt to index into it,
% conds/scans are just repmat(conds/scans,nTRs,1)

root = mr.root;

tpts = whichTRs;
%avg_over_TRs = 1;  % can increase this if we want temporal smoothing




data_orig=double(mr.data);

voxGood = ~isnan(nanmean(data_orig,1)); % indexes nonNaN voxels
data = mr.data;

runs = mr.r;

% build a design matrix for stims on the left and stims on the right,
% only based on all runs except the present run



startidx = 1;
cnt = 1;
for ii=runs
    
    
    
    
    pname = [root, mr.subName, '/', mr.rp, num2str(ii), '/', mr.subName, '_', mr.rp, num2str(ii), '.prt'];
    %fprintf('loading PRT %s...\n',pname);
    prt = BVQXfile(pname);
    curD = (data(cnt*mr.NofTRs-mr.NofTRs+1:cnt*mr.NofTRs, :));
    
    mname = pname;
    mname(end-3:end) = '.mat';
    load(mname);
    
    n_trials = size(stim.stimLocsDeg,1);
    
    if ii == 1
        
        nblank = n_trials * length(runs);
        bEst = nan(nblank,size(data,2));
        stimLocs = nan(nblank,2);
        targPresent = nan(nblank,1);

        scans = nan(nblank,1);

        clear nblank;
        
    end
    
    thisidx = startidx:(startidx+n_trials-1);
    
    try
        for t = 1:length(prt.Cond)
            %braw(t,:) = curD(ceil(prt.Cond(t).OnOffsets(1)/mr.TR),:);
            trial_start_idx = ceil(prt.Cond(t).OnOffsets(1)/mr.TR);
            %disp(sprintf('stim dur = %i, t_idx = :',stim(t).stimExposeDur));
            t_idx = whichTRs + trial_start_idx;%:(trial_start_idx+avg_over_TRs-1)];
            %braw(t,:) = curD(round(prt.Cond(t).OnOffsets(1)/mr.TR)+shift_TRs,:);
            braw(t,:) = mean(curD(t_idx,:),1);
            clear trial_start_idx;
        end
        
    catch
        disp('whoops');
    end
    
    
    
    stimLocs(thisidx,:) = stim.stimLocsDeg;

    
    bEst(thisidx,:) = braw;
    

    scans(thisidx) = ones(size(braw,1), 1)*ii ;
    targPresent(thisidx) = stim.targPresent;

    cnt = cnt + 1;
    startidx = thisidx(end)+1;
    clear stim;
    
end

bEst = bEst(:,voxGood);
end