function [bEst,conds, scans,targ_idx,attn_coord,dist_coord,ang_offset] = pri2stim_extractSignal_avg(mr,whichTRs)
% extract the signal in each voxel on each trial by averaging whichTRs relative to trial onset
% also return the conditions (attn task, side, contrast level, target?
% extracts signal from mr, at times of interest after trial start
% (returns one set of vox responses per trial)
% no need to compute design matrix here
%
%[tst, conds, tst_scans,tpts,targ_idx,tr_coord,tn_coord,ang_offset]


% bEst, etc is stacked from tpt 0:nTRs-1, use tpt to index into it,
% conds/scans are just repmat(conds/scans,nTRs,1)

root = mr.root;


data_orig=double(mr.data);

voxGood = ~isnan(nanmean(data_orig,1)); % indexes nonNaN voxels
data = mr.data;

runs = mr.r;


startidx = 1;
    cnt = 1;
    for ii=runs
        
        
        
        
        pname = [root, mr.subName, '/', mr.rp, num2str(ii), '/', mr.subName, '_', mr.rp, num2str(ii), '.prt'];

        prt = BVQXfile(pname);
        curD = (data(cnt*mr.NofTRs-mr.NofTRs+1:cnt*mr.NofTRs, :));
        
        mname = pname;
        mname(end-3:end) = '.mat';
        load(mname);
        
        n_trials = size(stim.attnCoordDeg,1);
        
        if ii == 1

            nblank = size(stim.attnCoordDeg,1) * length(runs);
            bEst = nan(nblank,size(data,2));
            conds = nan(nblank,size(stim.conditions,2)); % check dim2...
            scans = nan(nblank,1);
            targ_idx = nan(nblank,2);
            attn_coord = nan(nblank,2);
            dist_coord = nan(nblank,2);

            clear nblank;
            
        end
        
        thisidx = startidx:(startidx+n_trials-1);
        
        
        for t = 1:length(prt.Cond)
            %braw(t,:) = curD(ceil(prt.Cond(t).OnOffsets(1)/mr.TR),:);
            trial_start_idx = ceil(prt.Cond(t).OnOffsets(1)/mr.TR);
            %disp(sprintf('stim dur = %i, t_idx = :',stim(t).stimExposeDur));
            t_idx = whichTRs + trial_start_idx;%:(trial_start_idx+avg_over_TRs-1)];
            %braw(t,:) = curD(round(prt.Cond(t).OnOffsets(1)/mr.TR)+shift_TRs,:);
            braw(t,:) = mean(curD(t_idx,:),1);
            clear trial_start_idx;
        end

        
        

        condnew = stim.conditions;
        attnnew = stim.attnLocIdx;
        distnew = stim.distLocIdx;
        

        conds(thisidx,:) = condnew;
        bEst(thisidx,:) = braw;
        
        targ_idx(thisidx,:) = [stim.attnLocIdx stim.distLocIdx];
        attn_coord(thisidx,:) = stim.attnCoordDeg;
        dist_coord(thisidx,:) = stim.distCoordDeg;
        ang_offset(thisidx) = stim.allAngOffset;
        
        tr(thisidx) = attnnew;
        tn(thisidx) = distnew;
        
        scans(thisidx) = ones(size(braw,1), 1)*ii ;
        %tpt(thisidx) = ones(size(braw,1),1)*tpts(tt);
        cnt = cnt + 1;
        startidx = thisidx(end)+1;
        clear stim;
       
    end
    


bEst = bEst(:,voxGood);
end