
function pri2stim_concatMRs(subj,VOIs,trnPrefix,tstPrefix)
%
% 12/16/2015 (TCS)
% - only does trial-level signal extraction (right now, just avgs for test
% data); doesn't pull out timecourse, though we could add that back later
% on



if nargin < 1
    
    subj = {'AR151','AP151','AR151','AS151','BA151','BB151','BC151','BF151'}
    
end


if nargin < 2
   VOIs = {'V1','V2','V3','V3A','V4','IPS0','IPS1','IPS2','IPS3','sPCS'};

end

if nargin < 3
    trnPrefix = 'hexMap';
end

if nargin < 4

    tstPrefix = 'pri2stim';
end


locstr = 'wmLazLocAll';

mrstr = 'VT1_IntInf_Norm2_RmMo1_Vol_1_';

truncOutliers = 0;
normalize = 0;
vKeep = inf;
keep = inf;
rawData = 3;
thresh = '05';

whichTRs_trn = [3 4];
whichTRs_tst = [3 4];



hemis = {'LH','RH'};

root = '/usr/local/serenceslab/tommy/pri2stim/';

for ss = 1:size(subj,2)   % loop over subjects
    sn=char(subj{1,ss});
    
     
    %loop over VOIs
    vcnt = 1;

    
    for vv=1:length(VOIs)
        % build a struct array of mr's - one for each cond (for each VOI)
        mrfn = cell(2,1);
        for c = 1:2   % TRN & TEST
            
            
            for hh = 1:2
                
                % AP81_wmDropLoc_VT1_IntInf_Norm2_RmMo1_Vol_1_wmDrop_LH-V1.mat	
                if c == 2  % TESTING DATA
                    myfs = [root 'pri2stim_mrStructs/' subj{ss} '_' locstr '_' mrstr tstPrefix '_' hemis{hh} '-'  VOIs{vv} '.mat'];
                    myf = dir(myfs);
                    load([root 'pri2stim_mrStructs/' myf(1).name]);
                    clear myf myfs;
                elseif c == 1 % TRAINING DATA
                    
                    %load/supplement extracted VOI data for each task condition
                    
                    myfs = [root 'pri2stim_mrStructs/' subj{ss} '_' locstr '_' mrstr trnPrefix '_' hemis{hh} '-' VOIs{vv} '.mat'];
                    myf = dir(myfs);
                    try
                        fprintf('loading: %s\n',myf(1).name);
                    catch
                        disp('uh oh');
                    end
                    load([root 'pri2stim_mrStructs/' myf(1).name]);
                    clear myf myfs;
                end
                
                if c == 2
                    fprintf('Processing data from %s, %s numVox: %d\n', char(mr.voiNames{1}),[tstPrefix],size(mr.voiData(1).mod,2))
                elseif c == 1
                    fprintf('Processing data from %s, %s numVox: %d\n', char(mr.voiNames{1}),[trnPrefix],size(mr.voiData(1).mod,2))
                end
                % select the data type
                if rawData==1
                    tmpdata{hh} = mr.voiData(1).rd;
                    
                elseif rawData==2
                    tmpdata{hh} = mr.voiData(1).nd;
                    
                elseif rawData == 3
                    tmpdata{hh} = mr.voiData(1).mod;
                    
                end
                if hh == 1
                    clear mr;
                end
            end
             % now for RH...
            
            mr.data = [tmpdata{1} tmpdata{2}];
           
            
            mr.truncOutliers = truncOutliers;
            mr.rawData = rawData;
            mr.normalize = normalize;
            mr.thresh = thresh;
            %mr.numDirs = numDirs;
            
            
            mr.root = root;
            
            mr.locstr = locstr;
            
            mr.thisVOI = ['Bilat-' VOIs{vv}];

            if c == 2 % TEST
                mr.cond = tstPrefix;
                mr_tst = mr;
                
            elseif c == 1 % TRAIN
                mr.cond = trnPrefix;
                mr_trn = mr;
            end
            
            % save bilateral MR
            mrfn{c} = sprintf('%spri2stim_mrStructs/%s_%s_%s%s_Bilat-%s.mat',root,subj{ss},locstr,mrstr,mr.cond,VOIs{vv});
            fprintf('saving bilateral mr struct to %s\n\n',mrfn{c});
            save(mrfn{c},'mr');
            
            
             
            
        end % end cond loop
        
        
        
        %% extract timecourses from testing data
        
        
            
            % extract training data: avg of timepoints [3 4]
            fprintf('extracting training data (avg)\n');
            [a_trn, s_trn, stimLocs, targPresent] = hexMap_extractSignal_avg(mr_trn,whichTRs_trn);

            
            % also get betas
            fprintf('extracting training data (betas)\n');
            [b_trn] = hexMap_extractSignal_betas(mr_trn);
            
            

            fprintf('extracting testing data (beta timecourse)\n');
            % maybe make a separate function for loading conditions? sloppy
            % having these together
            [b_tst, conds, tst_scans,targ_idx,attn_coord,dist_coord,ang_offset] = pri2stim_extractSignal_betas(mr_tst);
            
            fprintf('extracting testing data (avg tpts)\n');
            [a_tst] = pri2stim_extractSignal_avg(mr_tst,whichTRs_tst);
            
            fn2s = sprintf('%spri2stim_trialData/%s_%s_pri2stim1.mat',root,subj{ss},VOIs{vv});
            fprintf('Saving extracted data to: %s...\n',fn2s);
            save(fn2s,'mrfn','locstr','a_trn','s_trn','stimLocs','targPresent','b_trn','b_tst','a_tst','conds','tst_scans', 'targ_idx','attn_coord','dist_coord','ang_offset','whichTRs_trn');
            clear fn2s a_trn a_tst s_trn stimLocs targPresent b_trn b_tst conds tst_scans targ_idx tr_coord_tn_coord ang_offset mrfn;

            
            
               
        
        
    end
end
return;