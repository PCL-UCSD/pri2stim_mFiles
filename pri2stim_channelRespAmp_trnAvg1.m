function pri2stim_channelRespAmp_trnAvg1(subj,VOIs,hex_size)
% TCS 12/16/2015 - adapted form wmLazSpace_channelRespAmp_trnAvg1.m
% - got rid of saving "tpts", just assumes one data point per trial right
% now

if nargin < 1
    
    subj = {'AI151','AP151','AR151','AS151','BA151','BB151','BC151','BF151'};
    
end


if nargin < 2
    VOIs = {'V1','V2','V3','V3A','V4','IPS0','IPS1','IPS2','IPS3','sPCS'};
end

if nargin < 3
    hex_size = 7;
end


locstr = 'wmLazLocAll';

trialDataStr = 'pri2stim1';

trnPrefix = 'IEM_hexMap';
tstPrefix = 'pri2stim';

% optimize scaling given stim parameters - horiz separation (for now)
% should be equal to radDeg

% FWHM:spacing should be ~1.25:1

filt_scale = 1.1; % FWHM: spacing ratio need to multiply by 1/rad2fwhm(1) to convert to cos7 size constant


root = '/usr/local/serenceslab/tommy/pri2stim/';

for ss = 1:size(subj,2)   % loop over subjects

    for vv = 1:length(VOIs)
        
        
        
        %% load trialData files
        trialData_fn = sprintf('%spri2stim_trialData/%s_%s_%s.mat',root,subj{ss},VOIs{vv},trialDataStr);
        fprintf('Loading trial data from: %s...\n',trialData_fn);
        load(trialData_fn);
        
        
        % tr_coord, tn_coord are in CARTESIAN coordinates (+,+ for quadrant
        % 1), ang_offset is in radians, was added to theta of stimLocsX/Y in
        % behavioral script BEFORE flipping Y, so when working in screen
        % coords, -ang_offset was added (so use +ang_offset to fix)
        
        %% generate design matrix for training data
        %
        % TODO: can do this just once, before loading each ROI's data for
        % each subj (design matrix should be the same for all ROIs)
        
        fprintf('generating design matrix\n');
        
        % load this subject's stimulus parameter file
        fn = sprintf('%sIEM_hexMap_params/%s_%s_params.mat',root,subj{ss},trnPrefix);
        load(fn);
        
        % stimulus was part of a make_hexagon(7) grid, made up of
        % equilateral triangles with wmDropHex_params.stepSize sides, so
        % the maximum distance from fixation the stimulus reached was:
        %   max_dist = 3 * wmMapHex_params.stepSize +
        %   wmMapHex_params.radDeg + wmMapHex_params.jitterRadDeg
        
        [rfX, rfY] = make_hex(hex_size); 
        max_dist = 3 * map_params.stepSize + map_params.radDeg + map_params.jitterRadDeg;
        scale_factor = 0.7*max_dist*(2/sqrt(3));
        %scale_factor = max_dist*(2/sqrt(3));
        rfX = rfX*scale_factor;  rfY = rfY*scale_factor;  
        
        basis_sep = rfX(2)-rfX(1); % how far, in DVA, the basis functions are spaced apart
        rfSize = basis_sep * filt_scale;  % size of basis functions, in FWHM units
        
        % TODO: rotate the points?wmMapSpace
        
        
        trnX = hexMap_makeX1(stimLocs,map_params.radDeg,rfX,rfY,rfSize,251);
        
        % we're only going to use target absent trials (for now)
        trnX = trnX(targPresent~=1,:);
        a_trn = a_trn(targPresent~=1,:);
        
        
        
        trnX = trnX/max(trnX(:));
         % should we add a constant term?
        
        clear map_params max_dist scale_factor;
        
        
        
        %% (maybe) optimize design matrix for training data
        
        
        
        
        %% compute channel weights using only training set (wmMapSpace)
        
        fprintf('computing channel weights\n');
        w = trnX\a_trn;   % a_trn: train w/ average
        
        
        %% use (optimized) design matrix to compute channel responses using testing data
        
        fprintf('computing channel responses\n');
        
        chan_resp = inv(w*w')*w*a_tst';
        
        
        %% save channel responses, etc
        
        fn2s = sprintf('%spri2stim_chanResp/%s_%s_pri2stim_hex%i_channelResp_trnAvg1.mat',root,subj{ss},VOIs{vv},hex_size);
        fprintf('saving to %s...\n\n',fn2s);
        save(fn2s,'chan_resp','w','rfX','rfY','rfSize','conds','tst_scans','attn_coord','dist_coord','targ_idx','ang_offset','locstr','whichTRs_trn');
        %clear fn2s w rfX rfY rfSize conds tst_scans tpts tr_coord tn_coord targ_idx ang_offset;
        clear a_trn s_trn stimLocs targPresent b_trn b_tst a_tst conds tst_scans  targ_idx attn_coord dist_coord ang_offset mrfn fn2s w rfX rfY rfSize;
        
    end % end VOI loop


end % end subj loop

end % end fcn


