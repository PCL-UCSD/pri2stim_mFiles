function pri2stim_computeReconstructions_posCoreg1(subj, VOIs,hex_size)
%
% TCS 12/16/2015 - removed the support for timecourses - just operates on
% single-trial data points (betas or average) (for now, can add back
% support for timecourse - see wmDrop_computeReconstructions_posCoreg1.m)
%


if nargin < 1
    subj = {'AI151','AP151','AR151','AS151','BA151','BB151','BC151','BF151'};
end

if nargin < 2 
    VOIs = {'V1','V2','V3','V3A','V4','IPS0','IPS1','IPS2','IPS3','sPCS'};
end

if nargin < 3
    hex_size = 7;
end

root = '/usr/local/serenceslab/tommy/pri2stim/';


ecc_to_align = 3.5; % align to this point

myTR = 2.00;  % repetition time, in sec

mycolors = [0 127 66; 56 9 124; 28 69 68]./255;


maxecc = 5.5; % dva from fixation
res = 101; % in x, y

[gridx,gridy] = meshgrid(linspace(-maxecc,maxecc,res),linspace(-maxecc,maxecc,res));
gridx = reshape(gridx,numel(gridx),1);gridy = reshape(gridy,numel(gridy),1);

for ss = 1:length(subj)
    
    for vv = 1:length(VOIs)
        
        % load data from /matData/
        % subj_ROI_channelRespAmp.mat
        chan_fn = sprintf('%s/pri2stim_chanResp/%s_%s_pri2stim_hex%i_channelResp_trnAvg1.mat',root,subj{ss},VOIs{vv},hex_size);
        fprintf('loading %s...\n',chan_fn);
        chan = load(chan_fn);
        
        [rfT, rfR] = cart2pol(chan.rfX,chan.rfY);
        
        % cartesian: x/y
        t_coord_cart{1} = chan.attn_coord;
        t_coord_cart{2} = chan.dist_coord;
        
        % theta/r  
        [t,r] = cart2pol(chan.attn_coord(:,1),chan.attn_coord(:,2));
        t_coord_pol{1} = [t r]; clear t r;
        [t,r] = cart2pol(chan.dist_coord(:,1),chan.dist_coord(:,2));
        t_coord_pol{2} = [t r]; clear t r;
        
        
        coreg_coord_pol{1} = nan(size(t_coord_pol{1}));
        coreg_coord_pol{2} = nan(size(t_coord_pol{2}));
        
        % each trial needs to be rotated by ang_offset so that stim discs
        % align, then rotated by (conds(:,1)-1)*60 deg so that relevant
        % target (TR/TN) is at left horizontal meridian
        % ADD ang_offset to t_coord_pol(:,1)
        %
        % trials for which displacement is -1 (4th column of conditions)
        % should be flipped up/down (rfY = rfY*-1) - REVERSE THIS SIGN FOR
        % ii==22 (TN)
        
        n_pos = length(unique(chan.conds(:,3))); % 5 here
        stim_pos = deg2rad((0:(n_pos-1)) * (360/n_pos)); % points along circle - 0:60:(360-60) (deg), converted to radians
        %chan.conds(chan.conds(:,3)==3,4) = 0; % oops - this didn't happen correctly when creating conditions matrix in scanner script
        

        % create a trial_num index so that we don't need to keep generating
        % basis sets

        n_trials = size(chan.conds,1);%/length(tu);
        trial_num = (1:n_trials);
        
        recons_vec = cell(2,1);
        
        for ii = 1:2   % TR, TN
            recons_vec{ii} = nan(size(chan.chan_resp,2),numel(gridx));
            for tt = 1:n_trials
                
                % all of this is in cartesian coords, so this is correct
                this_rfT = rfT-chan.ang_offset(tt);
                this_rfT = this_rfT-stim_pos(chan.targ_idx(tt,ii));
                
                this_rfR = rfR;
                
                coreg_coord_pol{ii}(tt,:) = [t_coord_pol{ii}(tt,1)-chan.ang_offset(tt)-stim_pos(chan.targ_idx(tt,ii)) t_coord_pol{ii}(tt,2)];
                
                % for plotting basis set
                
                [this_rfX, this_rfY] = pol2cart(this_rfT,this_rfR);
                
                
                if ii == 1 && chan.conds(tt,5)==-1
                    this_rfY = -1*this_rfY;
                elseif ii == 2 && chan.conds(tt,5)==1
                    this_rfY = -1*this_rfY;
                end
                
                
                % build basis
                basis_set = build_basis_pts(this_rfX,this_rfY,chan.rfSize/rad2fwhm(1),gridx,gridy);
                
                % get all timepoints from that trial
                thisidx = trial_num==tt;
                thisdata = chan.chan_resp(:,thisidx);
                
                % reconstruct (as in wmdelay_reconstruct*.m)
                recons_vec{ii}(thisidx,:) = thisdata'*basis_set';
                
            end
            [x,y] = pol2cart(coreg_coord_pol{ii}(:,1),coreg_coord_pol{ii}(:,2));
            coreg_coord_cart{ii} = [x y]; clear x y;
            
        end
        
        
        conds = chan.conds;
        
        % save file
        fn2s = sprintf('%spri2stim_recons/%s_%s_hex%i_pos_coreg1.mat',root,subj{ss},VOIs{vv},hex_size);
        fprintf('saving to %s...\n',fn2s);
        save(fn2s,'recons_vec','conds','t_coord_cart','chan_fn','res','maxecc');
        
        clear chan;
        
    end
    
end


return
