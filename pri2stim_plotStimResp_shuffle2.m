function pri2stim_plotStimResp_shuffle2(subj,VOIs,hex_size)
%
% adapted from pri2stim_plotStimResp_resample1.m, TCS 12/17/2015
%
% plots extracted signal from reconstruction at stimulus position, sorted
% by target & distractor contrast, and potentially by separation distance
% as well
%
% For attended stimulus, x-axis contrast is contrast of attended stim,
% averaged over trials w/ all distractor contrasts; similar for unattended
% stimulus (x-axis contrast is trials w/ that contrast for unattn stim,
% averaged over all attended stim contrasts)
%
% TCS 9/5/2016 - updated to use common RNG seed

rand_seed = load_rng_seed;
rng(rand_seed);

save_stats = 0;

if nargin < 1    
    subj = {'AI151','AP151','AR151','AS151','BA151','BB151','BC151','BF151'};
end

if nargin < 2 

    VOIs = {'V1','V2','V3','V3A','V4','IPS0','IPS1','IPS2','IPS3','sPCS'};

end

if nargin < 3
    hex_size = 7;
end


niter = 1000;

root = load_root;
load(sprintf('%spri2stim_colors.mat',root));


targ_pos = [3.5 0]; % reconstructions are aligned to this point
targ_window = 0.9; % extract data within this many dva of target position (stim size)

myTR = 2.00;  % repetition time, in sec


u_subj = unique(cellfun(@(s) s(1:end-1),subj,'uniformoutput',0));



maxecc = 6; % dva from fixation
res = 101; % in x, y

[gridx,gridy] = meshgrid(linspace(-maxecc,maxecc,res),linspace(-maxecc,maxecc,res));
gridx = reshape(gridx,numel(gridx),1);gridy = reshape(gridy,numel(gridy),1);

trials_per_superrun = 90;


nblank = length(subj) * length(VOIs) * 4 * 90 ; % 3 superruns, 54 trials per superrun, 12 tpts
all_conds = nan(nblank,9);

all_subj = nan(nblank,1);
all_vois = nan(nblank,1);
all_targ_resp = nan(nblank,2); % PT & NPT

% n_conds x n_VOIs
startidx = 1;
for ss = 1:length(subj)
    
    this_subj_id = find(strcmpi(u_subj,subj{ss}(1:end-1)));
    
    for vv = 1:length(VOIs)
        
        fn = sprintf('%spri2stim_recons/%s_%s_hex%i_pos_coreg1.mat',root,subj{ss},VOIs{vv},hex_size);
        load(fn);
        
        thisidx = startidx:(startidx+size(recons_vec{1},1)-1);
        
        if startidx == 1
            %res = sqrt(size(recons_vec{1},2));
            all_recons = nan(nblank,res^2);
            [gridx,gridy] = meshgrid(linspace(-maxecc,maxecc,res),linspace(-maxecc,maxecc,res));
            gridx = reshape(gridx,numel(gridx),1);gridy = reshape(gridy,numel(gridy),1);
            targ_idx = sqrt((gridx-targ_pos(1)).^2 + (gridy-targ_pos(2)).^2) < targ_window;
        end
        
        %all_recons(thisidx,:) = recons_vec{1};
        all_conds(thisidx,:) = conds;
        
        all_subj(thisidx) = this_subj_id*ones(size(conds,1),1);
        all_vois(thisidx) = vv*ones(size(conds,1),1);
                
        all_targ_resp(thisidx,1) = mean(recons_vec{1}(:,targ_idx'),2);
        all_targ_resp(thisidx,2) = mean(recons_vec{2}(:,targ_idx'),2);
        
        
        

        clear whichhalf
        
        startidx = thisidx(end)+1;
    end
end

valididx = 1:(startidx-1);
%all_recons = all_recons(valididx,:);
all_conds  = all_conds(valididx,:);
%all_tpts = all_tpts(valididx);
all_subj = all_subj(valididx);
all_vois = all_vois(valididx);
all_targ_resp = all_targ_resp(valididx,:);

% for x-axis
mycontrasts = 100*[0.2 0.4 0.8];
attnstr = {'Attended','Unattended'};


%% plot CRF - collapsed over irrelevant contrasts, attended target vs. unattended target
% (at first, just sorting by contrast, no sorting by position)


uc1 = unique(all_conds(:,1));
uc2 = unique(all_conds(:,2));

% put subj in last dimension, easy to average over
all_coeffs_log10 = cell(length(VOIs),2);
all_coeffs_lin   = cell(length(VOIs),2);
all_stats_crfs = cell(length(VOIs),1);

% main effect of attention, contrast, and interaction for each ROI
F2_real = nan(3,length(VOIs));

X_forstats = nan(2*length(u_subj)*length(uc1),length(VOIs));

overall_mean = nan(length(VOIs),2,length(u_subj));

figure; ax = [];
for vv = 1:length(VOIs)
    ax(end+1) = subplot(2,length(VOIs),vv);
    hold on;
    
    % save things in an easy format for rm_anova2.m
    s_forstats = nan(2*length(u_subj)*length(uc1),1);  % s, a, c are all identical across ROIs
    a_forstats = nan(2*length(u_subj)*length(uc1),1);
    c_forstats = nan(2*length(u_subj)*length(uc1),1);
    stats_idx = 1;
    
    % attended; unattended
    for ii = 1:2
        myd = nan(length(u_subj),length(uc1));
        all_coeffs_log10{vv,ii} = nan(length(u_subj),2);
        all_coeffs_lin{vv,ii} = nan(length(u_subj),2);
        
        for ss = 1:length(u_subj)
            
            % save mean across contrasts for future power analyses
            overall_mean(vv,ii,ss) = mean(all_targ_resp(all_subj==ss & all_vois==vv,ii));
            
            
            for cc1 = 1:length(uc1)
                if ii == 1
                    thisidx = all_subj==ss & all_vois == vv & all_conds(:,1)==uc1(cc1);
                else
                    thisidx = all_subj==ss & all_vois == vv & all_conds(:,2) == uc2(cc1);
                end
                myd(ss,cc1) = mean(all_targ_resp(thisidx,ii));
                
                %plot(mycontrasts(cc1),myd(ss,cc1),'o','color',mycolors(ii,:));
                %plot(mycontrasts(cc1),myd(ss,cc1),'-','color',mycolors(ii,:));
                
                X_forstats(stats_idx,vv) = myd(ss,cc1);
                s_forstats(stats_idx)= ss;
                c_forstats(stats_idx) = cc1;
                a_forstats(stats_idx) = ii;
                stats_idx = stats_idx+1;
            end
            
            % fit here for each subj
            all_coeffs_log10{vv,ii}(ss,:) = polyfit(log10(mycontrasts),myd(ss,:),1);
            all_coeffs_lin{vv,ii}(ss,:)   = polyfit(mycontrasts,myd(ss,:),1);
            

            
        end
        for cc1 = 1:length(uc1)
            mym = mean(myd(:,cc1),1);
            mye = std(myd(:,cc1),[],1)/sqrt(length(u_subj));
            plot(mycontrasts(cc1)*[1 1],mym + mye*[-1 1],'-','Color',mycolors(ii,:),'LineWidth',1.5);
        end
        
        % do fit
        %pred = []

        coeffs_log10 = polyfit(log10(mycontrasts),mean(myd,1),1);
        coeffs_lin   = polyfit(mycontrasts,mean(myd,1),1);
        mcoeffs_log10{vv,ii} = coeffs_log10;
        mcoeffs_lin{vv,ii} = coeffs_lin;
        %plot(mycontrasts,myd.','-','LineWidth',0.75,'color',mycolors(ii,:));
        
        % LOG-LINEAR
        %plot(mycontrasts,coeffs_log10(1)*log10(mycontrasts)+coeffs_log10(2),'-','LineWidth',1.5,'color',mycolors(ii,:));
        
        % LINEAR
        myc_fit = linspace(mycontrasts(1),mycontrasts(end),101);
        plot(myc_fit,coeffs_lin(1)*myc_fit+coeffs_lin(2),'-','LineWidth',1.5,'color',mycolors(ii,:));
        
        plot(mycontrasts,mean(myd,1),'o','LineWidth',1.5,'color',mycolors(ii,:),'MarkerFaceColor',[1 1 1]);
        
        % and error bar...
        
        clear myd coeffs;
    end
    
    hold off;
    
    if vv == 1
        ylabel('Reconstruction activation (BOLD Z-score)');
        
    else
        set(gca,'YTickLabel',[]);
    end
    
    
    title(VOIs{vv});
    
    
    % do the 2-way ANOVA for factors of contrast and attended/unattended
    
    all_stats_crfs{vv} = rm_anova2(X_forstats(:,vv),s_forstats,a_forstats,c_forstats,{'priority','salience'});
    
    
    F2_real(:,vv) = [all_stats_crfs{vv}{2:4,5}];%[all_stats_crfs{vv}(2,5); all_stats_crfs{vv}(3,5); all_stats_crfs{vv}(4,5)];
    
    
    % TODO: add differences here too
    ax(end+1) = subplot(2,length(VOIs),vv+length(VOIs));
    hold on;
    myd = nan(length(u_subj),length(uc1));
    
    for cc1 = 1:length(uc1)
        for ss = 1:length(u_subj)
            
            thisidx = all_subj==ss & all_vois == vv & all_conds(:,1)==uc1(cc1) & all_conds(:,2) == uc2(cc1);
            myd(ss,cc1) = mean(all_targ_resp(thisidx,1) - all_targ_resp(thisidx,2));
            
            %plot(mycontrasts(cc1),myd(ss,cc1),'o','color',mycolors(ii,:));
            %plot(mycontrasts(cc1),myd(ss,cc1),'-','color',mycolors(ii,:));
            
        end
        % error bar        
        mym = mean(myd(:,cc1),1);
        mye = std(myd(:,cc1),[],1)/sqrt(length(u_subj));
        plot(mycontrasts(cc1)*[1 1],mym + mye*[-1 1],'-','Color',[0.2 0.2 0.2],'LineWidth',1.5);

    end
    %plot(mycontrasts,myd.','-','LineWidth',0.75,'color',[0.2 0.2 0.2]);
    plot(mycontrasts,mean(myd,1),'o-','LineWidth',1.5,'color',[0.2 0.2 0.2],'MarkerFaceColor',[1 1 1]);
    
    % and error bar...
    
    clear myd;

    hold off;
    if vv == 1
        ylabel('Attn - unattn');
        xlabel('Contrast (%)');
    else
        set(gca,'YTickLabel',[]);
    end

    
end

match_ylim(ax);
set(ax,'XScale','log','XLim',[15 100],'XTick',[20 40 80]);
set(gcf,'Name','Contrast response functions');


%% plot best-fit coefficients for log-linear fit  
% to keep things tidy, removing this from plot list for now
% coeff_str = {'Slope','Constant'};
% 
% 
% x_offset = 0.08*[1 -1];
% figure;
% for ii = 1:2 % for each coefficient (NOT attention condition!!!!)
%     subplot(2,1,ii); hold on;
%     for vv = 1:length(VOIs)
%         this_coeffs = [all_coeffs_log10{vv,1}(:,ii) all_coeffs_log10{vv,2}(:,ii)];
%         
%         
%         % draw single-subj lines
%         %plot(vv + x_offset,this_coeffs.','-','LineWidth',0.5,'Color',[0.5 0.5 0.5]);
%         
%         % now draw mean, error bars for each condition
%         for aa = 1:2
%             mye = std(this_coeffs(:,aa))/sqrt(length(u_subj));
%             myx = x_offset(aa)+vv;
%             plot(myx*[1 1],mean(this_coeffs(:,aa))+mye*[-1 1],'-','LineWidth',1.5,'Color',mycolors(aa,:));
%             %plot(myx,mcoeffs{vv,aa}(ii),'+','MarkerSize',5,'LineWidth',1.5,'Color',mycolors(aa,:));
%             plot(myx,mean(this_coeffs(:,aa)),'o','MarkerSize',5,'MarkerFaceColor',[1 1 1],'Color',mycolors(aa,:),'LineWidth',1.5);
%         end
%         
%     end
%     set(gca,'XTick',1:length(VOIs),'XTickLabel',VOIs);
%     ylabel('Reconstruction activation');
%     title(['Log-Linear: ' coeff_str{ii}]);
% end
% set(gcf,'Name','Best-fit log-linear coefficients');

%% plot best-fit coefficients for LINEAR fit  

coeff_str = {'Slope','Constant'};

% top row: salience, bottom row: priority
x_offset = 0.08*[1 -1]; % unattn, then attn (left to right)
figure;
for ii = 1:2 % for each coefficient (NOT attention condition!!!!)
    subplot(2,1,ii); hold on;
    for vv = 1:length(VOIs)
        this_coeffs = [all_coeffs_lin{vv,1}(:,ii) all_coeffs_lin{vv,2}(:,ii)];
        
        
        % draw single-subj lines
        %plot(vv + x_offset,this_coeffs.','-','LineWidth',0.5,'Color',[0.5 0.5 0.5]);
        
        % now draw mean, error bars for each condition
        for aa = 1:2
            mye = std(this_coeffs(:,aa))/sqrt(length(u_subj));
            myx = x_offset(aa)+vv;
            plot(myx*[1 1],mean(this_coeffs(:,aa))+mye*[-1 1],'-','LineWidth',1.5,'Color',mycolors(aa,:));
            %plot(myx,mcoeffs{vv,aa}(ii),'+','MarkerSize',5,'LineWidth',1.5,'Color',mycolors(aa,:));
            plot(myx,mean(this_coeffs(:,aa)),'o','MarkerSize',5,'MarkerFaceColor',[1 1 1],'Color',mycolors(aa,:),'LineWidth',1.5);
        end
        
    end
    set(gca,'XTick',1:length(VOIs),'XTickLabel',VOIs);
    ylabel('Reconstruction activation');
    title(['Linear: ' coeff_str{ii}]);
end
set(gcf,'Name','Best-fit linear coefficients');


%% shuffled version of parametric stats

% use X_forstats, s_forstats, a_forstats, c_forstats from above; and
% all_stats_crfs{vv} for "real" F-scores
%

F2_shuf = nan(3,length(VOIs),niter);


for ii = 1:niter
    
    this_d = nan(size(X_forstats));
    
    % shuffle data across ROIs within each subj
    for ss = 1:length(u_subj)
        this_idx = find(s_forstats==ss);
        
        shuf_idx = this_idx(randperm(length(this_idx)));
        
        this_d(this_idx,:) = X_forstats(shuf_idx,:);
 
    end
    
    for vv = 1:length(VOIs)
        tmpstats = rm_anova2(this_d(:,vv),s_forstats,a_forstats,c_forstats,{'priority','salience'});
        F2_shuf(:,vv,ii) = [tmpstats{2,5};tmpstats{3,5};tmpstats{4,5}];
        clear tmpstats;
    end
    
    clear this_d;
end

%% compute p-values, FDR threshold, and print

comparison_str = {'stimulus','contrast','interaction'};

p2_shuf = nan(3,length(VOIs));
fprintf('\n\n2-WAY ANOVA FOR EACH ROI');
for vv = 1:length(VOIs)
    fprintf('\n%s:\n',VOIs{vv});
    for ii = 1:3
        p2_shuf(ii,vv) = mean(squeeze(F2_shuf(ii,vv,:))>=F2_real(ii,vv));
        fprintf('%s:\tF = %0.03f, p = %0.03f\n',comparison_str{ii},F2_real(ii,vv),p2_shuf(ii,vv));
    end
end
fdr_thresh2 = fdr(p2_shuf(:),0.05);
fprintf('\nFDR threshold (across all comparisons): p  %0.03f\n',fdr_thresh2);

%% plot bar graph of F-scores for prioritity (2,5), salience (3,5) and interaction (4, 5)
% p's are 2,6, 3, 6 and 4,6

which_factors = [2 3 4];
factor_names = {'priority', 'salience','interaction'};

factor_colors = [146 54 71; 103 94 169; 8 146 85]./255;

% collate F, p into comparison x ROI matrix
allF = nan(length(which_factors),length(VOIs));
allP = nan(length(which_factors),length(VOIs));

for cc = 1:length(which_factors)
    for vv = 1:length(VOIs)
        allF(cc,vv) = all_stats_crfs{vv}{which_factors(cc),5};
        allP(cc,vv) = all_stats_crfs{vv}{which_factors(cc),6};
    end
    
    fdr_thresh(cc) = fdr(allP(cc,:),0.05);
    bonf_thresh(cc) = 0.05/length(VOIs);
    
end

figure;
for cc = 1:length(which_factors)
    subplot(length(which_factors),1,cc); hold on; % could also plot bars goign up/down; or grouped?
    for vv = 1:length(VOIs)
        if allP(cc,vv)<=fdr_thresh(cc)
            bar(vv,allF(cc,vv),'EdgeColor',factor_colors(cc,:),'FaceColor',factor_colors(cc,:),'LineWidth',1.5);
        else
            bar(vv,allF(cc,vv),'EdgeColor',factor_colors(cc,:),'FaceColor',[1 1 1],'LineWidth',1.5);
        end
    end
    ylabel(factor_names{cc});
    set(gca,'XTick',1:length(VOIs),'XTickLabel',VOIs);
end
match_ylim(get(gcf,'Children'));
set(gcf,'Name','CRF F-scores');



%% 3-way ANOVA
X_forstats3 = reshape(X_forstats,numel(X_forstats),1);
g_forstats3 = repmat([a_forstats c_forstats nan(size(a_forstats)) s_forstats],length(VOIs),1); % stimulus, contrast, voi (subj)
% insert VOI label
g_forstats3(:,3) = reshape(repmat(1:length(VOIs),size(X_forstats,1),1),numel(X_forstats),1);

tic;
F3_real = RMAOV33_out([X_forstats3 g_forstats3],0.05); % adapted version that spits out F-values
toc;
% me of attn, me of contrast, interaction [3, 5, 6, 7 are me VOI, attn x VOI, contrast x VOI, 3-way interaction



%% shuffled 3-way ANOVA
%
% concatenate all data across ROIs, shuffle all datapoints within subj
% niter times (RMAOV33)


% do randomization before parfor
F3_shuf = nan(niter,size(F3_real,2));
X_forstats3_shuf = nan(size(X_forstats3,1),niter);
for ii = 1:niter
       
    for ss = 1:length(u_subj)
        this_idx = find(g_forstats3(:,4)==ss);
        
        shuf_idx = this_idx(randperm(length(this_idx)));
        
        X_forstats3_shuf(this_idx,ii) = X_forstats3(shuf_idx);
        
    end

end

mypool = parpool(12);
parfor ii = 1:niter
    F3_shuf(ii,:) = RMAOV33_out([X_forstats3_shuf(:,ii) g_forstats3],0.05);    
end
delete(mypool);


p3_shuf = mean(bsxfun(@ge,F3_shuf,F3_real),1);

% print out P-values, just going to do this manually
fprintf('\n\n\n3-WAY ANOVA w/ 3-WAY INTERACTION\n');
fprintf('ME stimulus ID:\t\tF = %0.03f,\tp = %0.03f\n',F3_real(1),p3_shuf(1));
fprintf('ME contrast:\t\tF = %0.03f,\tp = %0.03f\n',F3_real(2),p3_shuf(2));
fprintf('ME VOI:\t\t\tF = %0.03f,\tp = %0.03f\n',F3_real(3),p3_shuf(3));
fprintf('stim x contrast:\tF = %0.03f,\tp = %0.03f\n',F3_real(4),p3_shuf(4));
fprintf('stim x VOI:\t\tF = %0.03f,\tp = %0.03f\n',F3_real(5),p3_shuf(5));
fprintf('contrast x VOI:\t\tF = %0.03f,\tp = %0.03f\n',F3_real(6),p3_shuf(6));
fprintf('3-way interaction:\tF = %0.03f,\tp = %0.03f\n',F3_real(7),p3_shuf(7));

if save_stats==1 

    
    %% save
    
    fn = fullfile(root,'pri2stim_stats',sprintf('n%i_binnedContrasts_2way_and_3way_ANOVA_shuffled_%iIter_%s.mat',length(u_subj),niter,datestr(now,30)));
    fprintf('saving to %s\n',fn);
    save(fn,'F2_shuf','p2_shuf','F3_real','F3_shuf','p3_shuf','fdr_thresh2','VOIs','subj','rand_seed');
    
end

return

