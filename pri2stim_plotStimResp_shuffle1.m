function pri2stim_plotStimResp_shuffle1(subj,VOIs,hex_size)
%
% adapted from wmDrop_plotTargetResp_resample1.m, TCS 12/17/2015
%
% plots extracted signal from reconstruction at stimulus position, sorted
% by target & distractor contrast, and potentially by separation distance
% as well
%
% % TCS 9/5/2016, note: no resampling anywhere in here, should likely be
% renamed

rand_seed = load_rng_seed;
rng(rand_seed);

save_stats = 1;

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

targ_pos = [3.5 0]; % reconstructions are aligned to this point
targ_window = 0.9; % extract data within this many dva of target position (stim size)

myTR = 2.00;  % repetition time, in sec

load(fullfile(root,'pri2stim_colors.mat'));

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
        fprintf('loading %s\n',fn);
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


mycontrasts = 100*[0.2 0.4 0.8];
attnstr = {'Attended','Unattended'};


%% plot CRF - matched contrasts, attended target vs. unattended target
% (at first, just sorting by contrast, no sorting by position)

uc1 = unique(all_conds(:,1));
uc2 = unique(all_conds(:,2));

figure; ax = [];
for vv = 1:length(VOIs)
    ax(end+1) = subplot(2,length(VOIs),vv);
    hold on;
    
    % attended; unattended
    for ii = 1:2
        myd = nan(length(u_subj),length(uc1));
        
        for cc1 = 1:length(uc1)
            for ss = 1:length(u_subj)
                
                thisidx = all_subj==ss & all_vois == vv & all_conds(:,1)==uc1(cc1) & all_conds(:,2) == uc2(cc1);
                myd(ss,cc1) = mean(all_targ_resp(thisidx,ii));
                
                %plot(mycontrasts(cc1),myd(ss,cc1),'o','color',mycolors(ii,:));
                %plot(mycontrasts(cc1),myd(ss,cc1),'-','color',mycolors(ii,:));
                
                
                
            end
            
            mym = mean(myd(:,cc1),1);
            mye = std(myd(:,cc1),[],1)/sqrt(length(u_subj));
            plot(mycontrasts(cc1)*[1 1],mym + mye*[-1 1],'-','Color',mycolors(ii,:),'LineWidth',1.5);

            
        end
        
        
        % do fit
        %pred = []

        coeffs = polyfit(log10(mycontrasts),mean(myd,1),1);
        
        %plot(mycontrasts,myd.','-','LineWidth',0.75,'color',mycolors(ii,:));
        

        plot(mycontrasts,coeffs(1)*log10(mycontrasts)+coeffs(2),'-','LineWidth',1.5,'color',mycolors(ii,:));
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
set(gcf,'Name','Matched contrasts');


%% do stats - mean contrast 
%% all data (combine across separation conditions for now)
%
% n_ROIs columns, 2 rows (attended, unattended)
% within each subplot, 3 sets of lines next to one another, corresponding
% to the "un/attended contrast".
% each of those lines has 3 data points - corresponding to un/attended
% contrast
%
% data and labels for stats (attended, sorted by attended; unattended contrast)
d_forstats_attn = cell(length(VOIs),1);
g_forstats_attn = cell(length(VOIs),1); % 3 colums: attended contrast, unattended contrast, subj ID

F_real_attn = nan(3,length(VOIs)); % actual F-scores from 2-way ANOVA for each main effect & interaction (rows)  & ROI (col)
F_shuf_attn = nan(3,length(VOIs),niter); % shuffled F-scores

d_forstats_unattn = cell(length(VOIs),1);
g_forstats_unattn = cell(length(VOIs),1); % 3 colums: attended contrast, unattended contrast, subj ID

F_real_unattn = nan(3,length(VOIs)); % actual F-scores from 2-way ANOVA for each main effect & interaction (rows)  & ROI (col)
F_shuf_unattn = nan(3,length(VOIs),niter); % shuffled F-scores

figure; ax = [];
for vv = 1:length(VOIs)
    ax(end+1) = subplot(1,length(VOIs),vv);
    hold on;
    
    % attended; unattended
    %for ii = 1:2
    myd = nan(length(uc2),length(uc1),length(u_subj));

    
    d_forstats_attn{vv} = nan(length(uc2)*length(uc1)*length(u_subj),1);
    g_forstats_attn{vv} = nan(length(uc2)*length(uc1)*length(u_subj),3);
    statsidx = 1;
    
    % attended contrast
    for cc1 = 1:length(uc1)
        
        % unattended contrast
        for cc2 = 1:length(uc2)
            
            for ss = 1:length(u_subj)
                
                thisidx = all_subj==ss & all_vois == vv & all_conds(:,1)==uc1(cc1) & all_conds(:,2) == uc2(cc2);
                myd(cc2,cc1,ss) = mean(all_targ_resp(thisidx,1)); % attended
                
                %plot(mycontrasts(cc1),myd(ss,cc1),'o','color',mycolors(ii,:));
                
                % plug into stats variable
                d_forstats_attn{vv}(statsidx) = myd(cc2,cc1,ss);
                g_forstats_attn{vv}(statsidx,:) = [cc1 cc2 ss];
                statsidx = statsidx+1;
                
            end
            mym = mean(myd(cc2,cc1,:),3);
            mye = std(myd(cc2,cc1,:),[],3)/sqrt(length(u_subj));
            plot(cc2+(cc1-1)*length(uc2)*[1 1],mym + mye*[-1 1],'-','Color',mycolors(1,:),'LineWidth',1);
        end
        
        % plot the relevant line
        
        
        plot((1:length(uc2))+(cc1-1)*length(uc2),mean(myd(:,cc1,:),3),'o-','Color',mycolors(1,:),'LineWidth',1,'MarkerFaceColor',[1 1 1]);
        
    end
    %plot(mycontrasts,myd.','-','LineWidth',0.75,'color',mycolors(ii,:));
    %plot(mycontrasts,mean(myd,3),'-','LineWidth',1.5,'color',mycolors(ii,:));
    
    clear myd;
    %end
    
    hold off;
    
    
    set(gca,'XTick',[2 5 8],'XTickLabel',{'20%','40%','80%'});
    if vv == 1
        ylabel('Attended stimulus activation (BOLD Z-score)');
        xlabel('Attended contrast (%)');
    else
        set(gca,'YTickLabel',[]);
    end
    
    title(VOIs{vv});
    
    
    % compute real F-score (like pri2stim_plotBehav1.m)
    tmpstats = rm_anova2(d_forstats_attn{vv},g_forstats_attn{vv}(:,3),g_forstats_attn{vv}(:,1),g_forstats_attn{vv}(:,2),{'attn','unattn'});   % stats_acc_shuff = rm_anova2(this_d(:,1),my_subj,my_c1,my_c2,{'attn','unattn'});
    
    
    F_real_attn(1,vv) = tmpstats{2,5};
    F_real_attn(2,vv) = tmpstats{3,5};
    F_real_attn(3,vv) = tmpstats{4,5};
    
    clear tmpstats;
    
end


%% do shuffled stats - use same randomization for each ROI (so analysis does not depend on # ROIs)

% n_datapoints x n_ROIs
d_cat_attn = horzcat(d_forstats_attn{:});

% do shuffled stats
for ii = 1:niter
    this_d = nan(size(d_cat_attn));
    for ss = 1:length(u_subj)
        this_idx = find(g_forstats_attn{1}(:,3)==ss);
        
        shuf_idx = this_idx(randperm(length(this_idx)));
        
        this_d(this_idx,:) = d_cat_attn(shuf_idx,:);
        
    end
    
    % compute shuffled F for each ROI
    for vv = 1:length(VOIs)
        tmpstats = rm_anova2(this_d(:,vv),g_forstats_attn{vv}(:,3),g_forstats_attn{vv}(:,1),g_forstats_attn{vv}(:,2),{'attn','unattn'});   % stats_acc_shuff = rm_anova2(this_d(:,1),my_subj,my_c1,my_c2,{'attn','unattn'});
        
        F_shuf_attn(1,vv,ii) = tmpstats{2,5};
        F_shuf_attn(2,vv,ii) = tmpstats{3,5};
        F_shuf_attn(3,vv,ii) = tmpstats{4,5};
        
        
        clear tmpstats;
    end
end

%% compute and print out p-values

p_shuf_attn = nan(3,length(VOIs));
comparison_str = {'ME of attended contrast','ME of unattn contrast   ','interaction          '};
fprintf('\n\nAttended stimulus representation\n');
for vv = 1:length(VOIs)
    fprintf('\n%s\n',VOIs{vv});
    for ii = 1:3
        p_shuf_attn(ii,vv) = mean(squeeze(F_shuf_attn(ii,vv,:))>=F_real_attn(ii,vv)); 
        fprintf('%s\tF = %2.03f,\tp  %0.03f\n',comparison_str{ii},F_real_attn(ii,vv),p_shuf_attn(ii,vv));
    end
end



ax_a = ax;
%match_ylim(ax);

set(gcf,'Name','Attended response');





%% do stats - attended responses (2-way, factors of attended & unattended contrast
%  TCS 9/5/2016 - moved this above
% % FACTOR 1: attended contrast
% % FACTOR 2: unattended contrast
% for vv = 1:length(VOIs)
%     
%     % the data structure for RMAOV2
%     myX = nan(length(uc1)*length(uc2)*length(u_subj),4);
%     myidx = 1;
%     
%     
%     for cc1 = 1:length(uc1)
%         for cc2 = 1:length(uc2)
%             for ss = 1:length(u_subj)
%                 
%                 thisidx = all_subj==ss & all_vois == vv & all_conds(:,1)==uc1(cc1) & all_conds(:,2) == uc2(cc2);
%                 myX(myidx,:) = [mean(all_targ_resp(thisidx,1)) cc1 cc2 ss];
%                 myidx = myidx+1;
%             end
%         end
%         
%     end
%     
% 
%     fprintf('******************************\nROI: %s\n\n',VOIs{vv});
%     RMAOV2(myX);
%     clear myX;
% end
% 


%% same, but for unattended responses, sorted by attended contrast


figure; ax = [];
for vv = 1:length(VOIs)
    ax(end+1) = subplot(1,length(VOIs),vv);
    hold on;
    
    % attended; unattended
    %for ii = 1:2
    myd = nan(length(uc2),length(uc1),length(u_subj));

        
    d_forstats_unattn{vv} = nan(length(uc2)*length(uc1)*length(u_subj),1);
    g_forstats_unattn{vv} = nan(length(uc2)*length(uc1)*length(u_subj),3);
    statsidx = 1;
    
    % unattended contrast
    for cc2 = 1:length(uc2)
        
        
        % attended contrast
        for cc1 = 1:length(uc1)
            
            
            for ss = 1:length(u_subj)
                
                thisidx = all_subj==ss & all_vois == vv & all_conds(:,1)==uc1(cc1) & all_conds(:,2) == uc2(cc2);
                myd(cc1,cc2,ss) = mean(all_targ_resp(thisidx,2)); % unattended
                
                % plug into stats variable (note: now UNATTN, ATTN, SUBJ)
                d_forstats_unattn{vv}(statsidx) = myd(cc1,cc2,ss);
                g_forstats_unattn{vv}(statsidx,:) = [cc2 cc1 ss];
                statsidx = statsidx+1;
                
            end
            mym = mean(myd(cc1,cc2,:),3);
            mye = std(myd(cc1,cc2,:),[],3)/sqrt(length(u_subj));
            plot(cc1+(cc2-1)*length(uc1)*[1 1],mym + mye*[-1 1],'-','Color',mycolors(2,:),'LineWidth',1);
            

            
        end
        
        % plot the relevant line
        
        
        plot((1:length(uc1))+(cc2-1)*length(uc1),mean(myd(:,cc2,:),3),'o-','Color',mycolors(2,:),'LineWidth',1,'MarkerFaceColor',[1 1 1]);
        
    end
    %plot(mycontrasts,myd.','-','LineWidth',0.75,'color',mycolors(ii,:));
    %plot(mycontrasts,mean(myd,3),'-','LineWidth',1.5,'color',mycolors(ii,:));
    
    clear myd;
    %end
    
    hold off;
    
    
    set(gca,'XTick',[2 5 8],'XTickLabel',{'20%','40%','80%'});
    if vv == 1
        xlabel('Unattended contrast (%)');
        ylabel('Unattended stimulus activation (BOLD Z-score)');
    else
        set(gca,'YTickLabel',[]);
    end
    
    title(VOIs{vv});    % compute real F-score (like pri2stim_plotBehav1.m)
    tmpstats = rm_anova2(d_forstats_unattn{vv},g_forstats_unattn{vv}(:,3),g_forstats_unattn{vv}(:,1),g_forstats_unattn{vv}(:,2),{'unattn','attn'});   % stats_acc_shuff = rm_anova2(this_d(:,1),my_subj,my_c1,my_c2,{'attn','unattn'});
    
    
    F_real_unattn(1,vv) = tmpstats{2,5};
    F_real_unattn(2,vv) = tmpstats{3,5};
    F_real_unattn(3,vv) = tmpstats{4,5};
    
    clear tmpstats;
    
end

match_ylim([ax_a;ax]);

set(gcf,'Name','Unattended response');


%% do shuffled stats - use same randomization for each ROI (so analysis does not depend on # ROIs)

% n_datapoints x n_ROIs
d_cat_unattn = horzcat(d_forstats_unattn{:});

% do shuffled stats
for ii = 1:niter
    this_d = nan(size(d_cat_unattn));
    for ss = 1:length(u_subj)
        this_idx = find(g_forstats_unattn{1}(:,3)==ss);
        
        shuf_idx = this_idx(randperm(length(this_idx)));
        
        this_d(this_idx,:) = d_cat_unattn(shuf_idx,:);
        
    end
    
    % compute shuffled F for each ROI
    for vv = 1:length(VOIs)
        tmpstats = rm_anova2(this_d(:,vv),g_forstats_unattn{vv}(:,3),g_forstats_unattn{vv}(:,1),g_forstats_unattn{vv}(:,2),{'unattn','attn'});   % stats_acc_shuff = rm_anova2(this_d(:,1),my_subj,my_c1,my_c2,{'attn','unattn'});
        
        F_shuf_unattn(1,vv,ii) = tmpstats{2,5};
        F_shuf_unattn(2,vv,ii) = tmpstats{3,5};
        F_shuf_unattn(3,vv,ii) = tmpstats{4,5};
        
        
        clear tmpstats;
    end
end

%% compute and print out p-values

p_shuf_unattn = nan(3,length(VOIs));
comparison_str = {'ME of unattn contrast ','ME of attended contrast','interaction             '};
fprintf('\n\nUnattended stimulus representation\n');
for vv = 1:length(VOIs)
    fprintf('\n%s\n',VOIs{vv});
    for ii = 1:3
        p_shuf_unattn(ii,vv) = mean(squeeze(F_shuf_unattn(ii,vv,:))>=F_real_unattn(ii,vv)); 
        fprintf('%s\tF = %2.03f,\tp  %0.03f\n',comparison_str{ii},F_real_unattn(ii,vv),p_shuf_unattn(ii,vv));
    end
end


%% FDR threshold (concatenate interaction across attended & unattn)
fdr_thresh = fdr([p_shuf_attn(3,:) p_shuf_unattn(3,:)],0.05);
fprintf('FDR trheshold (interaction across both stimuli): %0.03f\n',fdr_thresh);




if save_stats == 1
    fn = fullfile(root,'pri2stim_stats',sprintf('n%i_fullyCrossed_2wayANOVA_shuffled_%iIter_%s.mat',length(u_subj),niter,datestr(now,30)));
    fprintf('saving to %s\n',fn);
    save(fn,'p_shuf_attn','p_shuf_unattn','F_real_attn','F_real_unattn','fdr_thresh','VOIs','subj','rand_seed');
end


%% do stats - unattended responses (2-way, factors of attended & unattended contrast

% FACTOR 1: unattended contrast
% FACTOR 2: aattended contrast
% for vv = 1:length(VOIs)
%     
%     % the data structure for RMAOV2
%     myX = nan(length(uc1)*length(uc2)*length(u_subj),4);
%     myidx = 1;
%     
%     
%     for cc1 = 1:length(uc1)
%         for cc2 = 1:length(uc2)
%             for ss = 1:length(u_subj)
%                 
%                 thisidx = all_subj==ss & all_vois == vv & all_conds(:,1)==uc1(cc1) & all_conds(:,2) == uc2(cc2);
%                 myX(myidx,:) = [mean(all_targ_resp(thisidx,2)) cc2 cc1 ss];
%                 myidx = myidx+1;
%             end
%         end
%         
%     end
%     
% 
%     fprintf('******************************\nROI: %s\n\n',VOIs{vv});
%     RMAOV2(myX);
%     clear myX;
% end







% TODO: stats w/ separation condition as a factor

%% plot ROI vs. ROI stim activation-o-grams
%
% is stimulus activation for attended/unattended stimulus in one ROI
% related to that in another ROI? to mean activation? to something else?
%
% take all activation across trials within a subj/ROI, correlate to same
% trials in another ROI, same subj - then z-score, average, and un-z-score
% for plotting
% (shoudl move this to a separate function soon...)

% for attended representation activation vs. attended representation
% activation, so it's symmetric
all_corr_attn = nan(length(VOIs),length(VOIs),length(u_subj));

for ss = 1:length(u_subj)
    
    % we want a trials x VOIs matrix of activation values
    
    this_activation = nan(sum(all_subj==ss & all_vois==1),length(VOIs));
    
    for vv = 1:length(VOIs)
        this_activation(:,vv) = all_targ_resp(all_subj==ss & all_vois==vv,1);
    end
    
    
    all_corr_attn(:,:,ss) = atanh(corr(this_activation));
end

figure;
imagesc(mean(all_corr_attn,3));axis square;colorbar;
set(gca,'XTick',1:length(VOIs),'XTickLabel',VOIs,'XTickLabelRotation',-45,'YTick',1:length(VOIs),'YTickLabel',VOIs);


return

