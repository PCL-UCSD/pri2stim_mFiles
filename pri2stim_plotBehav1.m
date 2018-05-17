% pri2stim_plotBehav1.m
% TCS 8/5/2016

function pri2stim_plotBehav1(subj)
% plots behavioral data (accuracy, RT) as a function of attended,
% unattended stim contrast
%
% also average coherence presented across a session?
% also acc through time? coh through time? 

if nargin < 1
    subj = {'AI151','AP151','AR151','AS151','BA151','BB151','BC151','BF151'};
end


% sort subj alphabetically - so all resampling, etc, is the same
subj = sort(subj);

root = load_root;

rng_seed = load_rng_seed;
niter = 1000; % # of iterations for shuffled 


% seed random # generator
rng(rng_seed);

% all subj had 360 trials
nblank = 360*length(subj);

all_conds = nan(nblank,9);
all_rt = nan(nblank,1);
all_resp = nan(nblank,1);
all_r = nan(nblank,2);
all_correct = nan(nblank,1);

all_targOnset = nan(nblank,2);


all_stimLocs{1} = nan(nblank,2);
all_stimLocs{2} = nan(nblank,2);

all_subj = nan(nblank,1);


all_m_acc = nan(length(subj),1); % mean accuracy overall


startidx = 1;

for ss = 1:length(subj)
    
    fn = fullfile(root,'pri2stim_behav',sprintf('%s_pri2stim_behav.mat',subj{ss}));
    fprintf('loading from %s...\n',fn);
    thisbehav = load(fn);
    
    thisidx = startidx:(startidx+size(thisbehav.all_conds,1)-1);
    
    all_conds(thisidx,:) = thisbehav.all_conds;
    all_rt(thisidx) = thisbehav.all_rt;
    all_resp(thisidx) = thisbehav.all_resp;
    all_correct(thisidx) = thisbehav.all_correct;
    
    all_r(thisidx,:) = thisbehav.all_r;
    all_targOnset(thisidx,:) = thisbehav.all_targOnset;
    
    all_stimLocs{1}(thisidx,:) = thisbehav.all_stimLocs{1};
    all_stimLocs{2}(thisidx,:) = thisbehav.all_stimLocs{2};
    
    all_subj(thisidx) = ss;
    
    all_m_acc(ss) = mean(thisbehav.all_correct);
    
    startidx = thisidx(end)+1;
    
    clear thisbehav fn thisidx;
    
end

fprintf('\n\nmean acc: %f, err: %f\n\n',mean(all_m_acc),std(all_m_acc)/sqrt(length(subj)));

mycontrasts = [20 40 80];



%% sort by contrast of target, distracter and average acc, RT


uc1 = unique(all_conds(:,1)); % attended stim contrast level 1:3
uc2 = unique(all_conds(:,2)); % unattended stim contrast level 1:3

all_acc_sorted = nan(length(mycontrasts),length(mycontrasts),length(subj));
all_rt_sorted = nan(length(mycontrasts),length(mycontrasts),length(subj));

all_coh_sorted = nan(length(mycontrasts),length(subj)); % just sort by target contrast

for ss = 1:length(subj)
    for cc1 = 1:length(uc1) % cols: attended
        
        
        thisidx_coh = all_subj==ss & all_conds(:,1)==uc1(cc1);
        all_coh_sorted(cc1,ss) = mean(all_conds(thisidx_coh,8));
        clear thisidx_coh;
        
        for cc2 = 1:length(uc2)  % rows: unattended
            
            thisidx = all_subj==ss & all_conds(:,1)==uc1(cc1) & all_conds(:,2)==uc2(cc2);
            
            all_acc_sorted(cc2,cc1,ss) = mean(all_correct(thisidx));
            all_rt_sorted(cc2,cc1,ss) = nanmean(all_rt(thisidx));

            
        end
    end
    
    
end


%% plot acc, RT as a function of attended, unattended contrast w/ error bars

figure;

x_offset = 0.2*[-1 0 1];

m_acc = 100*mean(all_acc_sorted,3);
m_rt = mean(all_rt_sorted,3);

e_acc = std(100*all_acc_sorted,[],3)/sqrt(length(subj));
e_rt  = std(all_rt_sorted,[],3)/sqrt(length(subj));


% first acc
subplot(1,5,[1 2]); hold on;
for cc1 = 1:length(uc1)
    
   % plot(cc1+x_offset,squeeze(100*all_acc_sorted(:,cc1,:)).','-','Color',[0.5 0.5 0.5]);
    
    for cc2 = 1:length(uc2)
        
        plot(cc1+x_offset(cc2)*[1 1],m_acc(cc2,cc1)+[-1 1]*e_acc(cc2,cc1),'k-');
        
    end
    
    
    
    plot(cc1+x_offset,m_acc(:,cc1),'ko-','MarkerSize',5,'MarkerFaceColor','w');
end
xlim([0 4]);
ylim([0 100]);
set(gca,'XTick',[x_offset + 1 x_offset + 2 x_offset + 3],'XTickLabel',{[] '20%',[],[],'40%',[],[],'80%',[]},...
    'YTick',0:25:100);

ylabel('Accuracy (%)');
xlabel('Attended contrast');

% then RT
subplot(1,5,[3 4]); hold on;
for cc1 = 1:length(uc1)
    
    %plot(cc1+x_offset,squeeze(all_rt_sorted(:,cc1,:)).','-','Color',[0.5 0.5 0.5]);
    
    for cc2 = 1:length(uc2)
        
        plot(cc1+x_offset(cc2)*[1 1],m_rt(cc2,cc1)+[-1 1]*e_rt(cc2,cc1),'k-');
        
    end
    
    
    
    plot(cc1+x_offset,m_rt(:,cc1),'ko-','MarkerSize',5,'MarkerFaceColor','w');
end

xlim([0 4]);
ylim([0 1.75]);
set(gca,'XTick',[x_offset + 1 x_offset + 2 x_offset + 3],'XTickLabel',{[] '20%',[],[],'40%',[],[],'80%',[]},...
    'YTick',0:.5:2);

ylabel('Response time (s)');
xlabel('Attended contrast');


subplot(1,5,5); hold on;
mc = 100*mean(all_coh_sorted,2);
ec = std(100*all_coh_sorted,[],2)/sqrt(length(subj));
hold on;

for cc1 = 1:length(uc1)
    plot(cc1*[1 1],mc(cc1)+[-1 1]*ec(cc1),'k-');
end
plot(uc1,mc,'ko-','MarkerSize',5,'MarkerFaceColor','w');

xlim([0 4]);
ylim([0 100]);
xlabel('Attended contrast');
ylabel('Coherence (%)');
set(gca,'XTick',uc1,'XTickLabel',{'20%','40%','80%'},'YTick',[0 25 50 75 100]);

hold off;



%% stats - 2way ANOVA (raw)
nblank = length(uc1)*length(uc2)*length(subj);
my_subj = nan(nblank,1);
my_c1 = nan(nblank,1);
my_c2 = nan(nblank,1);
my_d = nan(nblank,2); % acc & RT

idx = 1;
for ss = 1:length(subj)
    for cc1 = 1:length(uc1) 
        for cc2 = 1:length(uc2)
            
            my_subj(idx) = ss;
            my_c1(idx) = cc1;
            my_c2(idx) = cc2;
            my_d(idx,:) = [all_acc_sorted(cc2,cc1,ss) all_rt_sorted(cc2,cc1,ss)];
            
            idx = idx+1;
        end
    end
end

[tmp_s,tmp_c] = meshgrid(1:length(subj),1:length(uc1));
my_d_coh = all_coh_sorted(:);
my_subj_coh = tmp_s(:);
my_c1_coh = tmp_c(:); 
clear tmp_s tmp_c;

stats_acc = rm_anova2(my_d(:,1),my_subj,my_c1,my_c2,{'attn','unattn'})
stats_rt  = rm_anova2(my_d(:,2),my_subj,my_c1,my_c2,{'attn','unattn'})

F_real_coh = RMAOV1([my_d_coh, my_c1_coh, my_subj_coh],0.05);

% save the F-values for comparison to shuffled distribution below
F_real = [stats_acc{2,5} stats_rt{2,5}; stats_acc{3,5} stats_rt{3,5}; stats_acc{4,5} stats_rt{4,5}];

% shuffle data within subj niter times
all_F_behav = cell(3,1);
all_F_behav{1} =  nan(niter,2); % acc, RT - main effect of attended contrast
all_F_behav{2} =  nan(niter,2); % acc, RT - main effect of unattended contrast
all_F_behav{3} =  nan(niter,2); % acc, RT - interaction

all_F_coh = nan(niter,1); % only 1 main effect

for ii = 1:niter
    
    % this is where we'll put the shuffled data
    this_d = nan(size(my_d));
    this_d_coh = nan(size(my_d_coh));
    
    % on every iteration, shuffle datapoints within subj
    for ss = 1:length(subj)
        
        this_idx = find(my_subj==ss);
        
        shuf_idx = this_idx(randperm(length(this_idx)));
        
        this_d(this_idx,:) = my_d(shuf_idx,:);
        clear this_idx shuf_idx;
        
        % same for coherence
        this_idx = find(my_subj_coh==ss);
        shuf_idx = this_idx(randperm(length(this_idx)));
        
        this_d_coh(this_idx,:) = my_d_coh(shuf_idx,:);
    end
    
    stats_acc_shuff = rm_anova2(this_d(:,1),my_subj,my_c1,my_c2,{'attn','unattn'});
    stats_rt_shuff  = rm_anova2(this_d(:,2),my_subj,my_c1,my_c2,{'attn','unattn'});
    
    all_F_behav{1}(ii,:) = [stats_acc_shuff{2,5} stats_rt_shuff{2,5}]; 
    all_F_behav{2}(ii,:) = [stats_acc_shuff{3,5} stats_rt_shuff{3,5}]; 
    all_F_behav{3}(ii,:) = [stats_acc_shuff{4,5} stats_rt_shuff{4,5}]; 
    
    all_F_coh(ii) = RMAOV1_noOut([this_d_coh, my_c1_coh my_subj_coh],0.05);
    
    clear this_d stats_acc_shuff stats_rt_shuff;
    
end

data_str = {'Accuracy','RT'};
comparison_str = {'Main effect of attended contrast','Main effect of unattended contrast','Interaction'};

all_p_shuff = nan(3,2); % ME1, ME2, interaction x acc, RT
for jj = 1:2    % sorry - these are out-of-order because it prints more intuitively this way
    fprintf('\n%s\n',data_str{jj});
    for ii = 1:3

        all_p_shuff(ii,jj) = mean(all_F_behav{ii}(:,jj) >= F_real(ii,jj));
        fprintf('%s:\t F = %0.03f, p = %0.03f\n',comparison_str{ii},F_real(ii,jj),all_p_shuff(ii,jj));
    end
end
p_coh_shuf = mean(all_F_coh>=F_real_coh);
fprintf('\nCoherence:\t F = %0.03f, p  = %0.03f\n',F_real_coh,p_coh_shuf);

return