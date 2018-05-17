function pri2stim_plotReconstructions_posCoreg1(subj,VOIs,hex_size)
% adapted from wmDrop_plotReconstructions_posCoreg1
%
% TCS 12/16/2015 (no tpts)
close all;
if nargin < 1
    subj = {'AI151','AP151','AR151','AS151','BA151','BB151','BC151','BF151'};
    %subj = {'BA151'};
end

if nargin < 2
    VOIs = {'V1','V2','V3','V3A','V4','IPS0','IPS1','IPS2','IPS3','sPCS'}; % FEFnew
    %VOIs = {'V1','V3A','IPS0','IPS2','sPCS','SuperWMDrop'}; % FEFnew
    %VOIs = {'V1','IPS0'};
end


if nargin < 3
    hex_size = 7;
end

root = '/mnt/neurocube/local/serenceslab/tommy/pri2stim/';

clim_mode = 2; % % 2 = within-ROI, 1 = across-ROIs

myTR = 2.00;  % repetition time, in sec

%mycolors = [0 127 66; 56 9 124; 28 69 68]./255;
%load('/usr/local/serenceslab/tommy/wmDrop/wmDrop_colors.mat');



u_subj = unique(cellfun(@(s) s(1:end-1),subj,'uniformoutput',0));

save_figs = 0;
fig_path = '/mnt/neurocube/local/serenceslab/tommy/pri2stim/figs/manuscript_firstDraft/';

maxecc = 6; % dva from fixation
res = 101; % in x, y

[gridx,gridy] = meshgrid(linspace(-maxecc,maxecc,res),linspace(-maxecc,maxecc,res));
gridx = reshape(gridx,numel(gridx),1);gridy = reshape(gridy,numel(gridy),1);

trials_per_superrun = 90;

%tpts_of_interest = [2 3 4]; % to start with...

nblank = length(subj) * length(VOIs) * 4 * trials_per_superrun; % 3 superruns, 54 trials per superrun, 12 tpts
all_recons{1} = nan(nblank,res^2);
all_recons{2} = nan(nblank,res^2);
all_conds = nan(nblank,9);
all_subj = nan(nblank,1);
all_vois = nan(nblank,1);

% n_conds x n_VOIs
startidx = 1;
for ss = 1:length(subj)
    
    this_subj_id = find(strcmpi(u_subj,subj{ss}(1:end-1)));
    
    for vv = 1:length(VOIs)
        
        fn = sprintf('%spri2stim_recons/%s_%s_hex%i_pos_coreg1.mat',root,subj{ss},VOIs{vv},hex_size);
        fprintf('loading %s...\n',fn);
        load(fn);
        
        thisidx = startidx:(startidx+size(recons_vec{1},1)-1);
        
        
        all_recons{1}(thisidx,:) = recons_vec{1};
        all_recons{2}(thisidx,:) = recons_vec{2};
        all_conds(thisidx,:) = conds;
        %all_tpts(thisidx) = tpts;
        all_subj(thisidx) = this_subj_id*ones(size(conds,1),1);
        all_vois(thisidx) = vv*ones(size(conds,1),1);
        
        
        startidx = thisidx(end)+1;
        clear whichhalf
    end
end

valididx = 1:(startidx-1);
all_recons{1} = all_recons{1}(valididx,:);
all_recons{2} = all_recons{2}(valididx,:);
all_conds  = all_conds(valididx,:);
all_subj = all_subj(valididx);
all_vois = all_vois(valididx);

ax = [];
% for each ROI, 2 figures, 1 each for the 2 separation conditions
% each figure is 3x3 (column: attended contrast, 0 deg polar angle, row:
% unattended contrast, changes w/ polar angle pos)

su = unique(all_conds(:,4)); % 1, 2

c1u = unique(all_conds(:,1)); % attended
c2u = unique(all_conds(:,2)); % distractor

ROI_clim = nan(length(VOIs),2);

for vv = 1:length(VOIs)
    for sepidx = 1:length(su)
        rf(vv,sepidx) = figure;
        set(gcf,'Name',sprintf('%s - sep %i',VOIs{vv},sepidx),'NumberTitle','off');
        for ii = 2
            for cc1 = 1:length(c1u)
                for cc2 = 1:length(c2u)
                    
                    
                    ax(end+1) = subplot(length(c2u),length(c1u),cc1 + (cc2-1)*(length(c1u)));% sepidx+(vv-1)*length(su));
                    this_recon = nan(length(u_subj),size(all_recons{ii},2));
                    
                    % first average reconstructions for each subject, then average over
                    % subjects (so that different numbers of sessions don't contribute
                    % to differces in overall average) - this is how wmdelay plotting
                    % worked
                    for ss = 1:length(u_subj)
                        
                        thisidx = all_subj==ss & all_conds(:,1)==c1u(cc1) & all_conds(:,2)==c2u(cc2) & all_conds(:,4) == su(sepidx) & all_vois==vv ;
                        this_recon(ss,:) = mean(all_recons{ii}(thisidx,:),1);
                        
                    end
                    
                    mm = mean(this_recon,1);
                    imagesc(reshape(mm,res,res));colormap parula;
                    
%                     if vv == 1 && sepidx == ceil(length(su)/2)
%                         title(sprintf('%s - %s',condstr{cc},targstr{ii}));
%                     end
%                     
%                     if sepidx == 1
%                         ylabel(VOIs{vv});
%                     end
%                     
                     axis square xy;
                    
                    clear this_recon;
                    
                end
            end
        end
%         if length(VOIs)==9
             set(gcf,'Position',[680   621   508   474]);
%         end
    end
    
    if clim_mode == 2
        axs = get(rf(vv,:),'Children');
        axs = vertcat(axs{:});
        this_clim = get(axs,'CLim');
        this_clim = cell2mat(this_clim);
        ROI_clim(vv,:) = [min(this_clim(:,1)),max(this_clim(:,2))];
        set(axs,'CLim',ROI_clim(vv,:));
    end
    
end
if clim_mode == 1
    all_clim = cell2mat(get(ax,'CLim'));
    set(ax,'CLim',[min(all_clim(:,1)) max(all_clim(:,2))]);
end

set(ax,'Box','on','XTick',[],'YTick',[]);

% if clim_mode == 2, make a plot showing all the clims together
if clim_mode == 2
    figure; hold on;
    for vv = 1:length(VOIs)
        plot([vv vv],ROI_clim(vv,:),'-','LineWidth',3,'Color',[0 0 0]);
    end
    set(gca,'XTick',1:length(VOIs),'XTickLabel',VOIs);
    ylabel('BOLD Z-score');
    title('Color range');
end



if save_figs == 1
    for vv = 1:length(VOIs)
        for sepidx = 1:length(su)
            fn2s = sprintf('%sreconstructions/n%i_%s_sep%i',fig_path,length(u_subj),VOIs{vv},sepidx);
            if clim_mode == 2
                fn2s = [fn2s '_CLimROI'];
            end
            fprintf('saving figure to %s\n',fn2s);
            saveas(rf(vv,sepidx),[fn2s '.fig'],'fig');
            saveas(rf(vv,sepidx),[fn2s '.eps'],'eps');
        end
    end
end

return
