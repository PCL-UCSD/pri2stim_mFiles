
function pri2stim_compileData_hexMap(subs,vois)
% adapted from wmDrop_compileData_hexMap


if nargin<1
    subs = {'BB151'};
end

if nargin < 2
    vois = {'LH-V1','LH-V2','LH-V3','LH-V4','LH-V3A','LH-IPS0','LH-IPS1','LH-IPS2','LH-IPS3','RH-V1','RH-V2','RH-V3','RH-V4','RH-V3A','RH-IPS0','RH-IPS1','RH-IPS2','RH-IPS3','LH-sPCS','RH-sPCS'};
end

rootanat = '/usr/local/serenceslab/tommy/sot5/';
root = '/usr/local/serenceslab/tommy/pri2stim/'; 

loc_vmp = 'wmLazLocAll.vmp'; 
loc_name = 'wmLazLocAll'; % for output mat file


mr.rp=['hexMap'];

mr.tbtPRT = ['doPrt_pri2stim_map'];



mr.extendedTal = [0];
mr.funcExt = '_SCCAI_3DMCT_LTR_THP3c.vtc';
mr.funcType = '.vtc';    



ext = [''];%
%conds = 1:36;
mr.Resolution = 3;
cvoi=1;
chrf=0;
saveBestNVoxels=0;

%general input params, all of which should stay the same across subjects
mr.vNum = 2;

%timing stuff
mr.NofTRs = 171;        
mr.TR = 2000;

%processing params
mr.vthresh = 1;     %how many 1x1x1mm voxels must be in each 3x3x3mm VTC voxel
mr.intThresh = inf; %exclude voxels with a mean intensity across scans below this level
mr.stdThresh = inf; %exclude voxels with a std > this number (applied after timeseries normalization)
mr.meanSig = 2;     %0=none, 1=psc, 2=zscore % JOHN SUGGESTED TRYING ZSCORE INSTEAD OF PSC
mr.rmMo = 1;        %project out influence of motion based estimated motion from BV
mr.rmMia = 0;
mr.vOrMS = 1;
mr.saveRaw = 1;     %1=yes, 2=no
mr.saveBehav = 0;   %input rt data from special prt files? (typically this should be 0)


%if 1, then assumes positivie t-values correspond to RH regions, and
%negative values to LH regions, otherwise masking is carried out on JUST
%the positive values (so only 'red' voxels will be included)
mr.vmpLeftRight =0;
 
mr.plotFmrVOI = 0;

default_mr = mr;

%%% TRYING OUT individual trial, ERA - CHANGE BACK!
%%% TODO: make a PRT/run this where each location = separate condition, avg
%%% across runs

%fill in the hrf structure
hrf.vNum = 1;
hrf.notes = {};
hrf.fPlot = 0;
hrf.fit = 0;
hrf.dType = 4;      %4=mod, 5=mia  
hrf.estType = 3;    %1=separate, 2=separate, then average, 3=concat
hrf.hrfType = 1;    %1=decon, 2=era, 3=gamma model
hrf.avgVoxels = 0;
hrf.sHRF = 0;
hrf.eHRF = 12;

hrf.rThresh = 0;
hrf.doBehav = 0;

%fitting stuff
hrf.tr = mr.TR/1000;
hrf.dur = 20;
hrf.t = 0:hrf.tr:hrf.dur;
hrf.dt = 2;
hrf.tau = [1.2500 1.5];
hrf.n = [3 5];
hrf.a = [1 1];
hrf.minPeakTime = 2;
hrf.gamdur = 20;
hrf.evtDur = 1;
hrf.gam = zeros(size(hrf.t));%[0 0 7.6889 11.2610 2.0300 -4.2775 -5.5696 -4.4249 -2.8442 -1.6128 -0.8399 -0.4109];

%loop over subjects and compile the data into output mat files
for s=1:length(subs)      
    
    searchStr = sprintf('%s%s/%s*',root,subs{s},mr.rp);
    dd = dir(searchStr);
    if strfind(dd(end).name,'Loc')
        dd = dd(1:end-1);
    end

    
    runRange = 1:sum([dd.isdir]);
    
    mr.subName = char(subs{s});
    mr.subNameOrig = mr.subName;
    if length(mr.subName) < 3
        mr.subNameVOI = [mr.subName '2'];
    else 
  
        mr.subNameVOI = [mr.subName(1:2) '2'];
  
    end
    
    mr.expPath = [root, mr.subName];
    
    mr_reset = mr;
    
    for v = 1:length(vois)

        mr = mr_reset; % reset mr
        

        mr.outMatName = sprintf('%spri2stim_mrStructs/%s_%s_VT%i_Int%i_Norm%i_RmMo%i_Vol_%i_%s_%s.mat',root,mr.subName,loc_name,mr.vthresh,mr.intThresh,mr.meanSig,mr.rmMo,mr.vOrMS,mr.rp,vois{v});
        
        
        if cvoi
            
            % BELOW IS A HACK FOR NOW!!!!!
            mr.vmpName=[rootanat, mr.subNameVOI, '/Anat/', mr.subName(1:2), '_' loc_vmp]; % set up the vmp name  --- circle loc
          
            
            % for now, all subj are returning, so don't need this...
            % directories
            %if ismember(mr.subName,{'AF61','AE61','AD61','AC61'})
            %    mr.voiName = [root   , mr.subNameVOI, '/Anat/VOIs/', mr.subNameVOI, '_',vois{v}, '.voi']; % set up the voi from tommys
            %else
                mr.voiName = [rootanat, mr.subNameVOI, '/Anat/VOIs/', mr.subNameVOI, '_',vois{v}, '.voi']; % set up the voi from tommys
            %end
            
            MapInfo = BVQXfile(mr.vmpName);%readVMP(mr.vmpName, 0);
            mr.statThresh = MapInfo.Map.LowerThreshold;%ThreshCrit;
            mr.vmpResolution = MapInfo.Resolution;
            mr.VMPData = MapInfo.Map.VMPData;
            %figure out how many runs of this type were performed
            mr.r = runRange;

            mr.funcFile = [root, mr.subName, '/', mr.rp, num2str(mr.r(1)), '/', mr.subName, '_', mr.rp, num2str(mr.r(1)), mr.funcExt];
            
            % ADDED HERE on 7.12.2012 - read in a vtc file to get the
            % actual dimensions. this will accomodate your unusual file
            % sizes...
            vtc = BVQXfile(mr.funcFile);
            mr.vmpBoxStartX = vtc.XStart;
            mr.vmpBoxStartY = vtc.YStart;
            mr.vmpBoxStartZ = vtc.ZStart;

            mr.vtcBoxStartX = vtc.XStart;
            mr.vtcBoxStartY = vtc.YStart;
            mr.vtcBoxStartZ = vtc.ZStart;            
                   
            mr.DimX = (vtc.XEnd - vtc.XStart)/vtc.Resolution;
            mr.DimY = (vtc.YEnd - vtc.YStart)/vtc.Resolution;
            mr.DimZ = (vtc.ZEnd - vtc.ZStart)/vtc.Resolution;                
            vtc.ClearObject;    % clear the vtc object.
            mr.rmMean = 0;
            %process the data
            compVOI(mr);
        end
        
        % estimate the HRF for each region
        if chrf     
            %then compute the estimated HRFs 
            load(mr.outMatName);
            reWritePrt = 1;
            % check to make sure we have a prt with only 1 condition, where
            % that condition has ALL trials in it.
            
            prt = BVQXfile(mr.prtName);
            len = cellfun('length',mr.prtName);
            usedPRTs = find(len);
            numConditions = prt{usedPRTs(1)}.NrOfConditions;
            if numConditions~=1
%            if prt.NrOfConditions~=1
                %evalstr = [mr.hrfPRT, '({''', mr.subName, '''});']
                evalstr = [this_hrf_PRT, '({''', mr.subName, '''});']
                %eval([mr.hrfPRT, '({''', mr.subName, '''});']);
                eval(evalstr);
                for i=1:size(prt,2)
                    if ~isempty(mr.prtName{i})
                        disp(sprintf('clearing prt %i',i));
                        prt{i}.ClearObject;
                    end
                end
                [mr.sTimes, mr.seTimes, mr.avgEvtDur, mr.condNames] = evtTimes(mr);
                prt = BVQXfile(mr.prtName); % re-read them
                reWritePrt = 1;             % set a flag to switch prt files back to TBT version
            end

            hrf.selConds = 1:length(mr.condNames);%numConditions; 
            
            % by some god-knows-what way load the right fuckign prts...
            
           %for i=1:size(prt,2)
           %    if ~isempty(mr.prtName{i})
           %        prt{i}.ClearObject;
           %    end
           %end
            
            hrf.selVois = 1:length(mr.voiNames);            
            %compute hrf on only the best voxels
            if saveBestNVoxels
                for v=1:length(mr.voiNames);
                    %find the best n voxels based on their sensory response
                    if findstr(char(mr.voiNames{v}), 'LH')
                        [x,y]=sort(mr.voiData(v).statValue);
                    else
                        [x,y]=sort(mr.voiData(v).statValue, 'descend');
                    end                    
            %        y=y(mr.voiData(v).goodv_IntThresh);
                    %just keep the best sensory voxels
                    if length(y)<saveBestNVoxels
                        mr.voiData(v).mod=mr.voiData(v).mod(:,y);            
                    else
                        mr.voiData(v).mod=mr.voiData(v).mod(:,y(1:saveBestNVoxels));
                    end
                end        
            end
            hrf.voiNames = {''};
            for nam=1:length(hrf.selVois)
                hrf.voiNames{nam} = mr.voiNames{hrf.selVois(nam)};
            end
            cntn=1;
            for nam=hrf.selConds
                hrf.condNames{cntn} = mr.condNames{nam};
                cntn=cntn+1;
            end
            hrf.outMatName = [mr.outMatName(1:end-4), '_HRF.mat'];
            %computethe HRFs
            compHRF(mr,hrf); 
            
            if reWritePrt
                cmd = [mr.tbtPRT, '({''', mr.subName, '''});'];
                disp(sprintf('rewriting PRT with command: %s',cmd));
            	eval(cmd);%[mr.tbtPRT, '({''', mr.subName, '''});']);   
            end
            
        end
    end % end of VOIs
end









