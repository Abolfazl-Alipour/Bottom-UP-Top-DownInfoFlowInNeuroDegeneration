
function []= compConnFromSubjID(sourceFilePath,subjectID,destinationPath,locPath,locPathFrwrdMdel,chckModelConsistancy,drug,varargin)
%
%
%          []= connFromSubjID(subjectID,varargin)
%
% Computes connectivity estimates from a subjectID using the specified
% connectivity model and stores the result in ALLEEG.CAT.Conn and then
% saves it into your destination path
%
% Input:
%
%   sourceFilePath          path to the file where you have your
%   
%
%
%       subjectID           ID number of the subject from the excel file
% 
%
%
%       destinationPath     path in which you want the resulting EEG
%                           structure to be saved
%       
%       locPath             path to your default channel locations in
%                           EEGLab
%
%
%       locPathFrwrdMdel    path to your channel locations in
%                           EEGLab's forward model 
%
%       modelConsistancy    check model consistency, default is 0
%   
%       
%       Drug                if on processes drug on group if off takes drug off gourp
%       
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% opens and closes eeglab to ge some of the functions loaded
eeglab;
close all
if drug=='off'
    load([sourceFilePath,int2str(subjectID),'_2_PD_REST.mat'])
elseif drug=='onn'
    load([sourceFilePath,int2str(subjectID),'_1_PD_REST.mat'])
elseif drug=='Cnt'
    load([sourceFilePath,int2str(subjectID),'_1_PD_REST.mat'])
end

%% If you want, you could do some of this to make life easier:
                
% ---------- Remove X,Y,Z & VEOG
EEG.VEOG=squeeze(EEG.data(64,:,:));
EEG.X=squeeze(EEG.data(65,:,:));
EEG.Y=squeeze(EEG.data(66,:,:));
EEG.Z=squeeze(EEG.data(67,:,:));
EEG.data=EEG.data(1:63,:,:);
EEG.nbchan=63;
EEG.chanlocs(67)=[];  EEG.chanlocs(66)=[]; EEG.chanlocs(65)=[]; EEG.chanlocs(64)=[];
% ---------- Add CPz
EEG = pop_chanedit(EEG,'append',63,'changefield',{64 'labels' 'CPz'});
EEG = pop_chanedit(EEG,'lookup', locPath);
% ---------- Re-Ref to Average Ref and recover CPz
EEG = pop_reref(EEG,[],'refloc',struct('labels',{'CPz'},'type',{''},'theta',{180},'radius',{0.12662},'X',{-32.9279},'Y',{-4.0325e-15},'Z',{78.363},...
    'sph_theta',{-180},'sph_phi',{67.208},'sph_radius',{85},'urchan',{64},'ref',{''}),'keepref','on');
% ---------- Remove CONSISTENLY BAD channels now that CPz has been reconstructed from the total
EEG.MASTOIDS = squeeze(mean(EEG.data([10,21],:,:),1));
EEG.data = EEG.data([1:4,6:9,11:20,22:26,28:64],:,:);
EEG.nbchan=60;
EEG.chanlocs(27)=[];  EEG.chanlocs(21)=[];   EEG.chanlocs(10)=[];   EEG.chanlocs(5)=[];  % Have to be in this order!
% ---------- Re-ref to average again now that the contaminated channels are gone
EEG = pop_reref(EEG,[]);
% ---------- Remove mean
EEG = pop_rmbase(EEG,[],[]);


%%      preprocessing

EEG2=EEG;
rawData=EEG2.data;

% save(['E:\Abolfazl\OtherProjs\PDPredictiveCoding\Data\restingStatesLastWave(hopefully)\DownlaodedDataset\ProcessedData\Data_',int2str(subjectID),'.mat'],'rawData')
% eeglab
% EEG = pop_importdata('dataformat','matlab','nbchan',0,'data','E:\\Abolfazl\\OtherProjs\\PDPredictiveCoding\\Data\\restingStatesLastWave(hopefully)\\DownlaodedDataset\\ProcessedData\\Data_894.mat','srate',500,'pnts',0,'xmin',0);EEG.setname=int2str(subjectID);EEG = eeg_checkset( EEG );
% EEG=EEG2;

%importing 
% choosing the selecting Eyes closed resting state
EEG = pop_rmdat( EEG, {'S  3'},[-0.001 1.999] ,0);EEG.setname=int2str(subjectID);EEG = eeg_checkset( EEG );pop_eegplot( EEG, 1, 1, 1);
% filtering 1 to 58 Hz
EEG = pop_eegfiltnew(EEG, 'locutoff',1,'hicutoff',58,'plotfreqz',1);EEG.setname=int2str(subjectID);EEG = eeg_checkset( EEG );
% cleaning artifacts automatically
EEG = clean_artifacts(EEG, 'FlatlineCriterion',5,'ChannelCriterion',0.8,'LineNoiseCriterion',4,'Highpass','off','BurstCriterion',20,'WindowCriterion',0.25,'BurstRejection','on','Distance','Euclidian','WindowCriterionTolerances',[-Inf 7] );EEG = eeg_checkset( EEG );pop_eegplot( EEG, 1, 1, 1);

%running ICA
EEG = pop_runica(EEG, 'icatype', 'runica', 'extended',1,'interrupt','on');EEG = eeg_checkset( EEG );

%checking dataset and running IC label
EEG = eeg_checkset( EEG );EEG = pop_iclabel(EEG, 'default');EEG = eeg_checkset( EEG );

%removing ICs that were recognized as noise
tmp=EEG.etc.ic_classification.ICLabel.classifications;
[i,j]=find(tmp==max(tmp,[],2));

keepCell={};
ii=1;
for kk=1:size(tmp,1)
    if j(kk)==1
        keepCell{ii}=i(kk);
        ii=ii+1;
        
    end
end
keepVec=cell2mat(keepCell);
EEG = pop_subcomp( EEG,keepVec , 0);EEG = eeg_checkset( EEG );


if drug=='off'
    save([destinationPath,'EEGObj_offDrug',int2str(subjectID),'.mat'],'EEG')
elseif drug=='onn'
    save([destinationPath,'EEGObj_OnDrug',int2str(subjectID),'.mat'],'EEG')
elseif drug=='Cnt'
    save([destinationPath,'EEGObj_Cntrl',int2str(subjectID),'.mat'],'EEG')
end

% if drug=='off'
%     load([destinationPath,'EEGObj_offDrug',int2str(subjectID),'.mat'],'EEG')
% elseif drug=='onn'
%     load([destinationPath,'EEGObj_OnDrug',int2str(subjectID),'.mat'],'EEG')
% elseif drug=='Cnt'
%     load([destinationPath,'EEGObj_Cntrl',int2str(subjectID),'.mat'],'EEG')
% end



% selecting Sources for SIFTING
EEG = pop_pre_prepData(EEG,'nogui','SignalType','Components'); aEEG = eeg_checkset( EEG );

% estimating model order, rounding to the nearest integer, if not rounded,
% it will generate instability in connectivity calculation
winLength=round(EEG.xmax-5,0);
WindowStepSizeSec=0.03;
VERBOSITY_LEVEL=2;

EEG = pop_est_selModelOrder(EEG,'nogui','modelingApproach', ...  %% I have commented out the pop up window asking for model fitting inside the pop_est_selModelOder
        {'Segmentation VAR'     ...
        'algorithm' {'Vieira-Morf'} ...
        'winStartIdx' []    ...
        'winlen'  winLength   ...
        'winstep' WindowStepSizeSec  ...
        'taperfcn' 'blackmanharris'  ...
        'epochTimeLims' []      ...
        'prctWinToSample' 100   ...
        'normalize' {'method' {'time' 'ensemble'}} ...
        'detrend' {'method' 'constant'} ...
        'verb' VERBOSITY_LEVEL},      ...
        'morderRange',[1 10] ,  ...
        'downdate',true,        ...
        'runPll',[],            ...
        'icselector',{'hq'},  ...  %we are using hq to find the optimal model order
        'winStartIdx',[],       ...
        'epochTimeLims',[],     ...
        'prctWinToSample',100,   ...
        'plot', [], ...
        'verb', VERBOSITY_LEVEL);
    
    


%fitting the model
optimalModelOrder=EEG.CAT.IC.hq.popt(1);
EEG = pop_est_fitMVAR(EEG,'nogui','Algorithm','Vieira-Morf','ModelOrder',optimalModelOrder,'winlen' , winLength);



%checking model consistency
if chckModelConsistancy==1
    stats=est_checkMVARConsistency(EEG,EEG.CAT.MODEL);
    if stats<.85
        disp(['consitency is p' stats])
    else
        disp(['consistency is' stats])
    end
end

%estimating connectivity
EEG= pop_est_mvarConnectivity(EEG,'nogui','connmethods',{'RPDC'}); %{'DTF','dDTF','ffDTF','RPDC','Coh','pCoh'}


save([destinationPath,'EEGObj_Component',int2str(subjectID),'.mat'],'EEG')
ffDTF=EEG.CAT.Conn.ffDTF;
dDTF=EEG.CAT.Conn.dDTF;
RPDC=EEG.CAT.Conn.RPDC;
if drug=='off'
    save([destinationPath,'CompffDTFMatrx_',int2str(subjectID),'OFFDRUG','.mat'],'ffDTF')
    save([destinationPath,'CompdDTFMatrx_',int2str(subjectID),'OFFDRUG','.mat'],'dDTF')
    save([destinationPath,'CompRPDCMatrx_',int2str(subjectID),'OFFDRUG','.mat'],'RPDC')
elseif drug=='onn'
    save([destinationPath,'CompffDTFMatrx_',int2str(subjectID),'ONDRUG','.mat'],'ffDTF')
    save([destinationPath,'CompdDTFMatrx_',int2str(subjectID),'ONDRUG','.mat'],'dDTF')
    save([destinationPath,'CompRPDCMatrx_',int2str(subjectID),'ONDRUG','.mat'],'RPDC')
elseif drug=='Cnt'
    save([destinationPath,'CompffDTFMatrx_',int2str(subjectID),'.mat'],'ffDTF')
    save([destinationPath,'CompdDTFMatrx_',int2str(subjectID),'.mat'],'dDTF')
    save([destinationPath,'CompRPDCMatrx_',int2str(subjectID),'.mat'],'RPDC')
end

end















 
 
 %%%        graveyard %%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  tmp(:,2:end)=tmp(:,2:end)-rejectionThreshold;
% [i,j]=find(tmp()==max(tmp,[],2));
% 
% 
% keepCell={};
% ii=1;
% for kk=1:size(tmp,1)
%     if tmp(kk,1)>(sum(tmp(kk,2:7))-rejectionThreshold)
%         keepCell{ii}=kk;
%         ii=ii+1
%     end
% end
% keepVec=cell2mat(keepCell);
% EEG = pop_subcomp( EEG,keepVec , 0);EEG = eeg_checkset( EEG );
% 
% est_selModelOrder(EEG,'morder',[1 100])
% 'winStartIdx'
% 
% 
% 
% 
% 
% 
% 
% 
% [EEG,com]=pop_importdata( 'setname',int2str(subjectID),'data',['E:\Abolfazl\OtherProjs\PDPredictiveCoding\Data\restingStatesLastWave(hopefully)\DownlaodedDataset\ProcessedData\Data_',int2str(subjectID),'.mat'],...
%     'dataformat','matlab')
% pop_topoplot(EEG, 0, [1:13] ,'895 pruned with ICA',[4 4] ,0,'electrodes','on')
% 
% 
%  itpop_eegplot( EEG, 1, 1, 1);EEG = pop_rmdat( EEG, {'S  3'},[-0.1 1.9] ,0);EEG.setname='895';EEG = eeg_checkset( EEG );EEG = eeg_checkset( EEG );pop_eegplot( EEG, 1, 1, 1);pop_eegplot( EEG, 1, 1, 1);EEG = clean_artifacts(EEG, 'FlatlineCriterion',5,'ChannelCriterion',0.8,'LineNoiseCriterion',4,'Highpass','off','BurstCriterion',20,'WindowCriterion',0.25,'BurstRejection','on','Distance','Euclidian','WindowCriterionTolerances',[-Inf 7] );EEG = eeg_checkset( EEG );pop_eegplot( EEG, 1, 1, 1);EEG = eeg_checkset( EEG );EEG = pop_runica(EEG, 'icatype', 'runica', 'extended',1,'interrupt','on');EEG = eeg_checkset( EEG );EEG = pop_iclabel(EEG, default);EEG = eeg_checkset( EEG );EEG = pop_icflag(EEG, [NaN NaN;0.9 1;0.9 1;0.9 1;0.9 1;0.9 1;0.9 1]);EEG = eeg_checkset( EEG );EEG = pop_icflag(EEG, [NaN NaN;0.75 1;0.75 1;0.75 1;0.75 1;0.75 1;0.75 1]);EEG = eeg_checkset( EEG );EEG = pop_iclabel(EEG, default);EEG = eeg_checkset( EEG );pop_eegplot( EEG, 1, 1, 1);pop_topoplot(EEG, 0, [1:13] ,'895 pruned with ICA',[4 4] ,0,'electrodes','on');EEG = pop_subcomp( EEG, [1], 0);EEG = eeg_checkset( EEG );
% 
%  
%  
%      'icselector',{'hq'},'modelingApproach',         ...
%         {'Segmentation VAR'     ...
%         'algorithm' {'Vieira-Morf'} ...
%         'winStartIdx' []    ...
%         'winlen'  EEG.xmax-5    ...
%         'winstep' WindowStepSizeSec  ...
%         'taperfcn' 'blackmanharris'  ...
%         'epochTimeLims' []      ...
%         'prctWinToSample' 100   ...
%         'normalize' {'method' {'time' 'ensemble'}} ...
%         'detrend' {'method' 'constant'} ...
%         'verb' VERBOSITY_LEVEL},      ...
%         'morderRange',[1 30] ,  ...
%         'downdate',true,        ...
%         'runPll',[],            ...
%         'icselector',{'sbc' 'aic' 'fpe' 'hq'},  ...
%         'winStartIdx',[],       ...
%         'epochTimeLims',[],     ...
%         'prctWinToSample',100,   ...
%         'plot', [], ...
%         'verb', VERBOSITY_LEVEL);
%     
%     
%     
%     
%     
%     
%     
%     
%     
%     
%     
%     
% 
% 
% 
% 
% ,'icselector',{'hq'},'winlen',winLength);,'icselector','hq','winlen',(EEG.xmax-5));EEG = eeg_checkset( EEG );
% 
% EEG.etc.eeglabvers = '2019.1'; % this tracks which version of EEGLAB is being used, you may ignore itpop_eegplot( EEG, 1, 1, 1);EEG = pop_rmdat( EEG, {'S  3'},[-0.1 1.9] ,0);EEG.setname='895';EEG = eeg_checkset( EEG );EEG = eeg_checkset( EEG );pop_eegplot( EEG, 1, 1, 1);pop_eegplot( EEG, 1, 1, 1);EEG = clean_artifacts(EEG, 'FlatlineCriterion',5,'ChannelCriterion',0.8,'LineNoiseCriterion',4,'Highpass','off','BurstCriterion',20,'WindowCriterion',0.25,'BurstRejection','on','Distance','Euclidian','WindowCriterionTolerances',[-Inf 7] );EEG = eeg_checkset( EEG );pop_eegplot( EEG, 1, 1, 1);EEG = eeg_checkset( EEG );EEG = pop_runica(EEG, 'icatype', 'runica', 'extended',1,'interrupt','on');EEG = eeg_checkset( EEG );EEG = pop_iclabel(EEG, default);EEG = eeg_checkset( EEG );EEG = pop_icflag(EEG, [NaN NaN;0.9 1;0.9 1;0.9 1;0.9 1;0.9 1;0.9 1]);EEG = eeg_checkset( EEG );EEG = pop_icflag(EEG, [NaN NaN;0.75 1;0.75 1;0.75 1;0.75 1;0.75 1;0.75 1]);EEG = eeg_checkset( EEG );EEG = pop_iclabel(EEG, default);EEG = eeg_checkset( EEG );EEG = eeg_checkset( EEG );EEG = pop_iclabel(EEG, default);EEG = eeg_checkset( EEG );EEG = pop_icflag(EEG, [NaN NaN;0.9 1;0.9 1;0.9 1;0.9 1;0.9 1;0.9 1]);EEG = eeg_checkset( EEG );EEG = pop_iclabel(EEG, default);EEG = eeg_checkset( EEG );pop_eegplot( EEG, 1, 1, 1);EEG = pop_iclabel(EEG, default);EEG = eeg_checkset( EEG );EEG = pop_iclabel(EEG, default);EEG = eeg_checkset( EEG );EEG = pop_iclabel(EEG, default);EEG = eeg_checkset( EEG );pop_topoplot(EEG, 0, [1:38] ,'895 pruned with ICA',[6 7] ,0,'electrodes','on');pop_eegplot( EEG, 1, 1, 1);EEG = pop_forwardModel(EEG, 'C:\Users\aalipour\Documents\MATLAB\eeglab2019_1\plugins\dsi\headModel\resources\head_modelColin27_8003_Standard-10-5-Cap339.mat', [0.33       0.022        0.33], 0,1);EEG = pop_rsbl(EEG, 1, 1, 'power', 'bsbl', 1);EEG.etc.eeglabvers = '2019.1'; % this tracks which version of EEGLAB is being used, you may ignore itEEG = pop_pre_prepData(EEG);EEG = eeg_checkset( EEG );EEG = eeg_checkset( EEG );EEG = pop_pre_prepData(EEG);EEG = eeg_checkset( EEG );EEG = pop_est_selModelOrder(EEG,0);EEG = eeg_checkset( EEG );EEG = pop_pre_prepData(EEG);EEG = eeg_checkset( EEG );EEG = pop_est_selModelOrder(EEG,0);EEG = pop_est_selModelOrder(EEG,0);
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% EEG.etc.eeglabvers = '2019.1'; % this tracks which version of EEGLAB is being used, you may ignore itpop_eegplot( EEG, 1, 1, 1);EEG = pop_rmdat( EEG, {'S  3'},[-0.1 1.9] ,0);EEG.setname='895';EEG = eeg_checkset( EEG );EEG = eeg_checkset( EEG );pop_eegplot( EEG, 1, 1, 1);pop_eegplot( EEG, 1, 1, 1);EEG = clean_artifacts(EEG, 'FlatlineCriterion',5,'ChannelCriterion',0.8,'LineNoiseCriterion',4,'Highpass','off','BurstCriterion',20,'WindowCriterion',0.25,'BurstRejection','on','Distance','Euclidian','WindowCriterionTolerances',[-Inf 7] );EEG = eeg_checkset( EEG );pop_eegplot( EEG, 1, 1, 1);EEG = eeg_checkset( EEG );EEG = pop_runica(EEG, 'icatype', 'runica', 'extended',1,'interrupt','on');EEG = eeg_checkset( EEG );EEG = pop_iclabel(EEG, default);EEG = eeg_checkset( EEG );EEG = pop_icflag(EEG, [NaN NaN;0.9 1;0.9 1;0.9 1;0.9 1;0.9 1;0.9 1]);EEG = eeg_checkset( EEG );EEG = pop_icflag(EEG, [NaN NaN;0.75 1;0.75 1;0.75 1;0.75 1;0.75 1;0.75 1]);EEG = eeg_checkset( EEG );EEG = pop_iclabel(EEG, default);EEG = eeg_checkset( EEG );pop_eegplot( EEG, 1, 1, 1);pop_topoplot(EEG, 0, [1:13] ,'895 pruned with ICA',[4 4] ,0,'electrodes','on');EEG = pop_subcomp( EEG, [1], 0);EEG = eeg_checkset( EEG );EEG = eeg_checkset( EEG );
% 
% 
% 
% 
% [EEG,com]=pop_importdata( 'setname',int2str(subjectID),'data',['E:\Abolfazl\OtherProjs\PDPredictiveCoding\Data\restingStatesLastWave(hopefully)\DownlaodedDataset\ProcessedData\Data_',int2str(subjectID),'.mat'],...
%     'dataformat','matlab')
% pop_topoplot(EEG, 0, [1:13] ,'895 pruned with ICA',[4 4] ,0,'electrodes','on')
% 
% pop_eegplot( EEG, 1, 1, 1);EEG = pop_rmdat( EEG, {'S  3'},[-0.1 1.9] ,0);EEG.setname='895';EEG = eeg_checkset( EEG );EEG = eeg_checkset( EEG );pop_eegplot( EEG, 1, 1, 1);pop_eegplot( EEG, 1, 1, 1);EEG = clean_artifacts(EEG, 'FlatlineCriterion',5,'ChannelCriterion',0.8,'LineNoiseCriterion',4,'Highpass','off','BurstCriterion',20,'WindowCriterion',0.25,'BurstRejection','on','Distance','Euclidian','WindowCriterionTolerances',[-Inf 7] );EEG = eeg_checkset( EEG );pop_eegplot( EEG, 1, 1, 1);EEG = eeg_checkset( EEG );EEG = pop_runica(EEG, 'icatype', 'runica', 'extended',1,'interrupt','on');EEG = eeg_checkset( EEG );EEG = pop_iclabel(EEG, default);EEG = eeg_checkset( EEG );EEG = pop_icflag(EEG, [NaN NaN;0.9 1;0.9 1;0.9 1;0.9 1;0.9 1;0.9 1]);EEG = eeg_checkset( EEG );EEG = pop_icflag(EEG, [NaN NaN;0.75 1;0.75 1;0.75 1;0.75 1;0.75 1;0.75 1]);EEG = eeg_checkset( EEG );EEG = pop_iclabel(EEG, default);EEG = eeg_checkset( EEG );pop_eegplot( EEG, 1, 1, 1);pop_topoplot(EEG, 0, [1:13] ,'895 pruned with ICA',[4 4] ,0,'electrodes','on');EEG = pop_subcomp( EEG, [1], 0);EEG = eeg_checkset( EEG );
% 
% 
%  
%  
%  
%  
%  
%  
%  
%  
%  
%  
%  