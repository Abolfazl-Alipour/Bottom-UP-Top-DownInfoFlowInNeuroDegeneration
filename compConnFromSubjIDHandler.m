

subjectID=909;
destinationPath='PATH TO EEG OBJECTS';
sourceFilePath='PATH TO RAW DATA';
locPath=('PATH TO\MATLAB\eeglab2019_1\plugins\dipfit\standard_BESA\standard-10-5-cap385.elp');
locPathFrwrdMdel='PATH TO\MATLAB\eeglab2019_1\plugins\dsi\headModel\resources\head_modelColin27_8003_Standard-10-5-Cap339.mat';
chckModelConsistancy=0;



upperFrequencyBound=58;
% CntrlList=[894,908,8010,906,903,8060,893,909,911,895,913,900,896,899,914,910,890,891,912,905,904,892,902,901,898,897,907];
CntrlList=[908,8010,8060,893,909,900,914,910,891,892,902];


for i=1:length(CntrlList)
    subjectID=CntrlList(i);
    compConnFromSubjID(sourceFilePath,subjectID,destinationPath,locPath,locPathFrwrdMdel,chckModelConsistancy)
end


% PDList=[801	802	803	804	805	806	807	808	809	810	811	813	814	815	816	817	818	819	820	821	822	823	824	825	];%826	827	828	829];
PDListOffDrug=[802	803	806	807	808	813	816	817	819	823	824	827	828	829];%826	;
PDList=[801	802	803	804	805	806	807	808	809	810	811	813	815	816	817	818	819	820	821	822	823	824	825	827	828	829];

% PDRPDCMat=nan(length(PDList),upperFrequencyBound);
for i=1:length(PDListOffDrug)
    subjectID=PDListOffDrug(i);
    compConnFromSubjID(sourceFilePath,subjectID,destinationPath,locPath,locPathFrwrdMdel,chckModelConsistancy)
end



PDListOnDrug=[801	804	805	809	810	811	815	818	820	821	822	825]%	826];

for i=1:length(PDListOnDrug)
    subjectID=PDListOnDrug(i);
    compConnFromSubjID(sourceFilePath,subjectID,destinationPath,locPath,locPathFrwrdMdel,chckModelConsistancy)
end



PDListOnDrugCntrl=[894	906 903	911 895 913 899 890 912 905 904 901];%	826];

for i=1:length(PDListOnDrugCntrl)
    subjectID=PDListOnDrugCntrl(i);
    compConnFromSubjID(sourceFilePath,subjectID,destinationPath,locPath,locPathFrwrdMdel,chckModelConsistancy)
end







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

upperFrequencyBound=58;
CntrlList=[894,908,8010,906,903,8060,893,909,911,895,913,900,896,899,914,910,890,891,912,905,904,892,902,901,898,897,907,908,8010,8060,893,909,900,914,910,891,892,902];

CntrlRPDCMat=nan(length(CntrlList),upperFrequencyBound);
for i=1:length(CntrlList)
    subjectID=CntrlList(i);
    RPDCFileName=strcat('ffDTFMatrx_',int2str(subjectID),'.mat');
    load(strcat(destinationPath,RPDCFileName))
    
    tmpTime=mean(ffDTF(:,:,1:58,:),4);
    infOutflow=nan(1,size(tmpTime,3));
    
    for layerNum=1:size(tmpTime,3)
        tmp=tmpTime(:,:,layerNum);
        tmp(1:(size(tmp,1)+1):end)=0; %zeroing out the diagonal values
        infOutflow(1,layerNum)=sum(sum(tmp));
    end
    CntrlRPDCMat(i,:)=infOutflow;
end


PDListonDrug=[801	802	803	804	805	806	807	808	809	810	811	813	814	815	816	817	818	819	820	821	822	823	824	825	802	803	806	807	808	813	816	817	819	823	824	827	828	829];%826	;];%826	827	828	829];

PDRPDCMat=nan(length(PDListonDrug),upperFrequencyBound);
for i=1:length(PDListonDrug)
    subjectID=PDListonDrug(i);
    RPDCFileName=strcat('ffDTFMatrx_',int2str(subjectID),'.mat');
    load(strcat(destinationPath,RPDCFileName))
    
    tmpTime=mean(ffDTF(:,:,1:58,:),4);
    infOutflow=nan(1,size(tmpTime,3));
    
    for layerNum=1:size(tmpTime,3)
        tmp=tmpTime(:,:,layerNum);
        tmp(1:(size(tmp,1)+1):end)=0; %zeroing out the diagonal values
        infOutflow(1,layerNum)=sum(sum(tmp));
    end
    PDRPDCMat(i,:)=infOutflow;
end



















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%AFTERMATH

PDList=[801	802	803	804	805	806	807	808	809	810	811	813	815	816	817	818	819	820	821	822	823	824	825	827	828	829];
drug='off';
% PDRPDCMat=nan(length(PDList),upperFrequencyBound);
for i=1:length(PDList)
    subjectID=PDList(i);
    compConnFromSubjID(sourceFilePath,subjectID,destinationPath,locPath,locPathFrwrdMdel,chckModelConsistancy,drug)
end

drug='onn';
% PDRPDCMat=nan(length(PDList),upperFrequencyBound);
for i=1:length(PDList)
    subjectID=PDList(i);
    compConnFromSubjID(sourceFilePath,subjectID,destinationPath,locPath,locPathFrwrdMdel,chckModelConsistancy,drug)
end


CntrlList=[894	908	8010	906	903	8060	893	909	911	895	913	900		899	914	910	890	891	912	905	904	892	902	901		897		907 ];
% 8070 only 5 seconds, excluded

drug='Cnt';
% PDRPDCMat=nan(length(PDList),upperFrequencyBound);
for i=1:length(CntrlList)
    subjectID=CntrlList(i);
    compConnFromSubjID(sourceFilePath,subjectID,destinationPath,locPath,locPathFrwrdMdel,chckModelConsistancy,drug)
end


CntrlList=[894	908	8010	906	903	8060	893	909	911	895	913	900		899	914	910	890	891	912	905	904	892	902	901		897		907 ];
for i=1:length(CntrlList)
    subjectID=CntrlList(i);
    load([destinationPath,'EEGObj_Cntrl',int2str(subjectID),'.mat'],'EEG')
    

    % selecting Components for SIFTING
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




    %estimating connectivity
    EEG= pop_est_mvarConnectivity(EEG,'nogui','connmethods',{'RPDC'}); %{'DTF','dDTF','ffDTF','RPDC','Coh','pCoh'}


RPDC=EEG.CAT.Conn.RPDC;

    save([destinationPath,'CompRPDCMatrx_',int2str(subjectID),'.mat'],'RPDC')
end


