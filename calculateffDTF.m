clear
eeglab;
close all
destinationPath='E:\Abolfazl\OtherProjs\PDPredictiveCoding\Data\restingStatesLastWave(hopefully)\DownlaodedDataset\ProcessedData\';

CntrlList=[908,8010,8060,893,909,900,914,910,891,892,902];

for i=1:length(CntrlList)
    subjectID=CntrlList(i);
    load([destinationPath,'EEGObj_Component',int2str(subjectID),'.mat'],'EEG')
    EEG= pop_est_mvarConnectivity(EEG,'nogui','connmethods',{'ffDTF'}); %{'DTF','dDTF','ffDTF','RPDC','Coh','pCoh'}
    ffDTF=EEG.CAT.Conn.ffDTF;
    save([destinationPath,'CompffDTFMatrx_',int2str(subjectID),'.mat'],'ffDTF')
end



PDListOnDrug=[802	803	806	807	808	813	816	817	819	823	824	827	828	829];

for i=1:length(PDListOnDrug)
    subjectID=PDListOnDrug(i);
    load([destinationPath,'EEGObj_Component',int2str(subjectID),'.mat'],'EEG')
    EEG= pop_est_mvarConnectivity(EEG,'nogui','connmethods',{'ffDTF'}); %{'DTF','dDTF','ffDTF','RPDC','Coh','pCoh'}
    ffDTF=EEG.CAT.Conn.ffDTF;
    save([destinationPath,'CompffDTFMatrx_',int2str(subjectID),'.mat'],'ffDTF')
end



PDListOffDrug=[801	804	805	809	810	811	815	818	820	821	822	825]%	826];

for i=1:length(PDListOffDrug)
    subjectID=PDListOffDrug(i);
    load([destinationPath,'EEGObj_Component',int2str(subjectID),'.mat'],'EEG')
    EEG= pop_est_mvarConnectivity(EEG,'nogui','connmethods',{'ffDTF'}); %{'DTF','dDTF','ffDTF','RPDC','Coh','pCoh'}
    ffDTF=EEG.CAT.Conn.ffDTF;
    save([destinationPath,'CompffDTFMatrx_',int2str(subjectID),'.mat'],'ffDTF')
end
