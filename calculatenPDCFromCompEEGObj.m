clear
eeglab;
close all
destinationPath='PATH TO EEG OBJECTS';

CntrlList=[908,8010,8060,893,909,900,914,910,891,892,902];

for i=1:length(CntrlList)
    subjectID=CntrlList(i);
    load([destinationPath,'EEGObj_Component',int2str(subjectID),'.mat'],'EEG')
    EEG= pop_est_mvarConnectivity(EEG,'nogui','connmethods',{'nPDC'}); %{'DTF','dDTF','ffDTF','RPDC','Coh','pCoh'}
    nPDC=EEG.CAT.Conn.nPDC;
    save([destinationPath,'CompnPDCMatrx_',int2str(subjectID),'.mat'],'nPDC')
end



PDListOnDrug=[802	803	806	807	808	813	816	817	819	823	824	827	828	829];

for i=1:length(PDListOnDrug)
    subjectID=PDListOnDrug(i);
    load([destinationPath,'EEGObj_Component',int2str(subjectID),'.mat'],'EEG')
    EEG= pop_est_mvarConnectivity(EEG,'nogui','connmethods',{'nPDC'}); %{'DTF','dDTF','ffDTF','RPDC','Coh','pCoh'}
    nPDC=EEG.CAT.Conn.nPDC;
    save([destinationPath,'CompnPDCMatrx_',int2str(subjectID),'.mat'],'nPDC')
end



PDListOffDrug=[801	804	805	809	810	811	815	818	820	821	822	825];%	826];

for i=1:length(PDListOffDrug)
    subjectID=PDListOffDrug(i);
    load([destinationPath,'EEGObj_Component',int2str(subjectID),'.mat'],'EEG')
    EEG= pop_est_mvarConnectivity(EEG,'nogui','connmethods',{'nPDC'}); %{'DTF','dDTF','ffDTF','RPDC','Coh','pCoh'}
    nPDC=EEG.CAT.Conn.nPDC;
    save([destinationPath,'CompnPDCMatrx_',int2str(subjectID),'.mat'],'nPDC')
end


PDListOnDrugCntrl=[894	906 903	911 895 913 899 890 912 905 904 901];%	826];

for i=1:length(PDListOnDrugCntrl)
    subjectID=PDListOnDrugCntrl(i);
    load([destinationPath,'EEGObj_Component',int2str(subjectID),'.mat'],'EEG')
    EEG= pop_est_mvarConnectivity(EEG,'nogui','connmethods',{'nPDC'}); %{'DTF','dDTF','ffDTF','RPDC','Coh','pCoh'}
    nPDC=EEG.CAT.Conn.nPDC;
    save([destinationPath,'CompnPDCMatrx_',int2str(subjectID),'.mat'],'nPDC')
end

