clear
eeglab;
close all
destinationPath='E:\Abolfazl\OtherProjs\PDPredictiveCoding\Data\restingStatesLastWave(hopefully)\DownlaodedDataset\ProcessedData\';

CntrlList=[894	908	8010	906	903	8060	893	909	911	895	913	900		899	914	910	890	891	912	905	904	892	902	901		897		907];
CntrlListCumPFCComps={};
mm=1;
for i=1:length(CntrlList)
    subjectID=CntrlList(i);
    load([destinationPath,'EEGObj_Cntrl',int2str(subjectID),'.mat'],'EEG')
    EEG = pop_dipfit_settings( EEG, 'hdmfile','C:\\Users\\aalipour\\Documents\\MATLAB\\eeglab2019_1\\plugins\\dipfit\\standard_BEM\\standard_vol.mat','coordformat','MNI','mrifile','C:\\Users\\aalipour\\Documents\\MATLAB\\eeglab2019_1\\plugins\\dipfit\\standard_BEM\\standard_mri.mat','chanfile','C:\\Users\\aalipour\\Documents\\MATLAB\\eeglab2019_1\\plugins\\dipfit\\standard_BEM\\elec\\standard_1005.elc','chansel',[1:EEG.nbchan] );
    EEG = pop_multifit(EEG, [1:size(EEG.icaweights,1)] ,'threshold',100,'dipplot','on','plotopt',{'normlen' 'on'});
    
    PFCComps={};
    kk=1;
    for compNo=1:size(EEG.dipfit.model,2)
    tmp=EEG.dipfit.model(compNo).posxyz;
        if tmp(1)>-74 && tmp(1)<75 && tmp(2)>16 && tmp(2)<75 && tmp(3)>-40 && tmp(3)<75 ;
            PFCComps{kk}=compNo;
            kk=kk+1;
        end

    end
    
    pop_dipplot( EEG, [cell2mat(PFCComps)] ,'rvrange',[0 100] ,'mri','C:\\Users\\aalipour\\Documents\\MATLAB\\eeglab2019_1\\plugins\\dipfit\\standard_BEM\\standard_mri.mat','normlen','on');
    CntrlListCumPFCComps{mm}=cell2mat(PFCComps);
    mm=mm+1;
end

save([destinationPath,'CntrlListCumPFCCompsdDTF.mat'],'CntrlListCumPFCComps');

save([destinationPath,'CntrlListCumPFCCompsdDTF2.mat'],'CntrlListCumPFCComps');


PDListOFFDrug=[801	802	803	804	805	806	807	808	809	810	811	813	815	816	817	818	819	820	821	822	823	824	825	827	828	829];
PDListOFFDrugCumPFCComps={};
mm=1;
for i=1:length(PDListOFFDrug)
    subjectID=PDListOFFDrug(i);
    load([destinationPath,'EEGObj_offDrug',int2str(subjectID),'.mat'],'EEG')
    EEG = pop_dipfit_settings( EEG, 'hdmfile','C:\\Users\\aalipour\\Documents\\MATLAB\\eeglab2019_1\\plugins\\dipfit\\standard_BEM\\standard_vol.mat','coordformat','MNI','mrifile','C:\\Users\\aalipour\\Documents\\MATLAB\\eeglab2019_1\\plugins\\dipfit\\standard_BEM\\standard_mri.mat','chanfile','C:\\Users\\aalipour\\Documents\\MATLAB\\eeglab2019_1\\plugins\\dipfit\\standard_BEM\\elec\\standard_1005.elc','chansel',[1:EEG.nbchan] );
    EEG = pop_multifit(EEG, [1:size(EEG.icaweights,1)] ,'threshold',100,'dipplot','on','plotopt',{'normlen' 'on'});
    
    PFCComps={};
    kk=1;
    for compNo=1:size(EEG.dipfit.model,2)
    tmp=EEG.dipfit.model(compNo).posxyz;
        if tmp(1)>-74 && tmp(1)<75 && tmp(2)>16 && tmp(2)<75 && tmp(3)>-40 && tmp(3)<75 ;
            PFCComps{kk}=compNo;
            kk=kk+1;
        end

    end
    
    pop_dipplot( EEG, [cell2mat(PFCComps)] ,'rvrange',[0 100] ,'mri','C:\\Users\\aalipour\\Documents\\MATLAB\\eeglab2019_1\\plugins\\dipfit\\standard_BEM\\standard_mri.mat','normlen','on');
    PDListOFFDrugCumPFCComps{mm}=cell2mat(PFCComps);
    mm=mm+1;
    
end

save([destinationPath,'PDListOFFDrugCumPFCCompsdDTF_OFFDRUG.mat'],'PDListOFFDrugCumPFCComps');

PDListONDrug=[801	802	803	804	805	806	807	808	809	810	811	813	815	816	817	818	819	820	821	822	823	824	825	827	828	829];
PDListONDrugCumPFCComps={};
mm=1;
for i=1:length(PDListONDrug)
    subjectID=PDListONDrug(i);
    load([destinationPath,'EEGObj_onDrug',int2str(subjectID),'.mat'],'EEG')
    EEG = pop_dipfit_settings( EEG, 'hdmfile','C:\\Users\\aalipour\\Documents\\MATLAB\\eeglab2019_1\\plugins\\dipfit\\standard_BEM\\standard_vol.mat','coordformat','MNI','mrifile','C:\\Users\\aalipour\\Documents\\MATLAB\\eeglab2019_1\\plugins\\dipfit\\standard_BEM\\standard_mri.mat','chanfile','C:\\Users\\aalipour\\Documents\\MATLAB\\eeglab2019_1\\plugins\\dipfit\\standard_BEM\\elec\\standard_1005.elc','chansel',[1:EEG.nbchan] );
    EEG = pop_multifit(EEG, [1:size(EEG.icaweights,1)] ,'threshold',100,'dipplot','on','plotopt',{'normlen' 'on'});
    
    PFCComps={};
    kk=1;
    for compNo=1:size(EEG.dipfit.model,2)
    tmp=EEG.dipfit.model(compNo).posxyz;
        if tmp(1)>-74 && tmp(1)<75 && tmp(2)>16 && tmp(2)<75 && tmp(3)>-40 && tmp(3)<75 ;
            PFCComps{kk}=compNo;
            kk=kk+1;
        end

    end
    
    pop_dipplot( EEG, [cell2mat(PFCComps)] ,'rvrange',[0 100] ,'mri','C:\\Users\\aalipour\\Documents\\MATLAB\\eeglab2019_1\\plugins\\dipfit\\standard_BEM\\standard_mri.mat','normlen','on');
    PDListONDrugCumPFCComps{mm}=cell2mat(PFCComps);
    mm=mm+1;
    
end
save([destinationPath,'PDListONDrugCumPFCCompsdDTF_ONDRUG.mat'],'PDListONDrugCumPFCComps');

% PDListOnDrugCntrl=[894	908	8010	906	903	8060	893	909	911	895	913	900		899	914	910	890	891	912	905	904	892	902	901		897		907];
% PDListOnDrugCntrlCumPFCComps={};
% mm=1;
% for i=1:length(PDListOnDrugCntrl)
%     subjectID=PDListOnDrugCntrl(i);
%     load([destinationPath,'EEGObj_Component',int2str(subjectID),'.mat'],'EEG')
%     EEG = pop_dipfit_settings( EEG, 'hdmfile','C:\\Users\\aalipour\\Documents\\MATLAB\\eeglab2019_1\\plugins\\dipfit\\standard_BEM\\standard_vol.mat','coordformat','MNI','mrifile','C:\\Users\\aalipour\\Documents\\MATLAB\\eeglab2019_1\\plugins\\dipfit\\standard_BEM\\standard_mri.mat','chanfile','C:\\Users\\aalipour\\Documents\\MATLAB\\eeglab2019_1\\plugins\\dipfit\\standard_BEM\\elec\\standard_1005.elc','chansel',[1:EEG.nbchan] );
%     EEG = pop_multifit(EEG, [1:size(EEG.icaweights,1)] ,'threshold',100,'dipplot','on','plotopt',{'normlen' 'on'});
%     
%     PFCComps={};
%     kk=1;
%     for compNo=1:size(EEG.dipfit.model,2)
%     tmp=EEG.dipfit.model(compNo).posxyz;
%         if tmp(1)>-74 && tmp(1)<75 && tmp(2)>16 && tmp(2)<75 && tmp(3)>-40 && tmp(3)<75 ;
%             PFCComps{kk}=compNo;
%             kk=kk+1;
%         end
% 
%     end
%     
%     pop_dipplot( EEG, [cell2mat(PFCComps)] ,'rvrange',[0 100] ,'mri','C:\\Users\\aalipour\\Documents\\MATLAB\\eeglab2019_1\\plugins\\dipfit\\standard_BEM\\standard_mri.mat','normlen','on');
%     PDListOnDrugCntrlCumPFCComps{mm}=cell2mat(PFCComps);
%     mm=mm+1;
% end
% 
PDListOnDrugCntrlCumPFCComps=CntrlListCumPFCComps;
save([destinationPath,'PDListOnDrugCntrlCumPFCCompsdDTF.mat'],'PDListOnDrugCntrlCumPFCComps');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
