clear
eeglab;
close all
destinationPath='E:\Abolfazl\OtherProjs\PDPredictiveCoding\Data\restingStatesLastWave(hopefully)\DownlaodedDataset\ProcessedData\';

CntrlList=[908,8010,8060,893,909,900,914,910,891,892,902];
CntrlListCumPFCComps={};
mm=1;
for i=1:length(CntrlList)
    subjectID=CntrlList(2);
    load([destinationPath,'EEGObj_Component',int2str(subjectID),'.mat'],'EEG')
    EEG = pop_dipfit_settings( EEG, 'hdmfile','C:\\Users\\aalipour\\Documents\\MATLAB\\eeglab2019_1\\plugins\\dipfit\\standard_BEM\\standard_vol.mat','coordformat','MNI','mrifile','C:\\Users\\aalipour\\Documents\\MATLAB\\eeglab2019_1\\plugins\\dipfit\\standard_BEM\\standard_mri.mat','chanfile','C:\\Users\\aalipour\\Documents\\MATLAB\\eeglab2019_1\\plugins\\dipfit\\standard_BEM\\elec\\standard_1005.elc','chansel',[1:EEG.nbchan] );
    EEG = pop_multifit(EEG, [1:size(EEG.icaweights,1)] ,'threshold',100,'dipplot','on','plotopt',{'normlen' 'on'});
    
    PFCComps={};
    kk=1;
    for compNo=1:size(EEG.dipfit.model,2)
    tmp=EEG.dipfit.model(compNo).posxyz;
        if tmp(1)>-75 && tmp(1)<75 && tmp(2)>16 && tmp(2)<75 && tmp(3)>-40 && tmp(3)<75 ;
            PFCComps{kk}=compNo;
            kk=kk+1;
        end

    end
    pop_dipplot( EEG, [1:size(EEG.icaweights,1)] ,'rvrange',[0 100] ,'mri','C:\\Users\\aalipour\\Documents\\MATLAB\\eeglab2019_1\\plugins\\dipfit\\standard_BEM\\standard_mri.mat','normlen','on');%,'projlines','on'
   
    pop_dipplot( EEG, [cell2mat(PFCComps)] ,'rvrange',[0 100] ,'mri','C:\\Users\\aalipour\\Documents\\MATLAB\\eeglab2019_1\\plugins\\dipfit\\standard_BEM\\standard_mri.mat','normlen','on');
    CntrlListCumPFCComps{mm}=cell2mat(PFCComps);
    mm=mm+1;



end