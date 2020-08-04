
eeglab;
close all
destinationPath='E:\Abolfazl\OtherProjs\PDPredictiveCoding\Data\restingStatesLastWave(hopefully)\DownlaodedDataset\ProcessedData\';

CntrlList=[894	908	8010	906	903	8060	893	909	911	895	913	900		899	914	910	890	891	912	905	904	892	902	901		897		907];

for i=1:length(CntrlList)
    subjectID=CntrlList(i);
    load([destinationPath,'EEGObj_Component',int2str(subjectID),'.mat'],'EEG')
    EEG = pop_dipfit_settings( EEG, 'hdmfile','C:\\Users\\aalipour\\Documents\\MATLAB\\eeglab2019_1\\plugins\\dipfit\\standard_BESA\\standard_BESA.mat','coordformat','Spherical','mrifile','C:\\Users\\aalipour\\Documents\\MATLAB\\eeglab2019_1\\plugins\\dipfit\\standard_BESA\\avg152t1.mat','chanfile','C:\\Users\\aalipour\\Documents\\MATLAB\\eeglab2019_1\\plugins\\dipfit\\standard_BESA\\standard-10-5-cap385.elp','chansel',[1:EEG.nbchan] );
    EEG=pop_resample(EEG, 100, 0.8, 0.4);
    EEG=pop_saveset( EEG, 'filename',['ComponentSet',int2str(subjectID),'.set'],'filepath','E:\\Abolfazl\\OtherProjs\\PDPredictiveCoding\\Data\\restingStatesLastWave(hopefully)\\\DownlaodedDataset\\eeglabDsetfilesCntrl\\');
end



PDListOnDrug=[801	802	803	804	805	806	807	808	809	810	811	813	815	816	817	818	819	820	821	822	823	824	825	827	828	829];

for i=1:length(PDListOnDrug)
    subjectID=PDListOnDrug(i);
    load([destinationPath,'EEGObj_Component',int2str(subjectID),'.mat'],'EEG')
    load([destinationPath,'EEGObj_Component',int2str(subjectID),'.mat'],'EEG')
    EEG = pop_dipfit_settings( EEG, 'hdmfile','C:\\Users\\aalipour\\Documents\\MATLAB\\eeglab2019_1\\plugins\\dipfit\\standard_BESA\\standard_BESA.mat','coordformat','Spherical','mrifile','C:\\Users\\aalipour\\Documents\\MATLAB\\eeglab2019_1\\plugins\\dipfit\\standard_BESA\\avg152t1.mat','chanfile','C:\\Users\\aalipour\\Documents\\MATLAB\\eeglab2019_1\\plugins\\dipfit\\standard_BESA\\standard-10-5-cap385.elp','chansel',[1:EEG.nbchan] );
    EEG = pop_multifit(EEG, [1:size(EEG.icaact,1)] ,'threshold',100,'rmout','off','dipplot','off','plotopt',{'normlen' 'on'});
    EEG=pop_resample(EEG, 100, 0.8, 0.4);
    EEG=pop_saveset( EEG, 'filename',['ComponentSet',int2str(subjectID),'.set'],'filepath','E:\\Abolfazl\\OtherProjs\\PDPredictiveCoding\\Data\\restingStatesLastWave(hopefully)\\\DownlaodedDataset\\eeglabDsetfilesPD\\');
end



PDListOffDrug=[801	802	803	804	805	806	807	808	809	810	811	813	815	816	817	818	819	820	821	822	823	824	825	827	828	829];

for i=1:length(PDListOffDrug)
    subjectID=PDListOffDrug(i);
    load([destinationPath,'EEGObj_Component',int2str(subjectID),'.mat'],'EEG')
    load([destinationPath,'EEGObj_Component',int2str(subjectID),'.mat'],'EEG')
    EEG = pop_dipfit_settings( EEG, 'hdmfile','C:\\Users\\aalipour\\Documents\\MATLAB\\eeglab2019_1\\plugins\\dipfit\\standard_BESA\\standard_BESA.mat','coordformat','Spherical','mrifile','C:\\Users\\aalipour\\Documents\\MATLAB\\eeglab2019_1\\plugins\\dipfit\\standard_BESA\\avg152t1.mat','chanfile','C:\\Users\\aalipour\\Documents\\MATLAB\\eeglab2019_1\\plugins\\dipfit\\standard_BESA\\standard-10-5-cap385.elp','chansel',[1:EEG.nbchan] );
    EEG = pop_multifit(EEG, [1:size(EEG.icaact,1)] ,'threshold',100,'rmout','off','dipplot','off','plotopt',{'normlen' 'on'});
    EEG=pop_resample(EEG, 100, 0.8, 0.4);
    EEG=pop_saveset( EEG, 'filename',['ComponentSet',int2str(subjectID),'.set'],'filepath','E:\\Abolfazl\\OtherProjs\\PDPredictiveCoding\\Data\\restingStatesLastWave(hopefully)\\\DownlaodedDataset\\eeglabDsetfilesPDOff\\');
end

eeglab;
close all
destinationPath='E:\Abolfazl\OtherProjs\PDPredictiveCoding\Data\restingStatesLastWave(hopefully)\DownlaodedDataset\ProcessedData\';

CntrlList=[894	908	8010	906	903	8060	893	909	911	895	913	900		899	914	910	890	891	912	905	904	892	902	901		897		907];

for i=1:length(CntrlList)
    subjectID=CntrlList(i);
    
    EEG = pop_loadset('filename',['ComponentSet',int2str(subjectID),'.set'],'filepath','E:\\Abolfazl\\OtherProjs\\PDPredictiveCoding\\Data\\restingStatesLastWave(hopefully)\\DownlaodedDataset\\eeglabDsetfilesCntrl\\');
  
    EEG=pop_resample(EEG, 100, 0.8, 0.4);
    EEG=pop_saveset( EEG, 'filename',['ComponentSet',int2str(subjectID),'.set'],'filepath','E:\\Abolfazl\\OtherProjs\\PDPredictiveCoding\\Data\\restingStatesLastWave(hopefully)\\\DownlaodedDataset\\eeglabDsetfilesCntrl\\');
end



PDListOnDrug=[802	803	806	807	808	813	816	817	819	823	824	827	828	829];

for i=1:length(PDListOnDrug)
    subjectID=PDListOnDrug(i);
   EEG = pop_loadset('filename',['ComponentSet',int2str(subjectID),'.set'],'filepath','E:\\Abolfazl\\OtherProjs\\PDPredictiveCoding\\Data\\restingStatesLastWave(hopefully)\\DownlaodedDataset\\eeglabDsetfilesPD\\');
  
    EEG=pop_resample(EEG, 100, 0.8, 0.4);
     EEG=pop_saveset( EEG, 'filename',['ComponentSet',int2str(subjectID),'.set'],'filepath','E:\\Abolfazl\\OtherProjs\\PDPredictiveCoding\\Data\\restingStatesLastWave(hopefully)\\\DownlaodedDataset\\eeglabDsetfilesPD\\');
end



PDListOffDrug=[801	804	805	809	810	811	815	818	820	821	822	825	];

for i=1:length(PDListOffDrug)
    subjectID=PDListOffDrug(i);
    EEG = pop_loadset('filename',['ComponentSet',int2str(subjectID),'.set'],'filepath','E:\\Abolfazl\\OtherProjs\\PDPredictiveCoding\\Data\\restingStatesLastWave(hopefully)\\DownlaodedDataset\\eeglabDsetfilesPDOff\\');
  
    EEG=pop_resample(EEG, 100, 0.8, 0.4);
    EEG=pop_saveset( EEG, 'filename',['ComponentSet',int2str(subjectID),'.set'],'filepath','E:\\Abolfazl\\OtherProjs\\PDPredictiveCoding\\Data\\restingStatesLastWave(hopefully)\\\DownlaodedDataset\\eeglabDsetfilesPDOff\\');
end
