
destinationPath='E:\Abolfazl\OtherProjs\PDPredictiveCoding\Data\restingStatesLastWave(hopefully)\DownlaodedDataset\ProcessedData\';
sourceFilePath='E:\Abolfazl\OtherProjs\PDPredictiveCoding\Data\restingStatesLastWave(hopefully)\DownlaodedDataset\PD REST\';
locPath=('C:\Users\aalipour\Documents\MATLAB\eeglab2019_1\plugins\dipfit\standard_BESA\standard-10-5-cap385.elp');
locPathFrwrdMdel='C:\Users\aalipour\Documents\MATLAB\eeglab2019_1\plugins\dsi\headModel\resources\head_modelColin27_8003_Standard-10-5-Cap339.mat';
chckModelConsistancy=0;


PDList=[801	802	803	804	805	806	807	808	809	810	811	813	815	816	817	818	819	820	821	822	823	824	825	827	828	829];
drug='off';
% PDRPDCMat=nan(length(PDList),upperFrequencyBound);
for i=1:length(PDList)
    subjectID=PDList(i);
    compConnFromSubjIDEYESOPEN(sourceFilePath,subjectID,destinationPath,locPath,locPathFrwrdMdel,chckModelConsistancy,drug)
end

drug='onn';
% PDRPDCMat=nan(length(PDList),upperFrequencyBound);
for i=1:length(PDList)
    subjectID=PDList(i);
    compConnFromSubjIDEYESOPEN(sourceFilePath,subjectID,destinationPath,locPath,locPathFrwrdMdel,chckModelConsistancy,drug)
end


CntrlList=[894	908	8010	906	903	8060	893	909	911	895	913	900		899	914	910	890	891	912	905	904	892	902	901		897		907];
% 8070 only 5 seconds, excluded

drug='Cnt';
% PDRPDCMat=nan(length(PDList),upperFrequencyBound);
for i=1:length(CntrlList)
    subjectID=CntrlList(i);
    compConnFromSubjIDEYESOPEN(sourceFilePath,subjectID,destinationPath,locPath,locPathFrwrdMdel,chckModelConsistancy,drug)
end
