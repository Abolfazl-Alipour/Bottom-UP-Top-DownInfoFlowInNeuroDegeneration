
subjectID=894;
destinationPath='PATH TO EEG OBJECTS\ProcessedData\';
sourceFilePath='PATH TO RAW DATASET FILES';
locPath=('PATH TO LOCATION FILE\MATLAB\eeglab2019_1\plugins\dipfit\standard_BESA\standard-10-5-cap385.elp');
locPathFrwrdMdel='PATH TO HEAD MODEL\MATLAB\eeglab2019_1\plugins\dsi\headModel\resources\head_modelColin27_8003_Standard-10-5-Cap339.mat';
chckModelConsistancy=0;


CntrlList=[894,908,8010,906,903,8060,893,909,911,895,913,900,896,899,914,910,890,891,912,905,904,892,902,901,898,897,907];
%case num 8070, is not used due to noisyness
for i=1:length(CntrlList)
    subjectID=CntrlList(i);
    connFromSubjID(sourceFilePath,subjectID,destinationPath,locPath,locPathFrwrdMdel,chckModelConsistancy)
end

PDList=[801	802	803	804	805	806	807	808	809	810	811	813	814	815	816	817	818	819	820	821	822	823	824	825	826	827	828	829];
%subject 826 is too noisy
PDList=[827	828	829];
for i=1:length(PDList)
    subjectID=PDList(i);
    
    connFromSubjID(sourceFilePath,subjectID,destinationPath,locPath,locPathFrwrdMdel,chckModelConsistancy)
end
