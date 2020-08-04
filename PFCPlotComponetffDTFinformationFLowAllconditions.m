
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  now plots

clear
close all
destinationPath='E:\Abolfazl\OtherProjs\PDPredictiveCoding\Data\restingStatesLastWave(hopefully)\DownlaodedDataset\ProcessedData\';

load([destinationPath,'CntrlListCumPFCCompsdDTF.mat'],'CntrlListCumPFCComps');

upperFrequencyBound=58;
destinationPath='E:\Abolfazl\OtherProjs\PDPredictiveCoding\Data\restingStatesLastWave(hopefully)\DownlaodedDataset\ProcessedData\';

CntrlList=[894	908	8010	906	903	8060	893	909	911	895	913	900		899	914	910	890	891	912	905	904	892	902	901		897		907];

CntrlffDTFMat=nan(length(CntrlList),upperFrequencyBound);
for i=1:length(CntrlList)
    subjectID=CntrlList(i);
    RPDCFileName=strcat('CompffDTFMatrx_',int2str(subjectID),'.mat');
    load(strcat(destinationPath,RPDCFileName))
    PFCComps=CntrlListCumPFCComps{i};
    % Removes components that are not located in the PFC
    ffDTF=ffDTF(:,PFCComps,:,:);
    tmpTime=mean(ffDTF(:,:,1:58,:),4);
    infOutflow=nan(1,size(tmpTime,3));
    
    for layerNum=1:size(tmpTime,3)
        tmp=tmpTime(:,:,layerNum);
        tmp(1:(size(tmp,1)+1):end)=0; %zeroing out the diagonal values
        infOutflow(1,layerNum)=mean(mean(tmp));
    end
    CntrlffDTFMat(i,:)=infOutflow;
end


load([destinationPath,'PDListOFFDrugCumPFCCompsdDTF_OFFDRUG.mat'],'PDListOFFDrugCumPFCComps');

% PDList=[801	802	803	804	805	806	807	808	809	810	811	813	814	815	816	817	818	819	820	821	822	823	824	825	];%826	827	828	829];
PDOFFList=[801	802	803	804	805	806	807	808	809	810	811	813	815	816	817	818	819	820	821	822	823	824	825	827	828	829];
PDOFFListLED=[1275 600 520 550 1150 600 400 640 600 100 1175 1796 300 338];
PDOFFffDTFMat=nan(length(PDOFFList),upperFrequencyBound);
for i=1:length(PDOFFList)
    subjectID=PDOFFList(i);
    RPDCFileName=strcat('CompffDTFMatrx_',int2str(subjectID),'OFFDRUG','.mat');
    load(strcat(destinationPath,RPDCFileName))
    PFCComps=PDListOFFDrugCumPFCComps{i};
    
    
    ffDTFNONPFC=ffDTF;
    ffDTFNONPFC(:,PFCComps,:,:)=[];
    
    tmpTime2=mean(ffDTFNONPFC(:,:,1:58,:),4);
    infOutflow=nan(1,size(tmpTime2,3));
    
    for layerNum=1:size(tmpTime2,3)
        tmp2=tmpTime2(:,:,layerNum);
        tmp2(1:(size(tmp2,1)+1):end)=0; %zeroing out the diagonal values
        infOutflow2(1,layerNum)=mean(mean(tmp2));
    end
    PDOFFffDTFMatNONPFC(i,:)=infOutflow2;
    
    
    
    
    % Removes components that are not located in the PFC
    
    ffDTF=ffDTF(:,PFCComps,:,:);
    tmpTime=mean(ffDTF(:,:,1:58,:),4);
    infOutflow=nan(1,size(tmpTime,3));
    
    for layerNum=1:size(tmpTime,3)
        tmp=tmpTime(:,:,layerNum);
        tmp(1:(size(tmp,1)+1):end)=0; %zeroing out the diagonal values
        infOutflow(1,layerNum)=mean(mean(tmp));
    end
    PDOFFffDTFMat(i,:)=infOutflow;
end


figure; errorbar(mean(CntrlffDTFMat,1),(std(CntrlffDTFMat,1)./sqrt(size(CntrlffDTFMat,1))));
hold on; errorbar(mean(PDOFFffDTFMat,1),(std(PDOFFffDTFMat,1)./sqrt(size(PDOFFffDTFMat,1))))


figure; stdshade(CntrlffDTFMat,0.3,'b');
hold on; stdshade(PDOFFffDTFMat,0.3,'r');
title('Infromation outflow from PFC across different frequencies')
xlabel('Frequency') 
ylabel('Mean ff-DTF') 
pLabel="PD Patients Off Drug (PFC Components)";
h= findobj(gca,'Type','Line');
% legend([h(1),h(2)],{'Control',pLabel},'Location','northeast')
legend('boxoff')
ax = gca;
ax.FontSize = 32;
ax.YLim=[0 0.05];
set(get(get(h(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(h(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
legend({'Control (PFC Components)',pLabel},'Location','northeast')


load([destinationPath,'PDListONDrugCumPFCCompsdDTF_ONDRUG.mat'],'PDListONDrugCumPFCComps');

PDListONDrug=[801	802	803	804	805	806	807	808	809	810	811	813	815	816	817	818	819	820	821	822	823	824	825	827	828	829];
PDListONDrugLED=[697 600 250 400 300 200 500 600 700 60  600];
PDONffDTF4Mat=nan(length(PDListONDrug),upperFrequencyBound);
for i=1:length(PDListONDrug)
    subjectID=PDListONDrug(i);
    RPDCFileName=strcat('CompffDTFMatrx_',int2str(subjectID),'ONDRUG','.mat');
    load(strcat(destinationPath,RPDCFileName))
    PFCComps=PDListONDrugCumPFCComps{i};
    
    
    
    ffDTFNONPFC=ffDTF;
    ffDTFNONPFC(:,PFCComps,:,:)=[];
    
    tmpTime2=mean(ffDTFNONPFC(:,:,1:58,:),4);
    infOutflow=nan(1,size(tmpTime2,3));
    
    for layerNum=1:size(tmpTime2,3)
        tmp2=tmpTime2(:,:,layerNum);
        tmp2(1:(size(tmp2,1)+1):end)=0; %zeroing out the diagonal values
        infOutflow2(1,layerNum)=mean(mean(tmp2));
    end
    PDONffDTFMatNONPFC(i,:)=infOutflow2;
    
    
    
    % Removes components that are not located in the PFC

    ffDTF=ffDTF(:,PFCComps,:,:);
    
    tmpTime=mean(ffDTF(:,:,1:58,:),4);
    infOutflow=nan(1,size(tmpTime,3));
    
    for layerNum=1:size(tmpTime,3)
        tmp=tmpTime(:,:,layerNum);
        tmp(1:(size(tmp,1)+1):end)=0; %zeroing out the diagonal values
        infOutflow(1,layerNum)=mean(mean(tmp));
    end
    PDONffDTF4Mat(i,:)=infOutflow;
end

figure; errorbar(mean(PDONffDTFMatNONPFC,1),(std(PDONffDTFMatNONPFC,1)./sqrt(size(PDONffDTFMatNONPFC,1))));
hold on; errorbar(mean(PDONffDTF4Mat,1),(std(PDONffDTF4Mat,1)./sqrt(size(PDONffDTF4Mat,1))))


figure; stdshade(CntrlffDTFMat,0.3,'b');
hold on; stdshade(PDONffDTF4Mat,0.3,'green');
title('Infromation outflow from PFC across different frequencies')
xlabel('Frequency') 
ylabel('Mean ff-DTF') 
pLabel="PD Patients On Drug (PFC Components)";
h= findobj(gca,'Type','Line');
% legend([h(1),h(2)],{'Control',pLabel},'Location','northeast')
legend('boxoff')
ax = gca;
ax.FontSize = 32;
ax.YLim=[0 0.05];
set(get(get(h(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(h(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
legend({'Control (PFC Components)',pLabel},'Location','northeast')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  GAMMA
lowerBound=52;
higherBound=54;
PDOFFAnovaData=mean(PDOFFffDTFMat(:,lowerBound:higherBound),2)';
PDOFFAnovaLabels={};
for hh=1:length(PDOFFAnovaData)
    PDOFFAnovaLabels{hh}= 'PDOFF';
end

cntrlAnovaData=mean(CntrlffDTFMat(:,lowerBound:higherBound),2)';
cntrlAnovaLabels={};
for hh=1:length(cntrlAnovaData)
    cntrlAnovaLabels{hh}= 'Cnrl';
end

PDONAnovaData=mean(PDONffDTF4Mat(:,lowerBound:higherBound),2)';
PDONAnovaLabels={};
for hh=1:length(PDONAnovaData)
    PDONAnovaLabels{hh}= 'PDON';
end

[p,tbl,stats]=anova1([PDOFFAnovaData,cntrlAnovaData,PDONAnovaData],[PDOFFAnovaLabels,cntrlAnovaLabels,PDONAnovaLabels]);
[COMPARISON,MEANS,H,GNAMES] = multcompare(stats);
[GNAMES(COMPARISON(:,1)), GNAMES(COMPARISON(:,2)), num2cell(COMPARISON(:,3:6))];

FEXTRACT=[tbl{2,3},tbl{3,3},tbl{2,5},p];
TUKEYEXTRACT=[MEANS COMPARISON(:,6)];
% 
% [h,p,ci,stats]=ttest2(mean(CntrlffDTFMat(:,lowerBound:higherBound),2),mean(PDOFFffDTFMat(:,lowerBound:higherBound),2))

% computeCohen_d(mean(CntrlffDTFMat(:,lowerBound:higherBound),2),mean(PDOFFffDTFMat(:,lowerBound:higherBound),2))

x = 1:3;
Cntrl=mean(CntrlffDTFMat(:,lowerBound:higherBound),2);
PDOff=mean(PDOFFffDTFMat(:,lowerBound:higherBound),2);
PDOn=mean(PDONffDTF4Mat(:,lowerBound:higherBound),2);
data=[mean(PDOff),mean(Cntrl),mean(PDOn)];
SEMs=[sem(PDOff),sem(Cntrl),sem(PDOn)]
figure; 
b = bar(x,data,'LineWidth',2);
b.FaceColor = 'flat';
b.CData(1,:) = [1 0.5 .5];
b.CData(2,:) = [0.5 0.5 1];
b.CData(3,:) = [.5 1 .5];
%,(std(PDONdDTFCntrlMat,1)./sqrt(size(PDONdDTFCntrlMat,1))));
hold on
er = errorbar(x,data,[sem(PDOff),sem(Cntrl),sem(PDOn)],'CapSize',18,'LineWidth', 2);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
 
title('Infromation flow in gamma band')

ylabel('Average ff-DTF (52-54 Hz)') 
ax=gca;
ax.YLim=[0 0.02];
ax.XTickLabel={('Cntrl (PFC)'),('PD Off Drug (PFC)'),('PD On Drug (PFC)')};
ax.FontSize = 28;

hold off


%%%%%%%%%%%%%%%%%%%%%  LOW and MID BETA

lowerBound=13;
higherBound=20;
PDOFFAnovaData=mean(PDOFFffDTFMat(:,lowerBound:higherBound),2)';
PDOFFAnovaLabels={};
for hh=1:length(PDOFFAnovaData)
    PDOFFAnovaLabels{hh}= 'PDOFF';
end

cntrlAnovaData=mean(CntrlffDTFMat(:,lowerBound:higherBound),2)';
cntrlAnovaLabels={};
for hh=1:length(cntrlAnovaData)
    cntrlAnovaLabels{hh}= 'Cnrl';
end

PDONAnovaData=mean(PDONffDTF4Mat(:,lowerBound:higherBound),2)';
PDONAnovaLabels={};
for hh=1:length(PDONAnovaData)
    PDONAnovaLabels{hh}= 'PDON';
end

[p,tbl,stats]=anova1([PDOFFAnovaData,cntrlAnovaData,PDONAnovaData],[PDOFFAnovaLabels,cntrlAnovaLabels,PDONAnovaLabels]);
[COMPARISON,MEANS,H,GNAMES] = multcompare(stats);
[GNAMES(COMPARISON(:,1)), GNAMES(COMPARISON(:,2)), num2cell(COMPARISON(:,3:6))];

FEXTRACT=[tbl{2,3},tbl{3,3},tbl{2,5},p];
TUKEYEXTRACT=[MEANS COMPARISON(:,6)];
% 
% [h,p,ci,stats]=ttest2(mean(CntrlffDTFMat(:,lowerBound:higherBound),2),mean(PDOFFffDTFMat(:,lowerBound:higherBound),2))

% computeCohen_d(mean(CntrlffDTFMat(:,lowerBound:higherBound),2),mean(PDOFFffDTFMat(:,lowerBound:higherBound),2))

x = 1:3;
Cntrl=mean(CntrlffDTFMat(:,lowerBound:higherBound),2);
PDOff=mean(PDOFFffDTFMat(:,lowerBound:higherBound),2);
PDOn=mean(PDONffDTF4Mat(:,lowerBound:higherBound),2);
data=[mean(PDOff),mean(Cntrl),mean(PDOn)];
SEMs=[sem(PDOff),sem(Cntrl),sem(PDOn)]
figure; 
b = bar(x,data,'LineWidth',2);
b.FaceColor = 'flat';
b.CData(1,:) = [1 0.5 .5];
b.CData(2,:) = [0.5 0.5 1];
b.CData(3,:) = [.5 1 .5];
%,(std(PDONdDTFCntrlMat,1)./sqrt(size(PDONdDTFCntrlMat,1))));
hold on
er = errorbar(x,data,[sem(PDOff),sem(Cntrl),sem(PDOn)],'CapSize',18,'LineWidth', 2);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
 
title('Infromation flow in low and mid beta band')

ylabel('Average ff-DTF (13-20 Hz)') 
ax=gca;
ax.YLim=[0 0.02];
ax.XTickLabel={('Cntrl (PFC)'),('PD Off Drug (PFC)'),('PD On Drug (PFC)')};
ax.FontSize = 28;

hold off  
%%%%%%%%%%%%%%%%% now alpha


lowerBound=8;
higherBound=12;
[h,p,ci,stats]=ttest2(mean(CntrlffDTFMat(:,lowerBound:higherBound),2),mean(PDOFFffDTFMat(:,lowerBound:higherBound),2))

computeCohen_d(mean(CntrlffDTFMat(:,lowerBound:higherBound),2),mean(PDOFFffDTFMat(:,lowerBound:higherBound),2))

x = 1:2;
Cntrl=mean(CntrlffDTFMat(:,lowerBound:higherBound),2);
PD=mean(PDOFFffDTFMat(:,lowerBound:higherBound),2);
data=[mean(Cntrl),mean(PD)];
figure; 
b = bar(x,data,'LineWidth',2);
b.FaceColor = 'flat';
b.CData(1,:) = [0.5 0.5 1];
b.CData(2,:) = [1 0.5 .5];
%,(std(PDONdDTFCntrlMat,1)./sqrt(size(PDONdDTFCntrlMat,1))));
hold on
er = errorbar(x,data,[sem(Cntrl),sem(PD)],[sem(Cntrl),sem(PD)],'CapSize',18,'LineWidth', 2);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
title('Infromation outflow from PFC in alpha band')

ylabel('Average sum ff-DTF (8-12 Hz)') 
ax=gca;
ax.YLim=[0 17];
ax.XTickLabel={('Cntrl'),('PD Patients Off Drug')};

ax.FontSize = 32;

hold off


%%%%%%%%%%%%%%% now theta


lowerBound=4;
higherBound=7;
[h,p,ci,stats]=ttest2(mean(CntrlffDTFMat(:,lowerBound:higherBound),2),mean(PDOFFffDTFMat(:,lowerBound:higherBound),2))

computeCohen_d(mean(CntrlffDTFMat(:,lowerBound:higherBound),2),mean(PDOFFffDTFMat(:,lowerBound:higherBound),2))

x = 1:2;
Cntrl=mean(CntrlffDTFMat(:,lowerBound:higherBound),2);
PD=mean(PDOFFffDTFMat(:,lowerBound:higherBound),2);
data=[mean(Cntrl),mean(PD)];
figure; 
b = bar(x,data,'LineWidth',2);
b.FaceColor = 'flat';
b.CData(1,:) = [0.5 0.5 1];
b.CData(2,:) = [1 0.5 .5];
%,(std(PDONdDTFCntrlMat,1)./sqrt(size(PDONdDTFCntrlMat,1))));
hold on
er = errorbar(x,data,[sem(Cntrl),sem(PD)],[sem(Cntrl),sem(PD)],'CapSize',18,'LineWidth', 2);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
title('Infromation outflow from PFC in theta band')

ylabel('Average sum ff-DTF (4-7 Hz)') 
ax=gca;
ax.YLim=[0 20];
ax.XTickLabel={('Cntrl'),('PD Patients Off Drug')};

ax.FontSize = 32;

hold off


%%%%%%%%%%%%%%
% should ADD 814!!!!!!!!!!?? it is bad just like 826?

% 
% PDListOnDrugCntrl=[894	906 903	911 895 913 899 890 912 905 904 901];%	826];
% PDONCntrlffDTF4Mat=nan(length(PDListOnDrugCntrl),upperFrequencyBound);
% for i=1:length(PDListOnDrugCntrl)
%     subjectID=PDListOnDrugCntrl(i);
%     RPDCFileName=strcat('CompffDTFMatrx_',int2str(subjectID),'.mat');
%     load(strcat(destinationPath,RPDCFileName))
%     PFCComps=PDListOnDrugCntrlCumPFCComps{i};
%     % Removes components that are not located in the PFC
%     ffDTF=ffDTF(:,PFCComps,:,:);
%     
%     tmpTime=mean(ffDTF(:,:,1:58,:),4);
%     infOutflow=nan(1,size(tmpTime,3));
%     
%     for layerNum=1:size(tmpTime,3)
%         tmp=tmpTime(:,:,layerNum);
%         tmp(1:(size(tmp,1)+1):end)=0; %zeroing out the diagonal values
%         infOutflow(1,layerNum)=sum(sum(tmp));
%     end
%     PDONCntrlffDTF4Mat(i,:)=infOutflow;
% end


load([destinationPath,'PDListONDrugCumPFCCompsdDTF_ONDRUG.mat'],'PDListONDrugCumPFCComps');

PDListONDrug=[801	802	803	804	805	806	807	808	809	810	811	813	815	816	817	818	819	820	821	822	823	824	825	827	828	829];
PDListONDrugLED=[697 600 250 400 300 200 500 600 700 60  600];
PDONffDTF4Mat=nan(length(PDListONDrug),upperFrequencyBound);
for i=1:length(PDListONDrug)
    subjectID=PDListONDrug(i);
    RPDCFileName=strcat('CompffDTFMatrx_',int2str(subjectID),'ONDRUG','.mat');
    load(strcat(destinationPath,RPDCFileName))
    PFCComps=PDListONDrugCumPFCComps{i};
    
    
    
    ffDTFNONPFC=ffDTF;
    ffDTFNONPFC(:,PFCComps,:,:)=[];
    
    tmpTime2=mean(ffDTFNONPFC(:,:,1:58,:),4);
    infOutflow=nan(1,size(tmpTime2,3));
    
    for layerNum=1:size(tmpTime2,3)
        tmp2=tmpTime2(:,:,layerNum);
        tmp2(1:(size(tmp2,1)+1):end)=0; %zeroing out the diagonal values
        infOutflow2(1,layerNum)=mean(mean(tmp2));
    end
    PDONffDTFMatNONPFC(i,:)=infOutflow2;
    
    
    
    % Removes components that are not located in the PFC

    ffDTF=ffDTF(:,PFCComps,:,:);
    
    tmpTime=mean(ffDTF(:,:,1:58,:),4);
    infOutflow=nan(1,size(tmpTime,3));
    
    for layerNum=1:size(tmpTime,3)
        tmp=tmpTime(:,:,layerNum);
        tmp(1:(size(tmp,1)+1):end)=0; %zeroing out the diagonal values
        infOutflow(1,layerNum)=mean(mean(tmp));
    end
    PDONffDTF4Mat(i,:)=infOutflow;
end

figure; errorbar(mean(PDONffDTFMatNONPFC,1),(std(PDONffDTFMatNONPFC,1)./sqrt(size(PDONffDTFMatNONPFC,1))));
hold on; errorbar(mean(PDONffDTF4Mat,1),(std(PDONffDTF4Mat,1)./sqrt(size(PDONffDTF4Mat,1))))


figure; stdshade(CntrlffDTFMat,0.3,'b');
hold on; stdshade(PDONffDTF4Mat,0.3,'green');
title('Infromation outflow from PFC across different frequencies')
xlabel('Frequency') 
ylabel('Mean ff-DTF') 
pLabel="PD Patients On Drug (PFC Components)";
h= findobj(gca,'Type','Line');
% legend([h(1),h(2)],{'Control',pLabel},'Location','northeast')
legend('boxoff')
ax = gca;
ax.FontSize = 32;
ax.YLim=[0 0.05];
set(get(get(h(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(h(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
legend({'Control (PFC Components)',pLabel},'Location','northeast')



%%%%%%%%%%%%%%%%%%%%  GAMMA PEAK

lowerBound=52;
higherBound=54;
[h,p,ci,stats]=ttest2(mean(CntrlffDTFMat(:,lowerBound:higherBound),2),mean(PDONffDTF4Mat(:,lowerBound:higherBound),2))

computeCohen_d(mean(CntrlffDTFMat(:,lowerBound:higherBound),2),mean(PDONffDTF4Mat(:,lowerBound:higherBound),2))

x = 1:2;
Cntrl=mean(CntrlffDTFMat(:,lowerBound:higherBound),2);
PD=mean(PDONffDTF4Mat(:,lowerBound:higherBound),2);
data=[mean(Cntrl),mean(PD)];
SEMs=[sem(Cntrl),sem(PD)]
figure; 
b = bar(x,data,'LineWidth',2);
b.FaceColor = 'flat';
b.CData(1,:) = [0.5 0.5 1];
b.CData(2,:) =  [.5 1 .5];
%,(std(PDONdDTFCntrlMat,1)./sqrt(size(PDONdDTFCntrlMat,1))));
hold on
er = errorbar(x,data,[sem(Cntrl),sem(PD)],[sem(Cntrl),sem(PD)],'CapSize',18,'LineWidth', 2);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
title('Infromation outflow from PFC in gamma band')

ylabel('Average ff-DTF (52-54 Hz)') 
ax=gca;
ax.YLim=[0 0.02];
ax.XTickLabel={('Control (PFC)'),('PD Patients On Drug (PFC)')};

ax.FontSize = 32;

hold off


%%%%%%%%%%%%%%%%%%  LOW and MID beta

lowerBound=13;
higherBound=20;
[h,p,ci,stats]=ttest2(mean(CntrlffDTFMat(:,lowerBound:higherBound),2),mean(PDONffDTF4Mat(:,lowerBound:higherBound),2))

computeCohen_d(mean(CntrlffDTFMat(:,lowerBound:higherBound),2),mean(PDONffDTF4Mat(:,lowerBound:higherBound),2))

x = 1:2;
Cntrl=mean(CntrlffDTFMat(:,lowerBound:higherBound),2);
PD=mean(PDONffDTF4Mat(:,lowerBound:higherBound),2);
data=[mean(Cntrl),mean(PD)];
SEMs=[sem(Cntrl),sem(PD)]
figure; 
b = bar(x,data,'LineWidth',2);
b.FaceColor = 'flat';
b.CData(1,:) = [0.5 0.5 1];
b.CData(2,:) =  [.5 1 .5];
%,(std(PDONdDTFCntrlMat,1)./sqrt(size(PDONdDTFCntrlMat,1))));
hold on
er = errorbar(x,data,[sem(Cntrl),sem(PD)],[sem(Cntrl),sem(PD)],'CapSize',18,'LineWidth', 2);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
title('Infromation outflow from PFC in low and mid Beta band')

ylabel('Average ff-DTF (13-20 Hz)') 
ax=gca;
ax.YLim=[0 0.02];
ax.XTickLabel={('Control (PFC)'),('PD Patients On Drug (PFC)')};

ax.FontSize = 32;

hold off



lowerBound=13;
higherBound=30;
[h,p,ci,stats] =ttest2(mean(PDONCntrlffDTF4Mat(:,lowerBound:higherBound),2),mean(PDONffDTF4Mat(:,lowerBound:higherBound),2))

computeCohen_d(mean(PDONCntrlffDTF4Mat(:,lowerBound:higherBound),2),mean(PDONffDTF4Mat(:,lowerBound:higherBound),2))



x = 1:2;
Cntrl=mean(PDONCntrlffDTF4Mat(:,lowerBound:higherBound),2);
PD=mean(PDONffDTF4Mat(:,lowerBound:higherBound),2);
data=[mean(Cntrl),mean(PD)];
figure; 
b = bar(x,data,'LineWidth',2);
b.FaceColor = 'flat';
b.CData(1,:) = [0.5 0.5 1];
b.CData(2,:) = [.5 1 .5];
%,(std(PDONdDTFCntrlMat,1)./sqrt(size(PDONdDTFCntrlMat,1))));
hold on
er = errorbar(x,data,[sem(Cntrl),sem(PD)],[sem(Cntrl),sem(PD)],'CapSize',18,'LineWidth', 2);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
title('Infromation outflow from PFC in beta band')

ylabel('Average sum ff-DTF (13-30 Hz)') 
ax=gca;
ax.YLim=[0 12];
ax.XTickLabel={('Cntrl'),('PD Patients On Drug')};

ax.FontSize = 32;

hold off


%%%%%%%%%%%%%%%%% now alpha

lowerBound=8;
higherBound=12;
[h,p,ci,stats] =ttest2(mean(PDONCntrlffDTF4Mat(:,lowerBound:higherBound),2),mean(PDONffDTF4Mat(:,lowerBound:higherBound),2))

computeCohen_d(mean(PDONCntrlffDTF4Mat(:,lowerBound:higherBound),2),mean(PDONffDTF4Mat(:,lowerBound:higherBound),2))



x = 1:2;
Cntrl=mean(PDONCntrlffDTF4Mat(:,lowerBound:higherBound),2);
PD=mean(PDONffDTF4Mat(:,lowerBound:higherBound),2);
data=[mean(Cntrl),mean(PD)];
figure; 
b = bar(x,data,'LineWidth',2);
b.FaceColor = 'flat';
b.CData(1,:) = [0.5 0.5 1];
b.CData(2,:) = [.5 1 .5];
%,(std(PDONdDTFCntrlMat,1)./sqrt(size(PDONdDTFCntrlMat,1))));
hold on
er = errorbar(x,data,[sem(Cntrl),sem(PD)],[sem(Cntrl),sem(PD)],'CapSize',18,'LineWidth', 2);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
title('Infromation outflow from PFC in alpha band')

ylabel('Average sum ff-DTF (8-12 Hz)') 
ax=gca;
ax.YLim=[0 17];
ax.XTickLabel={('Cntrl'),('PD Patients On Drug')};

ax.FontSize = 32;

hold off


%%%%%%%%%%%%%%%%% now theta

lowerBound=4;
higherBound=7;
[h,p,ci,stats] =ttest2(mean(PDONCntrlffDTF4Mat(:,lowerBound:higherBound),2),mean(PDONffDTF4Mat(:,lowerBound:higherBound),2))

computeCohen_d(mean(PDONCntrlffDTF4Mat(:,lowerBound:higherBound),2),mean(PDONffDTF4Mat(:,lowerBound:higherBound),2))



x = 1:2;
Cntrl=mean(PDONCntrlffDTF4Mat(:,lowerBound:higherBound),2);
PD=mean(PDONffDTF4Mat(:,lowerBound:higherBound),2);
data=[mean(Cntrl),mean(PD)];
figure; 
b = bar(x,data,'LineWidth',2);
b.FaceColor = 'flat';
b.CData(1,:) = [0.5 0.5 1];
b.CData(2,:) = [.5 1 .5];
%,(std(PDONdDTFCntrlMat,1)./sqrt(size(PDONdDTFCntrlMat,1))));
hold on
er = errorbar(x,data,[sem(Cntrl),sem(PD)],[sem(Cntrl),sem(PD)],'CapSize',18,'LineWidth', 2);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
title('Infromation outflow from PFC in theta band')

ylabel('Average sum ff-DTF (4-7 Hz)') 
ax=gca;
ax.YLim=[0 20];
ax.XTickLabel={('Cntrl'),('PD Patients On Drug')};

ax.FontSize = 32;

hold off


%%



figure; errorbar(mean(PDOFFffDTFMat,1),(std(PDOFFffDTFMat,1)./sqrt(size(PDOFFffDTFMat,1))));
hold on; errorbar(mean(PDONffDTF4Mat,1),(std(PDONffDTF4Mat,1)./sqrt(size(PDONffDTF4Mat,1))))


figure; stdshade(PDOFFffDTFMat,0.3,'r');
hold on; stdshade(PDONffDTF4Mat,0.3,'green');
title('Infromation outflow from PFC across different frequencies')
xlabel('Frequency') 
ylabel('Mean ff-DTF') 
pLabel="PD Patients On Drug (PFC Components)";
h= findobj(gca,'Type','Line');
% legend([h(1),h(2)],{'Control',pLabel},'Location','northeast')
legend('boxoff')
ax = gca;
ax.FontSize = 32;
ax.YLim=[0 0.05];
set(get(get(h(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(h(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
legend({'PD Patients Off Drug (PFC Components)',pLabel},'Location','northeast')


%%%%%%%%%%%%%%%% both  gamma
lowerBound=52;
higherBound=54;
[h,p]=ttest2(mean(PDOFFffDTFMat(:,lowerBound:higherBound),2),mean(PDONffDTF4Mat(:,lowerBound:higherBound),2))

x = 1:2;
Cntrl=mean(PDOFFffDTFMat(:,lowerBound:higherBound),2);
PD=mean(PDONffDTF4Mat(:,lowerBound:higherBound),2);
data=[mean(Cntrl),mean(PD)];
figure; 
b = bar(x,data,'LineWidth',2);
b.FaceColor = 'flat';
b.CData(1,:) = [1 0.5 0.5];
b.CData(2,:) = [.5 1 0.5];
%,(std(PDONdDTFCntrlMat,1)./sqrt(size(PDONdDTFCntrlMat,1))));
hold on
er = errorbar(x,data,[sem(Cntrl),sem(PD)],[sem(Cntrl),sem(PD)],'CapSize',18,'LineWidth', 2);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
title('Infromation outflow from PFC in gamma band')

ylabel('Average ff-DTF (52-54 Hz)') 
ax=gca;
ax.YLim=[0 0.02];
ax.XTickLabel={('PD Patients Off Drug'),('PD Patients On Drug')};

ax.FontSize = 32;

hold off



%%%%%%%%%%%%%%%% both  low and mid beta
lowerBound=13;
higherBound=20;
[h,p]=ttest2(mean(PDOFFffDTFMat(:,lowerBound:higherBound),2),mean(PDONffDTF4Mat(:,lowerBound:higherBound),2))

x = 1:2;
Cntrl=mean(PDOFFffDTFMat(:,lowerBound:higherBound),2);
PD=mean(PDONffDTF4Mat(:,lowerBound:higherBound),2);
data=[mean(Cntrl),mean(PD)];
figure; 
b = bar(x,data,'LineWidth',2);
b.FaceColor = 'flat';
b.CData(1,:) = [1 0.5 0.5];
b.CData(2,:) = [.5 1 0.5];
%,(std(PDONdDTFCntrlMat,1)./sqrt(size(PDONdDTFCntrlMat,1))));
hold on
er = errorbar(x,data,[sem(Cntrl),sem(PD)],[sem(Cntrl),sem(PD)],'CapSize',18,'LineWidth', 2);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
title('Infromation outflow from PFC in low and mid beta band')

ylabel('Average ff-DTF (13-20 Hz)') 
ax=gca;
ax.YLim=[0 0.02];
ax.XTickLabel={('PD Patients Off Drug'),('PD Patients On Drug')};

ax.FontSize = 32;

hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% both alpha


lowerBound=8;
higherBound=12;
[h,p]=ttest2(mean(PDOFFffDTFMat(:,lowerBound:higherBound),2),mean(PDONffDTF4Mat(:,lowerBound:higherBound),2))

x = 1:2;
Cntrl=mean(PDOFFffDTFMat(:,lowerBound:higherBound),2);
PD=mean(PDONffDTF4Mat(:,lowerBound:higherBound),2);
data=[mean(Cntrl),mean(PD)];
figure; 
b = bar(x,data,'LineWidth',2);
b.FaceColor = 'flat';
b.CData(1,:) = [1 0.5 0.5];
b.CData(2,:) = [.5 1 0.5];
%,(std(PDONdDTFCntrlMat,1)./sqrt(size(PDONdDTFCntrlMat,1))));
hold on
er = errorbar(x,data,[sem(Cntrl),sem(PD)],[sem(Cntrl),sem(PD)],'CapSize',18,'LineWidth', 2);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
title('Infromation outflow from PFC in alpha band')

ylabel('Average sum ff-DTF (8-12 Hz)') 
ax=gca;
ax.YLim=[0 17];
ax.XTickLabel={('PD Patients Off Drug'),('PD Patients On Drug')};

ax.FontSize = 32;

hold off

%%%%%%%%%%%%%%%%%%%%%% both theta


lowerBound=4;
higherBound=7;
[h,p]=ttest2(mean(PDOFFffDTFMat(:,lowerBound:higherBound),2),mean(PDONffDTF4Mat(:,lowerBound:higherBound),2))

x = 1:2;
Cntrl=mean(PDOFFffDTFMat(:,lowerBound:higherBound),2);
PD=mean(PDONffDTF4Mat(:,lowerBound:higherBound),2);
data=[mean(Cntrl),mean(PD)];
figure; 
b = bar(x,data,'LineWidth',2);
b.FaceColor = 'flat';
b.CData(1,:) = [1 0.5 0.5];
b.CData(2,:) = [.5 1 0.5];
%,(std(PDONdDTFCntrlMat,1)./sqrt(size(PDONdDTFCntrlMat,1))));
hold on
er = errorbar(x,data,[sem(Cntrl),sem(PD)],[sem(Cntrl),sem(PD)],'CapSize',18,'LineWidth', 2);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
title('Infromation outflow from PFC in theta band')

ylabel('Average sum ff-DTF (4-7 Hz)') 
ax=gca;
ax.YLim=[0 20];
ax.XTickLabel={('PD Patients Off Drug'),('PD Patients On Drug')};

ax.FontSize = 32;

hold off
%%


PDOFFListLED=[1275 600 520 550 1150 600 400 640 600 100 1175 1796 300 338];
y=mean(PDOFFffDTFMat(:,lowerBound:higherBound),2)';


figure; s=scatter(PDOFFListLED,y,150,'filled'); 
s.MarkerEdgeColor = [1 0.5 0.5];
s.MarkerFaceColor = [1 0.5 0.5];hold on;  
coef = polyfit(PDOFFListLED,y, 1);
h = refline(coef(1), coef(2));
[r,p]=corrcoef(PDOFFListLED,y);
h.LineWidth=4;
h.Color='k';
title('Correlation Between Beta Band Information Flow and LED')
xlabel('LED (mg)') 
ylabel('Beta Band ff-DTF (sum)') 

legend('boxoff')
ax = gca;
ax.FontSize = 32;
ax.YLim=[0 35];
% set(get(get(h(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% set(get(get(h(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
legend({'PD Patients Off Drug','{\it R^2}=.00001, p=0.99'},'Location','northwest')



PDListONDrugLED=[697 600 250 400 300 200 500 600 700 60  600];
y2=mean(PDONffDTF4Mat(:,lowerBound:higherBound),2)';

figure; s=scatter(PDListONDrugLED,y2,150,'filled');
s.MarkerEdgeColor = [0.5 1 0.5];
s.MarkerFaceColor = [0.5 1 0.5];hold on;  
coef = polyfit(PDListONDrugLED,y2, 1);
h = refline(coef(1), coef(2));
[r,p]=corrcoef(PDListONDrugLED,y2);
h.LineWidth=4;
h.Color='k';
title('Correlation Between Beta Band Information Flow and LED')
xlabel('LED (mg)') 
ylabel('Beta Band ff-DTF (sum)') 

legend('boxoff')
ax = gca;
ax.FontSize = 32;
ax.YLim=[0 35];
% set(get(get(h(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% set(get(get(h(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
legend({'PD Patients On Drug','{\it R^2}=.12, p=0.28'},'Location','northwest')




