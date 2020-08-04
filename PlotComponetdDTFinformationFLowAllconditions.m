

%%%%%%%%%%%%%%%%%%      Fifth wave          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


upperFrequencyBound=58;
destinationPath='E:\Abolfazl\OtherProjs\PDPredictiveCoding\Data\restingStatesLastWave(hopefully)\DownlaodedDataset\ProcessedData\';
GammaYLim=0.003;
BetaYLim=.002;
AlphaYLim=0.0025;
ThetaYLim=0.003;

PDOFFList=[801	802	803	804	805	806	807	808	809	810	811	813	815	816	817	818	819	820	821	822	823	824	825	827	828	829];%826	;
PDOFFListLED=[697	1275	600	600	250	520	550	1150	400	300	200	600	500	400	640	600	600	700	60	1650	1000	1175	600	1796	300	338];
PDOFFdDTFMat=nan(length(PDOFFList),upperFrequencyBound);
for i=1:length(PDOFFList)
    subjectID=PDOFFList(i);
    RPDCFileName=strcat('CompdDTFMatrx_',int2str(subjectID),'OFFDRUG','.mat');
    load(strcat(destinationPath,RPDCFileName))
    
    tmpTime=mean(dDTF(:,:,1:58,:),4);
    infOutflow=nan(1,size(tmpTime,3));
    
    for layerNum=1:size(tmpTime,3)
        tmp=tmpTime(:,:,layerNum);
        tmp(1:(size(tmp,1)+1):end)=0; %zeroing out the diagonal values
        infOutflow(1,layerNum)=mean(mean(tmp));
    end
    PDOFFdDTFMat(i,:)=infOutflow;
end

CntrlList=[894	908	8010	906	903	8060	893	909	911	895	913	900		899	914	910	890	891	912	905	904	892	902	901		897		907];
%excluded 896 898 8070(too short)

CntrldDTFMat=nan(length(CntrlList),upperFrequencyBound);
for i=1:length(CntrlList)
    subjectID=CntrlList(i);
    RPDCFileName=strcat('CompdDTFMatrx_',int2str(subjectID),'.mat');
    load(strcat(destinationPath,RPDCFileName))
    
    tmpTime=mean(dDTF(:,:,1:58,:),4);
    infOutflow=nan(1,size(tmpTime,3));
    
    for layerNum=1:size(tmpTime,3)
        tmp=tmpTime(:,:,layerNum);
        tmp(1:(size(tmp,1)+1):end)=0; %zeroing out the diagonal values
        infOutflow(1,layerNum)=mean(mean(tmp));
    end
    CntrldDTFMat(i,:)=infOutflow;
end



figure; errorbar(mean(CntrldDTFMat,1),(std(CntrldDTFMat,1)./sqrt(size(CntrldDTFMat,1))));
hold on; errorbar(mean(PDOFFdDTFMat,1),(std(PDOFFdDTFMat,1)./sqrt(size(PDOFFdDTFMat,1))))


figure; stdshade(CntrldDTFMat,0.3,'b');
hold on; stdshade(PDOFFdDTFMat,0.3,'r');
title('Infromation flow across different frequencies')
xlabel('Frequency') 
ylabel('Mean dDTF') 
pLabel="PD Patients Off drug";
h= findobj(gca,'Type','Line');
legend([h(1),h(2)],{'Control',pLabel},'Location','northeast')
legend('boxoff')
ax = gca;
ax.FontSize = 32;
ax.YLim=[0 0.006];
set(get(get(h(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(h(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
legend({'Control',pLabel},'Location','northeast')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


PDListOnDrugCntrl=[894	908	8010	906	903	8060	893	909	911	895	913	900	899	914	910	890	891	912	905	904	892	902	901		897		907];
PDONdDTFCntrlMat=nan(length(PDListOnDrugCntrl),upperFrequencyBound);
for i=1:length(PDListOnDrugCntrl)
    subjectID=PDListOnDrugCntrl(i);
    RPDCFileName=strcat('CompdDTFMatrx_',int2str(subjectID),'.mat');
    load(strcat(destinationPath,RPDCFileName))
    
    tmpTime=mean(dDTF(:,:,1:58,:),4);
    infOutflow=nan(1,size(tmpTime,3));
    
    for layerNum=1:size(tmpTime,3)
        tmp=tmpTime(:,:,layerNum);
        tmp(1:(size(tmp,1)+1):end)=0; %zeroing out the diagonal values
        infOutflow(1,layerNum)=mean(mean(tmp));
    end
    PDONdDTFCntrlMat(i,:)=infOutflow;
end

PDListONDrug=[801	802	803	804	805	806	807	808	809	810	811	813	815	816	817	818	819	820	821	822	823	824	825	827	828	829];
PDListONDrugLED=[697	1275	600	600	250	520	550	1150	400	300	200	600	500	400	640	600	600	700	60	1650	1000	1175	600	1796	300	338];
PDONdDTFMat=nan(length(PDListONDrug),upperFrequencyBound);
for i=1:length(PDListONDrug)
    subjectID=PDListONDrug(i);
    RPDCFileName=strcat('CompdDTFMatrx_',int2str(subjectID),'ONDRUG','.mat');
    load(strcat(destinationPath,RPDCFileName))
    tmpTime=mean(dDTF(:,:,1:58,:),4);
    infOutflow=nan(1,size(tmpTime,3));
    
    for layerNum=1:size(tmpTime,3)
        tmp=tmpTime(:,:,layerNum);
        tmp(1:(size(tmp,1)+1):end)=0; %zeroing out the diagonal values
        infOutflow(1,layerNum)=mean(mean(tmp));
    end
    PDONdDTFMat(i,:)=infOutflow;
end



figure; errorbar(mean(PDONdDTFCntrlMat,1),(std(PDONdDTFCntrlMat,1)./sqrt(size(PDONdDTFCntrlMat,1))));
hold on; errorbar(mean(PDONdDTFMat,1),(std(PDONdDTFMat,1)./sqrt(size(PDONdDTFMat,1))))


figure; stdshade(PDONdDTFCntrlMat,0.3,'b');
hold on; stdshade(PDONdDTFMat,0.3,'green');
title('Infromation flow across different frequencies')
xlabel('Frequency') 
ylabel('Mean dDTF') 
pLabel="PD Patients On Drug";
h= findobj(gca,'Type','Line');
% legend([h(1),h(2)],{'Control',pLabel},'Location','northeast')
legend('boxoff')
ax = gca;
ax.FontSize = 32;
ax.YLim=[0 0.006];
set(get(get(h(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(h(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
legend({'Control',pLabel},'Location','northeast')





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GAMMA
lowerBound=52;
higherBound=54;
PDOFFAnovaData=mean(PDOFFdDTFMat(:,lowerBound:higherBound),2)';
PDOFFAnovaLabels={};
for hh=1:length(PDOFFAnovaData)
    PDOFFAnovaLabels{hh}= 'PDOFF';
end

cntrlAnovaData=mean(CntrldDTFMat(:,lowerBound:higherBound),2)';
cntrlAnovaLabels={};
for hh=1:length(cntrlAnovaData)
    cntrlAnovaLabels{hh}= 'Cnrl';
end

PDONAnovaData=mean(PDONdDTFMat(:,lowerBound:higherBound),2)';
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
Cntrl=mean(CntrldDTFMat(:,lowerBound:higherBound),2);
PDOff=mean(PDOFFdDTFMat(:,lowerBound:higherBound),2);
PDOn=mean(PDONdDTFMat(:,lowerBound:higherBound),2);
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

ylabel('Average dDTF (52-54 Hz)') 
ax=gca;
ax.YLim=[0 GammaYLim];
ax.XTickLabel={('PD Off Drug'),('Cntrl'),('PD On Drug')};

ax.FontSize = 32;

hold off



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  BETA
lowerBound=13;
higherBound=30;
PDOFFAnovaData=mean(PDOFFdDTFMat(:,lowerBound:higherBound),2)';
PDOFFAnovaLabels={};
for hh=1:length(PDOFFAnovaData)
    PDOFFAnovaLabels{hh}= 'PDOFF';
end

cntrlAnovaData=mean(CntrldDTFMat(:,lowerBound:higherBound),2)';
cntrlAnovaLabels={};
for hh=1:length(cntrlAnovaData)
    cntrlAnovaLabels{hh}= 'Cnrl';
end

PDONAnovaData=mean(PDONdDTFMat(:,lowerBound:higherBound),2)';
PDONAnovaLabels={};
for hh=1:length(PDONAnovaData)
    PDONAnovaLabels{hh}= 'PDON';
end

[p,tbl,stats]=anova1([PDOFFAnovaData,cntrlAnovaData,PDONAnovaData],[PDOFFAnovaLabels,cntrlAnovaLabels,PDONAnovaLabels]);
[COMPARISON,MEANS,H,GNAMES] = multcompare(stats);
[GNAMES(COMPARISON(:,1)), GNAMES(COMPARISON(:,2)), num2cell(COMPARISON(:,3:6))]

FEXTRACT=[tbl{2,3},tbl{3,3},tbl{2,5},p];
TUKEYEXTRACT=[MEANS COMPARISON(:,6)];
% 
% [h,p,ci,stats]=ttest2(mean(CntrlffDTFMat(:,lowerBound:higherBound),2),mean(PDOFFffDTFMat(:,lowerBound:higherBound),2))

% computeCohen_d(mean(CntrlffDTFMat(:,lowerBound:higherBound),2),mean(PDOFFffDTFMat(:,lowerBound:higherBound),2))

x = 1:3;
Cntrl=mean(CntrldDTFMat(:,lowerBound:higherBound),2);
PDOff=mean(PDOFFdDTFMat(:,lowerBound:higherBound),2);
PDOn=mean(PDONdDTFMat(:,lowerBound:higherBound),2);
data=[mean(PDOff),mean(Cntrl),mean(PDOn)];
% SEMs=[sem(PDOff),sem(Cntrl),sem(PDOn)]
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
title('Infromation flow in beta band')

ylabel('Average dDTF (13-30 Hz)') 
ax=gca;
ax.YLim=[0 BetaYLim];
ax.XTickLabel={('PD Off Drug'),('Cntrl'),('PD On Drug')};

ax.FontSize = 32;

hold off
 

%%%%%%%%%%%%%%%%%%%%%%% now alpha

lowerBound=8;
higherBound=12;
PDOFFAnovaData=mean(PDOFFdDTFMat(:,lowerBound:higherBound),2)';
PDOFFAnovaLabels={};
for hh=1:length(PDOFFAnovaData)
    PDOFFAnovaLabels{hh}= 'PDOFF';
end

cntrlAnovaData=mean(CntrldDTFMat(:,lowerBound:higherBound),2)';
cntrlAnovaLabels={};
for hh=1:length(cntrlAnovaData)
    cntrlAnovaLabels{hh}= 'Cnrl';
end

PDONAnovaData=mean(PDONdDTFMat(:,lowerBound:higherBound),2)';
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
Cntrl=mean(CntrldDTFMat(:,lowerBound:higherBound),2);
PDOff=mean(PDOFFdDTFMat(:,lowerBound:higherBound),2);
PDOn=mean(PDONdDTFMat(:,lowerBound:higherBound),2);
data=[mean(PDOff),mean(Cntrl),mean(PDOn)];
% SEMs=[sem(PDOff),sem(Cntrl),sem(PDOn)]
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

title('Infromation flow in alpha band')

ylabel('Average dDTF (8-12 Hz)') 
ax=gca;
ax.YLim=[0 AlphaYLim];
ax.XTickLabel={('PD Off Drug'),('Cntrl'),('PD On Drug')};

ax.FontSize = 32;

hold off
%%%%%%%%%%%%%%%%%%%%%%% now theta


lowerBound=4;
higherBound=7;
PDOFFAnovaData=mean(PDOFFdDTFMat(:,lowerBound:higherBound),2)';
PDOFFAnovaLabels={};
for hh=1:length(PDOFFAnovaData)
    PDOFFAnovaLabels{hh}= 'PDOFF';
end

cntrlAnovaData=mean(CntrldDTFMat(:,lowerBound:higherBound),2)';
cntrlAnovaLabels={};
for hh=1:length(cntrlAnovaData)
    cntrlAnovaLabels{hh}= 'Cnrl';
end

PDONAnovaData=mean(PDONdDTFMat(:,lowerBound:higherBound),2)';
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
Cntrl=mean(CntrldDTFMat(:,lowerBound:higherBound),2);
PDOff=mean(PDOFFdDTFMat(:,lowerBound:higherBound),2);
PDOn=mean(PDONdDTFMat(:,lowerBound:higherBound),2);
data=[mean(PDOff),mean(Cntrl),mean(PDOn)];
% SEMs=[sem(PDOff),sem(Cntrl),sem(PDOn)]
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

title('Infromation flow in theta band')

ylabel('Average dDTF (4-7 Hz)') 
ax=gca;
ax.YLim=[0 0.003];
ax.XTickLabel={('PD Off Drug'),('Cntrl'),('PD On Drug')};

ax.FontSize = 32;

hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GAMMA
lowerBound=52;
higherBound=54;
[h,p,ci,stats]=ttest2(mean(PDONdDTFCntrlMat(:,lowerBound:higherBound),2),mean(PDONdDTFMat(:,lowerBound:higherBound),2))

computeCohen_d(mean(PDONdDTFCntrlMat(:,lowerBound:higherBound),2),mean(PDONdDTFMat(:,lowerBound:higherBound),2))




x = 1:2;
Cntrl=mean(PDONdDTFCntrlMat(:,lowerBound:higherBound),2);
PD=mean(PDONdDTFMat(:,lowerBound:higherBound),2);
data=[mean(Cntrl),mean(PD)];
SEMs=[sem(Cntrl),sem(PD)]
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
title('Infromation flow in gamma band')

ylabel('Average dDTF (52-54 Hz)') 
ax=gca;
ax.YLim=[0 GammaYLim];
ax.XTickLabel={('Cntrl'),('PD Patients On Drug')};

ax.FontSize = 32;

hold off



%%%%%%%%%%%%%%%%%%%%%%%%%  BETA

lowerBound=13;
higherBound=30;
[h,p,ci,stats]=ttest2(mean(PDONdDTFCntrlMat(:,lowerBound:higherBound),2),mean(PDONdDTFMat(:,lowerBound:higherBound),2))

computeCohen_d(mean(PDONdDTFCntrlMat(:,lowerBound:higherBound),2),mean(PDONdDTFMat(:,lowerBound:higherBound),2))




x = 1:2;
Cntrl=mean(PDONdDTFCntrlMat(:,lowerBound:higherBound),2);
PD=mean(PDONdDTFMat(:,lowerBound:higherBound),2);
data=[mean(Cntrl),mean(PD)];
SEMs=[sem(Cntrl),sem(PD)]
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
title('Infromation flow in beta band')

ylabel('Average dDTF (13-30 Hz)') 
ax=gca;
ax.YLim=[0 BetaYLim];
ax.XTickLabel={('Cntrl'),('PD Patients On Drug')};

ax.FontSize = 32;

hold off

%%%%%%%%%%%%%%%%%%%  now alpha

lowerBound=8;
higherBound=12;
[h,p,ci,stats]=ttest2(mean(PDONdDTFCntrlMat(:,lowerBound:higherBound),2),mean(PDONdDTFMat(:,lowerBound:higherBound),2))

computeCohen_d(mean(PDONdDTFCntrlMat(:,lowerBound:higherBound),2),mean(PDONdDTFMat(:,lowerBound:higherBound),2))

x = 1:2;
Cntrl=mean(PDONdDTFCntrlMat(:,lowerBound:higherBound),2);
PD=mean(PDONdDTFMat(:,lowerBound:higherBound),2);
data=[mean(Cntrl),mean(PD)];
SEMs=[sem(Cntrl),sem(PD)]
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
title('Infromation flow in alpha band')

ylabel('Average dDTF (8-12 Hz)') 
ax=gca;
ax.YLim=[0 AlphaYLim];
ax.XTickLabel={('Cntrl'),('PD Patients On Drug')};

ax.FontSize = 32;

hold off
%%%%%%%%%%%%%%%%%%%%%%%% now theta

lowerBound=4;
higherBound=7;
[h,p,ci,stats]=ttest2(mean(PDONdDTFCntrlMat(:,lowerBound:higherBound),2),mean(PDONdDTFMat(:,lowerBound:higherBound),2))

computeCohen_d(mean(PDONdDTFCntrlMat(:,lowerBound:higherBound),2),mean(PDONdDTFMat(:,lowerBound:higherBound),2))

x = 1:2;
Cntrl=mean(PDONdDTFCntrlMat(:,lowerBound:higherBound),2);
PD=mean(PDONdDTFMat(:,lowerBound:higherBound),2);
data=[mean(Cntrl),mean(PD)];
SEMs=[sem(Cntrl),sem(PD)]
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
title('Infromation flow in theta band')

ylabel('Average dDTF (4-7 Hz)') 
ax=gca;
ax.YLim=[0 0.003];
ax.XTickLabel={('Cntrl'),('PD Patients On Drug')};

ax.FontSize = 32;

hold off

%%



figure; errorbar(mean(PDOFFdDTFMat,1),(std(PDOFFdDTFMat,1)./sqrt(size(PDOFFdDTFMat,1))));
hold on; errorbar(mean(PDONdDTFMat,1),(std(PDONdDTFMat,1)./sqrt(size(PDONdDTFMat,1))))


figure; stdshade(PDOFFdDTFMat,0.3,'r');
hold on; stdshade(PDONdDTFMat,0.3,'green');
title('Infromation flow across different frequencies')
xlabel('Frequency') 
ylabel('Mean dDTF') 
pLabel="PD Patients On Drug";
h= findobj(gca,'Type','Line');
% legend([h(1),h(2)],{'Control',pLabel},'Location','northeast')
legend('boxoff')
ax = gca;
ax.FontSize = 32;
ax.YLim=[0 0.006];
set(get(get(h(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(h(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
legend({'PD Patients Off Drug',pLabel},'Location','northeast')


%% plotting both
%%%%%%%%%%%%%%%%%%%%GAMMA
lowerBound=52;
higherBound=54;
[h,p]=ttest2(mean(PDOFFdDTFMat(:,lowerBound:higherBound),2),mean(PDONdDTFMat(:,lowerBound:higherBound),2))


x = 1:2;
Cntrl=mean(PDOFFdDTFMat(:,lowerBound:higherBound),2);
PD=mean(PDONdDTFMat(:,lowerBound:higherBound),2);
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
title('Infromation flow in gamma band')

ylabel('Average dDTF (52-54 Hz)') 
ax=gca;
ax.XTickLabel={('PD Patients Off Drug'),('PD Patients On Drug')};

ax.FontSize = 32;
ax.YLim=[0 GammaYLim];
hold off
%%%%%%%%%%%%%%%%%%%% BETA both
lowerBound=13;
higherBound=30;
[h,p]=ttest2(mean(PDOFFdDTFMat(:,lowerBound:higherBound),2),mean(PDONdDTFMat(:,lowerBound:higherBound),2))


x = 1:2;
Cntrl=mean(PDOFFdDTFMat(:,lowerBound:higherBound),2);
PD=mean(PDONdDTFMat(:,lowerBound:higherBound),2);
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
title('Infromation flow in beta band')

ylabel('Average dDTF (13-30 Hz)') 
ax=gca;
ax.XTickLabel={('PD Patients Off Drug'),('PD Patients On Drug')};

ax.FontSize = 32;
ax.YLim=[0 BetaYLim];
hold off

%%%%%%%%%%%%%%%%%%%%%% alpha both
lowerBound=8;
higherBound=12;
[h,p]=ttest2(mean(PDOFFdDTFMat(:,lowerBound:higherBound),2),mean(PDONdDTFMat(:,lowerBound:higherBound),2))

x = 1:2;
Cntrl=mean(PDOFFdDTFMat(:,lowerBound:higherBound),2);
PD=mean(PDONdDTFMat(:,lowerBound:higherBound),2);
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
title('Infromation flow in alpha band')

ylabel('Average dDTF (8-12 Hz)') 
ax=gca;
ax.XTickLabel={('PD Patients Off Drug'),('PD Patients On Drug')};

ax.FontSize = 32;
ax.YLim=[0 AlphaYLim];
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%% theta both
lowerBound=4;
higherBound=7;
[h,p]=ttest2(mean(PDOFFdDTFMat(:,lowerBound:higherBound),2),mean(PDONdDTFMat(:,lowerBound:higherBound),2))

x = 1:2;
Cntrl=mean(PDOFFdDTFMat(:,lowerBound:higherBound),2);
PD=mean(PDONdDTFMat(:,lowerBound:higherBound),2);
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
title('Infromation flow in theta band')

ylabel('Average dDTF (4-7 Hz)') 
ax=gca;
ax.XTickLabel={('PD Patients Off Drug'),('PD Patients On Drug')};

ax.FontSize = 32;
ax.YLim=[0 ThetaYLim];
hold off
%%%%%%%%%%%%%%%%%        CORRELATIONs




%%%%%%%%%%%Gamma

lowerBound=52;
higherBound=54;
PDOFFListLED=[697	1275	600	600	250	520	550	1150	400	300	200	600	500	400	640	600	600	700	60	1650	1000	1175	600	1796	300	338];

y=mean(PDOFFdDTFMat(:,lowerBound:higherBound),2)';


figure; s=scatter(PDOFFListLED,y,150,'filled'); 
s.MarkerEdgeColor = [1 0.5 0.5];
s.MarkerFaceColor = [1 0.5 0.5];hold on;  
coef = polyfit(PDOFFListLED,y, 1);
h = refline(coef(1), coef(2));
[r,p]=corrcoef(PDOFFListLED,y);
h.LineWidth=4;
h.Color='k';
title('Correlation Between Gamma Band Peak Information Flow and LED')
xlabel('LED (mg)') 
ylabel('Mean Gamma Band dDTF (52-54 Hz)') 

legend('boxoff')
ax = gca;
ax.FontSize = 27;
ax.YLim=[0 0.004];
% set(get(get(h(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% set(get(get(h(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
legend({'PD Patients Off Drug','{\it R^2}=.15, p=.051'},'Location','northwest')



PDListONDrugLED=[697	1275	600	600	250	520	550	1150	400	300	200	600	500	400	640	600	600	700	60	1650	1000	1175	600	1796	300	338];
y2=mean(PDONdDTFMat(:,lowerBound:higherBound),2)';

figure; s=scatter(PDListONDrugLED,y2,150,'filled');
s.MarkerEdgeColor = [0.5 1 0.5];
s.MarkerFaceColor = [0.5 1 0.5];hold on;  
coef = polyfit(PDListONDrugLED,y2, 1);
h = refline(coef(1), coef(2));
[r,p]=corrcoef(PDListONDrugLED,y2);
h.LineWidth=4;
h.Color='k';
title('Correlation Between Gamma Band Peak Information Flow and LED')
xlabel('LED (mg)') 
ylabel('Mean Gamma Band dDTF (52-54 Hz))') 

legend('boxoff')
ax = gca;
ax.FontSize = 27;
ax.YLim=[0 0.004];
% set(get(get(h(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% set(get(get(h(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
legend({'PD Patients On Drug','{\it R^2}=.0086, p=.65'},'Location','northwest')



%%%%%%%%%%%beta

lowerBound=13;
higherBound=30;
PDOFFListLED=[697	1275	600	600	250	520	550	1150	400	300	200	600	500	400	640	600	600	700	60	1650	1000	1175	600	1796	300	338];

y=mean(PDOFFdDTFMat(:,lowerBound:higherBound),2)';


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
ylabel('Beta Band dDTF (sum)') 

legend('boxoff')
ax = gca;
ax.FontSize = 32;
ax.YLim=[0 0.0014];
% set(get(get(h(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% set(get(get(h(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
legend({'PD Patients Off Drug','{\it R^2}=.0209, p=0.48'},'Location','northeast')



PDListONDrugLED=[697	1275	600	600	250	520	550	1150	400	300	200	600	500	400	640	600	600	700	60	1650	1000	1175	600	1796	300	338];
y2=mean(PDONdDTFMat(:,lowerBound:higherBound),2)';

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
ylabel('Beta Band dDTF (sum)') 

legend('boxoff')
ax = gca;
ax.FontSize = 28;
ax.YLim=[0 0.0014];
% set(get(get(h(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% set(get(get(h(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
legend({'PD Patients On Drug','{\it R^2}=.0004, p=.93'},'Location','northeast')



%%%%%%%%%%%%%% alpha
lowerBound=8;
higherBound=12;
PDOFFListLED=[697	1275	600	600	250	520	550	1150	400	300	200	600	500	400	640	600	600	700	60	1650	1000	1175	600	1796	300	338];

y=mean(PDOFFdDTFMat(:,lowerBound:higherBound),2)';


figure; s=scatter(PDOFFListLED,y,150,'filled'); 
s.MarkerEdgeColor = [1 0.5 0.5];
s.MarkerFaceColor = [1 0.5 0.5];hold on;  
coef = polyfit(PDOFFListLED,y, 1);
h = refline(coef(1), coef(2));
[r,p]=corrcoef(PDOFFListLED,y);
h.LineWidth=4;
h.Color='k';
title('Correlation Between Alpha Band Information Flow and LED')
xlabel('LED (mg)') 
ylabel('Alpha Band dDTF (Mean)') 

legend('boxoff')
ax = gca;
ax.FontSize = 28;
ax.YLim=[0 0.0016];
% set(get(get(h(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% set(get(get(h(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
legend({'PD Patients Off Drug','{\it R^2}=.0019, p=0.83'},'Location','northeast')



PDListONDrugLED=[697	1275	600	600	250	520	550	1150	400	300	200	600	500	400	640	600	600	700	60	1650	1000	1175	600	1796	300	338];
y2=mean(PDONdDTFMat(:,lowerBound:higherBound),2)';

figure; s=scatter(PDListONDrugLED,y2,150,'filled');
s.MarkerEdgeColor = [0.5 1 0.5];
s.MarkerFaceColor = [0.5 1 0.5];hold on;  
coef = polyfit(PDListONDrugLED,y2, 1);
h = refline(coef(1), coef(2));
[r,p]=corrcoef(PDListONDrugLED,y2);
h.LineWidth=4;
h.Color='k';
title('Correlation Between Alpha Band Information Flow and LED')
xlabel('LED (mg)') 
ylabel('Alpha Band dDTF (Mean)') 

legend('boxoff')
ax = gca;
ax.FontSize = 28;
ax.YLim=[0 0.0016];
% set(get(get(h(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% set(get(get(h(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
legend({'PD Patients On Drug','{\it R^2}=.0003, p=.94'},'Location','northeast')








%%%%%%%%%%%%%% theta
lowerBound=4
higherBound=7
PDOFFListLED=[697	1275	600	600	250	520	550	1150	400	300	200	600	500	400	640	600	600	700	60	1650	1000	1175	600	1796	300	338];

y=mean(PDOFFdDTFMat(:,lowerBound:higherBound),2)';


figure; s=scatter(PDOFFListLED,y,150,'filled'); 
s.MarkerEdgeColor = [1 0.5 0.5];
s.MarkerFaceColor = [1 0.5 0.5];hold on;  
coef = polyfit(PDOFFListLED,y, 1);
h = refline(coef(1), coef(2));
[r,p]=corrcoef(PDOFFListLED,y);
h.LineWidth=4;
h.Color='k';
title('Correlation Between Thata Band Information Flow and LED')
xlabel('LED (mg)') 
ylabel('Thata Band dDTF (Mean)') 

legend('boxoff')
ax = gca;
ax.FontSize = 28;
ax.YLim=[0 0.0025];
% set(get(get(h(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% set(get(get(h(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
legend({'PD Patients Off Drug','{\it R^2}=.0006, p=0.91'},'Location','northwest')



PDListONDrugLED=[697	1275	600	600	250	520	550	1150	400	300	200	600	500	400	640	600	600	700	60	1650	1000	1175	600	1796	300	338];
y2=mean(PDONdDTFMat(:,lowerBound:higherBound),2)';

figure; s=scatter(PDListONDrugLED,y2,150,'filled');
s.MarkerEdgeColor = [0.5 1 0.5];
s.MarkerFaceColor = [0.5 1 0.5];hold on;  
coef = polyfit(PDListONDrugLED,y2, 1);
h = refline(coef(1), coef(2));
[r,p]=corrcoef(PDListONDrugLED,y2);
h.LineWidth=4;
h.Color='k';
title('Correlation Between Theta Band Information Flow and LED')
xlabel('LED (mg)') 
ylabel('Theta Band dDTF (Mean)') 

legend('boxoff')
ax = gca;
ax.FontSize = 28;
ax.YLim=[0 0.0025];
% set(get(get(h(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% set(get(get(h(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
legend({'PD Patients On Drug','{\it R^2}=.0005, p=0.91'},'Location','northeAst')



