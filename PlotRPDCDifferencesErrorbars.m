


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%                         %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   Comp rPDC          %%%%%%%%%%%%%%%%%



upperFrequencyBound=58;

destinationPath='PATH TO EEG OBJECTS';

CntrlList=[894	908	8010	906	903	8060	893	909	911	895	913	900		899	914	910	890	891	912	905	904	892	902	901		897		907];
CntrlRPDCMat=nan(length(CntrlList),upperFrequencyBound);
for i=1:length(CntrlList)
    subjectID=CntrlList(i);
    RPDCFileName=strcat('CompRPDCMatrx_',int2str(subjectID),'.mat');
    load(strcat(destinationPath,RPDCFileName))
    
    tmpTime=mean(RPDC(:,:,1:58,:),4);
    infOutflow=nan(1,size(tmpTime,3));
    
    for layerNum=1:size(tmpTime,3)
        tmp=tmpTime(:,:,layerNum);
        tmp(1:(size(tmp,1)+1):end)=0; %zeroing out the diagonal values
        infOutflow(1,layerNum)=mean(mean(tmp));
    end
    CntrlRPDCMat(i,:)=infOutflow;
end

normalizedCntrlRPDCMat=CntrlRPDCMat;

for i=1:size(normalizedCntrlRPDCMat,1)
    normalizedCntrlRPDCMat(i,:)=normalizedCntrlRPDCMat(i,:)/mean(normalizedCntrlRPDCMat(i,:));
end

% PDList=[801	802	803	804	805	806	807	808	809	810	811	813	814	815	816	817	818	819	820	821	822	823	824	825	];%826	827	828	829];
PDList=[801	802	803	804	805	806	807	808	809	810	811	813	815	816	817	818	819	820	821	822	823	824	825	827	828	829];%826	;

PDOFFListLED=[1275 600 520 550 1150 600 400 640 600 100 1175 1796 300 338];

PDOFFRPDCMat=nan(length(PDList),upperFrequencyBound);
for i=1:length(PDList)
    subjectID=PDList(i);
    RPDCFileName=strcat('CompRPDCMatrx_',int2str(subjectID),'OFFDRUG','.mat');
    load(strcat(destinationPath,RPDCFileName))
    
    tmpTime=mean(RPDC(:,:,1:58,:),4);
    infOutflow=nan(1,size(tmpTime,3));
    % going through all frequencies now:
    for layerNum=1:size(tmpTime,3)
        tmp=tmpTime(:,:,layerNum);
        tmp(1:(size(tmp,1)+1):end)=0; %zeroing out the diagonal values
       %  summing up all values in a given frequency
        infOutflow(1,layerNum)=mean(mean(tmp));
    end
    PDOFFRPDCMat(i,:)=infOutflow;
end

normalizedPDOFFRPDCMat=PDOFFRPDCMat;
for i=1:size(normalizedPDOFFRPDCMat,1)
    normalizedPDOFFRPDCMat(i,:)=normalizedPDOFFRPDCMat(i,:)/mean(normalizedPDOFFRPDCMat(i,:));
end

figure; errorbar(mean(normalizedCntrlRPDCMat,1),(std(normalizedCntrlRPDCMat,1)./sqrt(size(normalizedCntrlRPDCMat,1))));
hold on; errorbar(mean(normalizedPDOFFRPDCMat,1),(std(normalizedPDOFFRPDCMat,1)./sqrt(size(normalizedPDOFFRPDCMat,1))))


figure; errorbar(mean(CntrlRPDCMat,1),(std(CntrlRPDCMat,1)./sqrt(size(CntrlRPDCMat,1))));
hold on; errorbar(mean(PDOFFRPDCMat,1),(std(PDOFFRPDCMat,1)./sqrt(size(PDOFFRPDCMat,1))))


figure; stdshade(normalizedCntrlRPDCMat,0.3,'b');
hold on; stdshade(normalizedPDOFFRPDCMat,0.3,'r');
title('Infromation flow across different frequencies')
xlabel('Frequency') 
ylabel('Mean rPDC (normalized)') 
pLabel="PD Patients Off Drug";
h= findobj(gca,'Type','Line');
% legend([h(1),h(2)],{'Control',pLabel},'Location','northeast')
legend('boxoff')
ax = gca;
ax.FontSize = 32;
set(get(get(h(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(h(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
legend({'Control',pLabel},'Location','northeast')
ax.YLim=[.4 1.6];







%%%%%%
PDListOnDrugCntrl=[894	908	8010	906	903	8060	893	909	911	895	913	900		899	914	910	890	891	912	905	904	892	902	901		897		907];
PDONCntrlRPDCMat=nan(length(PDListOnDrugCntrl),upperFrequencyBound);
for i=1:length(PDListOnDrugCntrl)
    subjectID=PDListOnDrugCntrl(i);
    RPDCFileName=strcat('CompRPDCMatrx_',int2str(subjectID),'.mat');
    load(strcat(destinationPath,RPDCFileName))
    
    tmpTime=mean(RPDC(:,:,1:58,:),4);
    infOutflow=nan(1,size(tmpTime,3));
    
    for layerNum=1:size(tmpTime,3)
        tmp=tmpTime(:,:,layerNum);
        tmp(1:(size(tmp,1)+1):end)=0; %zeroing out the diagonal values
        infOutflow(1,layerNum)=mean(mean(tmp));
    end
    PDONCntrlRPDCMat(i,:)=infOutflow;
end

normalizedPDONCntrlRPDCMat=PDONCntrlRPDCMat;

for i=1:size(normalizedPDONCntrlRPDCMat,1)
    normalizedPDONCntrlRPDCMat(i,:)=normalizedPDONCntrlRPDCMat(i,:)/mean(normalizedPDONCntrlRPDCMat(i,:));
end

PDListONDrug=[801	802	803	804	805	806	807	808	809	810	811	813	815	816	817	818	819	820	821	822	823	824	825	827	828	829];%826	;

PDListONDrugLED=[697 600 250 400 300 200 500 600 700 60  600];
PDONRPDCMat=nan(length(PDListONDrug),upperFrequencyBound);
for i=1:length(PDListONDrug)
    subjectID=PDListONDrug(i);
    RPDCFileName=strcat('CompRPDCMatrx_',int2str(subjectID),'ONDRUG','.mat');
    load(strcat(destinationPath,RPDCFileName))
    
    tmpTime=mean(RPDC(:,:,1:58,:),4);
    infOutflow=nan(1,size(tmpTime,3));
    
    for layerNum=1:size(tmpTime,3)
        tmp=tmpTime(:,:,layerNum);
        tmp(1:(size(tmp,1)+1):end)=0; %zeroing out the diagonal values
        infOutflow(1,layerNum)=mean(mean(tmp));
    end
    PDONRPDCMat(i,:)=infOutflow;
end


normalizedPDONRPDCMat=PDONRPDCMat;
for i=1:size(normalizedPDONRPDCMat,1)
    normalizedPDONRPDCMat(i,:)=normalizedPDONRPDCMat(i,:)/mean(normalizedPDONRPDCMat(i,:));
end





figure; errorbar(mean(normalizedPDONCntrlRPDCMat,1),(std(normalizedPDONCntrlRPDCMat,1)./sqrt(size(normalizedPDONCntrlRPDCMat,1))));
hold on; errorbar(mean(normalizedPDONRPDCMat,1),(std(normalizedPDONRPDCMat,1)./sqrt(size(normalizedPDONRPDCMat,1))))


figure; stdshade(normalizedPDONCntrlRPDCMat,0.3,'b');
hold on; stdshade(normalizedPDONRPDCMat,0.3,'green');
title('Infromation flow across different frequencies')
xlabel('Frequency') 
ylabel('Mean rPDC (normalized)') 
pLabel="PD Patients On Drug";
h= findobj(gca,'Type','Line');
% legend([h(1),h(2)],{'Control',pLabel},'Location','northeast')
legend('boxoff')
ax = gca;
ax.FontSize = 32;
set(get(get(h(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(h(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
legend({'Control',pLabel},'Location','northeast')
ax.YLim=[.4 1.6];




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  GAMMA
lowerBound=30;
higherBound=58;
PDOFFAnovaData=mean(normalizedPDOFFRPDCMat(:,lowerBound:higherBound),2)';
PDOFFAnovaLabels={};
for hh=1:length(PDOFFAnovaData)
    PDOFFAnovaLabels{hh}= 'PDOFF';
end

cntrlAnovaData=mean(normalizedPDONCntrlRPDCMat(:,lowerBound:higherBound),2)';
cntrlAnovaLabels={};
for hh=1:length(cntrlAnovaData)
    cntrlAnovaLabels{hh}= 'Cnrl';
end

PDONAnovaData=mean(normalizedPDONRPDCMat(:,lowerBound:higherBound),2)';
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
Cntrl=mean(normalizedPDONCntrlRPDCMat(:,lowerBound:higherBound),2);
PDOff=mean(normalizedPDOFFRPDCMat(:,lowerBound:higherBound),2);
PDOn=mean(normalizedPDONRPDCMat(:,lowerBound:higherBound),2);
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
ylabel('Average normalized rPDC (30-58 Hz)') 
ax=gca;
ax.XTickLabel={('Cntrl'),('PD Off Drug'),('PD On Drug')};
ax.YLim=[.5 1.2];
ax.FontSize = 32;

hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% high beta

lowerBound=20;
higherBound=30;
PDOFFAnovaData=mean(normalizedPDOFFRPDCMat(:,lowerBound:higherBound),2)';
PDOFFAnovaLabels={};
for hh=1:length(PDOFFAnovaData)
    PDOFFAnovaLabels{hh}= 'PDOFF';
end

cntrlAnovaData=mean(normalizedPDONCntrlRPDCMat(:,lowerBound:higherBound),2)';
cntrlAnovaLabels={};
for hh=1:length(cntrlAnovaData)
    cntrlAnovaLabels{hh}= 'Cnrl';
end

PDONAnovaData=mean(normalizedPDONRPDCMat(:,lowerBound:higherBound),2)';
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
Cntrl=mean(normalizedPDONCntrlRPDCMat(:,lowerBound:higherBound),2);
PDOff=mean(normalizedPDOFFRPDCMat(:,lowerBound:higherBound),2);
PDOn=mean(normalizedPDONRPDCMat(:,lowerBound:higherBound),2);
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
title('Infromation flow in high beta band')

ylabel('Average normalized rPDC (20-30 Hz)') 
ax=gca;
ax.XTickLabel={('Cntrl'),('PD Off Drug'),('PD On Drug')};
ax.YLim=[.8 1.4];
ax.FontSize = 32;

hold off
  



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% alpha

lowerBound=8;
higherBound=12;
PDOFFAnovaData=mean(normalizedPDOFFRPDCMat(:,lowerBound:higherBound),2)';
PDOFFAnovaLabels={};
for hh=1:length(PDOFFAnovaData)
    PDOFFAnovaLabels{hh}= 'PDOFF';
end

cntrlAnovaData=mean(normalizedPDONCntrlRPDCMat(:,lowerBound:higherBound),2)';
cntrlAnovaLabels={};
for hh=1:length(cntrlAnovaData)
    cntrlAnovaLabels{hh}= 'Cnrl';
end

PDONAnovaData=mean(normalizedPDONRPDCMat(:,lowerBound:higherBound),2)';
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
Cntrl=mean(normalizedPDONCntrlRPDCMat(:,lowerBound:higherBound),2);
PDOff=mean(normalizedPDOFFRPDCMat(:,lowerBound:higherBound),2);
PDOn=mean(normalizedPDONRPDCMat(:,lowerBound:higherBound),2);
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
title('Infromation flow in alpha band')

ylabel('Average normalized rPDC (8-12 Hz)') 
ax=gca;
ax.XTickLabel={('Cntrl'),('PD Off Drug'),('PD On Drug')};
ax.YLim=[.8 1.4];
ax.FontSize = 32;

hold off









%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% theta

lowerBound=4;
higherBound=7;
PDOFFAnovaData=mean(normalizedPDOFFRPDCMat(:,lowerBound:higherBound),2)';
PDOFFAnovaLabels={};
for hh=1:length(PDOFFAnovaData)
    PDOFFAnovaLabels{hh}= 'PDOFF';
end

cntrlAnovaData=mean(normalizedPDONCntrlRPDCMat(:,lowerBound:higherBound),2)';
cntrlAnovaLabels={};
for hh=1:length(cntrlAnovaData)
    cntrlAnovaLabels{hh}= 'Cnrl';
end

PDONAnovaData=mean(normalizedPDONRPDCMat(:,lowerBound:higherBound),2)';
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
Cntrl=mean(normalizedPDONCntrlRPDCMat(:,lowerBound:higherBound),2);
PDOff=mean(normalizedPDOFFRPDCMat(:,lowerBound:higherBound),2);
PDOn=mean(normalizedPDONRPDCMat(:,lowerBound:higherBound),2);
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
title('Infromation flow in theta band')

ylabel('Average normalized rPDC (4-7 Hz)') 
ax=gca;
ax.XTickLabel={('Cntrl'),('PD Off Drug'),('PD On Drug')};
ax.YLim=[0 1];
ax.FontSize = 32;

hold off








%%%%%%
PDListOnDrugCntrl=[894	908	8010	906	903	8060	893	909	911	895	913	900		899	914	910	890	891	912	905	904	892	902	901		897		907];
PDONCntrlRPDCMat=nan(length(PDListOnDrugCntrl),upperFrequencyBound);
for i=1:length(PDListOnDrugCntrl)
    subjectID=PDListOnDrugCntrl(i);
    RPDCFileName=strcat('CompRPDCMatrx_',int2str(subjectID),'.mat');
    load(strcat(destinationPath,RPDCFileName))
    
    tmpTime=mean(RPDC(:,:,1:58,:),4);
    infOutflow=nan(1,size(tmpTime,3));
    
    for layerNum=1:size(tmpTime,3)
        tmp=tmpTime(:,:,layerNum);
        tmp(1:(size(tmp,1)+1):end)=0; %zeroing out the diagonal values
        infOutflow(1,layerNum)=mean(mean(tmp));
    end
    PDONCntrlRPDCMat(i,:)=infOutflow;
end

normalizedPDONCntrlRPDCMat=PDONCntrlRPDCMat;

for i=1:size(normalizedPDONCntrlRPDCMat,1)
    normalizedPDONCntrlRPDCMat(i,:)=normalizedPDONCntrlRPDCMat(i,:)/mean(normalizedPDONCntrlRPDCMat(i,:));
end

PDListONDrug=[801	802	803	804	805	806	807	808	809	810	811	813	815	816	817	818	819	820	821	822	823	824	825	827	828	829];%826	;

PDListONDrugLED=[697 600 250 400 300 200 500 600 700 60  600];
PDONRPDCMat=nan(length(PDListONDrug),upperFrequencyBound);
for i=1:length(PDListONDrug)
    subjectID=PDListONDrug(i);
    RPDCFileName=strcat('CompRPDCMatrx_',int2str(subjectID),'ONDRUG','.mat');
    load(strcat(destinationPath,RPDCFileName))
    
    tmpTime=mean(RPDC(:,:,1:58,:),4);
    infOutflow=nan(1,size(tmpTime,3));
    
    for layerNum=1:size(tmpTime,3)
        tmp=tmpTime(:,:,layerNum);
        tmp(1:(size(tmp,1)+1):end)=0; %zeroing out the diagonal values
        infOutflow(1,layerNum)=mean(mean(tmp));
    end
    PDONRPDCMat(i,:)=infOutflow;
end


normalizedPDONRPDCMat=PDONRPDCMat;
for i=1:size(normalizedPDONRPDCMat,1)
    normalizedPDONRPDCMat(i,:)=normalizedPDONRPDCMat(i,:)/mean(normalizedPDONRPDCMat(i,:));
end





figure; errorbar(mean(normalizedPDONCntrlRPDCMat,1),(std(normalizedPDONCntrlRPDCMat,1)./sqrt(size(normalizedPDONCntrlRPDCMat,1))));
hold on; errorbar(mean(normalizedPDONRPDCMat,1),(std(normalizedPDONRPDCMat,1)./sqrt(size(normalizedPDONRPDCMat,1))))


figure; stdshade(normalizedPDONCntrlRPDCMat,0.3,'b');
hold on; stdshade(normalizedPDONRPDCMat,0.3,'green');
title('Infromation flow across different frequencies')
xlabel('Frequency') 
ylabel('Mean rPDC (normalized)') 
pLabel="PD Patients On Drug";
h= findobj(gca,'Type','Line');
% legend([h(1),h(2)],{'Control',pLabel},'Location','northeast')
legend('boxoff')
ax = gca;
ax.FontSize = 32;
set(get(get(h(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(h(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
legend({'Control',pLabel},'Location','northeast')
ax.YLim=[.4 1.6];


lowerBound=30;
higherBound=58;
[h,p,ci,stats]=ttest2(mean(normalizedPDONCntrlRPDCMat(:,lowerBound:higherBound),2),mean(normalizedPDONRPDCMat(:,lowerBound:higherBound),2))

computeCohen_d(mean(normalizedPDONCntrlRPDCMat(:,lowerBound:higherBound),2),mean(normalizedPDONRPDCMat(:,lowerBound:higherBound),2))


x = 1:2;
Cntrl=mean(normalizedPDONCntrlRPDCMat(:,lowerBound:higherBound),2);
PD=mean(normalizedPDONRPDCMat(:,lowerBound:higherBound),2);
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

ylabel('Average normalized rPDC (30-58 Hz)') 
ax=gca;
ax.XTickLabel={('Cntrl'),('PD Patients On Drug')};
ax.YLim=[.5 1.2];
ax.FontSize = 32;

hold off

%%%%%%%%%%%%%%%%%%%%%%% high beta



lowerBound=20;
higherBound=30;
[h,p,ci,stats]=ttest2(mean(normalizedPDONCntrlRPDCMat(:,lowerBound:higherBound),2),mean(normalizedPDONRPDCMat(:,lowerBound:higherBound),2))

computeCohen_d(mean(normalizedPDONCntrlRPDCMat(:,lowerBound:higherBound),2),mean(normalizedPDONRPDCMat(:,lowerBound:higherBound),2))


x = 1:2;
Cntrl=mean(normalizedPDONCntrlRPDCMat(:,lowerBound:higherBound),2);
PD=mean(normalizedPDONRPDCMat(:,lowerBound:higherBound),2);
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
title('Infromation flow in high beta band')

ylabel('Average normalized rPDC (20-30 Hz)') 
ax=gca;
ax.XTickLabel={('Cntrl'),('PD Patients On Drug')};
ax.YLim=[.8 1.4];
ax.FontSize = 32;

hold off


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   alpha

lowerBound=8;
higherBound=12;
[h,p,ci,stats]=ttest2(mean(normalizedPDONCntrlRPDCMat(:,lowerBound:higherBound),2),mean(normalizedPDONRPDCMat(:,lowerBound:higherBound),2))

computeCohen_d(mean(normalizedPDONCntrlRPDCMat(:,lowerBound:higherBound),2),mean(normalizedPDONRPDCMat(:,lowerBound:higherBound),2))


x = 1:2;
Cntrl=mean(normalizedPDONCntrlRPDCMat(:,lowerBound:higherBound),2);
PD=mean(normalizedPDONRPDCMat(:,lowerBound:higherBound),2);
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

ylabel('Average normalized rPDC (8-12 Hz)') 
ax=gca;
ax.XTickLabel={('Cntrl'),('PD Patients On Drug')};
ax.YLim=[.8 1.4];
ax.FontSize = 32;

hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% theta

lowerBound=4;
higherBound=7;
[h,p,ci,stats]=ttest2(mean(normalizedPDONCntrlRPDCMat(:,lowerBound:higherBound),2),mean(normalizedPDONRPDCMat(:,lowerBound:higherBound),2))

computeCohen_d(mean(normalizedPDONCntrlRPDCMat(:,lowerBound:higherBound),2),mean(normalizedPDONRPDCMat(:,lowerBound:higherBound),2))


x = 1:2;
Cntrl=mean(normalizedPDONCntrlRPDCMat(:,lowerBound:higherBound),2);
PD=mean(normalizedPDONRPDCMat(:,lowerBound:higherBound),2);
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

ylabel('Average normalized rPDC (4-7 Hz)') 
ax=gca;
ax.XTickLabel={('Cntrl'),('PD Patients On Drug')};
ax.YLim=[0 1];
ax.FontSize = 32;

hold off


%%%%%%%%%%%%%%%%%%%%%


figure; errorbar(mean(normalizedPDOFFRPDCMat,1),(std(normalizedPDOFFRPDCMat,1)./sqrt(size(normalizedPDOFFRPDCMat,1))));
hold on; errorbar(mean(normalizedPDONRPDCMat,1),(std(normalizedPDONRPDCMat,1)./sqrt(size(normalizedPDONRPDCMat,1))))


figure; stdshade(normalizedPDOFFRPDCMat,0.3,'r');
hold on; stdshade(normalizedPDONRPDCMat,0.3,'green');
title('Infromation flow across different frequencies')
xlabel('Frequency') 
ylabel('Mean rPDC (normalized)') 
pLabel="PD Patients On Drug";
h= findobj(gca,'Type','Line');
% legend([h(1),h(2)],{'Control',pLabel},'Location','northeast')
legend('boxoff')
ax = gca;
ax.FontSize = 32;
ax.YLim=[.4 1.6];
set(get(get(h(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(h(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
legend({'PD Patients Off Drug',pLabel},'Location','northeast')



lowerBound=30;
higherBound=40;
[h,p]=ttest2(mean(PDOFFRPDCMat(:,lowerBound:higherBound),2),mean(PDONRPDCMat(:,lowerBound:higherBound),2))






%%%%% plotting both PD bars
%%%%%%%%%%%%%%%%%%%%%  GAMMA
lowerBound=30;
higherBound=50;

x = 1:2;
Cntrl=mean(normalizedPDOFFRPDCMat(:,lowerBound:higherBound),2);
PD=mean(normalizedPDONRPDCMat(:,lowerBound:higherBound),2);
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

ylabel('Average normalized rPDC (30-58 Hz)') 
ax=gca;
ax.XTickLabel={('PD Patients Off Drug'),('PD Patients On Drug')};
ax.YLim=[.5 1.2];
ax.FontSize = 32;

hold off


%%%%%%%%%%%%%%%%%%%%%  high beta
lowerBound=20;
higherBound=30;

x = 1:2;
Cntrl=mean(normalizedPDOFFRPDCMat(:,lowerBound:higherBound),2);
PD=mean(normalizedPDONRPDCMat(:,lowerBound:higherBound),2);
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
title('Infromation flow in high beta band')

ylabel('Average normalized rPDC (20-30 Hz)') 
ax=gca;
ax.XTickLabel={('PD Patients Off Drug'),('PD Patients On Drug')};
ax.YLim=[.8 1.4];
ax.FontSize = 32;

hold off




%%%%%%%%%%%%%%%%%%%%%  alpha
lowerBound=8;
higherBound=12;

x = 1:2;
Cntrl=mean(normalizedPDOFFRPDCMat(:,lowerBound:higherBound),2);
PD=mean(normalizedPDONRPDCMat(:,lowerBound:higherBound),2);
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

ylabel('Average normalized rPDC (8-12 Hz)') 
ax=gca;
ax.XTickLabel={('PD Patients Off Drug'),('PD Patients On Drug')};

ax.YLim=[.8 1.4];
ax.FontSize = 32;

hold off





%%%%%%%%%%%%%%%%%%%%%  Theta
lowerBound=4;
higherBound=7;

x = 1:2;
Cntrl=mean(normalizedPDOFFRPDCMat(:,lowerBound:higherBound),2);
PD=mean(normalizedPDONRPDCMat(:,lowerBound:higherBound),2);
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

ylabel('Average normalized rPDC (4-7 Hz)') 
ax=gca;
ax.XTickLabel={('PD Patients Off Drug'),('PD Patients On Drug')};

ax.YLim=[0 1];
ax.FontSize = 32;

hold off





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lowerBound=30;
higherBound=50;

PDOFFListLED=[697	1275	600	600	250	520	550	1150	400	300	200	600	500	400	640	600	600	700	60	1650	1000	1175	600	1796	300	338];
y=mean(normalizedPDOFFRPDCMat(:,lowerBound:higherBound),2)';


figure; s=scatter(PDOFFListLED,y,150,'filled'); 
s.MarkerEdgeColor = [1 0.5 0.5];
s.MarkerFaceColor = [1 0.5 0.5];hold on;  
coef = polyfit(PDOFFListLED,y, 1);
h = refline(coef(1), coef(2));
[r,p]=corrcoef(PDOFFListLED,y);
h.LineWidth=4;
h.Color='k';
title('Correlation Between Gamma Band Information Flow and LED')
xlabel('LED (mg)') 
ylabel('Normalized gamma band rPDC (30-58 Hz)') 

legend('boxoff')
ax = gca;
ax.FontSize = 28;
ax.YLim=[0.8 1.1];
% set(get(get(h(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% set(get(get(h(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
legend({'PD Patients Off Drug','{\it R^2}=.0092, p=0.6'},'Location','northeast')




PDListONDrugLED=[697	1275	600	600	250	520	550	1150	400	300	200	600	500	400	640	600	600	700	60	1650	1000	1175	600	1796	300	338];
y2=mean(normalizedPDONRPDCMat(:,lowerBound:higherBound),2)';

figure; s=scatter(PDListONDrugLED,y2,150,'filled');
s.MarkerEdgeColor = [0.5 1 0.5];
s.MarkerFaceColor = [0.5 1 0.5];hold on;  
coef = polyfit(PDListONDrugLED,y2, 1);
h = refline(coef(1), coef(2));
[r,p]=corrcoef(PDListONDrugLED,y2);
h.LineWidth=4;
h.Color='k';
title('Correlation Between Gamma Band Information Flow and LED')
xlabel('LED (mg)') 
ylabel('Normalized gamma band rPDC (30-58 Hz)') 

legend('boxoff')
ax = gca;
ax.FontSize = 28;
ax.YLim=[0.8 1.1];
% set(get(get(h(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% set(get(get(h(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
legend({'PD Patients On Drug','{\it R^2}=.0376, p=0.34'},'Location','northeast')



