
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%                         %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   Comp rPDC          %%%%%%%%%%%%%%%%%



upperFrequencyBound=58;

destinationPath='E:\Abolfazl\OtherProjs\PDPredictiveCoding\Data\restingStatesLastWave(hopefully)\DownlaodedDataset\ProcessedData\';
CntrlList=[894	908	8010	906	903	8060	893	909	911	895	913	900		899	914	910	890	891	912	905	904	892	902	901		897		907];

load ([destinationPath,'CntrlListCumPFCCompsdDTF2.mat'],'CntrlListCumPFCComps');

CntrlRPDCMat=nan(length(CntrlList),upperFrequencyBound);
for i=1:length(CntrlList)
    subjectID=CntrlList(i);
    RPDCFileName=strcat('CompRPDCMatrx_',int2str(subjectID),'.mat');
                            
    load(strcat(destinationPath,RPDCFileName))
    PFCComps=CntrlListCumPFCComps{i};
    % Removes components that are not located in the PFC
    RPDC=RPDC(:,PFCComps,:,:);
%     RPDC(PFCComps,:,:,:)=[];
    tmpTime=mean(RPDC(:,:,1:58,:),4);
    infOutflow=nan(1,size(tmpTime,3));
    
    for layerNum=1:size(tmpTime,3)
        tmp=tmpTime(:,:,layerNum);
        tmp(1:(size(tmp,1)+1):end)=0; %zeroing out the diagonal values
        infOutflow(1,layerNum)=sum(sum(tmp));
    end
    CntrlRPDCMat(i,:)=infOutflow;
end


normalizedCntrlRPDCMat=CntrlRPDCMat;
for i=1:size(normalizedCntrlRPDCMat,1)
    normalizedCntrlRPDCMat(i,:)=normalizedCntrlRPDCMat(i,:)/mean(normalizedCntrlRPDCMat(i,:));
end

load([destinationPath,'PDListOFFDrugCumPFCCompsdDTF_OFFDRUG.mat'],'PDListOFFDrugCumPFCComps');

% PDList=[801	802	803	804	805	806	807	808	809	810	811	813	814	815	816	817	818	819	820	821	822	823	824	825	];%826	827	828	829];
PDList=[801	802	803	804	805	806	807	808	809	810	811	813	815	816	817	818	819	820	821	822	823	824	825	827	828	829];
PDOFFListLED=[1275 600 520 550 1150 600 400 640 600 100 1175 1796 300 338];

PDOFFRPDCMat=nan(length(PDList),upperFrequencyBound);
for i=1:length(PDList)
    subjectID=PDList(i);
    RPDCFileName=strcat('CompRPDCMatrx_',int2str(subjectID),'OFFDRUG','.mat');
    load(strcat(destinationPath,RPDCFileName))
    PFCComps=PDListOFFDrugCumPFCComps{i};
    
    
    
    
    RPDCNONPFC=RPDC;
    RPDCNONPFC(:,PFCComps,:,:)=[];
    
    tmpTime2=mean(RPDCNONPFC(:,:,1:58,:),4);
    infOutflow=nan(1,size(tmpTime2,3));
    
    for layerNum=1:size(tmpTime2,3)
        tmp2=tmpTime2(:,:,layerNum);
        tmp2(1:(size(tmp2,1)+1):end)=0; %zeroing out the diagonal values
        infOutflow2(1,layerNum)=mean(mean(tmp2));
    end
    PDOFFRPDCMatNONPFC(i,:)=infOutflow2;
    
    
    
    
    % Removes components that are not located in the PFC
    RPDC=RPDC(:,PFCComps,:,:);
%     RPDC(PFCComps,:,:,:)=[];
    tmpTime=mean(RPDC(:,:,1:58,:),4);
    infOutflow=nan(1,size(tmpTime,3));
    
    for layerNum=1:size(tmpTime,3)
        tmp=tmpTime(:,:,layerNum);
        tmp(1:(size(tmp,1)+1):end)=0; %zeroing out the diagonal values
        infOutflow(1,layerNum)=mean(mean(tmp));
    end
    PDOFFRPDCMat(i,:)=infOutflow;
end


normalizedPDOFFRPDCMatNONPFC=PDOFFRPDCMatNONPFC;
for i=1:size(normalizedPDOFFRPDCMatNONPFC,1)
    normalizedPDOFFRPDCMatNONPFC(i,:)=normalizedPDOFFRPDCMatNONPFC(i,:)/mean(normalizedPDOFFRPDCMatNONPFC(i,:));
end

normalizedPDOFFRPDCMat=PDOFFRPDCMat;
for i=1:size(normalizedPDOFFRPDCMat,1)
    normalizedPDOFFRPDCMat(i,:)=normalizedPDOFFRPDCMat(i,:)/mean(normalizedPDOFFRPDCMat(i,:));
end

figure; errorbar(mean(normalizedCntrlRPDCMat,1),(std(normalizedCntrlRPDCMat,1)./sqrt(size(normalizedCntrlRPDCMat,1))));
hold on; errorbar(mean(normalizedPDOFFRPDCMat,1),(std(normalizedPDOFFRPDCMat,1)./sqrt(size(normalizedPDOFFRPDCMat,1))))

figure; stdshade(normalizedCntrlRPDCMat ,0.3,'b'); 
hold on; stdshade(normalizedPDOFFRPDCMat,0.3,'r');
title('Infromation flow across different frequencies')
xlabel('Frequency') 
ylabel('Mean rPDC (normalized)') 
pLabel="PD Patients Off Drug (PFC Components)";
h= findobj(gca,'Type','Line');
% legend([h(1),h(2)],{'Control',pLabel},'Location','northeast')
legend('boxoff')
ax = gca;
ax.FontSize = 32;
set(get(get(h(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(h(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
legend({'Cntrl (PFC Components)',pLabel},'Location','northeast')
ax.YLim=[.4 2];



figure; errorbar(mean(normalizedPDOFFRPDCMatNONPFC,1),(std(normalizedPDOFFRPDCMatNONPFC,1)./sqrt(size(normalizedPDOFFRPDCMatNONPFC,1))));
hold on; errorbar(mean(normalizedPDOFFRPDCMat,1),(std(normalizedPDOFFRPDCMat,1)./sqrt(size(normalizedPDOFFRPDCMat,1))))


figure; stdshade(normalizedPDOFFRPDCMat,0.3,'r'); 
hold on; stdshade(normalizedPDOFFRPDCMatNONPFC,0.3,'c');
title('Infromation flow across different frequencies')
xlabel('Frequency') 
ylabel('Mean rPDC (normalized)') 
pLabel="PD Patients Off Drug (rest of the Components)";
h= findobj(gca,'Type','Line');
% legend([h(1),h(2)],{'Control',pLabel},'Location','northeast')
legend('boxoff')
ax = gca;
ax.FontSize = 32;
set(get(get(h(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(h(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
legend({'PD Patients Off Drug (PFC Components)',pLabel},'Location','northeast')
ax.YLim=[.4 2];




load([destinationPath,'PDListONDrugCumPFCCompsdDTF_ONDRUG.mat'],'PDListONDrugCumPFCComps');
PDListONDrug=[801	802	803	804	805	806	807	808	809	810	811	813	815	816	817	818	819	820	821	822	823	824	825	827	828	829];

PDListONDrugLED=[697 600 250 400 300 200 500 600 700 60  600];
PDONRPDCMat=nan(length(PDListONDrug),upperFrequencyBound);
for i=1:length(PDListONDrug)
    subjectID=PDListONDrug(i);
    RPDCFileName=strcat('CompRPDCMatrx_',int2str(subjectID),'.mat');
    load(strcat(destinationPath,RPDCFileName))
    PFCComps=PDListONDrugCumPFCComps{i};
    
    
    
    
    RPDCNONPFC=RPDC;
    RPDCNONPFC(:,PFCComps,:,:)=[];
    
    tmpTime2=mean(RPDCNONPFC(:,:,1:58,:),4);
    infOutflow=nan(1,size(tmpTime2,3));
    
    for layerNum=1:size(tmpTime2,3)
        tmp2=tmpTime2(:,:,layerNum);
        tmp2(1:(size(tmp2,1)+1):end)=0; %zeroing out the diagonal values
        infOutflow2(1,layerNum)=mean(mean(tmp2));
    end
    PDONRPDCMatNONPFC(i,:)=infOutflow2;
    
    
    
    % Removes components that are not located in the PFC
    RPDC=RPDC(PFCComps,:,:,:);
    tmpTime=mean(RPDC(:,:,1:58,:),4);
    infOutflow=nan(1,size(tmpTime,3));
    
    for layerNum=1:size(tmpTime,3)
        tmp=tmpTime(:,:,layerNum);
        tmp(1:(size(tmp,1)+1):end)=0; %zeroing out the diagonal values
        infOutflow(1,layerNum)=mean(mean(tmp));
    end
    PDONRPDCMat(i,:)=infOutflow;
end


normalizedPDONRPDCMatNONPFC=PDONRPDCMatNONPFC;
for i=1:size(normalizedPDONRPDCMatNONPFC,1)
    normalizedPDONRPDCMatNONPFC(i,:)=normalizedPDONRPDCMatNONPFC(i,:)/mean(normalizedPDONRPDCMatNONPFC(i,:));
end

normalizedPDONRPDCMat=PDONRPDCMat;
for i=1:size(normalizedPDONRPDCMat,1)
    normalizedPDONRPDCMat(i,:)=normalizedPDONRPDCMat(i,:)/mean(normalizedPDONRPDCMat(i,:));
end


figure; stdshade(normalizedCntrlRPDCMat ,0.3,'b'); 
hold on; stdshade(normalizedPDONRPDCMat,0.3,'g');
title('Infromation flow across different frequencies')
xlabel('Frequency') 
ylabel('Mean rPDC (normalized)') 
pLabel="PD Patients On Drug (PFC Components)";
h= findobj(gca,'Type','Line');
% legend([h(1),h(2)],{'Control',pLabel},'Location','northeast')
legend('boxoff')
ax = gca;
ax.FontSize = 32;
set(get(get(h(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(h(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
legend({'Cntrl (PFC Components)',pLabel},'Location','northeast')
ax.YLim=[.4 2];





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GAMMA
lowerBound=30;
higherBound=58;

PDOFFAnovaData=mean(normalizedPDOFFRPDCMat(:,lowerBound:higherBound),2)';
PDOFFAnovaLabels={};
for hh=1:length(PDOFFAnovaData)
    PDOFFAnovaLabels{hh}= 'PDOFF';
end

cntrlAnovaData=mean(normalizedCntrlRPDCMat(:,lowerBound:higherBound),2)';
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
Cntrl=mean(normalizedCntrlRPDCMat(:,lowerBound:higherBound),2);
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
ax.XTickLabel={('Cntrl (PFC)'),('PD Off Drug (PFC)'),('PD On Drug (PFC)')};
ax.YLim=[.5 1.1];
ax.FontSize = 28;

hold off




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% low and mid beta

lowerBound=13;
higherBound=20;

PDOFFAnovaData=mean(normalizedPDOFFRPDCMat(:,lowerBound:higherBound),2)';
PDOFFAnovaLabels={};
for hh=1:length(PDOFFAnovaData)
    PDOFFAnovaLabels{hh}= 'PDOFF';
end

cntrlAnovaData=mean(normalizedCntrlRPDCMat(:,lowerBound:higherBound),2)';
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
Cntrl=mean(normalizedCntrlRPDCMat(:,lowerBound:higherBound),2);
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
 
title('Infromation flow in low and mid beta band')

ylabel('Average normalized rPDC (13-20 Hz)') 
ax=gca;
ax.XTickLabel={('Cntrl (PFC)'),('PD Off Drug (PFC)'),('PD On Drug (PFC)')};
ax.YLim=[.8 2];
ax.FontSize = 28;


hold off

% 
% 
% lowerBound=31;
% higherBound=58;
% [h,p,ci,stats]=ttest2(mean(CntrlRPDCMat(:,lowerBound:higherBound),2),mean(PDOFFRPDCMat(:,lowerBound:higherBound),2))
% 
% 
% computeCohen_d(mean(CntrlRPDCMat(:,lowerBound:higherBound),2),mean(PDOFFRPDCMat(:,lowerBound:higherBound),2))
% 
% 
% x = 1:2;
% Cntrl=mean(CntrlRPDCMat(:,lowerBound:higherBound),2);
% PD=mean(PDOFFRPDCMat(:,lowerBound:higherBound),2);
% data=[mean(Cntrl),mean(PD)];
% figure; 
% b = bar(x,data,'LineWidth',2);
% b.FaceColor = 'flat';
% b.CData(1,:) = [0.5 0.5 1];
% b.CData(2,:) = [1 0.5 .5];
% %,(std(PDONdDTFCntrlMat,1)./sqrt(size(PDONdDTFCntrlMat,1))));
% hold on
% er = errorbar(x,data,[sem(Cntrl),sem(PD)],[sem(Cntrl),sem(PD)],'CapSize',18,'LineWidth', 2);    
% er.Color = [0 0 0];                            
% er.LineStyle = 'none';  
% title('Infromation flow in gamma band')
% 
% ylabel('Average sum rPDC (31-58 Hz)') 
% ax=gca;
% ax.XTickLabel={('Cntrl'),('PD Patients Off Drug')};
% ax.YLim=[0 .35];
% ax.FontSize = 32;
% 
% hold off

%%%%%%
% PDListOnDrugCntrl=[894	906 903	911 895 913 899 890 912 905 904 901];
% PDONCntrlRPDCMat=nan(length(PDListOnDrugCntrl),upperFrequencyBound);
% for i=1:length(PDListOnDrugCntrl)
%     subjectID=PDListOnDrugCntrl(i);
%     RPDCFileName=strcat('CompRPDCMatrx_',int2str(subjectID),'.mat');
%     load(strcat(destinationPath,RPDCFileName))
%     PFCComps=PDListOnDrugCntrlCumPFCComps{i};
%     % Removes components that are not located in the PFC
%     RPDC=RPDC(PFCComps,:,:,:);
%     tmpTime=mean(RPDC(:,:,1:58,:),4);
%     infOutflow=nan(1,size(tmpTime,3));
%     
%     for layerNum=1:size(tmpTime,3)
%         tmp=tmpTime(:,:,layerNum);
%         tmp(1:(size(tmp,1)+1):end)=0; %zeroing out the diagonal values
%         infOutflow(1,layerNum)=sum(sum(tmp));
%     end
%     PDONCntrlRPDCMat(i,:)=infOutflow;
% end
% 




figure; errorbar(mean(normalizedPDONRPDCMat,1),(std(normalizedPDONRPDCMat,1)./sqrt(size(normalizedPDONRPDCMat,1))));
hold on; errorbar(mean(normalizedPDONRPDCMatNONPFC,1),(std(normalizedPDONRPDCMatNONPFC,1)./sqrt(size(normalizedPDONRPDCMatNONPFC,1))))


figure; stdshade(normalizedPDONRPDCMat,0.3,'m');
hold on; stdshade(normalizedPDONRPDCMatNONPFC,0.3,'green');
title('Infromation flow across different frequencies')
xlabel('Frequency') 
ylabel('Mean rPDC (normalized') 
pLabel="PD Patients On Drug (rest of the Components)";
h= findobj(gca,'Type','Line');
% legend([h(1),h(2)],{'Control',pLabel},'Location','northeast')
legend('boxoff')
ax = gca;
ax.FontSize = 32;
set(get(get(h(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(h(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
legend({'PD Patients On Drug (PFC Components)',pLabel},'Location','northeast')
ax.YLim=[.4 2];



lowerBound=30;
higherBound=58;
[h,p,ci,stats]=ttest2(mean(normalizedCntrlRPDCMat(:,lowerBound:higherBound),2),mean(normalizedPDONRPDCMat(:,lowerBound:higherBound),2))

computeCohen_d(mean(normalizedCntrlRPDCMat(:,lowerBound:higherBound),2),mean(normalizedPDONRPDCMat(:,lowerBound:higherBound),2))


x = 1:2;
Cntrl=mean(normalizedCntrlRPDCMat(:,lowerBound:higherBound),2);
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
ax.XTickLabel={('Cntrl (PFC)'),('PD Patients On Drug (PFC)')};
ax.YLim=[.5 1.1];
ax.FontSize = 32;

hold off



%%%%%%%%%%%%%%%%%%%%% mid beta



lowerBound=13;
higherBound=20;
[h,p,ci,stats]=ttest2(mean(normalizedCntrlRPDCMat(:,lowerBound:higherBound),2),mean(normalizedPDONRPDCMat(:,lowerBound:higherBound),2))

computeCohen_d(mean(normalizedCntrlRPDCMat(:,lowerBound:higherBound),2),mean(normalizedPDONRPDCMat(:,lowerBound:higherBound),2))


x = 1:2;
Cntrl=mean(normalizedCntrlRPDCMat(:,lowerBound:higherBound),2);
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
title('Infromation flow in low and mid beta band')

ylabel('Average normalized rPDC (13-20 Hz)') 
ax=gca;
ax.XTickLabel={('Cntrl (PFC)'),('PD Patients On Drug (PFC)')};
ax.YLim=[.8 2];
ax.FontSize = 32;

hold off

%%%%%%%%%%%%%%%%%%%%
% lowerBound=30;
% higherBound=58;
% [h,p,ci,stats]=ttest2(mean(normalizedCntrlRPDCMat(:,lowerBound:higherBound),2),mean(normalizedPDOFFRPDCMat(:,lowerBound:higherBound),2))
% 
% 
% computeCohen_d(mean(normalizedCntrlRPDCMat(:,lowerBound:higherBound),2),mean(normalizedPDOFFRPDCMat(:,lowerBound:higherBound),2))
% 
% 
% x = 1:2;
% Cntrl=mean(normalizedCntrlRPDCMat(:,lowerBound:higherBound),2);
% PD=mean(normalizedPDOFFRPDCMat(:,lowerBound:higherBound),2);
% data=[mean(Cntrl),mean(PD)];
% figure; 
% b = bar(x,data,'LineWidth',2);
% b.FaceColor = 'flat';
% b.CData(1,:) = [0.5 0.5 1];
% b.CData(2,:) = [1 0.5 .5];
% %,(std(PDONdDTFCntrlMat,1)./sqrt(size(PDONdDTFCntrlMat,1))));
% hold on
% er = errorbar(x,data,[sem(Cntrl),sem(PD)],[sem(Cntrl),sem(PD)],'CapSize',18,'LineWidth', 2);    
% er.Color = [0 0 0];                            
% er.LineStyle = 'none';  
% title('Infromation flow in gamma band')
% 
% ylabel('Average normalized rPDC (30-58 Hz)') 
% ax=gca;
% ax.XTickLabel={('Cntrl'),('PD Patients Off Drug')};
% ax.YLim=[.5 1.2];
% ax.FontSize = 32;
% 
% hold off




% 
% 
% lowerBound=31;
% higherBound=58;
% [h,p,ci,stats]=ttest2(mean(PDONCntrlRPDCMat(:,lowerBound:higherBound),2),mean(PDONRPDCMat(:,lowerBound:higherBound),2))
% 
% computeCohen_d(mean(PDONCntrlRPDCMat(:,lowerBound:higherBound),2),mean(PDONRPDCMat(:,lowerBound:higherBound),2))
% 
% 
% x = 1:2;
% Cntrl=mean(PDONCntrlRPDCMat(:,lowerBound:higherBound),2);
% PD=mean(PDONRPDCMat(:,lowerBound:higherBound),2);
% data=[mean(Cntrl),mean(PD)];
% figure; 
% b = bar(x,data,'LineWidth',2);
% b.FaceColor = 'flat';
% b.CData(1,:) = [0.5 0.5 1];
% b.CData(2,:) = [.5 1 .5];
% %,(std(PDONdDTFCntrlMat,1)./sqrt(size(PDONdDTFCntrlMat,1))));
% hold on
% er = errorbar(x,data,[sem(Cntrl),sem(PD)],[sem(Cntrl),sem(PD)],'CapSize',18,'LineWidth', 2);    
% er.Color = [0 0 0];                            
% er.LineStyle = 'none';  
% title('Infromation flow in gamma band')
% 
% ylabel('Average sum rPDC (31-58 Hz)') 
% ax=gca;
% ax.XTickLabel={('Cntrl'),('PD Patients On Drug')};
% ax.YLim=[0 .35];
% ax.FontSize = 32;
% 
% hold off

%%%%%%%%%%%%%%%%%%%%%%% BOTH all freqs



figure; stdshade(normalizedPDOFFRPDCMat ,0.3,'r'); 
hold on; stdshade(normalizedPDONRPDCMat,0.3,'g');
title('Infromation flow across different frequencies')
xlabel('Frequency') 
ylabel('Mean rPDC (normalized)') 
pLabel="PD Patients Off Drug (PFC Components)";
h= findobj(gca,'Type','Line');
% legend([h(1),h(2)],{'Control',pLabel},'Location','northeast')
legend('boxoff')
ax = gca;
ax.FontSize = 32;
set(get(get(h(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(h(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
legend({'PD Patients On Drug (PFC Components)',pLabel},'Location','northeast')
ax.YLim=[.4 2];


%%%%% plotting both PD bars


lowerBound=30;
higherBound=58;
[h,p,ci,stats]=ttest2(mean(normalizedPDOFFRPDCMat(:,lowerBound:higherBound),2),mean(normalizedPDONRPDCMat(:,lowerBound:higherBound),2))

computeCohen_d(mean(normalizedPDOFFRPDCMat(:,lowerBound:higherBound),2),mean(normalizedPDONRPDCMat(:,lowerBound:higherBound),2))


x = 1:2;
Cntrl=mean(normalizedPDOFFRPDCMat(:,lowerBound:higherBound),2);
PD=mean(normalizedPDONRPDCMat(:,lowerBound:higherBound),2);
data=[mean(Cntrl),mean(PD)];
figure; 
b = bar(x,data,'LineWidth',2);
b.FaceColor = 'flat';
b.CData(1,:) = [1 0.5 0.5];
b.CData(2,:) = [.5 1 .5];
%,(std(PDONdDTFCntrlMat,1)./sqrt(size(PDONdDTFCntrlMat,1))));
hold on
er = errorbar(x,data,[sem(Cntrl),sem(PD)],[sem(Cntrl),sem(PD)],'CapSize',18,'LineWidth', 2);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
title('Infromation flow in gamma band')

ylabel('Average normalized rPDC (30-58 Hz)') 
ax=gca;
ax.XTickLabel={('PD Patients Off Drug (PFC)'),('PD Patients On Drug (PFC)')};
ax.YLim=[.5 1.1];
ax.FontSize = 26;

hold off



%%%%%%%%%%%%%%%%%%%%% mid beta



lowerBound=13;
higherBound=20;
[h,p,ci,stats]=ttest2(mean(normalizedPDOFFRPDCMat(:,lowerBound:higherBound),2),mean(normalizedPDONRPDCMat(:,lowerBound:higherBound),2))

computeCohen_d(mean(normalizedPDOFFRPDCMat(:,lowerBound:higherBound),2),mean(normalizedPDONRPDCMat(:,lowerBound:higherBound),2))


x = 1:2;
Cntrl=mean(normalizedPDOFFRPDCMat(:,lowerBound:higherBound),2);
PD=mean(normalizedPDONRPDCMat(:,lowerBound:higherBound),2);
data=[mean(Cntrl),mean(PD)];
figure; 
b = bar(x,data,'LineWidth',2);
b.FaceColor = 'flat';
b.CData(1,:) = [1 0.5 0.5];
b.CData(2,:) = [.5 1 .5];
%,(std(PDONdDTFCntrlMat,1)./sqrt(size(PDONdDTFCntrlMat,1))));
hold on
er = errorbar(x,data,[sem(Cntrl),sem(PD)],[sem(Cntrl),sem(PD)],'CapSize',18,'LineWidth', 2);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
title('Infromation flow in low and mid beta band')

ylabel('Average normalized rPDC (13-20 Hz)') 
ax=gca;
ax.XTickLabel={('PD Patients Off Drug (PFC)'),('PD Patients On Drug (PFC)')};
ax.YLim=[.8 2];
ax.FontSize = 26;

hold off

x = 1:2;
Cntrl=mean(PDOFFRPDCMat(:,lowerBound:higherBound),2);
PD=mean(PDONRPDCMat(:,lowerBound:higherBound),2);
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

ylabel('Average sum rPDC (31-58 Hz)') 
ax=gca;
ax.XTickLabel={('PD Patients Off Drug'),('PD Patients On Drug')};
ax.YLim=[0 .35];
ax.FontSize = 32;

hold off


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PDOFFListLED=[1275 600 520 550 1150 600 400 640 600 100 1175 1796 300 338];
y=mean(PDOFFRPDCMat(:,lowerBound:higherBound),2)';


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
ylabel('Gamma Band rPDC (sum)') 

legend('boxoff')
ax = gca;
ax.FontSize = 32;
% set(get(get(h(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% set(get(get(h(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
legend({'PD Patients Off Drug','{\it R^2}=.26, p=0.06'},'Location','northwest')



PDListONDrugLED=[697 600 250 400 300 200 500 600 700 60  600];
y2=mean(PDONRPDCMat(:,lowerBound:higherBound),2)';

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
ylabel('Gamma Band rPDC (sum)') 

legend('boxoff')
ax = gca;
ax.FontSize = 32;
% set(get(get(h(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% set(get(get(h(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
legend({'PD Patients On Drug','{\it R^2}=.44, p=0.02'},'Location','northwest')

%%%%%%%%%%%%%%%%%%%%%


figure; errorbar(mean(PDOFFRPDCMat,1),(std(PDOFFRPDCMat,1)./sqrt(size(PDOFFRPDCMat,1))));
hold on; errorbar(mean(PDONRPDCMat,1),(std(PDONRPDCMat,1)./sqrt(size(PDONRPDCMat,1))))


figure; stdshade(PDOFFRPDCMat,0.3,'r');
hold on; stdshade(PDONRPDCMat,0.3,'green');
title('Infromation flow across different frequencies')
xlabel('Frequency') 
ylabel('Sum rPDC') 
pLabel="PD Patients On Drug";
h= findobj(gca,'Type','Line');
% legend([h(1),h(2)],{'Control',pLabel},'Location','northeast')
legend('boxoff')
ax = gca;
ax.FontSize = 32;
ax.YLim=[.1 .45];
set(get(get(h(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(h(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
legend({'PD Patients Off Drug',pLabel},'Location','northeast')



lowerBound=30;
higherBound=40;
[h,p]=ttest2(mean(PDOFFRPDCMat(:,lowerBound:higherBound),2),mean(PDONRPDCMat(:,lowerBound:higherBound),2))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 
% 
% upperFrequencyBound=58;
% CntrlList=[894,908,8010,906,903,8060,893,909,911,895,913,900,896,899,914,910,890,891,912,905,904,892,902,901,898,897,907];
% CntrlList=[908,8010,8060,893,909,900,914,910,891,892,902];
% 
% CntrlRPDCMat=nan(length(CntrlList),upperFrequencyBound);
% for i=1:length(CntrlList)
%     subjectID=CntrlList(i);
%     RPDCFileName=strcat('RPDCMatrx_',int2str(subjectID),'.mat');
%     load(strcat(destinationPath,RPDCFileName))
%     
%     tmpTime=mean(RPDC(:,:,1:58,:),4);
%     infOutflow=nan(1,size(tmpTime,3));
%     
%     for layerNum=1:size(tmpTime,3)
%         tmp=tmpTime(:,:,layerNum);
%         tmp(1:(size(tmp,1)+1):end)=0; %zeroing out the diagonal values
%         infOutflow(1,layerNum)=sum(sum(tmp));
%     end
%     CntrlRPDCMat(i,:)=infOutflow;
% end
% 
% 
% PDList=[801	802	803	804	805	806	807	808	809	810	811	813	814	815	816	817	818	819	820	821	822	823	824	825	];%826	827	828	829];
% PDList=[802	803	806	807	808	813	816	817	819	823	824	827	828	829];%826	;
% 
% PDOFFRPDCMat=nan(length(PDList),upperFrequencyBound);
% for i=1:length(PDList)
%     subjectID=PDList(i);
%     RPDCFileName=strcat('RPDCMatrx_',int2str(subjectID),'.mat');
%     load(strcat(destinationPath,RPDCFileName))
%     
%     tmpTime=mean(RPDC(:,:,1:58,:),4);
%     infOutflow=nan(1,size(tmpTime,3));
%     
%     for layerNum=1:size(tmpTime,3)
%         tmp=tmpTime(:,:,layerNum);
%         tmp(1:(size(tmp,1)+1):end)=0; %zeroing out the diagonal values
%         infOutflow(1,layerNum)=sum(sum(tmp));
%     end
%     PDOFFRPDCMat(i,:)=infOutflow;
% end
% 
% 
% figure; errorbar(mean(CntrlRPDCMat,1),(std(CntrlRPDCMat,1)./sqrt(size(CntrlRPDCMat,1))));
% hold on; errorbar(mean(PDOFFRPDCMat,1),(std(PDOFFRPDCMat,1)./sqrt(size(PDOFFRPDCMat,1))))
% 
% destinationPath='E:\Abolfazl\OtherProjs\PDPredictiveCoding\Data\restingStatesLastWave(hopefully)\DownlaodedDataset\ProcessedData\';
% sourceFilePath='E:\Abolfazl\OtherProjs\PDPredictiveCoding\Data\restingStatesLastWave(hopefully)\DownlaodedDataset\PD REST\';
% locPath=('C:\Users\aalipour\Documents\MATLAB\eeglab2019_1\plugins\dipfit\standard_BESA\standard-10-5-cap385.elp');
% locPathFrwrdMdel='C:\Users\aalipour\Documents\MATLAB\eeglab2019_1\plugins\dsi\headModel\resources\head_modelColin27_8003_Standard-10-5-Cap339.mat';
% chckModelConsistancy=0;
% 
% 
% for i=1:length(PDList)
%     subjectID=PDList(i);
%     
%     connFromSubjID(sourceFilePath,subjectID,destinationPath,locPath,locPathFrwrdMdel,chckModelConsistancy)
% end
% 
% 
% 
% 
% upperFrequencyBound=58;
% CntrlList=[908,8010,8060,893,909,900,914,910,891,892,902];
% 
% CntrlRPDCMat=nan(length(CntrlList),upperFrequencyBound);
% for i=1:length(CntrlList)
%     subjectID=CntrlList(i);
%     RPDCFileName=strcat('sourcedDTFMatrx_',int2str(subjectID),'.mat');
%     load(strcat(destinationPath,RPDCFileName))
%     
%     tmpTime=mean(CompdDTF(:,:,1:58,:),4);
%     infOutflow=nan(1,size(tmpTime,3));
%     
%     for layerNum=1:size(tmpTime,3)
%         tmp=tmpTime(:,:,layerNum);
%         tmp(1:(size(tmp,1)+1):end)=0; %zeroing out the diagonal values
%         infOutflow(1,layerNum)=sum(sum(tmp));
%     end
%     CntrlRPDCMat(i,:)=infOutflow;
% end
% 
% 
% % PDList=[801	802	803	804	805	806	807	808	809	810	811	813	814	815	816	817	818	819	820	821	822	823	824	825	];%826	827	828	829];
% PDList=[802	803	806	807	808	813	816	817	819	823	824	827	828	829];%826	;
% 
% PDOFFRPDCMat=nan(length(PDList),upperFrequencyBound);
% for i=1:length(PDList)
%     subjectID=PDList(i);
%     RPDCFileName=strcat('sourcedDTFMatrx_',int2str(subjectID),'.mat');
%     load(strcat(destinationPath,RPDCFileName))
%     
%     tmpTime=mean(CompdDTF(:,:,1:58,:),4);
%     infOutflow=nan(1,size(tmpTime,3));
%     
%     for layerNum=1:size(tmpTime,3)
%         tmp=tmpTime(:,:,layerNum);
%         tmp(1:(size(tmp,1)+1):end)=0; %zeroing out the diagonal values
%         infOutflow(1,layerNum)=sum(sum(tmp));
%     end
%     PDOFFRPDCMat(i,:)=infOutflow;
% end
% 
% 
% figure; errorbar(mean(CntrlRPDCMat,1),(std(CntrlRPDCMat,1)./sqrt(size(CntrlRPDCMat,1))));
% hold on; errorbar(mean(PDOFFRPDCMat,1),(std(PDOFFRPDCMat,1)./sqrt(size(PDOFFRPDCMat,1))))
% 
% 
% [h,p]=ttest2(mean(CntrlRPDCMat(:,41:50),2),mean(PDOFFRPDCMat(:,41:50),2))
