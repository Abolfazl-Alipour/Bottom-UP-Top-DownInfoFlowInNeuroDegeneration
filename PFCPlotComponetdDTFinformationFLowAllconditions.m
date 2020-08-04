
clear
close all
destinationPath='E:\Abolfazl\OtherProjs\PDPredictiveCoding\Data\restingStatesLastWave(hopefully)\DownlaodedDataset\ProcessedData\';

load([destinationPath,'CntrlListCumPFCCompsdDTF.mat'],'CntrlListCumPFCComps');

upperFrequencyBound=58;
destinationPath='E:\Abolfazl\OtherProjs\PDPredictiveCoding\Data\restingStatesLastWave(hopefully)\DownlaodedDataset\ProcessedData\';

CntrlList=[894	908	8010	906	903	8060	893	909	911	895	913	900		899	914	910	890	891	912	905	904	892	902	901		897		907];

CntrldDTFMat=nan(length(CntrlList),upperFrequencyBound);
for i=1:length(CntrlList)
    subjectID=CntrlList(i);
    RPDCFileName=strcat('CompdDTFMatrx_',int2str(subjectID),'.mat');
    load(strcat(destinationPath,RPDCFileName))
    PFCComps=CntrlListCumPFCComps{i};
    % Removes components that are not located in the PFC
    dDTF=dDTF(:,PFCComps,:,:);
    tmpTime=mean(dDTF(:,:,1:58,:),4);
    infOutflow=nan(1,size(tmpTime,3));
    
    for layerNum=1:size(tmpTime,3)
        tmp=tmpTime(:,:,layerNum);
        tmp(1:(size(tmp,1)+1):end)=0; %zeroing out the diagonal values
        infOutflow(1,layerNum)=mean(mean(tmp));
    end
    CntrldDTFMat(i,:)=infOutflow;
end


load([destinationPath,'PDListOFFDrugCumPFCCompsdDTF_OFFDRUG.mat'],'PDListOFFDrugCumPFCComps');

% PDList=[801	802	803	804	805	806	807	808	809	810	811	813	814	815	816	817	818	819	820	821	822	823	824	825	];%826	827	828	829];
PDOFFList=[801	802	803	804	805	806	807	808	809	810	811	813	815	816	817	818	819	820	821	822	823	824	825	827	828	829];
PDOFFListLED=[1275 600 520 550 1150 600 400 640 600 100 1175 1796 300 338];
PDOFFdDTFMat=nan(length(PDOFFList),upperFrequencyBound);
for i=1:length(PDOFFList)
    subjectID=PDOFFList(i);
    RPDCFileName=strcat('CompdDTFMatrx_',int2str(subjectID),'OFFDRUG','.mat');
    load(strcat(destinationPath,RPDCFileName))
    PFCComps=PDListOFFDrugCumPFCComps{i};
    % Removes components that are not located in the PFC
    dDTFNONPFC=dDTF;
    dDTFNONPFC(:,PFCComps,:,:)=[];

    dDTF=dDTF(:,PFCComps,:,:);
    tmpTime=mean(dDTF(:,:,1:58,:),4);
    infOutflow=nan(1,size(tmpTime,3));
    
    for layerNum=1:size(tmpTime,3)
        tmp=tmpTime(:,:,layerNum);
        tmp(1:(size(tmp,1)+1):end)=0; %zeroing out the diagonal values
        infOutflow(1,layerNum)=mean(mean(tmp));
    end
    PDOFFdDTFMat(i,:)=infOutflow;
    
    
    
    tmpTime2=mean(dDTFNONPFC(:,:,1:58,:),4);
    infOutflow=nan(1,size(tmpTime2,3));
    
    for layerNum=1:size(tmpTime2,3)
        tmp2=tmpTime2(:,:,layerNum);
        tmp2(1:(size(tmp2,1)+1):end)=0; %zeroing out the diagonal values
        infOutflow2(1,layerNum)=mean(mean(tmp2));
    end
    PDOFFdDTFMatNONPFC(i,:)=infOutflow2;
end


figure; errorbar(mean(CntrldDTFMat,1),(std(CntrldDTFMat,1)./sqrt(size(CntrldDTFMat,1))));
hold on; errorbar(mean(PDOFFdDTFMat,1),(std(PDOFFdDTFMat,1)./sqrt(size(PDOFFdDTFMat,1))))


figure; stdshade(CntrldDTFMat,0.3,'b');
hold on; stdshade(PDOFFdDTFMat,0.3,'r');
title('Infromation outflow from PFC across different frequencies')
xlabel('Frequency') 
ylabel('Mean dDTF') 
pLabel="PD Patients Off drug (PFC Components)";
h= findobj(gca,'Type','Line');
legend([h(1),h(2)],{'PD Patients Off drug (rest of the Components)',pLabel},'Location','northeast')
legend('boxoff')
ax = gca;
ax.FontSize = 32;
ax.YLim=[0 0.013];
set(get(get(h(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(h(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
legend({'Control (PFC Components)',pLabel},'Location','northeast')



lowerBound=13;
higherBound=30;
[h,p,ci,stats]=ttest2(mean(CntrldDTFMat(:,lowerBound:higherBound),2),mean(PDOFFdDTFMat(:,lowerBound:higherBound),2))

computeCohen_d(mean(CntrldDTFMat(:,lowerBound:higherBound),2),mean(PDOFFdDTFMat(:,lowerBound:higherBound),2))


load([destinationPath,'PDListONDrugCumPFCCompsdDTF_ONDRUG.mat'],'PDListONDrugCumPFCComps');

PDListONDrug=[801	802	803	804	805	806	807	808	809	810	811	813	815	816	817	818	819	820	821	822	823	824	825	827	828	829];
PDListONDrugLED=[697 600 250 400 300 200 500 600 700 60  600];
PDONdDTFMat=nan(length(PDListONDrug),upperFrequencyBound);
for i=1:length(PDListONDrug)
    subjectID=PDListONDrug(i);
    RPDCFileName=strcat('CompdDTFMatrx_',int2str(subjectID),'ONDRUG','.mat');
    load(strcat(destinationPath,RPDCFileName))
    PFCComps=PDListONDrugCumPFCComps{i};
    
    dDTFNONPFC=dDTF;
    dDTFNONPFC(:,PFCComps,:,:)=[];
    
    
    
    tmpTime2=mean(dDTFNONPFC(:,:,1:58,:),4);
    infOutflow=nan(1,size(tmpTime2,3));
    
    for layerNum=1:size(tmpTime2,3)
        tmp2=tmpTime2(:,:,layerNum);
        tmp2(1:(size(tmp2,1)+1):end)=0; %zeroing out the diagonal values
        infOutflow2(1,layerNum)=mean(mean(tmp2));
    end
    PDOFFdDTFMatNONPFC(i,:)=infOutflow2;
    
    
    % Removes components that are not located in the PFC
    
    dDTF=dDTF(:,PFCComps,:,:);
    tmpTime=mean(dDTF(:,:,1:58,:),4);
    infOutflow=nan(1,size(tmpTime,3));
    
    for layerNum=1:size(tmpTime,3)
        tmp=tmpTime(:,:,layerNum);
        tmp(1:(size(tmp,1)+1):end)=0; %zeroing out the diagonal values
        infOutflow(1,layerNum)=mean(mean(tmp));
    end
    PDONdDTFMat(i,:)=infOutflow;
end

figure; errorbar(mean(CntrldDTFMat,1),(std(CntrldDTFMat,1)./sqrt(size(CntrldDTFMat,1))));
hold on; errorbar(mean(PDONdDTFMat,1),(std(PDONdDTFMat,1)./sqrt(size(PDONdDTFMat,1))))


figure; stdshade(CntrldDTFMat,0.3,'b');
hold on; stdshade(PDONdDTFMat,0.3,'green');
title('Infromation outflow from PFC across different frequencies')
xlabel('Frequency') 
ylabel('Sum dDTF') 
pLabel="PD Patients On Drug (PFC Components)";
h= findobj(gca,'Type','Line');
% legend([h(1),h(2)],{'Control',pLabel},'Location','northeast')
legend('boxoff')
ax = gca;
ax.FontSize = 32;
ax.YLim=[0 0.013];
set(get(get(h(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(h(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
legend({'Control (PFC Components)',pLabel},'Location','northeast')


x = 1:2;
Cntrl=mean(CntrldDTFMat(:,lowerBound:higherBound),2);
PD=mean(PDOFFdDTFMat(:,lowerBound:higherBound),2);
data=[mean(Cntrl),mean(PD)];
SEMs=[sem(Cntrl),sem(PD)]
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
title('Infromation outflow from PFC in beta band')

ylabel('Average dDTF (13-30 Hz)') 
ax=gca;
ax.YLim=[0 4];
ax.XTickLabel={('Cntrl'),('PD Patients Off Drug')};

ax.FontSize = 32;

hold off
%%%%%%%%%%%%%%%%%%%%%%% now alpha

lowerBound=8;
higherBound=12;
[h,p,ci,stats]=ttest2(mean(CntrldDTFMat(:,lowerBound:higherBound),2),mean(PDOFFdDTFMat(:,lowerBound:higherBound),2))

computeCohen_d(mean(CntrldDTFMat(:,lowerBound:higherBound),2),mean(PDOFFdDTFMat(:,lowerBound:higherBound),2))



x = 1:2;
Cntrl=mean(CntrldDTFMat(:,lowerBound:higherBound),2);
PD=mean(PDOFFdDTFMat(:,lowerBound:higherBound),2);
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

ylabel('Average dDTF (8-12 Hz)') 
ax=gca;
ax.YLim=[0 5];
ax.XTickLabel={('Cntrl'),('PD Patients Off Drug')};

ax.FontSize = 32;

hold off

%%%%%%%%%%%%%%%%%%%%%%% now theta


lowerBound=4;
higherBound=7;
[h,p,ci,stats]=ttest2(mean(CntrldDTFMat(:,lowerBound:higherBound),2),mean(PDOFFdDTFMat(:,lowerBound:higherBound),2))

computeCohen_d(mean(CntrldDTFMat(:,lowerBound:higherBound),2),mean(PDOFFdDTFMat(:,lowerBound:higherBound),2))



x = 1:2;
Cntrl=mean(CntrldDTFMat(:,lowerBound:higherBound),2);
PD=mean(PDOFFdDTFMat(:,lowerBound:higherBound),2);
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

ylabel('Average dDTF (4-7 Hz)') 
ax=gca;
ax.YLim=[0 6];
ax.XTickLabel={('Cntrl'),('PD Patients Off Drug')};

ax.FontSize = 32;

hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%should ADD 814!!!!!!!!!!?? it is bad just like 826?

% 
% PDListOnDrugCntrl=[894	908	8010	906	903	8060	893	909	911	895	913	900		899	914	910	890	891	912	905	904	892	902	901		897		907];
% PDONdDTFCntrlMat=nan(length(PDListOnDrugCntrl),upperFrequencyBound);
% for i=1:length(PDListOnDrugCntrl)
%     subjectID=PDListOnDrugCntrl(i);
%     RPDCFileName=strcat('CompdDTFMatrx_',int2str(subjectID),'.mat');
%     load(strcat(destinationPath,RPDCFileName))
%     PFCComps=PDListOnDrugCntrlCumPFCComps{i};
%     % Removes components that are not located in the PFC
%     dDTF=dDTF(:,PFCComps,:,:);
%     tmpTime=mean(dDTF(:,:,1:58,:),4);
%     infOutflow=nan(1,size(tmpTime,3));
%     
%     for layerNum=1:size(tmpTime,3)
%         tmp=tmpTime(:,:,layerNum);
%         tmp(1:(size(tmp,1)+1):end)=0; %zeroing out the diagonal values
%         infOutflow(1,layerNum)=mean(mean(tmp));
%     end
%     PDONdDTFCntrlMat(i,:)=infOutflow;
% end


load([destinationPath,'PDListONDrugCumPFCCompsdDTF_ONDRUG.mat'],'PDListONDrugCumPFCComps');

PDListONDrug=[801	802	803	804	805	806	807	808	809	810	811	813	815	816	817	818	819	820	821	822	823	824	825	827	828	829];
PDListONDrugLED=[697 600 250 400 300 200 500 600 700 60  600];
PDONdDTFMat=nan(length(PDListONDrug),upperFrequencyBound);
for i=1:length(PDListONDrug)
    subjectID=PDListONDrug(i);
    RPDCFileName=strcat('CompdDTFMatrx_',int2str(subjectID),'ONDRUG','.mat');
    load(strcat(destinationPath,RPDCFileName))
    PFCComps=PDListONDrugCumPFCComps{i};
    
    dDTFNONPFC=dDTF;
    dDTFNONPFC(:,PFCComps,:,:)=[];
    
    
    
    tmpTime2=mean(dDTFNONPFC(:,:,1:58,:),4);
    infOutflow=nan(1,size(tmpTime2,3));
    
    for layerNum=1:size(tmpTime2,3)
        tmp2=tmpTime2(:,:,layerNum);
        tmp2(1:(size(tmp2,1)+1):end)=0; %zeroing out the diagonal values
        infOutflow2(1,layerNum)=mean(mean(tmp2));
    end
    PDOFFdDTFMatNONPFC(i,:)=infOutflow2;
    
    
    % Removes components that are not located in the PFC
    
    dDTF=dDTF(:,PFCComps,:,:);
    tmpTime=mean(dDTF(:,:,1:58,:),4);
    infOutflow=nan(1,size(tmpTime,3));
    
    for layerNum=1:size(tmpTime,3)
        tmp=tmpTime(:,:,layerNum);
        tmp(1:(size(tmp,1)+1):end)=0; %zeroing out the diagonal values
        infOutflow(1,layerNum)=mean(mean(tmp));
    end
    PDONdDTFMat(i,:)=infOutflow;
end

figure; errorbar(mean(CntrldDTFMat,1),(std(CntrldDTFMat,1)./sqrt(size(CntrldDTFMat,1))));
hold on; errorbar(mean(PDONdDTFMat,1),(std(PDONdDTFMat,1)./sqrt(size(PDONdDTFMat,1))))


figure; stdshade(CntrldDTFMat,0.3,'b');
hold on; stdshade(PDONdDTFMat,0.3,'green');
title('Infromation outflow from PFC across different frequencies')
xlabel('Frequency') 
ylabel('Sum dDTF') 
pLabel="PD Patients On Drug (PFC Components)";
h= findobj(gca,'Type','Line');
% legend([h(1),h(2)],{'Control',pLabel},'Location','northeast')
legend('boxoff')
ax = gca;
ax.FontSize = 32;
ax.YLim=[0 0.013];
set(get(get(h(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(h(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
legend({'Control (PFC Components)',pLabel},'Location','northeast')

lowerBound=13;
higherBound=20;
[h,p,ci,stats]=ttest2(mean(CntrldDTFMat(:,lowerBound:higherBound),2),mean(PDONdDTFMat(:,lowerBound:higherBound),2))

computeCohen_d(mean(CntrldDTFMat(:,lowerBound:higherBound),2),mean(PDONdDTFMat(:,lowerBound:higherBound),2))

x = 1:2;
Cntrl=mean(CntrldDTFMat(:,lowerBound:higherBound),2);
PD=mean(PDONdDTFMat(:,lowerBound:higherBound),2);
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

ylabel('Average dDTF (13-20 Hz)') 
ax=gca;
ax.YLim=[0 0.0035];
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

ylabel('Average sum dDTF (8-12 Hz)') 
ax=gca;
ax.YLim=[0 5];
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

ylabel('Average sum dDTF (4-7 Hz)') 
ax=gca;
ax.YLim=[0 6];
ax.XTickLabel={('Cntrl'),('PD Patients On Drug')};

ax.FontSize = 32;

hold off

%%



figure; errorbar(mean(PDOFFdDTFMat,1),(std(PDOFFdDTFMat,1)./sqrt(size(PDOFFdDTFMat,1))));
hold on; errorbar(mean(PDONdDTFMat,1),(std(PDONdDTFMat,1)./sqrt(size(PDONdDTFMat,1))))


figure; stdshade(PDOFFdDTFMat,0.3,'r');
hold on; stdshade(PDONdDTFMat,0.3,'green');
title('Infromation outflow from PFC across different frequencies')
xlabel('Frequency') 
ylabel('Sum dDTF') 
pLabel="PD Patients On Drug (PFC Components)";
h= findobj(gca,'Type','Line');
% legend([h(1),h(2)],{'Control',pLabel},'Location','northeast')
legend('boxoff')
ax = gca;
ax.FontSize = 32;
ax.YLim=[0 0.013];
set(get(get(h(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(h(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
legend({'PD Patients Off Drug (PFC Components)',pLabel},'Location','northeast')


%% plotting both
%%%%%%%%%%%%%%%%%%%%beta
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
title('Infromation outflow from PFC in beta band')

ylabel('Average sum dDTF (13-30 Hz)') 
ax=gca;
ax.XTickLabel={('PD Patients Off Drug'),('PD Patients On Drug')};

ax.FontSize = 32;
ax.YLim=[0 4];
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
title('Infromation outflow from PFC in alpha band')

ylabel('Average sum dDTF (8-12 Hz)') 
ax=gca;
ax.XTickLabel={('PD Patients Off Drug'),('PD Patients On Drug')};

ax.FontSize = 32;
ax.YLim=[0 5];
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
title('Infromation outflow from PFC in theta band')

ylabel('Average sum dDTF (4-7 Hz)') 
ax=gca;
ax.XTickLabel={('PD Patients Off Drug'),('PD Patients On Drug')};

ax.FontSize = 32;
ax.YLim=[0 6];
hold off
%%%%%%%%%%%%%%%%%        CORRELATIONs


PDOFFListLED=[1275 600 520 550 1150 600 400 640 600 100 1175 1796 300 338];
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
ax.YLim=[0 10];
% set(get(get(h(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% set(get(get(h(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
legend({'PD Patients Off Drug','{\it R^2}=.0025, p=0.87'},'Location','northwest')



PDListONDrugLED=[697 600 250 400 300 200 500 600 700 60  600];
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
ax.FontSize = 32;
ax.YLim=[0 10];
% set(get(get(h(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% set(get(get(h(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
legend({'PD Patients On Drug','{\it R^2}=.06, p=0.45'},'Location','northwest')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%pop_dipplot( EEG, [1:size(EEG.icaweights,1)] ,'rvrange',[0 100] ,'mri','C:\\Users\\aalipour\\Documents\\MATLAB\\eeglab2019_1\\plugins\\dipfit\\standard_BEM\\standard_mri.mat','normlen','on');



newa=a(:,cell2mat(PFCComps),:,:);
CntrlffDTFMat=nan(length(CntrlList),upperFrequencyBound);
for i=1:length(CntrlList)
    subjectID=CntrlList(i);
    RPDCFileName=strcat('CompffDTFMatrx_',int2str(subjectID),'.mat');
    load(strcat(destinationPath,RPDCFileName))
    
    tmpTime=mean(ffDTF(:,:,1:58,:),4);
    infOutflow=nan(1,size(tmpTime,3));
    
    for layerNum=1:size(tmpTime,3)
        tmp=tmpTime(:,:,layerNum);
        tmp(1:(size(tmp,1)+1):end)=0; %zeroing out the diagonal values
        infOutflow(1,layerNum)=sum(sum(tmp));
    end
    CntrlffDTFMat(i,:)=infOutflow;
end





