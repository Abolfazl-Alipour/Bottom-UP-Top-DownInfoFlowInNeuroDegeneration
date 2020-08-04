EEG = pop_dipfit_settings( EEG, 'hdmfile','C:\\Users\\aalipour\\Documents\\MATLAB\\eeglab2019_1\\plugins\\dipfit\\standard_BEM\\standard_vol.mat','coordformat','MNI','mrifile','C:\\Users\\aalipour\\Documents\\MATLAB\\eeglab2019_1\\plugins\\dipfit\\standard_BEM\\standard_mri.mat','chanfile','C:\\Users\\aalipour\\Documents\\MATLAB\\eeglab2019_1\\plugins\\dipfit\\standard_BEM\\elec\\standard_1005.elc','chansel',[1:EEG.nbchan] );
EEG = pop_multifit(EEG, [1:size(EEG.icaweights,1)] ,'threshold',100,'dipplot','on','plotopt',{'normlen' 'on'});
%pop_dipplot( EEG, [1:size(EEG.icaweights,1)] ,'rvrange',[0 100] ,'mri','C:\\Users\\aalipour\\Documents\\MATLAB\\eeglab2019_1\\plugins\\dipfit\\standard_BEM\\standard_mri.mat','normlen','on');



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
