% calculate SC from NeMo outputs
% predict observed (partial) FC from estimated SC (from NeMo Tool) and see if there is a shift over time in the signaling mechanism that explains the most variance
% September 14th 2020

studydir = pwd;

% Load average healthy connectome & individual healthy connectomes.
allref=load(strcat(pwd,'/subject_data/fs86_avg/allref_denom.mat'))
allref=allref.allref_denom;
for i=1:420
    neg_idx{i}=find(allref{i}<0);
end


full3d=[];
for i=1:420
    full3d=cat(3, full3d, full(allref{i}));
end

avg_connectome=squeeze(mean(full3d,3));

% Load subject chacoconn scores
for i=1:23
    tmp=load(strcat(studydir,'/subject_data/SUB', num2str(i), '_lesion_1mmMNI_fs86subj_mean_chacoconn.mat'));
    structname=strcat('SUB', num2str(i), 'chacoconn');
    all_patients_chacoconn{i}=tmp.(structname);
end


%% Multiply subject chacoconn scores by avg healhy connectome to get
% 'disrupted' connectome 
for i=1:23
    all_patients_disrupted{i}=all_patients_chacoconn{i}.*avg_connectome;
end

% Subtract 'disrupted' connectome from average healthy connectome to get
% 'spared' connectome
for i=1:23
    all_patients_spared{i}=avg_connectome-all_patients_disrupted{i};
end

for i=1:23
    sc=all_patients_spared{i};
    save(strcat(pwd, '/subject_data/SUB',num2str(i), '_sc.mat'), 'sc')
end

