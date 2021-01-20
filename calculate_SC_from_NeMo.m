% calculate SC from NeMo outputs
% predict observed (partial) FC from estimated SC (from NeMo Tool) and see if there is a shift over time in the signaling mechanism that explains the most variance
% Jan 19th 2021

% Load subject chacoconn scores
for i=1:23
    tmp=load(strcat(studydir,'/subject_data/fs86_subj/SUB', num2str(i), '_lesion_1mmMNI_fs86subj_mean_chacoconn.mat'));
    structname=strcat('SUB', num2str(i), 'chacoconn');
    all_patients_chacoconn{i}=full(tmp.(structname));
end

for i=1:23
    A=load(strcat(studydir,'/subject_data/fs86_subj/SUB', num2str(i), '_lesion1mm_nemo_output_chacoconn_fs86subj_allref_denom.mat'));
    avg_connectome{i}=mean(A.allref,3)
end
%% Multiply subject chacoconn scores by avg healhy connectome to get
% 'disrupted' connectome 
for i=1:23
    all_patients_disrupted{i}=all_patients_chacoconn{i}.*avg_connectome{i};
end

% Subtract 'disrupted' connectome from average healthy connectome to get
% 'spared' connectome
for i=1:23
    all_patients_spared{i}=avg_connectome{i}-all_patients_disrupted{i};
end

for i=1:23
    sc=all_patients_spared{i};
    save(strcat(pwd, '/subject_data/SUB',num2str(i), '_sc.mat'), 'sc')
end

