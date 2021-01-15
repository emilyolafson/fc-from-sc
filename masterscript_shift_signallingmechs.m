% predict observed (partial) FC from estimated SC (from NeMo Tool) and see if there is a shift over time in the signaling mechanism that explains the most variance
% September 14th 2020

studydir = '/Users/emilyolafson/Documents/Thesis/shortest_paths/';

% Load average healthy connectome & individual healthy connectomes.
load('allref.mat')
avg_connectome=squeeze(mean(allref,1));
all_controls=load('allref.mat');
all_controls=all_controls.allref;
numscans1_11=[4 3 3 3 3 3 3 3 3 3 3;4 3 3 3 3 3 3 4 3 3 3;3 3 3 3 3 3 3 3 3 3 3;4 4 3 3 3 3 3 3 3 3 3 ; 3 3 3 3 3 0 3 3 3 3 3];
numscans12_23=[3 3 3 3 3 2 2 2 2 2 2 2; 3 3 3 4 3 2 2 2 2 2 2 2 ;3 3 3 3 3 2 2 2 0 2 2 2;0 3 3 3 3 2 2 2 0 2 2 2;0 3 4 3 3 2 2 2 0 2 2 2];
nscans = [numscans1_11, numscans12_23];
load('num_visits.mat') %nsess

% Load subject chacoconn scores
for i=1:23
    tmp=load(strcat(studydir,'NeMo2_shen_outputs/SUB', num2str(i), '_lesion_1mmMNI_shen268_mean_chacoconn.mat'));
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

%% Threshold using consistency-based thresholding.
% consistency-based thresholding on the healthy controls to identify the
% edge weights which are observed most consistently across controls.
thr=0.3;

all_controls_symm=[];

for i=1:420
    subj=squeeze(all_controls(i,:,:));
    bottomleft=rot90(fliplr(subj));
    subj=subj+bottomleft;
    all_controls_symm=cat(3, all_controls_symm,subj);
end

%determine weights which are most consistent across 420 healthy controls
all_controls_connectome=all_controls_symm; %input needs to be NxNxm

thresholded_ctl = threshold_consistency(all_controls_connectome,thr);
controls_consistency_thresholded=logical(thresholded_ctl);
imagesc(controls_consistency_thresholded)
sum(sum(controls_consistency_thresholded))

% threshold patients
for i=1:23 
    % apply control binary matrix to all individual patients SC (weighted)
    bottomleft=rot90(fliplr(all_patients_spared{i}));
    all_patients_spared{i}=all_patients_spared{i}+bottomleft;
    all_patients_spared_thresholded{i}=all_patients_spared{i}.*controls_consistency_thresholded;
end

%% Calculate the shortest paths between each brain region (betweenness centrality)

% Convert weight matrices to distance matrices.
% patients
for i=1:23
    patients_distance{i}=all_patients_spared_thresholded{i}.^(-1);
end

% Calculate the betweenness centrality.
clear bc
bc_all=[];
for i=1:23
    bc{i}=distance_wei(patients_distance{i});
    bc_all=cat(3,bc_all, bc{i});
end

for i=1:23
    bc_std{i}=normalize(bc{i});
    bc_all=cat(3,mfpt_all, mfpt_std_symm{i});
end


%% Calculate diffusion distance (Mean First Passage Time)
% patients
for i=1:23
    mfpt{i}=mean_first_passage_time(all_patients_spared_thresholded{i});
end
mfpt_all=[]
for i=1:23
    mfpt_std{i}=normalize(mfpt{i});
    mfpt_std_triu{i}=triu(mfpt_std{i});
    bottomleft=rot90(fliplr(mfpt_std_triu{i}));
    mfpt_std_symm{i}=mfpt_std_triu{i}+bottomleft;
    mfpt_all=cat(3,mfpt_all, mfpt_std_symm{i});
end


%% calculate FC (FC was overwritten by z-score FC)
allts=load('/Users/emilyolafson/Documents/Thesis/SUB1_23_data/ts_GSR_shen268_allsub_allsess.mat')
allts=allts.ts_shen268;
for i=1:23
    for j=1:nsess(i)
       a=cell2mat(allts{i,j});
       imagesc(a)
    end
end


%% Load partial correlation data.
session1=session1.C;
for i=1:23
    pFC{i}=squeeze(session1(i,:,:));
    pFC_xfm{i}=atanh(abs(pFC{i}))+1.0*10^(-4);
end

%% GLM with FC
fc=load('/Users/emilyolafson/Documents/Thesis/SUB1_23_data/FC.mat')
fc=fc.fc;

% take the absolute value of FC correlation coefficients and perform fishers r-to-z
% transformation.
for i=1:23
    for j=1:nsess(i)
        fc_xfm{i,j}=atanh(abs(fc{i,j})); 
       % fc_xfm{i,j}(logical(eye(268)))=0;
    end
end

% Set number of iterations for glmfit.
opts = statset('glmfit');
opts.MaxIter = 1000; % default value for glmfit is 100.

indexes=logical(triu(ones(268,268)));
e=eye(268);
indexes=indexes-e;
indexes=logical(indexes)

clear pval_diff
clear pval_sp
clear beta_diff
clear beta_sp
for i=1:23
    X=[mfpt_std_symm{i}(indexes),bc_std{i}(indexes)];
    for j=1:nsess(i)
        fc_scores=cell2mat(fc_xfm(i,j));
        Y=fc_scores(indexes);
        Y=Y+1.0*10^(-8);
        [b, ~, stats{i,j}]=glmfit(X, Y, 'gamma','link', 'log', 'options', opts);
        pval_diff{i,j}=stats{i,j}.p(2);
        pval_sp{i,j}=stats{i,j}.p(3);

        beta_diff{i,j}=stats{i,j}.beta(2);
        beta_sp{i,j}=stats{i,j}.beta(3);
    end
end


beta_diff{6,5}=single(0)
beta_diff{12,5}=single(0)
beta_diff{12,4}=single(0)
beta_diff{20,5}=single(0)
beta_diff{20,4}=single(0)
beta_diff{20,3}=single(0)
beta_diff=cell2mat(beta_diff);

beta_sp{6,5}=single(0)
beta_sp{12,5}=single(0)
beta_sp{12,4}=single(0)
beta_sp{20,5}=single(0)
beta_sp{20,4}=single(0)
beta_sp{20,3}=single(0)
beta_sp=cell2mat(beta_sp);

pval_diff{6,5}=single(0)
pval_diff{12,5}=single(0)
pval_diff{12,4}=single(0)
pval_diff{20,5}=single(0)
pval_diff{20,4}=single(0)
pval_diff{20,3}=single(0)
pval_diff=cell2mat(pval_diff);

pval_sp{6,5}=single(0)
pval_sp{12,5}=single(0)
pval_sp{12,4}=single(0)
pval_sp{20,5}=single(0)
pval_sp{20,4}=single(0)
pval_sp{20,3}=single(0)
pval_sp=cell2mat(pval_sp);
pvals_FDR=fdr_bh([pval_diff,pval_sp], 0.05)

mean_beta_diff_S1=mean(beta_diff(:,1))
mean_beta_sp_S1=mean(beta_sp(:,1))

mean_beta_diff_S2=mean(beta_diff(:,1))
mean_beta_sp_S2=mean(beta_sp(:,1))

mean_beta_diff_S3=mean(beta_diff(:,1))
mean_beta_sp_S3=mean(beta_sp(:,1))

mean_beta_diff_S4=mean(beta_diff(:,1))
mean_beta_sp_S4=mean(beta_sp(:,1))

figure(1)
for i=1:20
    subplot(4,5,i)
    plot(beta_diff(i,:), 'r')
    hold on
    %plot(beta_sp(i,:), 'b')
    %hold on
   % ylim([-0.4 0.3])
    title(['SUB', num2str(i)])
    ylabel('\beta coefficient')
    xlabel('time point')
  %  legend('Diffusion', 'Shortest paths')
end


figure(2)
for i=21:23
    subplot(1,4,i-20)
    plot(beta_diff(i,:), 'r')
    hold on
    plot(beta_sp(i,:), 'b')
    hold on
    ylim([-0.4 0.3])
     title(['SUB', num2str(i)])
    ylabel('\beta coefficient')
    xlabel('time point')
    legend('Diffusion', 'Shortest paths')
end

subplot(1,4,4)
axis off
subplot(4,5,10)

load('/Users/emilyolafson/Documents/Thesis/SUB1_23_data/realizedrecovery.mat')
plot(basline, beta_diff(:,1), '*r')
[rho,p]=corr(basline, beta_diff(:,1))

%plot realtive ratio between betas
for i=1:23
    subtr(i,:)=abs(beta_diff(i,:))-abs(beta_sp(i,:))
    diff(i,:)=beta_diff(i,:)./beta_sp(i,:)

end
figure(1)
for i=1:20
    subplot(4,5,i)
     plot(subtr(i,:), 'b')
    hold on
    title(['SUB', num2str(i)])
    ylabel('\beta diff - \beta sp')
    xlabel('time point')

  %  ylim([-0.4 0.3])

  %  legend('Diffusion', 'Shortest paths')
  xlim([1 5])
end

for i=1:23
    a=cell2mat(ts_shen268{i,5}); %time x rois
    save(['/Users/emilyolafson/Documents/Thesis/SUB1_23_data/timeseriesSUB1_23/session5/SUB',num2str(i),'.mat'], 'a')
end


a=double(mean(subtr, 'omitnan'));
b=double(std(subtr, 'omitnan'));
boundedline([7 14 30 90 180],a,b, 'r')


figure(3)
a=double(mean(beta_diff, 'omitnan'));
b=double(std(beta_diff, 'omitnan'));
boundedline([7 14 30 90 180],a,b, 'r')
hold on;
a=double(mean(beta_sp, 'omitnan'));
b=double(std(beta_sp, 'omitnan'));
boundedline([7 14 30 90 180],a,b, 'b')
legend({'+/-  \sigma \beta_1 (diffusion (MFPT))',' \mu \beta_1 (diffusion (MFPT))', '+/- \sigma \beta_2 (shortest path length)','\mu \beta_2 (shortest path length)'})
set(gca,'XTick', [7 14 30 90 180], 'FontSize', 14)
xlabel('Days post stroke')
ylabel('\beta coefficient value')

b=std(beta_diff)

ratio=beta_diff./beta_sp
figure(4)
for i=1:20
    subplot(4,5,i)
    bar(ratio(i,:))
end

a=mean(ratio, 'omitnan')
b=std(ratio, 'omitnan')

figure(1)
boundedline(double([1:5]), double(a), double(b))


figure(1)
for i=1:23
    plot([1:5],beta_diff(i,:), '-r')
    hold on
    plot([1:5], beta_sp(i,:), '-b')
end

figure(2)
for i=1:20
    j=1
    subplot(4,5,i)
    fc_scores=cell2mat(fc_xfm(i,j));
    Y=fc_scores(indexes);
    X=mfpt_std_symm{i}(indexes);
    scatter(X,Y, 3,'red','filled', 'MarkerFaceAlpha', 0.05)
    title(['SUB', num2str(i)])
    ylabel('z-FC')
    xlabel('z-scored MFPT (diffusion)')
end

for i=21:23
    j=1
    subplot(1,3,i-20)
    fc_scores=cell2mat(fc_xfm(i,j));
    Y=fc_scores(indexes);
    X=mfpt_std_symm{i}(indexes);
    scatter(X,log(Y), 3,'red','filled', 'MarkerFaceAlpha', 0.05)
    title(['SUB', num2str(i)])
    ylabel('log(FC)')
    xlabel('z-scored MFPT (diffusion)')
    
end
figure(3)
for i=1:20
    j=1
    subplot(4,5,i)
    fc_scores=cell2mat(fc_xfm(i,j));
    Y=fc_scores(indexes);
    X=bc_std{i}(indexes);
    scatter(X,log(Y), 3,'filled', 'MarkerFaceAlpha', 0.05)
      title(['SUB', num2str(i)])
    ylabel('log(FC)')
    xlabel('z-scored shortest path length')
    
end
figure(4)
for i=21:23
    j=1
    subplot(1,3,i-20)
    fc_scores=cell2mat(fc_xfm(i,j));
    Y=fc_scores(indexes);
    X=bc_std{i}(indexes);
    scatter(X,log(Y), 3,'filled', 'MarkerFaceAlpha', 0.05)
      title(['SUB', num2str(i)])
    ylabel('log(FC)')
    xlabel('z-scored shortest path length')
   
end

bc_std{i}(indexes)


%% GLM with partial FC.


indexes=logical(triu(ones(268,268)));
e=eye(268);
indexes=indexes-e;
indexes=logical(indexes)

clear pval_diff
clear pval_sp
clear beta_diff
clear beta_sp
for i=1:23
    X=[mfpt_std_symm{i}(indexes),bc_std{i}(indexes)];
    for j=1:nsess(i)
        fc_scores=cell2mat(pFC_xfm(i,j));
        Y=fc_scores(indexes);
        Y=Y+1.0*10^(-8);
        [b, ~, stats{i,j}]=glmfit(X, Y, 'gamma','link', 'log', 'options', opts);
        pval_diff{i,j}=stats{i,j}.p(2);
        pval_sp{i,j}=stats{i,j}.p(3);

        beta_diff{i,j}=stats{i,j}.beta(2);
        beta_sp{i,j}=stats{i,j}.beta(3);
    end
end


beta_diff{6,5}=single(0)
beta_diff{12,5}=single(0)
beta_diff{12,4}=single(0)
beta_diff{20,5}=single(0)
beta_diff{20,4}=single(0)
beta_diff{20,3}=single(0)
beta_diff=cell2mat(beta_diff);

beta_sp{6,5}=single(0)
beta_sp{12,5}=single(0)
beta_sp{12,4}=single(0)
beta_sp{20,5}=single(0)
beta_sp{20,4}=single(0)
beta_sp{20,3}=single(0)
beta_sp=cell2mat(beta_sp);

pval_diff{6,5}=single(0)
pval_diff{12,5}=single(0)
pval_diff{12,4}=single(0)
pval_diff{20,5}=single(0)
pval_diff{20,4}=single(0)
pval_diff{20,3}=single(0)
pval_diff=cell2mat(pval_diff);

pval_sp{6,5}=single(0)
pval_sp{12,5}=single(0)
pval_sp{12,4}=single(0)
pval_sp{20,5}=single(0)
pval_sp{20,4}=single(0)
pval_sp{20,3}=single(0)
pval_sp=cell2mat(pval_sp);
pvals_FDR=fdr_bh([pval_diff,pval_sp], 0.05)

mean_beta_diff_S1=mean(beta_diff(:,1))
mean_beta_sp_S1=mean(beta_sp(:,1))

mean_beta_diff_S2=mean(beta_diff(:,1))
mean_beta_sp_S2=mean(beta_sp(:,1))

mean_beta_diff_S3=mean(beta_diff(:,1))
mean_beta_sp_S3=mean(beta_sp(:,1))

mean_beta_diff_S4=mean(beta_diff(:,1))
mean_beta_sp_S4=mean(beta_sp(:,1))

