% allref = 420x268x268 single 

avg_connectome=squeeze(mean(allref,1));
diagonal=eye(268);
imagesc(diagonal)
imagesc(avg_connectome.*~diagonal)
kmeans