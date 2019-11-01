% rmpath('scgeapp/');
% datatag='smpl0_GM12878_scRNA_seq_original';
% datatag='smpl1_GSM3204304_5_lung_airway_epithelial_cells';
datatag='smpl2_E-MTAB-5989_dermal_fibroblasts';


load(datatag,'X0','core_idx0','genelist','s_phate0g1','cellcycletag0');
X0g1=X0(:,cellcycletag0=="G1");
s=s_phate0g1;

figure;
scatter3(s(:,1),s(:,2),s(:,3),10,'filled');
hold on
scatter3(s(core_idx0,1),s(core_idx0,2),s(core_idx0,3),10,'filled');

X=X0g1(:,core_idx0);
% X=X0g1;

% i=randperm(size(X0,2));
% X0=X0(:,i);
% X=X0(:,1:1000);
if strcmp(datatag,'smpl2_E-MTAB-5989_dermal_fibroblasts')
    [X,genelist]=sc_selectg(X,genelist);
else
    [X,genelist]=sc_selectg(X,genelist,1,4);
end

figure;
[Xnorm]=sc_norm(X,'type','deseq');
T=sc_hvg(Xnorm,genelist,true,true,false,true);

return;


figure;
loglog(T.u, T.residualcv2,'o')
xlabel('Mean expression, log')
ylabel('Residual CV^2, log')


% T=sc_veg(Xnorm,genelist,true,true,false);
%grid on
%box on
%[Xnorm]=sc_norm(X,'type','libsize');
% T=sc_veg(Xnorm,genelist,true,true);

% T=sc_splinefit(Xnorm,genelist,true,true);
T=T(~isnan(T.fitratio),:);
try
T=T(T.u>0.01,:);
catch
T=T(T.lgu>log(0.01),:);    
end
T2=T(T.fdr<0.01 & T.fitratio>1.5,:);
T3=T(T.fdr==0,:);

%writetable(T2,'smpl0_topXXXgene.txt','Delimiter','\t');
%writetable(T3,'smpl0_topXXgene.txt','Delimiter','\t');
run_enrichr(T2.genes)

% 
% % run_gorilla(T2.genes)
% run_gorilla(T.genes(1:5000))

