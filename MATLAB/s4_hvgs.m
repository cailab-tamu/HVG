%%%%%%% DERMAL FIBROBLASTS %%%%%%%

% Data loading and filtering
datatag='smpl2_E-MTAB-5989_dermal_fibroblasts';
load(datatag,'X0','core_idx0','genelist','s_phate0g1','cellcycletag0');
X0g1=X0(:,cellcycletag0=="G1");
s=s_phate0g1;

% Core Selection
% figure;
% scatter3(s(:,1),s(:,2),s(:,3),10,'filled');
% hold on
% scatter3(s(core_idx0,1),s(core_idx0,2),s(core_idx0,3),10,'filled');
X=X0g1(:,core_idx0);

% Filtering genes with low expression
if strcmp(datatag,'smpl2_E-MTAB-5989_dermal_fibroblasts')
    [X,genelist]=sc_selectg(X,genelist);
else
    [X,genelist]=sc_selectg(X,genelist,1,4);
end

% Regression Plot
% figure;
% [Xnorm]=sc_norm(X,'type','deseq');
% T=sc_hvg(Xnorm,genelist,true,true,false,true);
% return;
% figure;
% loglog(T.u, T.residualcv2,'o')
% xlabel('Mean expression, log')
% ylabel('Residual CV^2, log')

% Filtering genes with FDR < 0.01 and FC > 0.15
T=T(~isnan(T.fitratio),:);
try
T=T(T.u>0.01,:);
catch
T=T(T.lgu>log(0.01),:);    
end
T2=T(T.fdr<0.01 & T.fitratio>0.15,:);

% Writing output
writetable(T2,'HVG_DF.txt','Delimiter','\t');
% run_enrichr(T2.genes)
% run_gorilla(T2.genes)

%%%%%%% GM12878 %%%%%%%

% Data loading and filtering
datatag='smpl0_GM12878_scRNA_seq_original';
load(datatag,'X0','core_idx0','genelist','s_phate0g1','cellcycletag0');
X0g1=X0(:,cellcycletag0=="G1");
s=s_phate0g1;

% Core Selection
% figure;
% scatter3(s(:,1),s(:,2),s(:,3),10,'filled');
% hold on
% scatter3(s(core_idx0,1),s(core_idx0,2),s(core_idx0,3),10,'filled');
X=X0g1(:,core_idx0);

% Filtering genes with low expression
if strcmp(datatag,'smpl2_E-MTAB-5989_dermal_fibroblasts')
    [X,genelist]=sc_selectg(X,genelist);
else
    [X,genelist]=sc_selectg(X,genelist,1,4);
end

% Regression Plot
% figure;
% [Xnorm]=sc_norm(X,'type','deseq');
% T=sc_hvg(Xnorm,genelist,true,true,false,true);
% return;
% figure;
% loglog(T.u, T.residualcv2,'o')
% xlabel('Mean expression, log')
% ylabel('Residual CV^2, log')

% Filtering genes with FDR < 0.01 and FC > 0.15
T=T(~isnan(T.fitratio),:);
try
T=T(T.u>0.01,:);
catch
T=T(T.lgu>log(0.01),:);    
end
T2=T(T.fdr<0.01 & T.fitratio>0.15,:);

% Writing output
writetable(T2,'HVG_GM12878.txt','Delimiter','\t');
% run_enrichr(T2.genes)
% run_gorilla(T2.genes)

%%%%%%% LUNG AIRWAY EPITHELIUM %%%%%%%

% Data loading and filtering
datatag='smpl1_GSM3204304_5_lung_airway_epithelial_cells';
load(datatag,'X0','core_idx0','genelist','s_phate0g1','cellcycletag0');
X0g1=X0(:,cellcycletag0=="G1");
s=s_phate0g1;

% Core Selection
% figure;
% scatter3(s(:,1),s(:,2),s(:,3),10,'filled');
% hold on
% scatter3(s(core_idx0,1),s(core_idx0,2),s(core_idx0,3),10,'filled');
X=X0g1(:,core_idx0);

% Filtering genes with low expression
if strcmp(datatag,'smpl2_E-MTAB-5989_dermal_fibroblasts')
    [X,genelist]=sc_selectg(X,genelist);
else
    [X,genelist]=sc_selectg(X,genelist,1,4);
end

% Regression Plot
% figure;
% [Xnorm]=sc_norm(X,'type','deseq');
% T=sc_hvg(Xnorm,genelist,true,true,false,true);
% return;
% figure;
% loglog(T.u, T.residualcv2,'o')
% xlabel('Mean expression, log')
% ylabel('Residual CV^2, log')

% Filtering genes with FDR < 0.01 and FC > 0.15
T=T(~isnan(T.fitratio),:);
try
T=T(T.u>0.01,:);
catch
T=T(T.lgu>log(0.01),:);    
end
T2=T(T.fdr<0.01 & T.fitratio>0.15,:);

% Writing output
writetable(T2,'HVG_LAEC.txt','Delimiter','\t');
% run_enrichr(T2.genes)
% run_gorilla(T2.genes)
