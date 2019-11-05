% s4_hvgs;
% addpath('scGEApp\src');
datatag='smpl0_GM12878_scRNA_seq_original';

% datatag='smpl1_GSM3204304_5_lung_airway_epithelial_cells';
% datatag='smpl2_E-MTAB-5989_dermal_fibroblasts';

load(datatag,'X0','core_idx0','genelist','s_phate0g1','cellcycletag0');
X=X0(:,cellcycletag0=="G1"); withcycle=false;
% X=X0(:,core_idx0); withcycle=false;
% X=X0; withcycle=true;
%%
X=run_magic(X,true);

%%
figure; i_myscatter(s_phate0g1,X(genelist=="IRF4",:));
% figure; i_myscatter(s_phate0g1,X(genelist=="PRDM1",:));
%             xlabel 'PHATE1'
%             ylabel 'PHATE2'
%             zlabel 'PHATE3'
%             title 'PHATE 3D'
view(27,13)

%%
i1=find(genelist=="AICDA");
i2=find(genelist=="BACH2");
i3=find(genelist=="BCL6");
i4=find(genelist=="IRF4");
i5=find(genelist=="PAX5");
i6=find(genelist=="PRDM1");
i7=find(genelist=="REL");
i8=find(genelist=="RELA");

x1=X(i1,:); x2=X(i2,:); x3=X(i3,:); x4=X(i4,:);
x5=X(i5,:); x6=X(i6,:); x7=X(i7,:); x8=X(i8,:);

% figure;
% subplot(2,2,1), stem(x1,'marker','none')
% subplot(2,2,2), stem(x2,'marker','none')
% subplot(2,2,3), stem(x3,'marker','none')
% subplot(2,2,4), stem(x4,'marker','none')
% 
% figure;
% subplot(2,2,1), stem(x5,'marker','none')
% subplot(2,2,2), stem(x6,'marker','none')
% subplot(2,2,3), stem(x7,'marker','none')
% subplot(2,2,4), stem(x8,'marker','none')

DataTable=table(x1',x2',x3',x4',x5',x6',x7',x8');
DataTable.Properties.VariableNames={'AICDA','BACH2','BCL6','IRF4','PAX5','PRDM1','REL','RELA'};

% DataTable=table(x1',x4',x5',x6',x7',x8');
% DataTable.Properties.VariableNames={'AICDA','IRF4','PAX5','PRDM1','REL','RELA'};
% corrplot(DataTable)
figure; corrplot(DataTable,'testR','on','alpha',0.01./nchoosek(8,2))
%%

figure; scatter(x7,x6,'.'); xlabel('REL (cRel)'); ylabel('PRDM1 (Blimp-1)');

if withcycle
    figure; 
    gscatter(x6,x7,cellcycletag0,'',''); xlabel('PRDM1 (Blimp-1)'); ylabel('REL (cRel)');
end
% figure; 
% hold on
% i=cellcycletag0=='G1';
% scatter(x6(i),x7(i)); 
% i=cellcycletag0=='S';
% scatter(x6(i),x7(i)); 
% i=cellcycletag0=='G2/M';
% scatter(x6(i),x7(i)); xlabel('REL (cRel)'); ylabel('PRDM1 (Blimp-1)');
% legend({'G1','S','G2/M'})

figure; 
scatter(x4,x1,15,x6,'filled'); xlabel('IRF4'); ylabel('AICDA (AID)'); title('PRDM1 (Blimp-1)'); colorbar;
box on

%figure; scatter(x7,x1,'.'); xlabel('REL (cRel)'); ylabel('AICDA (AID)');
%figure; scatter(x6,x1,'.'); xlabel('PRDM1 (Blimp-1)'); ylabel('AICDA (AID)');


%%
if withcycle
    figure;
    subplot(2,2,1)
    scatter(x4,x1,20,x6,'filled'); xlabel('IRF4'); ylabel('AICDA (AID)'); title('PRDM1 (Blimp-1) All');
    colorbar
    subplot(2,2,2)
    i=cellcycletag0=="G1";
    scatter(x4(i),x1(i),20,x6(i),'filled'); xlabel('IRF4'); ylabel('AICDA (AID)'); title('PRDM1 (Blimp-1) G1');
    colorbar
    subplot(2,2,3)
    i=cellcycletag0=="S";
    scatter(x4(i),x1(i),20,x6(i),'filled'); xlabel('IRF4'); ylabel('AICDA (AID)'); title('PRDM1 (Blimp-1) S');
    colorbar
    subplot(2,2,4)
    i=cellcycletag0=="G2M";
    scatter(x4(i),x1(i),20,x6(i),'filled'); xlabel('IRF4'); ylabel('AICDA (AID)'); title('PRDM1 (Blimp-1) G2/M');
    colorbar
end
% figure; scatter3(x4,x1,x6,'filled'); xlabel('IRF4'); ylabel('AICDA (AID)'); zlabel('PRDM1 (Blimp-1)')
% figure; gscatter(x6,x7,cellcycletag0); xlabel('PRDM1 (Blimp)'); ylabel('REL (cRel)'); 
% figure; gscatter(x7,x1,cellcycletag0); xlabel('REL (cRel)'); ylabel('AICDA (AID)');
% figure; gscatter(x6,x1,cellcycletag0); xlabel('PRDM1 (Blimp-1)'); ylabel('AICDA (AID)');
% figure; scatter3(x6,x1,x4,'.'); xlabel('PRDM1 (Blimp-1)'); ylabel('AICDA (AID)'); zlabel('IRF4');
%i=cellcycletag0=="G1";
%figure; scatter3(x1(i),x3(i),x4(i),'.')

