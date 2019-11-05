s4_hvgs;

%%
Xmagic=run_magic(X);
%%
topg=T.genes(1:50);
%  i=startsWith(topg,'IG');
%  topg(i)=[];

% topg=T3.genes;
[~,idx]=ismember(topg,genelist);

A=Xmagic(idx,:);
B=corr(A');
C=B>0.65;
C=C-diag(diag(C));
% C=logical(C);
% B(~C)=0;

i=sum(C)==0;
C(i,:)=[]; C(:,i)=[];
g2=topg; g2(i)=[];
disp('ok')
%
% [x,y,v]=find(triu(B));
% for k=1:length(x)
%     
% end

% bg1 = biograph(C>0,g2);
% figure;
% view(bg1);


addpath('C:\Users\jcai\Documents\GitHub\SBEToolbox_lite');
addpath('C:\Users\jcai\Documents\GitHub\SBEToolbox_lite\bgl');

%
%kamada_kawai_spring_layout
%gursoy_atun_layout
%circle_graph_layout
%fruchterman_reingold_force_directed_layout
    xy=fruchterman_reingold_force_directed_layout(double(sparse(C)),...
        'iterations',1000);
    figure;
    plotnet(C,xy,g2)

    try
    xy=kamada_kawai_spring_layout(double(sparse(C)),...
        'iterations',1000);
    figure;
    plotnet(C,xy,g2)
    catch
    end
    try
    xy=gursoy_atun_layout(double(sparse(C)),...
        'iterations',1000);
    figure;
    plotnet(C,xy,g2)
    catch
    end
    

   dgree=sum(C,2);
   btwness=betweenness_centrality(double(sparse(C)));
   T4=table(g2,dgree,btwness);
   T4 = sortrows(T4,'dgree','descend');
   
   % viewnetprotovis(C,g2)
   
writesbe2sif(C,g2,'forDF.sif');
% system('java -Xmx512M -jar cytoscape.jar --network input.sif -p plugins &');