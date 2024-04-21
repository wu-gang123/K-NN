clear
close all
fileID = fopen('Nk.txt', 'w'); %
%% 根据k-NN生成图G
for m = 1:25
    N = 50 * m;  %节点数
    NN(m) = N;
    M = 1;   %图信号类型数
    rng(1)   %保持rand随机分布不变
    fea = rand(N,M); % 图信号
    coords = rand( N, 2 ); %节点坐标
    arrange = 0.2*N;% k的取值范围
    for i = 1:arrange
        options = [];
        options.NeighborMode = 'KNN';
        options.k = i;
        options.WeightMode = 'HeatKernel';
        options.t = 1;
        W = constructW( fea, options );

        G = gsp_graph( W, coords );   %得到图G
        %%
        L = G.L;   %拉普拉斯矩阵
        L_n = diag(diag(L).^(-0.5))*L*diag(diag(L).^(-0.5));%对称归一化拉普拉斯矩阵
        [eigenvectors,eigenvalues,~] = svd(full(L_n));
        [E, inds] = sort(diag(eigenvalues),'ascend');  %%将特征值以升序重新排列
        % eigenvectors=eigenvectors(:,inds');  %%特征向量和排列后的特征值位置对应
        % %
        % U = eigenvectors;  %eigenvectors matricx
        lambda = E;
        sm(i) = fea' * L_n * fea;
        %sm(i)=trae(fea'*L_n*fea)
    end
    %% 最小平滑性k的取值
    %将输出写入空白文档
    [minsm,posi] = min(sm);
    rrate(m) = (posi/N)*100;
    rate = (posi/N)*100;
    fprintf(fileID,'N=%d, k=%d, rate=%.2f%% \n',N,posi,rate);

end
fclose(fileID);
%% 图示
% k占N的比率
figure(1)
plot(NN,rrate,'-*')
xlim([0 50*(m+1)])
ylim([0 20])
xlabel('the number of vertics')
ylabel('the rate(%)')
title('rate for different N')
% 拉普拉斯矩阵的特征值
figure(2)
stem(1:size(lambda),lambda,'r');
xlabel('the order of lambda');
ylabel('the value of lambda');

% k的值对平滑性的影响
figure(3)
plot(sm,'-o');
xlabel('k');
ylabel('the index of smooth');
title('The influence of the choice of k on the smoothness of graph signals');
