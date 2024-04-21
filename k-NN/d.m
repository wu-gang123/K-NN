clear
close all
fileID = fopen('NEpsilon.txt','w');% 创建文本
%%
for m = 1:25
    % 设置参数
    N = 50 * m; % 节点数量
    NN(m) = N;
    M = 1; % 图信号类型
    rng(1)
    fea = rand( N, M); % 随机图信号
    coords = rand(N, 2);% 生成随机节点坐标
    %%
    % 计算节点之间的距离
    dis = zeros(N, N);
    for i = 1:N
        for j = 1:N
            dis(i, j) = norm(coords(i, :) - coords(j, :));
        end
    end
    % 设置比率delta
    delta = 0.1:0.05:0.9;
    %%
    % 构建最近邻图
    for n = 1:length(delta)
        Epsilon(n) = delta(n) * (max(dis(:)) - min(dis(:)));  % 设置阈值Epsilon
        ddis = dis;
        for i = 1:N
            for j = 1:N
                if ddis(i,j) > Epsilon(n)
                    ddis(i,j) = 0;
                end
            end
        end
        W = ddis;
        G = gsp_graph( W, coords ); %get graph
        L = G.L;   %拉普拉斯矩阵
        L_n = diag(diag(L).^(-0.5))*L*diag(diag(L).^(-0.5));%对称归一化拉普拉斯矩阵
        [eigenvectors,eigenvalues,~] = svd(full(L_n));
        [E, inds] = sort(diag(eigenvalues),'ascend');  %%将特征值以升序重新排列
        eigenvectors=eigenvectors(:,inds');  %%特征向量和排列后的特征值位置对应
        %
        U = eigenvectors;  %eigenvectors matricx
        lambda = E;
        sm(n)= fea' * L_n * fea; % 平滑性
    end
    %% 平滑性最优时
    [minsm,posi] = min(sm);
    rrate(m) = delta(posi)*100;
    rate =  delta(posi)*100;
    fprintf(fileID,'N=%d, rate=%.2f%% \n',N,rate);
end
fclose(fileID);
%% 图示
% Epsilon占距离差之比
figure(1)
plot(NN,rrate,'-*',LineWidth',2);
ylim([0 100])
xlabel('the number of vertics');
ylabel('rate(%)');
title('chose rate for N ');
% 图信号的特征值
figure(2)
stem(1:size(lambda),lambda,'r');
xlabel('the order of lambda');
ylabel('the value of lambda');
% Epsilon的选择对图信号光滑性的影响
figure(3)
plot(Epsilon,sm,'-o');
xlabel('Epsilon');
ylabel('the index of smooth');
title('Epsilon的选择对图信号平滑性的影响');



