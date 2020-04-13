function [y]=test(N,Alpha,Beta,q0,Rho,Generation,DIM,filename)
global ant;
global pheromone;
global len; %% 相邻结点之间的距离
global bestpath;
global least_cost;
global position; %%城市坐标
tao=0.0;
unvisited=cell(N,1);% 建立未访问城市列表
nei_dis=zeros(DIM,DIM); % 邻居距离列表
nei_index=zeros(DIM,DIM); % 邻居的索引
%-------------------------------------------------------------------------%
init(filename,DIM);%% 生成地图

%可视化地图
figure(1)
plot(position(:,1),position(:,2),'.','MarkerSize',15,'color','r');
title('城市地图');
hold off


for i = 1:Generation% 迭代
    for j = 1:N % 更新未访问集合
        unvisited{j}=colon(1,DIM); % 所有城市
        %unvisited=[unvisited;temp];
    end
    
    CreateStartingNode(N,DIM); %每只蚂蚁选择一个初始位置
    for j = 1:N     %在未访问列表中去掉初始位置
        tmp=unvisited{j};
        temp=find(tmp==ant(j,1));
        tmp(temp)=[];
        unvisited{j}=[];
        unvisited{j}=tmp;
    end
    
    LNN=zeros(1,N); %最邻近启发式得到的路长
    for j = 1:DIM % 建立最近邻居列表
        nei_dis(j,:)=len(j,:);
        nei_dis(j,:)=sort(nei_dis(j,:));
        y=[];
        tmp=len(j,:);
        for k = 1:DIM
            if length(y)==0
                y=find(tmp==nei_dis(j,k));
            end
            if length(y)==1
                nei_index(j,k)=y;
                y=[];
            else
                nei_index(j,k)=y(1);
                y(1)=[];
            end
        end
    end
    for j = 1:N %开始添加路径
        visited=[]; % 访问列表
        visited=[visited,ant(j,1)]; % 将起点加入访问列表
        nei_tmp=cell(DIM,1);
        for k = 1:DIM % 复制邻居索引集
            nei_tmp{k}=nei_index(k,:); 
        end
        for l = 1:DIM-1
            for k = 1:DIM % 在邻居集删除当前城市
                tmp=nei_tmp{k};
                tmp(find(tmp==visited(l)))=[];
                nei_tmp{k}=tmp;
            end
            visited=[visited,nei_tmp{visited(l)}(1)]; % 当前城市的第一近邻居城市
        end
        for k = 1:DIM % 计算路长
            if k~=DIM
                LNN(j)=LNN(j)+len(visited(k),visited(k+1));
            else
                LNN(j)=LNN(j)+len(visited(k),visited(1));
            end            
        end
    end
    
    %开始行走
    for j = 2:DIM %行走的每一步
        for k = 1:N %每只蚂蚁
            rng('shuffle');
            q=rand(1,1); % 产生随机数q判断下一步状态
            if q<=q0    % 开发状态
                city=-1;% 被选中的城市
                tmp=[];
                for iter = 1:length(unvisited{k})  %遍历未访问列表
                    tmp=[tmp,pheromone(ant(k,j-1),unvisited{k}(iter))/power(len(ant(k,j-1),unvisited{k}(iter)),Beta)];
                end
                max_index=find(tmp==max(tmp)); % 找到最大值
                rng('shuffle');
                city=round(1+rand(1,1)*(length(max_index)-1));
                city=max_index(city);
                city=unvisited{k}(city);
                
                ant(k,j)=city; %选中目标城市
                
                tmp=[];
                tmp=unvisited{k}; %从未访问列表中移除
                temp=find(tmp==city);
                tmp(temp)=[];
                unvisited{k}=tmp;
                
            %-------------------------------------------------------------------%
            else  % 探索状态
                pk=zeros(1,length(unvisited{k})); %轮盘赌占比
                for iter = 1:length(unvisited{k})
                    pk(iter)=pheromone(ant(k,j-1),unvisited{k}(iter))/power(len(ant(k,j-1),unvisited{k}(iter)),Beta);
                end
                SUM=sum(pk); % 总和
                if SUM~=0
                    pk(1)=pk(1)/SUM;
                    for iter = 2:length(unvisited{k})  % 累积概率
                        pk(iter)=pk(iter)/SUM;
                        pk(iter)=pk(iter)+pk(iter-1);
                    end
                else % 全部为0  直接随机选取
                    pk=colon(1.0/length(unvisited{k}),1.0,length(unvisited{k}));
                end
                    
                rng('shuffle');
                p=rand(1,1); % 根据概率选择所属区间选择目标城市
                for iter = 2:length(unvisited{k})
                    if p<=pk(iter) && p>pk(iter-1)
                        ant(k,j)=unvisited{k}(iter); % 选定目标
                        break;
                    end
                end
                if p<=pk(1)
                    ant(k,j)=unvisited{k}(1); % 选定目标
                end
                    
                tmp=unvisited{k}; %从未访问列表中移除
                temp=find(tmp==ant(k,j));
                tmp(temp)=[];
                unvisited{k}=tmp;
            end
            %------------------------------------------------------------------------%
            
            % 本地更新信息素
            tao=power(DIM*LNN(k),-1); %计算tao
            pheromone(ant(k,j-1),ant(k,j))=(1-Rho)*pheromone(ant(k,j-1),ant(k,j))+Rho*tao;
            pheromone(ant(k,j),ant(k,j-1))=pheromone(ant(k,j-1),ant(k,j));
            if j==DIM %回到原点的本地更新
                pheromone(ant(k,j),1)=(1-Rho)*pheromone(ant(k,j),1)+Rho*tao;
                pheromone(1,ant(k,j))=pheromone(ant(k,j),1);
            end
        end
    end
    
    shortest_index=-1;% 找全局最佳蚂蚁
    %计算路径
    for j = 1:N
       SUM=0;
       for k = 1:DIM
           if k~=DIM
               SUM=SUM+len(ant(j,k),ant(j,k+1));
           else
               SUM=SUM+len(ant(j,k),ant(j,1));
           end           
       end
       if least_cost>SUM
           shortest_index=j;
           least_cost=SUM;
       end
    end
    if shortest_index~=-1 %更新最佳路径
        bestpath=ant(shortest_index,:);
    end
    
    % 画图
    %close all;
    figure(1)
    plot(1,1);
    plot(position(:,1),position(:,2),'.','MarkerSize',15,'color','r');
    text(position(bestpath(1),1),position(bestpath(1),2),'*','color','g','FontSize',40);
    text(position(bestpath(1),1),position(bestpath(1),2),'起点','color','k');
    for j = 1:DIM
        figure(1)
        if j~=DIM
            line([position(bestpath(j),1),position(bestpath(j+1),1)],[position(bestpath(j),2),position(bestpath(j+1),2)],'linestyle','-','color','k');
            title('城市地图');
        else
            line([position(bestpath(j),1),position(bestpath(1),1)],[position(bestpath(j),2),position(bestpath(1),2)],'linestyle','-','color','k');
            title('城市地图');
        end
    end
    hold off
    
    for j = 1:DIM %信息素全局更新
        for k = 1:DIM
            pheromone(j,k)=(1-Alpha)*pheromone(j,k);
            pheromone(k,j)=pheromone(j,k);
        end
    end
    for j = 1:DIM
        if j~=DIM
            pheromone(bestpath(j),bestpath(j+1)) = pheromone(bestpath(j),bestpath(j+1))+Alpha/least_cost;
            pheromone(bestpath(j+1),bestpath(j)) = pheromone(bestpath(j),bestpath(j+1));
        else
            pheromone(bestpath(j),bestpath(1)) = pheromone(bestpath(j),bestpath(1))+Alpha/least_cost;
            pheromone(bestpath(1),bestpath(j)) = pheromone(bestpath(j),bestpath(1));
        end
    end
    fprintf('第%d代\n',i);
    disp(bestpath);
    fprintf('花销为%d\n',least_cost);
    
    all_cost(i)=least_cost;
end
end
