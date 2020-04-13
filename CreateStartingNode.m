function [y]=CreateStartingNode(N,DIM)
global ant;
%-------------------------------------------------------------------------%
temp=randperm(DIM,N); % 返回长度为N的向量，范围为1-DIM
for i = 1:N %给每只蚂蚁出发点
    ant(i,1)=temp(i);
end
end
