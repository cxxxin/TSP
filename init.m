function [y]=init(filename,DIM)
global len;
global position;
%-------------------------------------------------------------------------%
% 开文件
data=load(filename);
x=data(:,1);
y=data(:,2);

for i = 1:DIM % 城市坐标赋值
    position(i,1)=x(i); % x坐标
    position(i,2)=y(i); % y坐标
end

% 计算距离
for i = 1:DIM
    for j = i:DIM
        len(i,j)=sqrt((x(i)-x(j))^2+(y(i)-y(j))^2);
        len(i,j)=round(len(i,j));
        if i~=j
            len(j,i)=len(i,j);
        end
    end
end
end

