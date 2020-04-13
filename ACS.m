%%%ACS解决TSP问题%%%%%%%  
%close all; %清图  
%clc; %清屏  

N=10; %% 蚂蚁个数 
Alpha=0.1;%% Alpha 信息素衰减参数 
Beta=2; %% Beta 信息素和距离的权重控制
q0 = 0.9; %% 状态选择概率
Rho = 0.1;%% 边缘探测参数
Generation=500;%%最大迭代次数 
DIM=15;%% 城市个数
filename='./data/p01.txt';% 选择文件

global len;
len=zeros(DIM,DIM);%% 相邻结点之间的距离
global position;
position=zeros(DIM,2); %%城市坐标
global ant;
ant=zeros(N,DIM);%% 蚂蚁的行走路径
global pheromone;
pheromone=zeros(DIM,DIM);%% 路径信息素积累量
global bestpath;
bestpath=zeros(1,DIM);%% 最优路径
global least_cost;
least_cost=intmax;%% 最小开销

test(N,Alpha,Beta,q0,Rho,Generation,DIM,filename);
