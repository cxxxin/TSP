// TSP-ACS.cpp

#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <ctime>
#include <stdio.h>
#include <Windows.h>
#include <cmath>
#include <algorithm>
#include <vector>
#include <map>
using namespace std;

//参数设置
filename='gr17.txt';//文件名
constexpr auto b = 2;//信息素和距离的权重控制
constexpr auto q0 = 0.9;//
constexpr auto a = 0.1;//信息素衰减参数
constexpr auto ρ = 0.1;//边缘探测参数
double r0 = 0;
constexpr auto DIM = 17;//城市个数
constexpr auto GENERATION = 1000; //迭代代数
constexpr auto N = 10;//蚂蚁数量

/***/
int ant[N][DIM] = { 0 };//蚂蚁的行走路径
double pheromone[DIM][DIM] = {0.0};//路径信息素积累量
int len[DIM][DIM];//相邻结点之间的距离
int bestpath[DIM] = { 0 };//最优路径
vector<int> unvisited[N];//对于每只蚂蚁都建立未访问城市集合
double f[N] = { 0 };//路径总长度
double least_cost = MAXINT;//最小开销
multimap<int,int> neighbour[DIM];//最近邻居
/**/

/**step_func*/
void test();
void init();
void cal_best();
/**/

/**辅助函数*/
void CreateStartingNode(int min=0, int max=DIM-1, int num=N)//生成不重复随机数
{
	vector<int> res;
	res.clear();
	if (max - min + 1 < num)
	{
		return;
	}
	srand(time(0));
	for (auto i{ 0 }; i < num; i++)
	{
		while (true)
		{
			auto temp{ rand() % (max + 1 - min) + min };
			auto iter{ find(res.begin(),res.end(),temp) };
			if (res.end() == iter)
			{
				res.push_back(temp);
				break;
			}
		}
	}
	vector<int> ::iterator iter;
	for (int j = 0; j < res.size(); j++)
	{
		ant[j][0] = res[j];//赋值给蚂蚁的起点
		for (iter = unvisited[j].begin(); iter != unvisited[j].end(); iter++)//更新未访问城市集合
		{
			if (*iter == ant[j][0])
			{
				unvisited[j].erase(iter);
				break;
			}
		}
	}
}
void nearestNeighbour()
{
	for (int i = 0; i < DIM; i++)
	{
		neighbour[i].clear();
		for (int j = 0; j < DIM; j++)
		{
			if(j!=i)
				neighbour[i].insert(pair<int, int>(len[i][j],j));//路长-目标城市
		}
	}
}
/**/

int main()
{
	test();
	system("pause");
	return 0;
}

void test()
{
	LARGE_INTEGER seed;
	QueryPerformanceFrequency(&seed);
	QueryPerformanceCounter(&seed);
	srand(seed.QuadPart);//初始化一个以微秒为单位的时间种子

	init();//地图
	for (int i = 0; i < GENERATION; i++)//circle0 迭代
	{	
		for (int j = 0; j < N; j++)//更新未访问集合
		{
			unvisited[j].clear();//清空
			for (int k = 0; k < DIM; k++)
				unvisited[j].insert(unvisited[j].end(),k);
		}
		CreateStartingNode();//每只蚂蚁选择一个初始位置（各不相同
		double LNN[N] = { 0 };//最邻近启发式得到的路长
		for (int j = 0; j < N; j++)
		{
			nearestNeighbour();
			vector<int> visited;
			visited.insert(visited.end(), ant[j][0]);//将起点加入访问城市的列表
			for (int l = 0; visited.size()!=DIM; l++)
			{
				multimap<int, int>::iterator iter;
				for (int k = 0; k < DIM; k++)//将该城市从所有邻居中删除
				{
					for (iter = neighbour[k].begin(); iter != neighbour[k].end(); iter++)
					{
						if (iter->second == visited[l])
						{
							neighbour[k].erase(iter);
							break;
						}
					}
				}
				iter = neighbour[visited[l]].begin();//选择下一个城市
				visited.insert(visited.end(), iter->second);//加入访问城市列表
			}
			for (int k = 0; k < DIM; k++)//计算路长
			{
				LNN[j] += len[visited[k]][visited[(k+1)%DIM]];
			}
			visited.clear();
		}
		
		//开始行走
		for (int j = 1; j < DIM; j++)//circle1 每一步
		{
			for (int k = 0; k < N; k++)//每只蚂蚁
			{
				double q = rand() / double(RAND_MAX);//判断下一步的状态
				if (q <= q0)	//开发
				{
					double max_phe = 0.0;//最多的信息素
					int city=-1;//被选中的城市
					vector<int>::iterator iter;
					vector<int> equ;
					for (iter = unvisited[k].begin(); iter != unvisited[k].end(); iter++)//遍历未访问列表，寻找最大值
					{
						if (pheromone[ant[k][j-1]][*iter] / pow(len[ant[k][j - 1]][*iter], b) > max_phe)
						{
							equ.clear();
							equ.insert(equ.end(), *iter);
							max_phe = pheromone[ant[k][j - 1]][*iter] / pow(len[ant[k][j - 1]][*iter], b);
							city = *iter;
						}
						else if (pheromone[ant[k][j - 1]][*iter] / pow(len[ant[k][j - 1]][*iter], b) == max_phe)
						{
							equ.insert(equ.end(), *iter); //一样都是最大的列表
						}
					}
					if (equ.size() >= 1) //不止一个最大选项
					{
						srand((unsigned)time(NULL));
						city = equ[rand()%equ.size()];
					}
					equ.clear();
					ant[k][j] = city;
					for (iter = unvisited[k].begin(); iter != unvisited[k].end(); iter++)//选定城市，从未访问列表中删除
					{
						if (*iter == city)
						{
							unvisited[k].erase(iter);
							break;
						}
					}
				}
				else	//探索
				{
					vector<int>::iterator iter;
					double *pk = new double[unvisited[k].size()];
					int l = 0;
					double sum = 0;
					for (iter = unvisited[k].begin(); iter != unvisited[k].end(); iter++,l++)//遍历未访问城市列表，计算概率区域
					{
						pk[l] = pheromone[ant[k][j - 1]][*iter] / pow(len[ant[k][j - 1]][*iter], b);
						sum += pk[l];
					}
					if (sum != 0) //总和不为0
					{
						pk[0] = pk[0] / sum;
						for (l = 1; l < unvisited[k].size(); l++)
						{
							pk[l] = pk[l] / sum;
							pk[l] += pk[l-1];
						}
					}
					else //总和为0 平均概率
					{
						for (l = 0; l < unvisited[k].size(); l++)
						{
							pk[l] = (1.0 / unvisited[k].size())*(l + 1);
						}
					}
					double p = rand() / double(RAND_MAX);//根据概率所属区间选择城市
					for (l = 1; l < unvisited[k].size(); l++)
					{
						if (p <= *(pk + l) && p > *(pk + l-1))
						{
							ant[k][j] = unvisited[k][l];
							break;
						}
					}
					if (p <= *(pk + 0))
						ant[k][j] = unvisited[k][0];
										
					for (iter = unvisited[k].begin(); iter != unvisited[k].end(); iter++)//从未访问列表中删除
					{
						if (*iter == ant[k][j])
						{
							unvisited[k].erase(iter);
							break;
						}
					}
					delete []pk;
					pk = NULL;
				}
				
				r0 =  1 / (DIM * LNN[k]);//计算r0
				pheromone[ant[k][j - 1]][ant[k][j]] = (1 - ρ)*pheromone[ant[k][j - 1]][ant[k][j]] + ρ * r0;//局部更新
				pheromone[ant[k][j]][ant[k][j - 1]] = (1 - ρ)*pheromone[ant[k][j]][ant[k][j - 1]] + ρ * r0;
				if (j == DIM - 1)//回到原点
				{
					pheromone[ant[k][j]][0] = (1 - ρ)*pheromone[ant[k][j]][0] + ρ * r0;//局部更新
					pheromone[0][ant[k][j]] = (1 - ρ)*pheromone[0][ant[k][j]] + ρ * r0;
				}
			}
		}
		
		cal_best();//找全局最佳蚂蚁
		for (int j = 0; j < DIM; j++)
		{
			for (int k = 0; k < DIM ; k++)
			{
				pheromone[j][k] = (1 - a)*pheromone[j][k];
				pheromone[k][j] = (1 - a)*pheromone[k][j];
			}
		}
		for (int j = 0; j < DIM; j++)//信息素全局更新
		{
			pheromone[bestpath[j]][bestpath[(j + 1) % DIM]] += a / least_cost;
			pheromone[bestpath[(j + 1) % DIM]][bestpath[j]] += a / least_cost;
		}
	}

	cout << "最短路长:" << least_cost << endl;//输出结果
	cout << "最短路径:";
	for (int i = 0; i < DIM; i++)
		cout << bestpath[i] << " ";
	cout << endl;
}

void init()
{
	//初始化相邻路径长度
	FILE *fp;//文件指针
	fopen_s(&fp, filename, "r");//以文本方式打开文件
	if (fp == NULL) //打开文件出错。
		return;
	for (int i = 0; i < DIM; i++)
	{
		for (int j = 0; j < DIM; j++)
		{
			if (!feof(fp))
				fscanf_s(fp, "%d", &len[i][j]); //读取数据到数组，直到文件结尾(返回EOF)
		}
	}
	fclose(fp);//关闭文件
}

void cal_best()
{
	int shortest_index = -1;
	//计算路径
	for (int i = 0; i < N; i++)
	{
		f[i] = 0;
		for (int j = 0; j < DIM; j++)
		{
			f[i] += len[ant[i][j]][ant[i][(j + 1) % DIM]];
		}
		if (least_cost > f[i])
		{
			shortest_index = i;
			least_cost=f[i];
		}
	}
	if (shortest_index != -1)//更新最优值
	{
		for (int i = 0; i < DIM; i++)
			bestpath[i] = ant[shortest_index][i];
	}
}
