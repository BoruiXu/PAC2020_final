#include<cstdio>
#include<iostream>
#include<map>
#include<vector>
#include<algorithm>
#include<queue>
using namespace std;
const int X=5,Y=5,Z=5;
int main(){
    map<int,int> d_2={{0,1},{1,2},{3,3},{2,1}};
    map<int,int> d_3={{0,1},{1,2},{3,3},{2,4},{6,5},{7,6},{5,7},{4,8}};
    pair<int,int> task_2d[X][Y];
    pair<pair<int,int>,int> task_3d[X][Y][Y];
    cout<<"starting encode gray"<<endl;
    //starting encode gray
    task_2d[0][0].first=0;task_2d[0][0].second=0;
    task_3d[0][0][0].first.first=0;task_3d[0][0][0].first.second=0;task_3d[0][0][0].second=0;
    for(int i=1;i<X;i++)//raw first
    {
        task_2d[i][0].second=1-task_2d[i-1][0].second;
        task_2d[i][0].first=task_2d[i-1][0].first;
        task_3d[i][0][0].first.first=1-task_3d[i-1][0][0].first.first;
        task_3d[i][0][0].first.second=task_3d[i-1][0][0].first.second;
        task_3d[i][0][0].second=task_3d[i-1][0][0].second;
    }

    for(int i=0;i<X;i++){
        for(int j=1;j<Y;j++){
            task_2d[i][j].second=task_2d[i][j-1].second;
            task_2d[i][j].first=1-task_2d[i][j-1].first;
            task_3d[i][j][0].first.first=task_3d[i][j-1][0].first.first;
            task_3d[i][j][0].first.second=1-task_3d[i][j-1][0].first.second;
            task_3d[i][j][0].second=task_3d[i][j-1][0].second;
        }
    }

    for(int i=0;i<X;i++)
    {
        for(int j=0;j<Y;j++)
        {
            for(int k=1;k<Z;k++)
            {
                task_3d[i][j][k].first.first=task_3d[i][j][k-1].first.first;
                task_3d[i][j][k].first.second=task_3d[i][j][k-1].first.second;
                task_3d[i][j][k].second=1-task_3d[i][j][k-1].second;
            }
        }
    }
    cout<<"finished encode gray"<<endl;
    //finished encode gray

    //starting get TDG and it contains two parts
    //first get Adjacency Matrix and indegree
    cout<<"first part:"<<endl;
    vector<pair<int,int>>adj_2d[X][Y];
    int indegree_2d[X][Y];
    vector<pair<pair<int,int>,int>>adj_3d[X][Y][Z];
    int indegree_3d[X][Y][Z];
    fill(indegree_2d[0],indegree_2d[0]+X*Y,0);
    fill(indegree_3d[0][0],indegree_3d[0][0]+X*Y*Z,0);
    cout<<"start 2d:"<<endl;
    //for 2d
    for(int i=0;i<X;i+=2)
        for(int j=0;j<Y;j+=2){
            // (i,j)
            if(i-1>=0){
                adj_2d[i][j].push_back(make_pair(i-1,j));
                indegree_2d[i-1][j]++;
            }
            
            if(i+1<X){
                adj_2d[i][j].push_back(make_pair(i+1,j));
                indegree_2d[i+1][j]++;
            }
            
            //(i+1,j)
            if(i+1<X){
                if(j-1>=0){
                    adj_2d[i+1][j].push_back(make_pair(i+1,j-1));
                    indegree_2d[i+1][j-1]++;
                }
                
                if(j+1<Y){
                    adj_2d[i+1][j].push_back(make_pair(i+1,j+1));
                    indegree_2d[i+1][j+1]++;
                }
            }
            //(i+1,j+1)
            if((j+1<Y)&&(i+1<X)){
                adj_2d[i+1][j+1].push_back(make_pair(i,j+1));
                indegree_2d[i][j+1]++;
                if(i+2<Y){
                    adj_2d[i+1][j+1].push_back(make_pair(i+2,j+1));
                    indegree_2d[i+2][j+1]++;
                }
            }

        }
    cout<<"start 3d:"<<endl;
    //for 3d
    for(int i=0;i<X;i+=2){
        for(int j=0;j<Y;j+=2){
            for(int k=0;k<Z;k+=2){
                //(i,j,k)
                cout<<"1"<<endl;
                if(k-1>=0){
                    adj_3d[i][j][k].push_back(make_pair(make_pair(i,j),k-1));
                    indegree_3d[i][j][k-1]++;
                }
                if(k+1<Z){
                    adj_3d[i][j][k].push_back(make_pair(make_pair(i,j),k+1));
                    indegree_3d[i][j][k+1]++;
                }
                //(i,j,k+1)
                cout<<"2"<<endl;
                if(k+1<Z){
                    if(j-1>=0){
                        adj_3d[i][j][k+1].push_back(make_pair(make_pair(i,j-1),k+1));
                        indegree_3d[i][j-1][k+1]++;
                    }
                    if(j+1<Y){
                        adj_3d[i][j][k+1].push_back(make_pair(make_pair(i,j+1),k+1));
                        indegree_3d[i][j+1][k+1]++;
                    }
                }
                //(i,j+1,k+1)
                cout<<"3"<<endl;
                if((j+1<Y)&&(k+1<Z)){
                    if(k+1+1<Z){
                        adj_3d[i][j+1][k+1].push_back(make_pair(make_pair(i,j+1),k+2));
                        indegree_3d[i][j+1][k+2]++;
                    }
                    if(k+1-1>=0){
                        adj_3d[i][j+1][k+1].push_back(make_pair(make_pair(i,j+1),k));
                        indegree_3d[i][j+1][k]++;
                    }
                }
                //(i,j+1,k)
                cout<<"4"<<endl;
                if(j+1<Y){
                    if(i-1>=0){
                        adj_3d[i][j+1][k].push_back(make_pair(make_pair(i-1,j+1),k));
                        indegree_3d[i-1][j+1][k]++;
                    }
                    if(i+1<X){
                        adj_3d[i][j+1][k].push_back(make_pair(make_pair(i+1,j+1),k));
                        indegree_3d[i+1][j+1][k]++;
                    }
                }

                //(i+1,j+1,k)
                cout<<"5"<<endl;
                if((i+1<X)&&(j+1<Y)){
                    if(k-1>=0){
                        adj_3d[i+1][j+1][k].push_back(make_pair(make_pair(i+1,j+1),k-1));
                        indegree_3d[i+1][j+1][k-1]++;
                    }
                    if(k+1<Z){
                        adj_3d[i+1][j+1][k].push_back(make_pair(make_pair(i+1,j+1),k+1));
                        indegree_3d[i+1][j+1][k+1]++;
                    }
                }
                //(i+1,j+1,k+1)
                cout<<"6"<<endl;
                if((i+1<X)&&(j+1<Y)&&(k+1)<Z){
                    if(j+1-1>=0){
                        adj_3d[i+1][j+1][k+1].push_back(make_pair(make_pair(i+1,j),k+1));
                        indegree_3d[i+1][j][k+1]++;
                    }
                    if(j+1+1<Y){
                        adj_3d[i+1][j+1][k+1].push_back(make_pair(make_pair(i+1,j+2),k+1));
                        indegree_3d[i+1][j+2][k+1]++;
                    }
                }
                //(i+1,j,k+1)
                cout<<"7"<<endl;
                if((i+1<X)&&(k+1<Z)){
                    if(k+1-1>=0){
                        adj_3d[i+1][j][k+1].push_back(make_pair(make_pair(i+1,j),k));
                        indegree_3d[i+1][j][k]++;
                    }
                    if(k+1+1<Z){
                        adj_3d[i+1][j][k+1].push_back(make_pair(make_pair(i+1,j),k+2));
                        indegree_3d[i+1][j][k+2]++;
                    }
                }
            }
            
        }
    }
    cout<<"second part"<<endl;
    //second part  Topological Order
    //2d:
    queue<pair<int,int> > q_2d;
    for(int i=0;i<X;i++)
        for(int j=0;j<Y;j++){
            if(indegree_2d[i][j]==0)
                q_2d.push(make_pair(i,j));
        }
    while(!q_2d.empty()){
        pair<int ,int> temp=q_2d.front();
        q_2d.pop();
        cout<<"point ( "<<temp.first<<" , "<<temp.second<<" ) With gray code: "
        <<task_2d[temp.first][temp.second].first<<" "
        <<task_2d[temp.first][temp.second].second<<endl;
        for(int i=0;i<adj_2d[temp.first][temp.second].size();i++){
            pair<int,int> t=adj_2d[temp.first][temp.second][i];
            indegree_2d[t.first][t.second]--;
            if(indegree_2d[t.first][t.second]==0)
                q_2d.push(t);
        }
    }

    //3d
    queue<pair<pair<int,int>,int> > q_3d;
    for(int i=0;i<X;i++)
        for(int j=0;j<Y;j++)
            for(int k=0;k<Z;k++){
                if(indegree_3d[i][j][k]==0)
                    q_3d.push(make_pair(make_pair(i,j),k));
            }
    while(!q_3d.empty()){
        pair<pair<int ,int >,int> temp=q_3d.front();
        q_3d.pop();
        cout<<"point ( "<<temp.first.first<<" , "<<temp.first.second<<" , "<<temp.second<<" ) With gray code: "
        <<task_3d[temp.first.first][temp.first.second][temp.second].first.first<<" "
        <<task_3d[temp.first.first][temp.first.second][temp.second].first.second<<" "
        <<task_3d[temp.first.first][temp.first.second][temp.second].second<<endl;
        for(int i=0;i<adj_3d[temp.first.first][temp.first.second][temp.second].size();i++){
            pair<pair<int ,int >,int >t=adj_3d[temp.first.first][temp.first.second][temp.second][i];
            indegree_3d[t.first.first][t.first.second][t.second]--;
            if(indegree_3d[t.first.first][t.first.second][t.second]==0)
                q_3d.push(t);
        }
    }
    return 0;
}
