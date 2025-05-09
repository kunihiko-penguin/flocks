//2d vicsek model
//这版还不是topological interaction，是半径r内的interaction
#include<iostream>
#include<cmath>
#include<random>
#include<fstream>
#include<vector>
#include<iomanip>
#include<array>
#include<sstream>
#include<omp.h>

using namespace std;

double eta = 0.0;
double L = 3.1;  // 移除 constexpr
int TOT_particles = 40;  // 移除 constexpr
constexpr int STEPS = 1000;
constexpr double SPEED=0.03;
constexpr double RADIUS=1.0; //this need to be revised into topological interactions,between the first n neighbours
constexpr double RADIUS2=RADIUS*RADIUS;
constexpr double COUPLING=1.0;//这个没用到
constexpr float dt=1;
constexpr int n_topos=7;//range of topological interactions
constexpr double PI = 3.14159265358979323846;

// Output folder
const string folder="eta";

//random number generator
static std::mt19937 rng(std::random_device{}());
static uniform_real_distribution<> dis(0,L);
static std::uniform_real_distribution<> dis_theta(-PI, PI); 

// function prototypes
struct particle;
void initiate(std::vector<particle>& particles);
vector<vector<int>> CellList(std::vector<particle>& particles);
vector<particle> local_neigh(int id, vector<particle>& particles, vector<vector<int>>& cell_list);
void UpdateState(std::vector<particle>& particles);
double periodic_dis(double x1, double x2);

//store position and orientation of particles
struct particle 
{
    double x, y,theta;
};

void initiate(std::vector<particle>& particles)
{
    particles.resize(TOT_particles);
    for (int i=0;i < TOT_particles;++i)
    {
        particle &p = particles[i];
        p.x = dis(rng);
        p.y = dis(rng);
        p.theta = dis_theta(rng);
    }
}

//shortest distance under PBC
double periodic_dis(double x1, double x2)
{
    double dx = x1-x2;
    dx -= L * std::round(dx/L);
    return dx;
}

//update state of particles
void UpdateState(std::vector<particle>& particles)
{
    vector<particle> copy_particles = particles;
    vector<vector<int>> cell_list = CellList(particles);
    #pragma omp parallel
    {
        // 每个线程创建自己的随机数生成器
        unsigned int seed = std::random_device{}() + omp_get_thread_num();
        std::mt19937 local_rng(seed);
        std::uniform_real_distribution<> local_dis_eta(-1, 1);
        
        #pragma omp for
        for (int i = 0; i < TOT_particles; ++i) 
        {
            particle &p = particles[i];
            double new_theta = 0.0;
            int n_neigh = 0; //这里还不是拓扑相互作用，是r以内的align
            double tan_neigh = 0.0;
            double noise = 0.0;

            //update position
            double dx = SPEED * cos(p.theta) * dt;
            double dy = SPEED * sin(p.theta) * dt;
            particle &copy_p = copy_particles[i];
            copy_p.x += dx;
            copy_p.y += dy;

            //periodic boundary conditions
            if (copy_p.x < 0) copy_p.x += L;
            if (copy_p.x >= L) copy_p.x -= L;
            if (copy_p.y < 0) copy_p.y += L;
            if (copy_p.y >= L) copy_p.y -= L;

            //update orientation
            vector<particle> neighbours = local_neigh(i, particles,cell_list); 

            for (particle j_p : neighbours)
            {
                double dist_x = periodic_dis(p.x, j_p.x);
                double dist_y = periodic_dis(p.y, j_p.y);
                if (dist_x*dist_x + dist_y*dist_y < RADIUS2) 
                {
                    double tan_j = tan(j_p.theta);
                    tan_neigh += tan_j;
                    n_neigh++;
                }
            }

            tan_neigh += tan(p.theta); //include particle i
            tan_neigh /= (n_neigh+1);
            new_theta = std::atan(tan_neigh);
            noise = local_dis_eta(local_rng)*eta;
            new_theta += noise;

            copy_p.theta = new_theta;
        }
    }

    particles = std::move(copy_particles);
}

vector<vector<int>> CellList(vector<particle>& particles)
{
    vector<vector<int>> cell_list;
    int N_cell = floor(L/RADIUS);
    double L_cell = L/N_cell;
    //binning the particles into cells
    cell_list.resize(N_cell*N_cell);
    for (int i=0; i< TOT_particles; ++i)
    {
        int cx=0,cy=0,id_cell=0;
        particle p = particles[i];
        cx = floor(p.x/L_cell);
        cy = floor(p.y/L_cell);
        if (cx>=N_cell) cx -= N_cell;
        if (cy>=N_cell) cy -= N_cell;
        id_cell = cx + cy*N_cell;
        cell_list[id_cell].push_back(i);
    }
    return cell_list;
}

vector<particle> local_neigh(int id, vector<particle>& particles, vector<vector<int>>& cell_list)
{
    vector<particle> neighbours;
    int N_cell = floor(L/RADIUS);
    double L_cell = L/N_cell;
    //find the cell of particle id
    particle p = particles[id];
    int cx = floor(p.x/L_cell);
    int cy = floor(p.y/L_cell);
    if (cx>=N_cell) cx -= N_cell;
    if (cy>=N_cell) cy -= N_cell;
    int id_cell = cx + cy*N_cell;
    //find the neighbours of cell of id
    array<int,9> cell_neigh_list = {0,0,0,0,0,0,0,0,0};
    int x_left=cx-1<0?N_cell+cx-1:cx-1;
    int x_right=cx+1>=N_cell?cx+1-N_cell:cx+1;
    int y_up=cy-1<0?N_cell+cy-1:cy-1;
    int y_down=cy+1>=N_cell?cy+1-N_cell:cy+1;
    array<int,3> cell_x_list={x_left,cx,x_right};
    array<int,3> cell_y_list={y_up,cy,y_down};
    for (int i=0;i<3;++i)
    {
        for (int j=0;j<3;++j)
        {
            int cn = cell_x_list[i] + cell_y_list[j]*N_cell;
            cell_neigh_list[i+j*3] = cn;
        }
    }
    for (int j : cell_neigh_list)
    {
        for (int k : cell_list[j])
        {
            if (k != id)//exclude itself
            {
                neighbours.push_back(particles[k]);
            }
        }
    }
    return neighbours;
}

int main() 
{
    // 定义不同的参数组合
    struct SimParams {
        int particles;
        double size;
        string foldername;
    };
    
    vector<SimParams> params = {
        {40, 3.1, "eta/N40"},
        {100, 5.0, "eta/N100"},
        {400, 10.0, "eta/N400"},
        {4000, 31.6, "eta/N4000"},
        {10000, 50.0, "eta/N10000"}
    };
    
    // 设定 eta 的范围和步长
    const double eta_start = 0.0;
    const double eta_end = 5.0;
    const double eta_step = 0.1;
    
    // 对每组参数进行模拟
    for(const auto& param : params) 
    {
        std::cout << "Starting simulation for N = " << param.particles 
             << ", L = " << param.size << endl;
        
        // 更新全局变量
        TOT_particles = param.particles;
        L = param.size;
        
        // 对每个 eta 值执行模拟
        for(eta = eta_start; eta <= eta_end; eta += eta_step) 
        {
            // 构建文件名
            ostringstream oss;
            oss << fixed << setprecision(1) << eta;
            string filename = param.foldername + "/particles_eta_" + oss.str() + ".csv";
            std::ofstream file(filename);
            
            if (!file.is_open()) {
                std::cerr << "无法打开文件: " << filename << endl;
                continue;
            }
            
            // 初始化粒子系统
            std::vector<particle> particles;
            initiate(particles);
            
            // 写入CSV文件头
            file << "step,id,x,y,theta\n";
            file << std::fixed << std::setprecision(6);
            std::ostringstream buffer;
            // 运行模拟
            for (int step = 0; step < STEPS; ++step) 
            {
                for (int i = 0; i < TOT_particles; ++i) 
                {
                    const particle& par = particles[i];
                    buffer << step << "," << i << "," << par.x << "," << par.y << "," << par.theta << "\n";
                }

                if (step % 50 == 0 || step == STEPS - 1) {
                    file << buffer.str();
                    buffer.str(""); // 清空内容
                    buffer.clear(); // 重置状态
                }

                UpdateState(particles);
            }
            
            file.close();
            std::cout << "eta = " << eta << endl;
        }
        std::cout << " N = " << param.particles  << endl;
    }
    
    std::cout << "completed" << std::endl;
    return 0;
}