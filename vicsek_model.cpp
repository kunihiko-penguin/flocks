//2d vicsek model
//这版还不是topological interaction，是半径r内的interaction
#include<iostream>
#include<cmath>
#include<random>
#include<fstream>
#include<vector>
#include<iomanip>

using namespace std;

const double eta = 0.1;
constexpr int TOT_particles = 512;
constexpr int STEPS = 1000;
constexpr double SPEED=1.0;
constexpr double RADIUS=1.0; //this need to be revised into topological interactions,between the first n neighbours
constexpr double RADIUS2=RADIUS*RADIUS;
constexpr double COUPLING=1.0;//这个没用到
constexpr float dt=0.1;
constexpr float L = 1.0f;
constexpr int n_topos=7;//range of topological interactions

//random number generator
static std::mt19937 rng(std::random_device{}());
static uniform_real_distribution<> dis(0,L);
static std::uniform_real_distribution<> dis_theta(-eta, eta); 

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

//update state of particles
void UpdateState(std::vector<particle>& particles)
{
    vector<particle> copy_particles = particles;
    for (int i = 0; i < TOT_particles; ++i) 
    {
        particle &p = particles[i];
        double new_theta = 0.0;
        int neighbours = 0; //这里还不是拓扑相互作用，是r以内的align
        double tan_neigh = 0.0;
        double noise = 0.0;

        //update position
        double dx = SPEED * cos(p.theta) * dt;
        double dy = SPEED * sin(p.theta) * dt;
        particle copy_p = copy_particles[i];
        copy_p.x += dx;
        copy_p.y += dy;

        //periodic boundary conditions
        if (copy_p.x < 0) copy_p.x += L;
        if (copy_p.x >= L) copy_p.x -= L;
        if (copy_p.y < 0) copy_p.y += L;
        if (copy_p.y >= L) copy_p.y -= L;

        //update orientation(注意这里和文章不一样，文章是平均的向量，这里是tan_theta)
        for (int j = 0; j < TOT_particles; ++j) 
        {
            //judge if the particle is within range r
            if (i != j) 
            {
                particle j_p = particles[j];
                double dist_x = p.x - j_p.x;
                double dist_y = p.y - j_p.y;
                if (dist_x*dist_x + dist_y*dist_y < RADIUS2) 
                {
                    double tan_j = tan(j_p.theta);
                    tan_neigh += tan_j;
                    neighbours++;
                }
            }
            if (i==j) continue;
        }

        tan_neigh += tan(p.theta); //include particle i
        tan_neigh /= (neighbours+1);
        new_theta = std::atan(tan_neigh);
        noise = dis_theta(rng);
        new_theta += noise;

        copy_p.theta = new_theta;
    }

    particles = std::move(copy_particles);
}

int main()
{
    std::vector<particle> particles;

    //initialization
    initiate(particles);

    std::ofstream file("particles.csv");
    file<<"step,id,x,y,theta\n";
    file<< std::fixed << std::setprecision(6);

    for (int step = 0; step<STEPS;++step)
    {
        particle par;
        for (int i =0; i<TOT_particles;++i)
        {
            par=particles[i];
            file<<step<<","<<i<<","<<par.x<<","<<par.y<<","<<par.theta<<"\n";
        }
        UpdateState(particles);
    }
    file.close();
    cout<<"completed"<<endl;
    return 0;
}