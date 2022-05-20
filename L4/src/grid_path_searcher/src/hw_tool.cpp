#include <hw_tool.h>

using namespace std;
using namespace Eigen;

void Homeworktool::initGridMap(double _resolution, Vector3d global_xyz_l, Vector3d global_xyz_u, int max_x_id, int max_y_id, int max_z_id)
{   
    gl_xl = global_xyz_l(0);
    gl_yl = global_xyz_l(1);
    gl_zl = global_xyz_l(2);

    gl_xu = global_xyz_u(0);
    gl_yu = global_xyz_u(1);
    gl_zu = global_xyz_u(2);
    
    GLX_SIZE = max_x_id;
    GLY_SIZE = max_y_id;
    GLZ_SIZE = max_z_id;
    GLYZ_SIZE  = GLY_SIZE * GLZ_SIZE;
    GLXYZ_SIZE = GLX_SIZE * GLYZ_SIZE;

    resolution = _resolution;
    inv_resolution = 1.0 / _resolution;    

    data = new uint8_t[GLXYZ_SIZE];
    memset(data, 0, GLXYZ_SIZE * sizeof(uint8_t));
}

void Homeworktool::setObs(const double coord_x, const double coord_y, const double coord_z)
{   
    if( coord_x < gl_xl  || coord_y < gl_yl  || coord_z <  gl_zl || 
        coord_x >= gl_xu || coord_y >= gl_yu || coord_z >= gl_zu )
        return;

    int idx_x = static_cast<int>( (coord_x - gl_xl) * inv_resolution);
    int idx_y = static_cast<int>( (coord_y - gl_yl) * inv_resolution);
    int idx_z = static_cast<int>( (coord_z - gl_zl) * inv_resolution);      
    
    data[idx_x * GLYZ_SIZE + idx_y * GLZ_SIZE + idx_z] = 1;
}

bool Homeworktool::isObsFree(const double coord_x, const double coord_y, const double coord_z)
{
    Vector3d pt;
    Vector3i idx;
    
    pt(0) = coord_x;
    pt(1) = coord_y;
    pt(2) = coord_z;
    idx = coord2gridIndex(pt);

    int idx_x = idx(0);
    int idx_y = idx(1);
    int idx_z = idx(2);

    return (idx_x >= 0 && idx_x < GLX_SIZE && idx_y >= 0 && idx_y < GLY_SIZE && idx_z >= 0 && idx_z < GLZ_SIZE && 
           (data[idx_x * GLYZ_SIZE + idx_y * GLZ_SIZE + idx_z] < 1));
}

Vector3d Homeworktool::gridIndex2coord(const Vector3i & index) 
{
    Vector3d pt;

    pt(0) = ((double)index(0) + 0.5) * resolution + gl_xl;
    pt(1) = ((double)index(1) + 0.5) * resolution + gl_yl;
    pt(2) = ((double)index(2) + 0.5) * resolution + gl_zl;

    return pt;
}

Vector3i Homeworktool::coord2gridIndex(const Vector3d & pt) 
{
    Vector3i idx;
    idx <<  min( max( int( (pt(0) - gl_xl) * inv_resolution), 0), GLX_SIZE - 1),
            min( max( int( (pt(1) - gl_yl) * inv_resolution), 0), GLY_SIZE - 1),
            min( max( int( (pt(2) - gl_zl) * inv_resolution), 0), GLZ_SIZE - 1);                  
  
    return idx;
}

Eigen::Vector3d Homeworktool::coordRounding(const Eigen::Vector3d & coord)
{
    return gridIndex2coord(coord2gridIndex(coord));
}

double Homeworktool::OptimalBVP(Eigen::Vector3d _start_position,Eigen::Vector3d _start_velocity,Eigen::Vector3d _target_position)
{
    double optimal_cost = 100000; // this just to initial the optimal_cost, you can delete it 
    /*
                    



    STEP 2: go to the hw_tool.cpp and finish the function Homeworktool::OptimalBVP
    the solving process has been given in the document

    because the final point of trajectory is the start point of OBVP, so we input the pos,vel to the OBVP

    after finish Homeworktool::OptimalBVP, the Trajctory_Cost will record the optimal cost of this trajectory


    */

    // 速度V是未指定的
    // delta_px = pxf-px0-vx*T; delta_py = pyf-py0-vy*T; delta_pz = pzf-pz0-vz*T;
    // a1 = -3 * delta_px / T^3; a2 = -3 * delta_py / T^3; a3 = -3 * delta_pz / T^3;
    // J = T + 1/3 * (a1^2) * T^3 + 1/3 * (a2^2) * T^3 + 1/3 * (a3^2) * T^3
    // 展开：J = T + 3(dx^2 + dy^2 + dz^2)*T^-3 - 6(dx*vx + dy*vy + dz*vz)*T^-2 + (vx^2 + vy^2 + vz^2)*T^-1;
    // dx = pxf-px0; dy = pyf-py0; dz = pzf-pz0;
    // 求导后： J' = 1 - 9(dx^2 + dy^2 + dz^2)*T^-4 + 12(dx*vx + dy*vy + dz*vz)*T^-3 - (vx^2 + vy^2 + vz^2)*T^-2;

    double dx = _target_position(0) - _start_position(0);
    double dy = _target_position(1) - _start_position(1);
    double dz = _target_position(2) - _start_position(2);
    double vx = _start_velocity(0);
    double vy = _start_velocity(1);
    double vz = _start_velocity(2);


    Eigen::PolynomialSolver<double, Eigen::Dynamic> solver; //多项式求根
    Eigen::VectorXd coeff(5);
    coeff[0] = -9.0 * (dx*dx + dy*dy + dz*dz);
    coeff[1] = 12.0 * (dx*vx + dy*vy + dz*vz);
    coeff[2] = -(vx*vx + vy*vy + vz*vz);
    coeff[3] = 0.0;
    coeff[4] = 1.0;
    solver.compute(coeff);    // compute()计算多项式复数根
    const Eigen::PolynomialSolver<double, Eigen::Dynamic>::RootsType &r = solver.roots(); // roots()返回多项式复数根
    std::vector<double> real_root;
    for (int i = 0; i < r.rows(); ++i)
    {
        if (r[i].imag() != 0.0)
        {
            continue;
        }
        real_root.push_back(r[i].real());
    }

    std::multimap<double,double> traj_cost;
    for (int i = 0; i < real_root.size(); ++i)
    {
        // std::cout << "======================" << std::endl;
        // cout << "多项式方程的实根:" << real_root[i] << "  J = " << real_root[i] + 3.0 * (dx * dx + dy * dy + dz * dz) * pow(real_root[i], -3) - 6.0 * (dx * vx + dy * vy + dz * vz) * pow(real_root[i], -2) + (vx * vx + vy * vy + vz * vz) * pow(real_root[i], -1) << std::endl;
        double cost = -1.0;
        if (real_root[i] < 0.0)
            continue;
        cost = real_root[i] + 3.0 * (dx * dx + dy * dy + dz * dz) * pow(real_root[i], -3) - 
                    6.0 * (dx * vx + dy * vy + dz * vz) * pow(real_root[i], -2) + 
                    (vx * vx + vy * vy + vz * vz) * pow(real_root[i], -1);
        if(cost > 0)
            traj_cost.insert({cost, real_root[i]});
    }

    if (!traj_cost.empty())
        optimal_cost = traj_cost.begin()->first;

    return optimal_cost;
}
