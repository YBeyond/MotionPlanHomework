#include "lec5_hw/visualizer.hpp"
#include "lec5_hw/trajectory.hpp"

#include <ros/ros.h>
#include <geometry_msgs/Point.h>
#include <geometry_msgs/PoseStamped.h>

#include <cmath>
#include <iostream>
#include <vector>

struct Config
{
    std::string targetTopic;
    double clickHeight;
    std::vector<double> initialVel;
    std::vector<double> initialAcc;
    std::vector<double> terminalVel;
    std::vector<double> terminalAcc;
    double allocationSpeed;
    double allocationAcc;
    int maxPieceNum;

    Config(const ros::NodeHandle &nh_priv)
    {
        nh_priv.getParam("TargetTopic", targetTopic);
        nh_priv.getParam("ClickHeight", clickHeight);
        nh_priv.getParam("InitialVel", initialVel);
        nh_priv.getParam("InitialAcc", initialAcc);
        nh_priv.getParam("TerminalVel", terminalVel);
        nh_priv.getParam("TerminalAcc", terminalAcc);
        nh_priv.getParam("AllocationSpeed", allocationSpeed);
        nh_priv.getParam("AllocationAcc", allocationAcc);
        nh_priv.getParam("MaxPieceNum", maxPieceNum);
    }
};


// 时间间隔离散
double timeTrapzVel(const double dist,
                    const double vel,
                    const double acc)
{
    const double t = vel / acc;
    const double d = 0.5 * acc * t * t;

    if (dist < d + d)
    {
        return 2.0 * sqrt(dist / acc);
    }
    else
    {
        return 2.0 * t + (dist - 2.0 * d) / vel;
    }
}

void minimumJerkTrajGen(
    // Inputs:
    const int pieceNum,
    const Eigen::Vector3d &initialPos,
    const Eigen::Vector3d &initialVel,
    const Eigen::Vector3d &initialAcc,
    const Eigen::Vector3d &terminalPos,
    const Eigen::Vector3d &terminalVel,
    const Eigen::Vector3d &terminalAcc,
    const Eigen::Matrix3Xd &intermediatePositions,
    const Eigen::VectorXd &timeAllocationVector,
    // Outputs:
    Eigen::MatrixX3d &coefficientMatrix)
{
    // coefficientMatrix is a matrix with 6*piece num rows and 3 columes
    // As for a polynomial c0+c1*t+c2*t^2+c3*t^3+c4*t^4+c5*t^5,
    // each 6*3 sub-block of coefficientMatrix is
    // --              --
    // | c0_x c0_y c0_z |
    // | c1_x c1_y c1_z |
    // | c2_x c2_y c2_z |
    // | c3_x c3_y c3_z |
    // | c4_x c4_y c4_z |
    // | c5_x c5_y c5_z |
    // --              --
    // Please computed coefficientMatrix of the minimum-jerk trajectory
    // in this function

    // ------------------------ Put your solution below ------------------------

    // 使用eigen求逆，直接求解解析解
    int order = 5, segment = pieceNum;
    int var_dim = 3; // 考虑p,v,a,非x,y,z分量
    // 定义求解矩阵
    Eigen::MatrixXd Q, M, C, df;

    // 构建M矩阵
    // p = [1 t t^2 t^3 t^4 t^5]
    // v = [0 1 2t 3t^2 4t^3 5t^4]
    // a = [0 0 2 6t 12t^2 20t^3]
    M.resize(2 * segment * var_dim, (order + 1) * segment);
    M = Eigen::MatrixXd::Zero(2 * segment * var_dim, (order + 1) * segment);
    double t_s = 0.0, t_e = 0.0;
    Eigen::MatrixXd M_j = Eigen::MatrixXd::Zero(2 * var_dim, order + 1);
    for (int i = 0; i < (int)timeAllocationVector.rows(); i++)
    {
        
        // int row_inc = 0;
        // t_s = t_e;
        t_e = t_s + timeAllocationVector(i);
        M_j.row(0) << 1.0, t_s, pow(t_s, 2), pow(t_s, 3), pow(t_s, 4), pow(t_s, 5);
        M_j.row(1) << 0.0, 1.0, 2.0 * t_s, 3.0 * pow(t_s, 2), 4.0 * pow(t_s, 3), 5.0 * pow(t_s, 5);
        M_j.row(2) << 0.0, 0.0, 2.0, 6.0 * t_s, 12.0 * pow(t_s, 2), 20.0 * pow(t_s, 3);
        M_j.row(3) << 1.0, t_e, pow(t_e, 2), pow(t_e, 3), pow(t_e, 4), pow(t_e, 5);
        M_j.row(4) << 0.0, 1.0, 2.0 * t_e, 3.0 * pow(t_e, 2), 4.0 * pow(t_e, 3), 5.0 * pow(t_e, 5);
        M_j.row(5) << 0.0, 0.0, 2.0, 6.0 * t_e, 12.0 * pow(t_e, 2), 20.0 * pow(t_e, 3);
        // M.block(i * 2 * var_dim, i * (order + 1), 2 * var_dim, order + 1) = M_j;
        M.block<6, 6>(i * 6, i * 6) = M_j;
    }
    // ROS_INFO("完成M矩阵计算");
    // debug message
    #define M_DEBUG 1
    #if M_DEBUG
    std::cout << "time rows :" << timeAllocationVector.rows() << " cols :" << timeAllocationVector.cols() << std::endl;
    std::cout << "time allocate :\n"
              << timeAllocationVector << std::endl;
    std::cout << "Matrix M :\n" << M << std::endl;
    #endif

    // 构建Q矩阵
    Q.resize((order+1) * segment, (order+1) * segment);
    Q = Eigen::MatrixXd::Zero((order + 1) * segment, (order + 1) * segment);
    // std::cout << "Matrix Q resize:\n" << Q << std::endl;
    double t_n = 0.0;
    for (int i = 0; i < (int)timeAllocationVector.rows(); i++)
    {
        Eigen::MatrixXd Q_n = Eigen::MatrixXd::Zero(order + 1, order + 1);
        t_n = timeAllocationVector(i);
        for (int j = var_dim; j < order + 1; j++)
        {
            for (int k = var_dim; k < order + 1; k++) {
                int iden = j + k - order;
                Q_n(j, k) = (double)j * (j - 1) * (j - 2) * k * (k - 1) * (k - 2) / iden * pow(t_n, iden);
            }
        }
        // Q.block(i * (order + 1), i * (order + 1), order + 1, order + 1) = Q_n;
        Q.block<6, 6>(i * 6, i * 6) = Q_n;
    }
    // ROS_INFO("完成Q矩阵计算");
    // debug message
    #define Q_DEBUG 1
    #if Q_DEBUG
    std::cout << "Matrix Q :\n" << Q << std::endl;
    #endif

    // 构建df矩阵
    df.resize(var_dim * 2 + segment - 1, 3);

    df.row(0) = initialPos.transpose();
    df.row(1) = initialVel.transpose();
    df.row(2) = initialAcc.transpose();
    for (int i = 0; i < segment - 1; i++)
    {
        df.row(i + var_dim) = intermediatePositions.col(i).transpose();
    }
    df.row(var_dim + segment - 1) = terminalPos.transpose();
    df.row(var_dim + segment ) = terminalVel.transpose();
    df.row(var_dim + segment + 1) = terminalAcc.transpose();
    // ROS_INFO("完成df矩阵计算");
    // debug message
    #define DF_DEBUG 1
    #if DF_DEBUG
    std::cout << "Matrix df :\n" << df << std::endl;
    std::cout << "Matrix terminal pos :\n" << terminalPos << std::endl;
    std::cout << "Matrix terminal vel :\n" << terminalVel << std::endl;
    std::cout << "Matrix terminal acc :\n" << terminalAcc << std::endl;
    #endif

    // 构建C矩阵
    C.resize(2 * var_dim * segment, var_dim * (segment + 1));
    C = Eigen::MatrixXd::Zero(2 * var_dim * segment, var_dim * (segment + 1));
    // p0,v0,a0
    for (int i = 0; i < var_dim; i++)
    {
        C(i, i) = 1;
    }

    // p1,p2,...,p(n-1)
    for (int i = 1; i < segment; i++)
    {
        C(var_dim * (i * 2 - 1), var_dim + i -1) = 1;
        C(var_dim * i * 2, var_dim + i - 1) = 1;
    }

    // pn,vn,an
    for (int i = 0; i < var_dim; i++)
    {
        C(2 * var_dim * segment - var_dim + i, var_dim + segment - 1 + i) = 1;
    }

    // v1,a1,...,v(n-1),a(n-1)
    for (int i = 1; i < segment; i++)
    {
        C(var_dim * (i * 2 - 1) + 1, 2 * var_dim + segment - 2 + i * 2 - 1) = 1;
        C(var_dim * i * 2 + 1, 2 * var_dim + segment - 2 + i * 2 - 1) = 1;
        C(var_dim * (i * 2 - 1) + 2, 2 * var_dim + segment - 2 + i * 2) = 1;
        C(var_dim * i * 2 + 2, 2 * var_dim + segment - 2 + i * 2) = 1;
    }
    // ROS_INFO("完成C矩阵计算");
    // debug message
    #define C_DEBUG 1
    #if C_DEBUG
    std::cout << "Matrix C :\n" << C << std::endl;
    #endif

    // 求解矩阵dp
    Eigen::MatrixXd dp = Eigen::MatrixXd::Zero((var_dim-1) * (segment -1) ,3);
    Eigen::MatrixXd R = Eigen::MatrixXd::Zero(2 * var_dim * segment, 2 * var_dim * segment);
    // std::cout << "C size << " << C.rows() << " x " << C.cols() << std::endl;
    // std::cout << "M size << " << M.rows() << " x " << M.cols() << std::endl;
    // std::cout << "Q size << " << Q.rows() << " x " << Q.cols() << std::endl;
    R = C.transpose() * (M.inverse().transpose()) * Q * M.inverse() * C;
    // ROS_INFO("完成R矩阵计算");
    dp = -R.block(df.rows(), df.rows(), dp.rows(), dp.rows()).inverse() * R.block(0, df.rows(), df.rows(), dp.rows()).transpose() * df;
    // ROS_INFO("完成dp矩阵计算");
    // 求解每段曲线系数矩阵 coefficientMatrix
    Eigen::MatrixX3d dd = Eigen::MatrixXd::Zero(var_dim * (segment+1),3);
    dd.topRows(df.rows()) = df;
    dd.bottomRows(dp.rows()) = dp;
    // debug message
    #define dd_DEBUG 1
    #if dd_DEBUG
    std::cout << "Matrix dd :\n" <<  dd << std::endl;
    std::cout << "Matrix d :\n" << C * dd << std::endl;
    #endif
    // ROS_INFO("完成d矩阵计算");
    coefficientMatrix = M.inverse() * C * dd;
    // ROS_INFO("完成coeff矩阵计算");
    // debug message
#define coeff_DEBUG 1
#if coeff_DEBUG
    std::cout << "Matrix coeff :\n"
              << coefficientMatrix << std::endl;
#endif
    // ------------------------ Put your solution above ------------------------
}

class ClickGen
{
private:
    Config config;

    ros::NodeHandle nh;
    ros::Subscriber targetSub;

    Visualizer visualizer;

    Eigen::Matrix3Xd positions;
    Eigen::VectorXd times;
    int positionNum;
    Trajectory<5> traj;

public:
    ClickGen(const Config &conf,
             ros::NodeHandle &nh_)
        : config(conf),
          nh(nh_),
          visualizer(nh),
          positions(3, config.maxPieceNum + 1),
          times(config.maxPieceNum),
          positionNum(0)
    {
        targetSub = nh.subscribe(config.targetTopic, 1,
                                 &ClickGen::targetCallBack, this,
                                 ros::TransportHints().tcpNoDelay());
    }

    void targetCallBack(const geometry_msgs::PoseStamped::ConstPtr &msg)
    {
        if (positionNum > config.maxPieceNum)
        {
            positionNum = 0;
            traj.clear();
        }

        positions(0, positionNum) = msg->pose.position.x;
        positions(1, positionNum) = msg->pose.position.y;
        positions(2, positionNum) = std::fabs(msg->pose.orientation.z) * config.clickHeight;

        if (positionNum > 0)
        {
            const double dist = (positions.col(positionNum) - positions.col(positionNum - 1)).norm();
            times(positionNum - 1) = timeTrapzVel(dist, config.allocationSpeed, config.allocationAcc);
        }

        ++positionNum;

        if (positionNum > 1)
        {
            const int pieceNum = positionNum - 1;
            const Eigen::Vector3d initialPos = positions.col(0);
            const Eigen::Vector3d initialVel(config.initialVel[0], config.initialVel[1], config.initialVel[2]);
            const Eigen::Vector3d initialAcc(config.initialAcc[0], config.initialAcc[1], config.initialAcc[2]);
            const Eigen::Vector3d terminalPos = positions.col(pieceNum);
            const Eigen::Vector3d terminalVel(config.terminalVel[0], config.terminalVel[1], config.terminalVel[2]);
            const Eigen::Vector3d terminalAcc(config.terminalAcc[0], config.terminalAcc[1], config.terminalAcc[2]);
            const Eigen::Matrix3Xd intermediatePositions = positions.middleCols(1, pieceNum - 1);
            const Eigen::VectorXd timeAllocationVector = times.head(pieceNum);

            Eigen::MatrixX3d coefficientMatrix = Eigen::MatrixXd::Zero(6 * pieceNum, 3);

            minimumJerkTrajGen(pieceNum,
                               initialPos, initialVel, initialAcc,
                               terminalPos, terminalVel, terminalAcc,
                               intermediatePositions,
                               timeAllocationVector,
                               coefficientMatrix);

            traj.clear();
            traj.reserve(pieceNum);
            for (int i = 0; i < pieceNum; i++)
            {
                traj.emplace_back(timeAllocationVector(i),
                                  coefficientMatrix.block<6, 3>(6 * i, 0).transpose().rowwise().reverse());
            }
        }

        visualizer.visualize(traj, positions.leftCols(positionNum));

        return;
    }
};

int main(int argc, char **argv)
{
    ros::init(argc, argv, "click_gen_node");
    ros::NodeHandle nh_;
    ClickGen clickGen(Config(ros::NodeHandle("~")), nh_);
    ros::spin();
    return 0;
}
