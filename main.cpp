#include <iostream>
#include <iomanip>
#include "robotKinematics.hpp"
#include "Timer.h"

int main()
{
    Timer timer;
    myRoboKine::RobotKinematics kine(0.0838, 0.2, 0.2);
    Eigen::Vector3d f_pos{-0.05, -0.01, 0.2};
    std::vector<Eigen::Vector3d> joint_angle;

    timer.start();
//    auto r = kine.calIk_Left(f_pos, joint_angle);
    auto r = kine.calIk_Right(f_pos, joint_angle);
    auto end_t9 = timer.getNs();
    std::cout << "time9: " << end_t9 << " ns." << std::endl;

    for (auto &elem : joint_angle)
    {
        std::cout << std::setprecision(15) << elem.transpose() << std::endl;
    }

    return 0;
}
