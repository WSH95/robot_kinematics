//
// Created by wsh on 2021/12/3.
//

#ifndef ROBOT_KINEMATICS_ROBOTKINEMATICS_HPP
#define ROBOT_KINEMATICS_ROBOTKINEMATICS_HPP

#include <cmath>
#include <Eigen/Core>
#include <vector>

#define ANY_ANGLE 123

namespace myRoboKine
{
    constexpr double PI = 3.1415926535897932384626;
    constexpr double N_PI = -3.1415926535897932384626;
    constexpr double PI2 = 2 * PI;

    class RobotKinematics
    {
    public:
        RobotKinematics()
        {
            l1 = 0.0;
            l2 = 0.0;
            l3 = 0.0;
            low1 = N_PI;
            up1 = PI;
            low2 = N_PI;
            up2 = PI;
            low3 = N_PI;
            up3 = PI;

            bound_s_r_min = sqrt(pow(l1, 2) + pow((l2 - l3), 2));
            bound_s_r_max = sqrt(pow(l1, 2) + pow((l2 + l3), 2));
            bound_c_r = l1;

            x_abs_limit = sqrt(pow(bound_s_r_max, 2) - pow(bound_c_r, 2));
        }

        RobotKinematics(double hip_len_, double thigh_len_, double shank_len_) : l1(hip_len_), l2(thigh_len_),
                                                                                 l3(shank_len_)
        {
            low1 = N_PI;
            up1 = PI;
            low2 = N_PI;
            up2 = PI;
            low3 = N_PI;
            up3 = PI;

            bound_s_r_min = sqrt(pow(l1, 2) + pow((l2 - l3), 2));
            bound_s_r_max = sqrt(pow(l1, 2) + pow((l2 + l3), 2));
            bound_c_r = l1;

            x_abs_limit = sqrt(pow(bound_s_r_max, 2) - pow(bound_c_r, 2));
        }

        ~RobotKinematics() = default;

        /// unit: m
        void setLinkLength(double hip_len_, double thigh_len_, double shank_len_);

        /// For left leg. unit: rad. Coordinate: x->forward, y->left.
        void setBound_Left(double low1_, double up1_, double low2_, double up2_, double low3_, double up3_);

//        /// unit: rad
//        void setJointAngle(double theta1_, double theta2_, double theta3_);

        void calFk_Left(Eigen::Vector3d &foot_pos_, double theta1_, double theta2_, double theta3_);

        Eigen::Vector3d calFk_Left(double theta1_, double theta2_, double theta3_);

        void calFk_Right(Eigen::Vector3d &foot_pos_, double theta1_, double theta2_, double theta3_);

        Eigen::Vector3d calFk_Right(double theta1_, double theta2_, double theta3_);

        void calJ_Left(Eigen::Matrix3d &J_, double theta1_, double theta2_, double theta3_);

        Eigen::Matrix3d calJ_Left(double theta1_, double theta2_, double theta3_);

        void calJ_Right(Eigen::Matrix3d &J_, double theta1_, double theta2_, double theta3_);

        Eigen::Matrix3d calJ_Right(double theta1_, double theta2_, double theta3_);

        /*!
         * @brief inverse kinematics of left leg
         * @param foot_pos_ (x, y, z) (unit: m)
         * @param joint_angles_ all solutions. find the nearest solution.
         * @return negative if no solution
         */
        double calIk_Left(const Eigen::Vector3d &foot_pos_, std::vector<Eigen::Vector3d> &joint_angles_);

        double calIkOptimal_Left(const Eigen::Vector3d &foot_pos_, const Eigen::Vector3d &current_angles,
                                 Eigen::Vector3d &joint_angles_);

        double normalizeAngleToPi(double angle_)
        {
            if ((angle_ >= N_PI) && (angle_ <= PI))
                return angle_;
            else
            {
                double norm_angle = fmod(angle_, PI2);
                if (norm_angle < N_PI)
                    norm_angle += PI2;
                else if (norm_angle > PI)
                    norm_angle -= PI2;

                return norm_angle;
            }
        }

        void normalizeAngleToPi(Eigen::Vector3d &angle_)
        {
            angle_[0] = normalizeAngleToPi(angle_[0]);
            angle_[1] = normalizeAngleToPi(angle_[1]);
            angle_[2] = normalizeAngleToPi(angle_[2]);
        }

    private:
        double l1, l2, l3; /// unit: m
//        double theta1, theta2, theta3; /// unit: rad
        double low1, up1, low2, up2, low3, up3; /// unit: rad
        double bound_s_r_min, bound_s_r_max, bound_c_r; /// s: sphere  c: cylinder
        double x_abs_limit;
    };

    double RobotKinematics::calIk_Left(const Eigen::Vector3d &foot_pos_, std::vector<Eigen::Vector3d> &joint_angles_)
    {
        joint_angles_.clear();
        double reward = 0;
        Eigen::Vector3d tmp_angle;
        Eigen::Vector3d new_pos;
        new_pos = foot_pos_;

        double distance_origin2foot;
        double distance2_origin2projectedPointYZ;
        double distance_tanPoint2projectedPointYZ;
        double y_tanPoint, z_tanPoint;

        distance2_origin2projectedPointYZ = new_pos.tail(2).squaredNorm();
        if (sqrt(distance2_origin2projectedPointYZ) - bound_c_r < 1.0e-10) /// inner small cylinder
        {
            if (new_pos[0] + x_abs_limit < 1.0e-10)
            {
                reward += pow(100 * (new_pos[0] + x_abs_limit), 2);
                new_pos[0] = -1.0 * x_abs_limit;
            }
            else if (new_pos[0] - x_abs_limit > -1.0e-10)
            {
                reward += pow(100 * (new_pos[0] - x_abs_limit), 2);
                new_pos[0] = x_abs_limit;
            }
            reward += (pow(100 * new_pos[1], 2) + pow(100 * (new_pos[2] + bound_c_r), 2));
            new_pos[1] = 0.0;
            new_pos[2] = -1.0 * bound_c_r;
        }

        distance_origin2foot = new_pos.norm();

        if (distance_origin2foot - bound_s_r_min < 1.0e-10) /// inner small sphere
        {
            reward = reward + 1.0 * pow(100.0 * (distance_origin2foot - bound_s_r_min), 2);
            if (bound_s_r_min - l1 < 1.0e-10) /// l2 = l3
            {
                tmp_angle[0] = atan2(foot_pos_[2], foot_pos_[1]);
                tmp_angle[1] = ANY_ANGLE;
                tmp_angle[2] = N_PI;
                normalizeAngleToPi(tmp_angle);
                joint_angles_.push_back(tmp_angle);

                return -sqrt(reward);
            }
            else
            {
                /// projection of the points onto the small sphere.
                new_pos = new_pos * ((bound_s_r_min + 1.0e-7) / distance_origin2foot);
            }
        }
        else if (distance_origin2foot - bound_s_r_max > -1.0e-10) /// outer max sphere
        {
            reward = reward + 1.0 * pow(100.0 * (distance_origin2foot - bound_s_r_max), 2);
            /// projection of the points onto the big sphere.
            new_pos = new_pos * ((bound_s_r_max - 1.0e-7) / distance_origin2foot);
        }

        distance2_origin2projectedPointYZ = new_pos.tail(2).squaredNorm();
        distance_tanPoint2projectedPointYZ = sqrt(distance2_origin2projectedPointYZ - pow(l1, 2));

        /// situation 1
        // cancel l1 in numerator and 'distance2_origin2projectedPointYZ' in denominator.
        y_tanPoint = l1 * new_pos[1] + new_pos[2] * distance_tanPoint2projectedPointYZ;
        z_tanPoint = l1 * new_pos[2] - new_pos[1] * distance_tanPoint2projectedPointYZ;
        tmp_angle[0] = atan2(z_tanPoint, y_tanPoint);

        // two link
        double beata, alpha, gamma;
        double x2addz2 = pow(new_pos[0], 2) + pow(distance_tanPoint2projectedPointYZ, 2);
        beata = acos((l2 * l2 + l3 * l3 - x2addz2) / (2 * l2 * l3));
        alpha = acos((l2 * l2 + x2addz2 - l3 * l3) / (2 * l2 * sqrt(x2addz2)));
        gamma = atan2(-1.0 * new_pos[0], -1.0 * distance_tanPoint2projectedPointYZ);

        // situation 1.1
        tmp_angle[1] = gamma - alpha;
        tmp_angle[2] = PI - beata;
        normalizeAngleToPi(tmp_angle);
        joint_angles_.push_back(tmp_angle);

        // situation 1.2
        tmp_angle[1] = gamma + alpha;
        tmp_angle[2] = beata - PI;
        normalizeAngleToPi(tmp_angle);
        joint_angles_.push_back(tmp_angle);

        /// situation 2
        // cancel l1 in numerator and 'distance2_origin2projectedPointYZ' in denominator.
        y_tanPoint = l1 * new_pos[1] - new_pos[2] * distance_tanPoint2projectedPointYZ;
        z_tanPoint = l1 * new_pos[2] + new_pos[1] * distance_tanPoint2projectedPointYZ;
        tmp_angle[0] = atan2(z_tanPoint, y_tanPoint);

        gamma = atan2(-1.0 * new_pos[0], distance_tanPoint2projectedPointYZ);

        // situation 2.1
        tmp_angle[1] = gamma - alpha;
        tmp_angle[2] = PI - beata;
        normalizeAngleToPi(tmp_angle);
        joint_angles_.push_back(tmp_angle);

        // situation 2.2
        tmp_angle[1] = gamma + alpha;
        tmp_angle[2] = beata - PI;
        normalizeAngleToPi(tmp_angle);
        joint_angles_.push_back(tmp_angle);

        return -sqrt(reward);
    }
}


#endif //ROBOT_KINEMATICS_ROBOTKINEMATICS_HPP
