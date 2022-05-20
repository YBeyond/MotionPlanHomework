/*
Copyright (C) 2022 Hongkai Ye (kyle_yeh@163.com)
Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO
EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY
OF SUCH DAMAGE.
*/
#ifndef _BIAS_SAMPLER_
#define _BIAS_SAMPLER_

#include <ros/ros.h>
#include <Eigen/Eigen>
#include <random>

class BiasSampler
{
public:
  BiasSampler()
  {
    std::random_device rd;
    gen_ = std::mt19937_64(rd());
    uniform_rand_ = std::uniform_real_distribution<double>(0.0, 1.0);
    normal_rand_ = std::normal_distribution<double>(0.0, 1.0);
    range_.setZero();
    origin_.setZero();
  };

  void setSamplingRange(const Eigen::Vector3d origin, const Eigen::Vector3d range)
  {
    origin_ = origin;
    range_ = range;
  }

  void samplingOnce(Eigen::Vector3d &sample)
  {
    sample[0] = uniform_rand_(gen_);
    sample[1] = uniform_rand_(gen_);
    sample[2] = uniform_rand_(gen_);
    sample.array() *= range_.array();
    sample += origin_;
  };

  void setInformedSamplingRange(const Eigen::Vector3d scale, const Eigen::Vector3d center, 
                                const Eigen::Matrix3d rotation)
  {
    scale_ = scale;
    center_ = center;
    rotation_ = rotation;
    
  }

  void samplingByEllipsoid(Eigen::Vector3d &sample){
    // x = r * theta, x为椭球内任意一点，r为x的模，theta为方位向量
    Eigen::Vector3d theta;
    theta[0] = normal_rand_(gen_);
    theta[1] = normal_rand_(gen_);
    theta[2] = normal_rand_(gen_);
    
    double r = pow(uniform_rand_(gen_), 1 / 3);  // r = y^(1/n),y为[0,1]上均匀分布的随机变量
    sample = r * theta.normalized() * uniform_rand_(gen_);

    // 转化为实际椭球空间
    sample.array() *= scale_.array();
    sample = rotation_ * sample;
    sample += center_;
  };

  // (0.0 - 1.0)
  double getUniRandNum()
  {
    return uniform_rand_(gen_);
  }

private:
  Eigen::Vector3d range_, origin_;
  std::mt19937_64 gen_;
  std::uniform_real_distribution<double> uniform_rand_;
  std::normal_distribution<double> normal_rand_;

  Eigen::Vector3d scale_,center_;
  Eigen::Matrix3d rotation_;
};

#endif
