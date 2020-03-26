// Copyright (C) 2008 Wim Meeussen <meeussen at willowgarage com>
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation; either version 2.1 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
//

#include <bfl/wrappers/rng/rng.h>  // Wrapper around several rng libraries
#include <robot_pose_ekf/nonlinearanalyticconditionalgaussianodo.h>
#define NUMCONDARGUMENTS_MOBILE 2

namespace BFL {
using namespace MatrixWrapper;

/*
AnalyticConditionalGaussianAdditiveNoise:这个类表示的概率密度模型是
               P(A|B,C,D,...)----> P(x_t|u_t,x_{t-1})
					其中，mu_A=f(B,C,D,...)+mu_additiveNoise
					sigma_A=sigma_additiveNoise
					A=N(mu_A,sigma_A)
需要两个参数，一个是系统不确定性（BFL::Gaussian，即x_t），
第二个是条件参数的数量（num_conditional_arguments，即u_t,x_{t-1},2个）默认值是1
*/

/*
EKF 本质上将非线性模型线性化，这样满足误差还是满足高斯分布
    而将非线性线性化，就得根据实际情况进行重新定义，重新定义也就是重载的两个虚函数，预测值 和  雅克比矩阵F
    预测值，直接使用定义好函数（即f = ExpectedValueGet），即预测不确定性的miu
    雅克比矩阵，本质上用来对不确定度高斯分布进行变换， 即预测不确定的sigma
*/

/*
  更新部分因为只有长度，所以还是用线性分析条件概率高斯分布，即本身是线性的，沿用旧的类接口即可
*/

NonLinearAnalyticConditionalGaussianOdo::NonLinearAnalyticConditionalGaussianOdo(const Gaussian& additiveNoise)
    : AnalyticConditionalGaussianAdditiveNoise(additiveNoise, NUMCONDARGUMENTS_MOBILE),
      df(6, 6)  //f的导数
{
  // initialize df matrix,现在是单位矩阵了
  for (unsigned int i = 1; i <= 6; i++) {
    for (unsigned int j = 1; j <= 6; j++) {
      if (i == j)
        df(i, j) = 1;
      else
        df(i, j) = 0;
    }
  }
}

NonLinearAnalyticConditionalGaussianOdo::~NonLinearAnalyticConditionalGaussianOdo() {}

/*
   x(t)=f(x(t-1),u(t))+Q(t),Q白噪声均值; Q 即为addtiveNoiseGaussian
  根据系统模型和k-1时刻数据得到k时刻的预测值：
  x_k=x_{k-1}+u1_k*cos(yaw_{k-1})
  y_k=y_{k-1}+u1_k*sin(yaw_{k-1})
  yaw_k=yaw_{k-1}+u2_k
  
  矩阵表示f(x,u)= [
					x_k   =x_{k-1}     +u1_k*cos(yaw_{k-1});  //f(1)
					y_k   =y_{k-1}     +u1_k*sin(yaw_{k-1})   //f(2)
					yaw_k =yaw_{k-1}   +u2_k                  //f(3)  用来求f函数的一阶雅克比矩阵
				]
  
  */
ColumnVector NonLinearAnalyticConditionalGaussianOdo::ExpectedValueGet() const {
  //系统参数[x,u],ConditionalArgumentGet(第几个参数)
  //state：[x,y,z,pitch,roll,yaw]
  ColumnVector state = ConditionalArgumentGet(0);  //这个很重要，filter.update， 系统模型和input 信息，和这部分关联
  //u=[u1,u2]=[v*T,w*T]
  ColumnVector vel = ConditionalArgumentGet(1);  //这个很重要
  state(1) += cos(state(6)) * vel(1);
  state(2) += sin(state(6)) * vel(1);
  state(6) += vel(2);
  return state + AdditiveNoiseMuGet();
}

//Get the i-th matrix of the system.
/*
  x(t)=f(x(t-1),u(t))+Q(t)
  那么求偏导df/dx=
     [1, 0, 0, 0, 0, -v*T*sin(yaw), 
			0, 1, 0, 0, 0, v*T*cos(yaw),
			0, 0, 1, 0, 0, 0,
			0, 0, 0, 1, 0, 0,
			0, 0, 0, 0, 1, 0,
			0, 0, 0, 0, 0, 1]
			
			不应该是这样吗？
    */
Matrix NonLinearAnalyticConditionalGaussianOdo::dfGet(unsigned int i) const {
  if (i == 0)  //derivative to the first conditional argument (x)
  {
    //平移距离,也就是v*T
    double vel_trans = ConditionalArgumentGet(1)(1);
    //w*T
    double yaw = ConditionalArgumentGet(0)(6);

    df(1, 3) = -vel_trans * sin(yaw);
    df(2, 3) = vel_trans * cos(yaw);

    return df;
  } else {
    if (i >= NumConditionalArgumentsGet()) {
      cerr << "This pdf Only has " << NumConditionalArgumentsGet() << " conditional arguments\n";
      exit(-BFL_ERRMISUSE);
    } else {
      cerr << "The df is not implemented for the" << i << "th conditional argument\n";
      exit(-BFL_ERRMISUSE);
    }
  }
}

}  //namespace BFL
