# 第五章作业
## 参考公式
假设以五次多项式表示距离$x$，可以使得jerk连续可导，即$x(t) = C_0 + C_1t + C_2 t^2 + C_3 t^3 + C_4 t^4 + C_5 t^5$，在空间中的路径可表示为：
$$
\begin{cases}p_x= C_{0_x}+C_{1_x}t+C_{2_x}t^2 + C_{3_x}t^3 + C_{4_x}t^4 + C_{5_x}t^5\\ v_x= C_{1_x}+2C_{2_x}t + 3C_{3_x}t^2 + 4C_{4_x}t^3 + 5C_{5_x}t^4\\ a_x= 2C_{2_x}+ 6C_{3_x}t + 12C_{4_x}t^2 + 20C_{5_x}t^3\\ j_x= 6C_{3_x} + 24C_{4_x}t + 60C_{5_x}t^2\end{cases}
$$

$$
\begin{cases}p_y= C_{0_y}+C_{1_y}t+C_{2_y}t^2 + C_{3_y}t^3 + C_{4_y}t^4 + C_{5_y}t^5\\ v_y= C_{1_y}+2C_{2_y}t + 3C_{3_y}t^2 + 4C_{4_y}t^3 + 5C_{5_y}t^4\\ a_y= 2C_{2_y}+ 6C_{3_y}t + 12C_{4_y}t^2 + 20C_{5_y}t^3\\ j_y= 6C_{3_y} + 24C_{4_y}t + 60C_{5_y}t^2\end{cases}
$$


$$
\begin{cases}p_z= C_{0_z}+C_{1_z}t+C_{2_z}t^2 + C_{3_z}t^3 + C_{4_z}t^4 + C_{5_z}t^5\\ v_z= C_{1_z}+2C_{2_z}t + 3C_{3_z}t^2 + 4C_{4_z}t^3 + 5C_{5_z}t^4\\ a_z= 2C_{2_z}+ 6C_{3_z}t + 12C_{4_z}t^2 + 20C_{5_z}t^3\\ j_z= 6C_{3_z} + 24C_{4_z}t + 60C_{5_z}t^2\end{cases}
$$

jerk可表示为$J(T) = \int_{T_{j-1}}^{T_j}(j_x(t)^2 + j_y(t)^2+j_z(t)^2 )dt$

根据$J = \left[ \begin{matrix} p_1 \\ \vdots \\ p_M \end{matrix} \right]^T  \left[ \begin{matrix} Q_1 & 0 & 0 \\ 0 & \ddots & 0 \\ 0 & 0 & Q_M \end{matrix} \right]  \left[ \begin{matrix} p_1 \\ \vdots \\ p_M \end{matrix} \right]$ 和 $ M_j p_j = d_j $，对$J$进行展开，得到：
$$
J = \left[ \begin{matrix} d_1 \\ \vdots \\ d_M \end{matrix} \right]^T \left[ \begin{matrix} M_1 & 0 & 0 \\ 0 & \ddots & 0 \\ 0 & 0 & M_M \end{matrix} \right]^{-T} \left[ \begin{matrix} Q_1 & 0 & 0 \\ 0 & \ddots & 0 \\ 0 & 0 & Q_M \end{matrix} \right] \left[ \begin{matrix} M_1 & 0 & 0 \\ 0 & \ddots & 0 \\ 0 & 0 & M_M \end{matrix} \right]^{-1} \left[ \begin{matrix} d_1 \\ \vdots \\ d_M \end{matrix} \right]
$$
其中$C^T \left[ \begin{matrix} d_F \\ d_P \end{matrix} \right] = \left[ \begin{matrix} d_1 \\ \vdots\\ d_M\end{matrix} \right]$，分解展开得到$J = d_F^TR_{FF}d_F + d_F^TR_{FP}d_P + d_P^TR_{PF}d_F + d_P^TR_{PP}d_P$

对$d_P$求导，得到$d_P^* = -R_{PP}^{-1}R_{FP}^Td_F$

其中多项式系数矩阵$Coeff$，映射矩阵$M$，jerk中心矩阵$Q$，以及选择矩阵$C$分别表示为
$$
Coeff = \left[ \begin{matrix} C_{0_x} & C_{0_y} & C_{0_z} \\ C_{1_x} & C_{1_y} & C_{1_z}  \\ C_{2_x} & C_{2_y} & C_{2_z} \\ C_{3_x} & C_{3_y} & C_{3_z} \\ C_{4_x} & C_{4_y} & C_{4_z} \\ C_{5_x} & C_{5_y} & C_{5_z}  \end{matrix} \right]
$$

$$
M_j = \left[ \begin{matrix} 1 & t_j & t_j^2 & t_j^3 & t_j^4 & t_j^5  \\ 0 & 1 & 2t_j & 3t_j^2 & 4t_j^3 & 5t_j^4   \\ 0 & 0 & 2 & 6t_j & 12t_j^2 & 20t_j^3  \\ 1 & t_{j+1} & t_{j+1}^2 & t_{j+1}^3 & t_{j+1}^4 & t_{j+1}^5  \\ 0 & 1 & 2t_{j+1} & 3t_{j+1}^2 & 4t_{j+1}^3 & 5t_{j+1}^4   \\ 0 & 0 & 2 & 6t_{j+1} & 12t_{j+1}^2 & 20t_{j+1}^3   \end{matrix} \right]
$$

$$
Q_j = \left[ \begin{matrix} \vdots & \vdots & \vdots \\ \cdots & \frac{i(i-1)(i-2)l(l-1)(l-2)}{i+l-5}t_j^{i+l-5} & \cdots \\ \vdots & \vdots & \vdots \end{matrix} \right]
$$

选择矩阵C示例，包含起点$p_0$，中间点$p_1$,终点$p_k$：  
$$
\left[ \begin{matrix}p_0 \\ v_0 \\ a_0 \\ p_1 \\ v_1 \\ a_1 \\ p_1 \\ v_1 \\ a_1 \\ p_k \\ v_k \\ a_k \end{matrix}\right]=\left[ \begin{matrix}1&0&0&0&0&0&0&0&0 \\ 0&1&0&0&0&0&0&0&0 \\ 0&0&1&0&0&0&0&0&0 \\ 0&0&0&1&0&0&0&0&0 \\ 0&0&0&0&0&0&0&1&0 \\ 0&0&0&0&0&0&0&0&1 \\ 0&0&0&1&0&0&0&0&0 \\ 0&0&0&0&0&0&0&1&0 \\ 0&0&0&0&0&0&0&0&1 \\ 0&0&0&0&1&0&0&0&0 \\ 0&0&0&0&0&1&0&0&0 \\ 0&0&0&0&0&0&1&0&0\end{matrix}\right] \left[ \begin{matrix} p_0  \\ v_0  \\ a_0 \\ p_1 \\ p_k \\ v_k \\ a_k \\ v_1 \\ a_1  \end{matrix} \right]
$$

当求得$d_P$时，根据$M_jp_j = d_j$即可求得每一段曲线的系数。


## 实现效果
目前代码还存在bug...
![运行结果](https://gitee.com/lxyclara/motion-plan-homework/raw/lxy/L5/pic/result.jpg "算法运行结果")



## 所遇问题
1.使用求解器获得数值解：使用了[OSQP](https://osqp.org/)  
2.Closed-form获得解析解：依赖于[Eigen](https://blog.csdn.net/s12k39/article/details/108381018)  
3.使用Eigen库中矩阵resize()函数后，如果矩阵中元素没有被重新复制，会出现异常值，因此建议resize后对矩阵赋0操作(zero())

