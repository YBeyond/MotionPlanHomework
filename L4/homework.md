# 第四章作业  
## 两个主要概念
- 控制空间采样，也是一种前向积分的过程。通过给定（并离散）控制量，通过运动模型，计算得到采样点位置  
- 状态空间采样，是一种后向计算的过程。通过一定规则或要求（目标导向）选择需要采样的状态或位置，计算初始状态到下一步状态的连接。当采样的层数较多时，会形成指数级的关联关系  


## 实现效果  
- 作业1： 推导  
![推导结果](https://gitee.com/lxyclara/motion-plan-homework/raw/lxy/L4/pic/formula.jpg "推导结果")


- 作业2： coding   
![运行结果](https://gitee.com/lxyclara/motion-plan-homework/raw/lxy/L4/pic/result.jpg "算法运行结果")


## 所遇问题
1.对J求极值：使用了Eigen多项式求解器

