#计算临界指数（参考关于vicsek model的文章，和2012年的那篇）
#先复现vicsek model的结果

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# 参数：粒子数和速度（与 C++ 中保持一致）
TOT_particles = 512  # 例如 L=50 时
SPEED = 0.03  # 你的设定速度，根据需要修改

# 读取 CSV 文件
df = pd.read_csv("particles.csv")

# 初始化结果列表
order_params = []

# 获取所有 step 的范围
steps = df["step"].unique()
steps.sort()

# 遍历每一个时间步，计算对应的序参量
for step in steps:
    df_step = df[df["step"] == step]

    # 将每个粒子的朝向转换为单位速度矢量
    vx = np.cos(df_step["theta"].values)
    vy = np.sin(df_step["theta"].values)

    # 所有粒子的速度矢量求和
    Vx = np.sum(vx)
    Vy = np.sum(vy)

    # 求合速度矢量的模长并归一化
    v_avg = np.sqrt(Vx**2 + Vy**2) / (TOT_particles * SPEED)
    order_params.append(v_avg)

# 绘图：序参量随时间演化
plt.figure(figsize=(8, 5))
plt.plot(steps, order_params, label='Average Normalized Velocity')
plt.xlabel("Step")
plt.ylabel("Order Parameter")
plt.title("Time Evolution of Order Parameter")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
