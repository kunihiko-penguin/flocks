#计算临界指数（参考关于vicsek model的文章，和2012年的那篇）
#先复现vicsek model的结果

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# ==== 读取数据 ====
df = pd.read_csv("particles.csv")

# ==== 第一张图：最后一步的静态箭头图 ====
last_step = df["step"].max()
df_last = df[df["step"] == last_step]

x = df_last["x"].values
y = df_last["y"].values
theta = df_last["theta"].values
u = np.cos(theta)
v = np.sin(theta)

plt.figure(figsize=(6, 6))
plt.quiver(x, y, u, v, angles='xy', scale_units='xy', scale=1, width=0.003)
plt.title(f"Particle Directions at Step {last_step}")
plt.xlabel("x")
plt.ylabel("y")
plt.axis("equal")
plt.grid(True)
plt.tight_layout()
plt.show()  # 显示第一张图


# ==== 第二张图：动画动图（保存 + 可选展示） ====
steps = sorted(df["step"].unique())
fig, ax = plt.subplots(figsize=(6, 6))
quiver = ax.quiver([], [], [], [], angles='xy', scale_units='xy', scale=1, width=0.003)

def init():
    ax.set_xlim(0, 5)  # 若模拟区域不是 [0,1]，修改为实际范围
    ax.set_ylim(0, 5)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title("Particle Directions Over Time")
    ax.grid(True)
    return quiver,

def update(frame):
    step = steps[frame]
    df_step = df[df["step"] == step]
    x = df_step["x"].values
    y = df_step["y"].values
    theta = df_step["theta"].values
    u = np.cos(theta)
    v = np.sin(theta)
    quiver.set_offsets(np.c_[x, y])
    quiver.set_UVC(u, v)
    ax.set_title(f"Step: {step}")
    return quiver,

ani = animation.FuncAnimation(fig, update, frames=len(steps), init_func=init, blit=True)

# 保存动画（可选：显示动画）
ani.save("particles_animation.gif", fps=10)
plt.show()  # 显示第二张图（只在某些环境下能看到动图）


# ==== 第三张图：序参量随时间演化 ====
TOT_particles = len(df[df["step"] == 0])  # 自动获取粒子总数
SPEED = 0.03  # 你的模拟速度

order_params = []

for step in steps:
    df_step = df[df["step"] == step]
    vx = np.cos(df_step["theta"].values)
    vy = np.sin(df_step["theta"].values)
    Vx = np.sum(vx)
    Vy = np.sum(vy)
    v_avg = np.sqrt(Vx**2 + Vy**2) / (TOT_particles * SPEED)
    order_params.append(v_avg)

plt.figure(figsize=(8, 5))
plt.plot(steps, order_params, label='Average Normalized Velocity')
plt.xlabel("Step")
plt.ylabel("Order Parameter")
plt.title("Time Evolution of Order Parameter")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()  # 显示第三张图
