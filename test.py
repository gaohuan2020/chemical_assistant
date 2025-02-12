import numpy as np
import matplotlib.pyplot as plt
from library.chemistry import Species, Equilibrium
from library.equilibria import EqSystem
from matplotlib.font_manager import FontManager
import matplotlib.font_manager as fm

# 获取系统中文字体
fm = FontManager()
font_names = sorted([f.name for f in fm.ttflist])
# 设置中文字体
plt.rcParams['font.sans-serif'] = ['SimHei']  # 用来正常显示中文标签
plt.rcParams['axes.unicode_minus'] = False  # 用来正常显示负号

# 1. 定义参与反应的物质
substance_names = ["H+", "OH-", "NH4+", "NH3", "H2O"]
subst = {n: Species.from_formula(n) for n in substance_names}

# 2. 设定初始浓度
init_conc = {
    "H+": 1e-7,  # 0.0000001 M
    "OH-": 1e-7,  # 0.0000001 M
    "NH4+": 1e-7,  # 0.0000001 M
    "NH3": 1.0,  # 1 M
    "H2O": 55.5,  # 55.5 M (纯水的浓度)
}

# 3. 定义化学平衡
# 水的电离平衡：H2O ⇌ H+ + OH-
w_autop = Equilibrium(
    {"H2O": 1}, {"H+": 1, "OH-": 1}, 10**-14 / 55.5  # 反应物  # 生成物  # 平衡常数
)

# 铵离子的解离平衡：NH4+ ⇌ H+ + NH3
NH4p_pr = Equilibrium(
    {"NH4+": 1}, {"H+": 1, "NH3": 1}, 10**-9.26  # 反应物  # 生成物  # 平衡常数
)

equilibria = (w_autop, NH4p_pr)

# 4. 创建平衡系统
eqsys = EqSystem(equilibria, subst)

# 5. 计算不同H+浓度下的平衡
nc = 60  # 取60个点
Hp_0 = np.logspace(-3, 0, nc)  # H+浓度从0.001到1

# 6. 绘图
plt.figure(figsize=(10, 6))
plt.grid(True)

# 计算并绘制结果
results = eqsys.roots(init_conc, Hp_0, "H+")
x, info, sane = results

# 设置最小浓度阈值，避免对数刻度的问题
min_conc = 1e-20  # 设置一个极小的浓度值作为下限
x = np.maximum(x, min_conc)  # 将所有小于min_conc的值替换为min_conc

# 在绘图之前添加调试信息
print("数据范围：")
for i, name in enumerate(substance_names):
    print(f"{name}: {x[:, i].min():.2e} to {x[:, i].max():.2e}")

# 在添加pH刻度之前
print("pH范围：", min(-np.log10(Hp_0)), "到", max(-np.log10(Hp_0)))

# 为每个物质画一条线
for i, name in enumerate(substance_names):
    plt.plot(Hp_0, x[:, i], label=name, linewidth=2)

# 设置图表格式
plt.xscale("log")
plt.yscale("log")
plt.xlabel("H+ 浓度 (M)", fontsize=12)
plt.ylabel("物质浓度 (M)", fontsize=12)
plt.title("氨水体系中各物质浓度随pH变化的关系", fontsize=14)
plt.legend(fontsize=10, loc='best')
plt.grid(True, which="both", ls="-", alpha=0.2)

# 添加pH刻度
ax2 = plt.gca().twiny()
pH_values = [-np.log10(h) for h in Hp_0]
ax2.set_xlim(plt.gca().get_xlim())

plt.tight_layout()
plt.show()
