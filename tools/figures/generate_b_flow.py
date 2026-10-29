"""Generate the adopted SNS-FEM matrix-to-internal-force workflow."""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = ROOT / "site/images/notes/axisymmetric/02_B_to_internal_force_current.png"

# Prefer the project's Windows host font; allow a CJK font on other platforms.
for candidate in (Path("C:/Windows/Fonts/msyh.ttc"), Path("/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc")):
    if candidate.exists():
        font_manager.fontManager.addfont(str(candidate))
        FONT = font_manager.FontProperties(fname=str(candidate)).get_name()
        break
else:
    raise RuntimeError("Install Microsoft YaHei or Noto Sans CJK to render this figure.")

plt.rcParams.update({
    "font.family": FONT,
    "mathtext.fontset": "stix",
    "axes.unicode_minus": False,
    "savefig.facecolor": "white",
})
fig, ax = plt.subplots(figsize=(18, 7.2), dpi=120)
fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
ax.set(xlim=(0, 1800), ylim=(720, 0))
ax.axis("off")

INK = "#21354F"
BLUE = "#397CBD"
MUTED = "#536A82"
GREEN = "#16816E"


def label(x, y, value, size=20, color=INK, weight="normal"):
    return ax.text(x, y, value, ha="center", va="center", fontsize=size,
                   color=color, fontweight=weight)


def panel(x, number, title):
    ax.add_patch(FancyBboxPatch((x, 130), 390, 510,
                               boxstyle="round,pad=0,rounding_size=18",
                               facecolor="#F5F8FC", edgecolor="#D7E2EF", linewidth=1.3))
    ax.add_patch(Circle((x + 42, 176), 23, color=BLUE))
    label(x + 42, 176, str(number), size=20, color="white", weight="bold")
    label(x + 225, 176, title, size=21, weight="bold")


def arrow(start, end, color=BLUE):
    ax.add_patch(FancyArrowPatch(start, end, arrowstyle="-|>",
                                mutation_scale=16, linewidth=1.8, color=color))


label(900, 52, "由平均矩阵到骨架内力", size=28, weight="bold")
xs = [30, 480, 930, 1380]
titles = ["精确体积平均", "构造四点矩阵", "本构更新", "组装骨架内力"]
for n, (x, title) in enumerate(zip(xs, titles), 1):
    panel(x, n, title)
for x in xs[:-1]:
    arrow((x + 398, 390), (x + 442, 390))

# 1. The adopted exact volume average, with no unresolved alternative.
label(225, 267, r"原始 $\mathbf{B}$", size=21)
label(225, 311, "按实际旋转体积平均", size=21)
label(225, 393, r"$\overline{\mathbf{B}}_k^V=\dfrac{1}{V_k}\int_{D_k} \mathbf{B}\,\mathrm{d}V$", size=25)
label(225, 478, "环向行精确积分", size=21, color=BLUE, weight="bold")
ax.add_patch(FancyBboxPatch((83, 537), 284, 52,
                           boxstyle="round,pad=0,rounding_size=12",
                           facecolor="#E5EFFB", edgecolor="#BED4EF"))
label(225, 563, "方法三 · 本文采用", size=20, color=BLUE, weight="bold")

# 2. Four stored matrix states, obtained directly by adding/subtracting gradients.
label(675, 267, "平均矩阵 ± 梯度项", size=21)
cx, cy, radius = 675, 415, 90
ax.add_patch(Circle((cx, cy), radius, fill=False, edgecolor="#9CB1C9",
                    linewidth=1.2, linestyle=(0, (4, 4))))
ax.plot([cx-radius, cx+radius], [cy, cy], color="#C9D5E3", linewidth=1.3)
ax.plot([cx, cx], [cy-radius, cy+radius], color="#C9D5E3", linewidth=1.3)
ax.add_patch(Circle((cx, cy), 5, color=MUTED))
for n, dx, dy, tx, ty in (
    (1, -radius, 0, -25, 0), (2, radius, 0, 25, 0),
    (3, 0, radius, 0, 27), (4, 0, -radius, 0, -27),
):
    ax.add_patch(Circle((cx+dx, cy+dy), 8, color=BLUE))
    label(cx+dx+tx, cy+dy+ty, str(n), size=19, color=BLUE, weight="bold")
label(675, 578, r"得到四个虚拟点 $\mathbf{B}$", size=21)

# 3. Constitutive updates use total increments from the converged state.
label(1125, 267, r"四点 $\mathbf{B}$ 与热应变矩阵", size=20)
label(1125, 311, "本步位移、温度总增量", size=20)
label(1125, 355, "上一步已收敛状态", size=20)
arrow((1125, 391), (1125, 428), color=MUTED)
ax.add_patch(FancyBboxPatch((982, 442), 286, 63,
                           boxstyle="round,pad=0,rounding_size=14",
                           facecolor="#E3F2ED", edgecolor="#B3DCCD"))
label(1125, 474, "弹性 / J2 / DP", size=22, color=GREEN, weight="bold")
arrow((1125, 520), (1125, 550), color=MUTED)
label(1125, 585, "得到四点有效应力", size=21, color=GREEN, weight="bold")

# 4. Only the solid-skeleton contribution is shown.
label(1575, 267, "每个虚拟点分别计算", size=20)
ax.add_patch(FancyBboxPatch((1415, 313), 320, 130,
                           boxstyle="round,pad=0,rounding_size=14",
                           facecolor="#FFF4E2", edgecolor="#ECCB94"))
label(1575, 377, r"$\widetilde{\mathbf{B}}_{kg}^{\,T}\,\boldsymbol{\sigma}'_{kg}\,\dfrac{V_k}{4}$", size=27, color="#8C560D")
arrow((1575, 462), (1575, 499), color=MUTED)
label(1575, 532, "四点求和", size=21, weight="bold")
label(1575, 572, "得到光滑域骨架内力", size=20, weight="bold")
label(1575, 611, "孔压内力另行组装", size=16, color=MUTED)

label(900, 684, "位移、温度总增量均相对上一步已收敛状态计算。", size=19, color=MUTED)
OUTPUT.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(OUTPUT, dpi=120)
plt.close(fig)
print(OUTPUT)
