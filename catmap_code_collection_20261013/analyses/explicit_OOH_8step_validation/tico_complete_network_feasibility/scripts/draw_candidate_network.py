"""Draw an evidence-bounded diagram for the Ti-Co NEB/cycle feasibility audit."""
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

OUT = Path(__file__).resolve().parents[1] / "figures" / "TICO_CANDIDATE_REACTION_NETWORK.png"


def node(ax, x, y, title, lines, *, width=1.8, height=1.05, color="#35536b",
         face="#ffffff", ts=False):
    box = FancyBboxPatch(
        (x - width / 2, y - height / 2), width, height,
        boxstyle="round,pad=0.035,rounding_size=0.06",
        linewidth=1.5, edgecolor=color, facecolor=face,
    )
    ax.add_patch(box)
    ax.text(x, y + 0.22, title, ha="center", va="center", fontsize=9.4,
            fontweight="bold", color="#17324a")
    ax.text(x, y - 0.12, lines, ha="center", va="center", fontsize=7.4,
            color="#263746", linespacing=1.22)
    if ts:
        ax.text(x, y - height / 2 - 0.13, "真实 NEB TS", ha="center",
                va="top", fontsize=7.1, color="#26805a")


def link(ax, x1, x2, label, *, validated=False, label_y=4.05):
    color = "#26805a" if validated else "#c65d20"
    style = "-" if validated else (0, (4, 3))
    ax.add_patch(FancyArrowPatch(
        (x1, 3.0), (x2, 3.0), arrowstyle="-|>", mutation_scale=12,
        linewidth=1.8, linestyle=style, color=color, shrinkA=2, shrinkB=2,
    ))
    if label:
        ax.text((x1 + x2) / 2, label_y, label, ha="center", va="center",
                fontsize=7.4, color=color,
                bbox=dict(facecolor="white", edgecolor="none", alpha=0.96, pad=1.5))


def main():
    plt.rcParams.update({"font.family": "Microsoft YaHei", "axes.unicode_minus": False})
    fig, ax = plt.subplots(figsize=(17.4, 7.4), dpi=180)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")
    ax.set_xlim(-0.2, 17.2)
    ax.set_ylim(-0.8, 5.6)
    ax.axis("off")

    xs = [0.9, 3.0, 5.0, 7.0, 9.4, 11.4, 13.4, 15.6]
    node(ax, xs[0], 3.0, "独立前体", "Ha* + OOH* + C6H12*",
         width=1.8, face="#f4f7f9")
    node(ax, xs[1], 3.0, "NEB1-IS", "Ha + OOH + C6H12\n共吸附态")
    node(ax, xs[2], 3.0, "NEB1-TS", "ΔG‡f = 0.556344 eV\nΔG‡r = 0.796030 eV",
         width=1.75, color="#26805a", face="#f0faf4", ts=True)
    node(ax, xs[3], 3.0, "NEB1-FS", "Ha + C6H12O* + O* + Hb*")
    node(ax, xs[4], 3.0, "NEB2-IS", "Ha + Hb* + O*")
    node(ax, xs[5], 3.0, "NEB2-TS", "ΔG‡f = 0.225994 eV\nΔG‡r = 0.706574 eV",
         width=1.75, color="#26805a", face="#f0faf4", ts=True)
    node(ax, xs[6], 3.0, "NEB2-FS", "Ha + OH*")
    node(ax, xs[7], 3.0, "独立水态", "H2O*", width=1.55, face="#f4f7f9")

    link(ax, 1.82, 2.08, "ΔG候选 −2.663642 eV\n连接 TS/势垒缺失",
         validated=False, label_y=4.15)
    link(ax, 3.92, 4.08, "NEB1局部路径", validated=True, label_y=4.46)
    link(ax, 5.92, 6.08, "NEB1局部路径", validated=True, label_y=4.46)
    link(ax, 7.92, 8.48, "ΔG候选 +0.811237 eV\n释放 C6H12O(g)；连接 TS/势垒缺失",
         validated=False, label_y=4.15)
    link(ax, 10.32, 10.48, "NEB2局部路径", validated=True, label_y=4.46)
    link(ax, 12.32, 12.48, "NEB2局部路径", validated=True, label_y=4.46)
    link(ax, 14.32, 14.82, "ΔG候选 −0.601109 eV\n成水/重排 TS 或势垒缺失",
         validated=False, label_y=4.15)

    # Product gas released in connection 2; this is a formula-level bookkeeping branch.
    node(ax, 8.25, 1.05, "气相产物", "C6H12O(g)", width=1.9, height=0.8,
         color="#71828f", face="#f4f7f9")
    ax.annotate("", xy=(8.2, 1.48), xytext=(8.15, 2.43),
                arrowprops=dict(arrowstyle="-|>", color="#71828f", lw=1.2,
                                linestyle=(0, (3, 3))))

    ax.text(8.5, 5.35, "Ti–Co：共吸附 NEB 与完整催化循环的最小可行性检查",
            ha="center", va="center", fontsize=15, fontweight="bold", color="#17324a")
    ax.text(8.5, 4.92,
            "绿色实线：已通过 CatMAP 0.3.1 单步复现　|　橙色虚线：端点 ΔG 候选已记录，但连接动力学缺失",
            ha="center", va="center", fontsize=9.0, color="#405463")
    ax.text(0.0, -0.1,
            "Ha 在 NEB 状态公式中按既有定义保留；具体构型连续性/原子路径与共吸附态真实位点占据仍未验证。",
            ha="left", va="center", fontsize=8.2, color="#4d5c66")
    ax.text(0.0, -0.43,
            "单步测试中的 n_sites=1 仅为形式记账，不作为真实几何占位证据；本图不假设连接势垒，也不是可计算完整循环的模型。",
            ha="left", va="center", fontsize=8.2, color="#4d5c66")

    fig.tight_layout(pad=0.6)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=240, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(OUT)


if __name__ == "__main__":
    main()
