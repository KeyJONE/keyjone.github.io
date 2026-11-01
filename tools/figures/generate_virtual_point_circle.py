"""Regenerate Figure 4 for the axisymmetric notes with legible point labels."""

from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import Circle, Arc

ROOT = Path(__file__).resolve().parents[2]
font_manager.fontManager.addfont(r'C:\Windows\Fonts\msyh.ttc')
plt.rcParams.update({
    'font.family': 'Microsoft YaHei', 'font.size': 13,
    'mathtext.fontset': 'stix', 'axes.unicode_minus': False,
    'savefig.facecolor': 'white',
})
navy, blue, orange, muted = '#21334d', '#397bbb', '#d88726', '#60748d'
fig = plt.figure(figsize=(9, 4.5), dpi=200, facecolor='white')
fig.text(.5, .945, '四个等权虚拟点：同一圆周，相邻夹角为 90°',
         ha='center', color=navy, fontsize=17, weight='bold')
ax = fig.add_axes([.025, .065, .535, .82])
ax.set_aspect('equal'); ax.set_xlim(-1.55, 1.55); ax.set_ylim(-1.30, 1.40)
ax.axis('off')
ax.add_patch(Circle((0, 0), 1, fc='#edf4fc', ec=blue, lw=1.8))
for end in [(1.40, 0), (0, 1.31)]:
    start = (-1.35, 0) if end[0] else (0, -1.23)
    ax.annotate('', xy=end, xytext=start,
                arrowprops={'arrowstyle':'->', 'lw':1.0, 'color':muted})
ax.text(1.44, -.04, r'$\delta r$', fontsize=16, color=navy)
ax.text(.07, 1.27, r'$\delta z$', fontsize=16, color=navy)
ax.scatter([0], [0], s=16, color=navy, zorder=5)
ax.text(-.11, -.16, r'$0$', fontsize=14, color=navy)
# Keep the horizontal labels clear of the axis and the top label above the circle.
points=[(-1,0,'1',orange,-.16,.18,'right'),
        (1,0,'2',orange,.14,.18,'left'),
        (0,-1,'3',blue,.12,-.10,'left'),
        (0,1,'4',blue,.20,.14,'left')]
for x,y,label,color,dx,dy,align in points:
    ax.scatter([x],[y],s=65,c=color,edgecolors='white',linewidths=1.2,zorder=6)
    ax.text(x+dx,y+dy,r'$g='+label+'$',ha=align,va='center',color=color,fontsize=15,
            zorder=7)
ax.add_patch(Arc((0,0), .56, .56, theta1=0, theta2=90, lw=1.4, color=orange))
ax.text(.31,.31,r'$90^\circ$',color=orange,fontsize=14)
ax.annotate('', xy=(-np.sqrt(.5),np.sqrt(.5)), xytext=(0,0),
            arrowprops={'arrowstyle':'<->','lw':1.1,'color':muted,'shrinkA':0,'shrinkB':0})
ax.text(-.50,.61,r'$\ell_k$',fontsize=17,color=navy)
ax.text(0,-1.27,'均匀二维圆盘',ha='center',va='top',fontsize=13,color=muted)

right=fig.add_axes([.60,.12,.36,.73]); right.axis('off')
right.text(.5,.96,'各点的径向偏移',ha='center',va='center',color=navy,fontsize=15,weight='bold')
right.text(.22,.80,r'$g$',ha='center',fontsize=16,color=navy)
right.text(.72,.80,r'$\delta r_g$',ha='center',fontsize=16,color=navy)
right.plot([.06,.94],[.75,.75],color='#c8d7e6',lw=1)
for row,(g,value) in enumerate([(1,r'-\ell_k'),(2,r'+\ell_k'),(3,'0'),(4,'0')]):
    y=.65-row*.13
    color=orange if g in (1,2) else blue
    right.text(.22,y,f'${g}$',ha='center',va='center',fontsize=17,color=color)
    right.text(.72,y,f'${value}$',ha='center',va='center',fontsize=18,color=color)
right.plot([.06,.94],[.18,.18],color='#c8d7e6',lw=1)
right.text(.5,.08,r'每点归一化权重 $1/4$',ha='center',va='center',color=navy,fontsize=14)
right.set_xlim(0,1);right.set_ylim(0,1)
out=ROOT/'site/images/notes/axisymmetric/04_virtual_point_circle.png'
fig.savefig(out,dpi=200)
plt.close(fig)
print(out)
