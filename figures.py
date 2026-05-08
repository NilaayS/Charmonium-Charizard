"""Generate the three figures for the charmonium poster."""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams

rcParams['font.family'] = 'serif'
rcParams['font.size'] = 11
rcParams['mathtext.fontset'] = 'dejavuserif'
rcParams['axes.linewidth'] = 0.8
rcParams['xtick.major.width'] = 0.8
rcParams['ytick.major.width'] = 0.8

# ---------------------------------------------------------------
# Figure 1: Cornell potential and its three pieces
# ---------------------------------------------------------------
# units: r in fm, V in GeV.  GeV * fm conversion: hbar c = 0.1973 GeV*fm
hbarc = 0.1973  # GeV * fm

# Parameters (Ali et al. 2016)
alpha_s = 0.4827
b_GeV2  = 0.1488          # GeV^2  (string tension)
m_q     = 1.4399          # GeV
sigma   = 1.2819          # GeV    -> 1/sigma in fm = hbarc/sigma
S       = 1               # spin-1 contribution

# Convert r-axis from fm to "natural" GeV^-1: r[GeV^-1] = r[fm]/hbarc
r_fm = np.linspace(0.001, 2.0, 1500)
r_nat = r_fm / hbarc                       # in GeV^-1

# Coulomb-like piece: -4/3 * alpha_s / r  (in GeV when r in GeV^-1)
V_coul = -4.0 / 3.0 * alpha_s / r_nat

# Linear piece: b * r   (b in GeV^2, r in GeV^-1 -> GeV)
V_lin = b_GeV2 * r_nat

# Spin-spin smeared delta:
#   delta_sigma(r) = (sigma^3 / pi^{3/2}) * exp(-sigma^2 r^2)    (in GeV^3 when r in GeV^-1)
#   coefficient   = 16 pi alpha_s / (9 m_q^2) * (S(S+1) - 3/2)
delta_sig = (sigma**3 / np.pi**1.5) * np.exp(-(sigma * r_nat)**2)
spin_factor = (S * (S + 1) - 1.5)
V_ss = (16.0 * np.pi * alpha_s) / (9.0 * m_q**2) * delta_sig * spin_factor

V_total = V_coul + V_lin + V_ss

fig, ax = plt.subplots(figsize=(6.0, 4.2))
ax.plot(r_fm, V_total, color='black', lw=2.0, label=r'Total $V(r)$')
ax.plot(r_fm, V_coul,  color='#0067A5', lw=1.4, ls='--', label=r'$-\dfrac{4}{3}\alpha_s/r$')
ax.plot(r_fm, V_lin,   color='#3A8C3A', lw=1.4, ls='-.', label=r'$b\,r$')
ax.plot(r_fm, V_ss,    color='#C1272D', lw=1.4, ls=':',  label=r'Spin-spin $(S\!=\!1)$')

ax.axhline(0, color='gray', lw=0.5)
ax.set_xlim(0, 2.0)
ax.set_ylim(-1.2, 1.2)
ax.set_xlabel(r'$r$ (fm)')
ax.set_ylabel(r'$V(r)$ (GeV)')
ax.legend(loc='lower right', frameon=False, fontsize=10)
ax.grid(alpha=0.2, lw=0.4)
plt.tight_layout()
plt.savefig('/home/claude/cornell.pdf', bbox_inches='tight')
plt.close()
print("cornell.pdf done")

# ---------------------------------------------------------------
# Figure 2: Charmonium spectrum (bar diagram, three columns S/P/D)
# ---------------------------------------------------------------
# Data: state, expt, ours, label
S_data = [
    ('1S0', 'eta_c(1S)',  2.9803, 3.0336),
    ('1S1', 'J/psi',      3.0969, 3.0843),
    ('2S0', 'eta_c(2S)',  3.6370, 3.6360),
    ('2S1', 'psi(2S)',    3.6861, 3.6598),
    ('3S1', 'psi(3S)',    4.0390, 4.0661),
    ('4S1', 'psi(4S)',    4.4210, 4.4088),
]
P_data = [
    ('1P3', 'chi_c(1P)',  3.5251, 3.5048),
    ('1P1', 'h_c(1P)',    3.5254, 3.5016),
]
D_data = [
    ('1D2', 'psi_2(1D)',  3.7729, 3.7809),
    ('2D2', 'psi_2(2D)',  4.1530, 4.1488),
]

fig, ax = plt.subplots(figsize=(6.6, 5.2))

labels = {
    'eta_c(1S)': r'$\eta_c(1S)$',
    'J/psi':     r'$J/\psi$',
    'eta_c(2S)': r'$\eta_c(2S)$',
    'psi(2S)':   r'$\psi(2S)$',
    'psi(3S)':   r'$\psi(3S)$',
    'psi(4S)':   r'$\psi(4S)$',
    'chi_c(1P)': r'$\langle\chi_c(1P)\rangle$',
    'h_c(1P)':   r'$h_c(1P)$',
    'psi_2(1D)': r'$\psi_2(1D)$',
    'psi_2(2D)': r'$\psi_2(2D)$',
}

def draw_column(x_center, half_w, data, color, label_offsets):
    for (tag, name, expt, ours), (dy, ha_side) in zip(data, label_offsets):
        ax.hlines(expt, x_center - half_w, x_center + half_w,
                  colors='black', lw=2.4, zorder=3)
        ax.hlines(ours, x_center - half_w, x_center + half_w,
                  colors=color, lw=2.4, zorder=2, alpha=0.85)
        y_lab = (expt + ours) / 2 + dy
        if ha_side == 'right':
            xt, ha = x_center + half_w + 0.04, 'left'
        else:
            xt, ha = x_center - half_w - 0.04, 'right'
        ax.text(xt, y_lab, labels.get(name, name), ha=ha, va='center',
                fontsize=9.5, color='black')

# Stagger labels so close-lying states do not overlap
S_offsets = [(-0.05, 'right'),  # eta_c(1S)
             ( 0.05, 'right'),  # J/psi
             (-0.05, 'right'),  # eta_c(2S)
             ( 0.05, 'right'),  # psi(2S)
             ( 0.00, 'right'),  # psi(3S)
             ( 0.00, 'right')]  # psi(4S)
P_offsets = [( 0.07, 'right'),  # chi_c(1P)
             (-0.07, 'right')]  # h_c(1P)
D_offsets = [( 0.00, 'right'),
             ( 0.00, 'right')]

draw_column(0.0, 0.35, S_data, '#C1272D', S_offsets)
draw_column(1.0, 0.35, P_data, '#0067A5', P_offsets)
draw_column(2.0, 0.35, D_data, '#3A8C3A', D_offsets)

ax.set_xticks([0, 1, 2])
ax.set_xticklabels([r'$S$ states', r'$P$ states', r'$D$ states'], fontsize=11)
ax.set_xlim(-0.7, 2.95)
ax.set_ylim(2.85, 4.55)
ax.set_ylabel(r'Mass (GeV)')
ax.set_title(r'Black bars = experiment $\;|\;$ Coloured bars = this work', fontsize=10)
ax.grid(axis='y', alpha=0.25, lw=0.4)
ax.set_axisbelow(True)
plt.tight_layout()
plt.savefig('/home/claude/spectrum.pdf', bbox_inches='tight')
plt.close()
print("spectrum.pdf done")

# ---------------------------------------------------------------
# Figure 3: J/psi mass convergence with N_max for several a
# ---------------------------------------------------------------
# Data taken from Belhouari Table III
N = np.array([300, 400, 500, 600, 700, 800, 900, 1000])
a10 = np.array([3.09558, 3.09557, 3.09556, 3.09556, 3.09556, 3.09556, 3.09556, 3.09556])
a30 = np.array([3.09600, 3.09575, 3.09565, 3.09561, 3.09559, 3.09558, 3.09558, 3.09557])
a50 = np.array([3.09738, 3.09640, 3.09600, 3.09582, 3.09572, 3.09567, 3.09563, 3.09561])

fig, ax = plt.subplots(figsize=(6.0, 4.2))
ax.plot(N, a10, 'o-', color='#C1272D', lw=1.6, ms=5, label=r'$a = 10$ fm')
ax.plot(N, a30, 's-', color='#0067A5', lw=1.6, ms=5, label=r'$a = 30$ fm')
ax.plot(N, a50, '^-', color='#3A8C3A', lw=1.6, ms=5, label=r'$a = 50$ fm')
ax.axhline(3.0969, color='black', ls=':', lw=1.0, label=r'expt $J/\psi$')

ax.set_xlabel(r'$N_{\max}$')
ax.set_ylabel(r'$M_{J/\psi}$ (GeV)')
ax.set_xlim(250, 1050)
ax.set_ylim(3.0950, 3.0980)
ax.legend(loc='upper right', frameon=False, fontsize=10)
ax.grid(alpha=0.25, lw=0.4)
plt.tight_layout()
plt.savefig('/home/claude/convergence.pdf', bbox_inches='tight')
plt.close()
print("convergence.pdf done")
