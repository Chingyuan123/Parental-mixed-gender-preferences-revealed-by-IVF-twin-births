import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# ======== Set font ========
all_fontsize = 7
plt.rcParams['font.family'] = 'Times new roman'
plt.rcParams['font.size'] = all_fontsize

# ======== read alpha data ========
alpha_path = 'alpha10_0304.xlsx'
alpha_data = pd.read_excel(alpha_path)
alpha_data_sorted = alpha_data.sort_values(by=['age_group', 'twr'])

age_group_order = ['16-29', '30-32', '33-34', '35-36', '37-50'] # label: top = 37-50, bottom = 16-29
color_map = plt.cm.get_cmap('tab10', len(age_group_order))

# ======== read beta data ========
beta_path = 'beta10_0304.xlsx'
beta_data = pd.read_excel(beta_path)
wealth_group_numbers = beta_data.index + 1
coefficients_b = beta_data['Coef.']
lower_error_b = coefficients_b - beta_data['95Lower']
upper_error_b = beta_data['95Upper'] - coefficients_b

# ======== Set picutre size ========
fig, axs = plt.subplots(1, 2, figsize=(17.8/2.54, 12.7/2.54))  # 17.8cm x 12.7cm

## ===== alpha plot (left) =====
for i, age_group in enumerate(age_group_order):
    group_data = alpha_data_sorted[alpha_data_sorted['age_group'] == age_group]
    coefficients = group_data['coef'].values
    ci_lengths = group_data['95Upper'].values - group_data['95Lower'].values

    axs[0].errorbar(coefficients, np.arange(len(coefficients)) + i * (len(coefficients) + 1),
                    xerr=ci_lengths / 2, fmt='o',
                    color=color_map(len(age_group_order) - i - 1),
                    label=age_group, markersize=2.5, elinewidth=1)

axs[0].set_xticks(np.arange(-0.06, 0.07, 0.01)) # tick= -0.06 to 0.07
y_positions = [i * (len(coefficients) + 1) + len(coefficients)/2 for i in range(len(age_group_order))]
axs[0].set_yticks(y_positions)
axs[0].set_yticklabels(age_group_order)
axs[0].set_ylabel('age')
axs[0].set_title(r'$\hat{\alpha}_{aw}$ in each age group', fontsize = all_fontsize)
axs[0].axvline(x=0, color='black', linestyle='--', linewidth=1)

for i in range(len(age_group_order) - 1):
    axs[0].axhline(y=i * (len(coefficients) + 1) + len(coefficients), color='grey', linestyle='-', linewidth=1)

axs[0].set_ylim(-1, y_positions[-1] + 5)
axs[0].set_facecolor('white')
for spine in axs[0].spines.values():
    spine.set_edgecolor('black')

## ===== beta plot (right) =====
axs[1].errorbar(coefficients_b, wealth_group_numbers,
                xerr=[lower_error_b, upper_error_b],
                fmt='o', color='navy', ecolor='navy', capsize=0, 
                elinewidth=1, markersize=2.5)

axs[1].set_yticks(ticks=range(1, 11))
axs[1].set_yticklabels(range(10, 0, -1))  # label: top = 10, bottom = 1
axs[1].set_xticks(ticks=[i / 100 for i in range(10, 41, 5)])
axs[1].set_ylabel('wealth group')
axs[1].set_title(r'$\hat{\beta}_w$ in wealth group', fontsize = all_fontsize)
axs[1].set_facecolor('white')
axs[1].grid(False)
axs[1].set_ylim(10.5, 0.5) # wealth group value max = 10.5, min = 0.5

for spine in axs[1].spines.values():
    spine.set_edgecolor('black')

## ===== save file =====
plt.tight_layout()
fig.savefig('Figure5.png', dpi=800)
plt.show()
