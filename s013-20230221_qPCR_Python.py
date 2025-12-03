```python

import pandas as pd
import os

# Use relative path or current directory instead of hardcoded absolute path
# os.chdir("C:/Users/sxl/Desktop")  # Commented out for portability

# 读取 Excel 文件
data = pd.read_excel('qPCR_demo.xlsx')

# 获取每组内参基因的 Ct 均值
ck_tub2_ct = data.loc[(data['Group'] == 'ck') & (data['Gene'] == 'TUB2'), 'Ct'].mean()
treat_tub2_ct = data.loc[(data['Group'] == 'treat') & (data['Gene'] == 'TUB2'), 'Ct'].mean()

# 计算每组的 ∆Ct 值 - Optimized using vectorized operations instead of apply()
# Create a mapping dictionary and use vectorized map operation
ct_reference = {'ck': ck_tub2_ct, 'treat': treat_tub2_ct}
data['∆Ct'] = data['Ct'] - data['Group'].map(ct_reference)

data['∆Ct']

# Optimized: Combine groupby operations to process both groups at once
# For large datasets, this reduces memory overhead and improves performance
# For small datasets, the difference may be negligible
delta_ct_means = data.groupby(['Group', 'Gene'], as_index=False)['∆Ct'].mean()
ck_ΔCt_means = delta_ct_means[delta_ct_means['Group'] == 'ck'][['Gene', '∆Ct']].reset_index(drop=True)
treat_ΔCt_means = delta_ct_means[delta_ct_means['Group'] == 'treat'][['Gene', '∆Ct']].reset_index(drop=True)

# 输出 ck_ΔCt_means
print(ck_ΔCt_means)

# 输出 treat_ΔCt_means
print(treat_ΔCt_means)


# 将每个基因在ck组的ΔCt均值保存到字典中
ck_dict = dict(zip(ck_ΔCt_means['Gene'], ck_ΔCt_means['∆Ct']))

ck_dict

# 计算每个样本的 ∆∆Ct 值 - Optimized using vectorized map operation instead of apply()
data['∆∆Ct'] = data['∆Ct'] - data['Gene'].map(ck_dict)
data

# 计算每个样本的相对表达量
data['Rel Exp'] = 2 ** (-data['∆∆Ct'])
data

#relative_expression = data[['Rel Exp']].copy()

# 显示 DataFrame
# print(relative_expression)

# 获取每个基因在treat组的相对表达量均值
#treat_rel_exp_means = data.loc[data['Group'] == 'treat', ['Gene', 'Rel Exp']].groupby('Gene').mean().reset_index()
#treat_rel_exp_means

# 计算每个样本相对表达量的均值 - Use as_index=False consistently
mean_rel_exp = data.groupby(['Gene', 'Group'], as_index=False)['Rel Exp'].mean()

# 更改列名
mean_rel_exp = mean_rel_exp.rename(columns={'Rel Exp': 'mean_Rel_Exp'})

mean_rel_exp

# 保存结果到Excel文件
with pd.ExcelWriter('output4.xlsx') as writer:
    ck_ΔCt_means.to_excel(writer, sheet_name='ck_ΔCt_means', index=False)
    treat_ΔCt_means.to_excel(writer, sheet_name='treat_ΔCt_means', index=False)
    data.to_excel(writer, sheet_name='data', index=False)
    # 保存 mean_rel_exp，并将其添加到 data 的 sheet 中
    mean_rel_exp.to_excel(writer, sheet_name='data', index=False, startcol=6)


# In[18]:


### 绘制 log2meanRelExp的柱形图

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Read in data
data = pd.read_excel('output4.xlsx', sheet_name='Sheet1')

# Set plot style and font scale
sns.set(style='whitegrid', font_scale=1.5)

# Set figure size
fig, ax = plt.subplots(figsize=(10, 7))

# Create bar plot
sns.barplot(x='Gene', y='log2mean_Rel_Exp', data=data, palette='viridis', ax=ax)

# Set x-axis label rotation
plt.xticks(rotation=90)

# Set y-axis label
ax.set_ylabel('log2mean_Rel_Exp')

# Set title
plt.title('lncRNA expression')

# Save figure
plt.savefig('gene_expression.png', bbox_inches='tight')

# Show plot
plt.show()


# In[ ]:

```


