import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import sys

def load_data_from_folder(folder_path):
    all_data = []

    # 遍历文件夹中的所有CSV文件
    for file_name in os.listdir(folder_path):
        if file_name.endswith('.csv'):
            # 提取软件名称，假设文件名格式为 softwareName_xxx.csv
            software_name = file_name.split('_')[0].capitalize()  # 首字母大写
            file_path = os.path.join(folder_path, file_name)

            # 读取CSV文件，添加Software列
            df = pd.read_csv(file_path)
            df['Software'] = software_name

            # 选择相关的列，保持一致性
            df = df[['Software', 'Database', 'Taxonomic Rank', 'Accuracy', 'Precision', 'Recall', 'F1 Score']]

            all_data.append(df)

    # 将所有数据合并为一个DataFrame
    return pd.concat(all_data, ignore_index=True)

def plot_performance(data, folder_path):
    # 确定数据库、分类级别（Taxonomic Rank）和性能指标
    levels = ['Superkingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
    databases = data['Database'].unique()
    metrics = ['Accuracy', 'Precision', 'Recall', 'F1 Score']

    # 设置图形风格
    sns.set_theme(style="ticks")
    sns.set_context("paper")  # 根据论文要求调整字体大小
    sns.set_palette(sns.color_palette("deep"))  # 使用专业的调色板

    # 设置字体
    plt.rcParams.update({'font.family':'Arial', 'font.size':14})

    # 设置图形网格布局
    fig, axes = plt.subplots(len(levels), len(databases), figsize=(20, 28), sharey=False)

    # 绘制每个子图
    for j, database in enumerate(databases):
        for i, level in enumerate(levels):
            # 筛选数据
            subset = data[(data['Taxonomic Rank'] == level) & (data['Database'] == database)]

            # 计算每个软件在各指标下的平均值
            mean_values = subset.groupby('Software')[metrics].mean().reset_index()

            # 如果数据为空，跳过此子图
            if mean_values.empty:
                axes[i, j].axis('off')
                continue

            # 将 F1 Score 四舍五入到两位小数用于排序
            mean_values['F1 Score Rounded'] = mean_values['F1 Score'].round(2)

            # 按照 F1 Score（两位小数）和软件名称排序
            mean_values = mean_values.sort_values(by=['F1 Score Rounded', 'Software'], ascending=[False, True])

            # 更新软件顺序
            software_order = mean_values['Software']

            # 将数据转换为长格式
            melted_data = pd.melt(mean_values, id_vars=['Software'], value_vars=metrics,
                                  var_name='Metric', value_name='Value')

            ax = axes[i, j]
            sns.barplot(x='Software', y='Value', hue='Metric', data=melted_data, ax=ax,
                        order=software_order, hue_order=metrics, edgecolor='black')

            ax.set_ylim(0.0, 1.0)  # 固定纵坐标范围

            if i == len(levels) - 1:
                # 在每列最下方添加数据库名称
                ax.set_xlabel(database, fontsize=16, fontweight='bold')
                ax.xaxis.set_label_coords(0.5, -0.35)
            else:
                ax.set_xlabel('')

            if j == 0:
                ax.set_ylabel(level, fontsize=16, fontweight='bold')
            else:
                ax.set_ylabel('')

            # 增大x轴标签（软件名称）字体大小，并旋转以避免重叠
            ax.tick_params(axis='x', rotation=45, labelsize=14)
            ax.tick_params(axis='y', labelsize=14)
            ax.get_legend().remove()  # 去除每个子图的图例

            # 去除顶部和右侧边框
            sns.despine(ax=ax)

            # 添加数值标签，过滤掉高度过小的柱子
            for p in ax.patches:
                height = p.get_height()
                if not pd.isna(height) and height >= 0.01:
                    ax.text(p.get_x() + p.get_width() / 2.,
                            height + 0.02,
                            f'{height:.2f}',
                            ha="center", fontsize=12)
                else:
                    continue

    # 调整图例字体大小，并放在图表顶部
    handles, labels = ax.get_legend_handles_labels()
    legend = fig.legend(handles, labels, loc='upper center', ncol=4, fontsize=16, frameon=True, framealpha=1, edgecolor='black')
    legend.get_frame().set_linewidth(1.0)

    # 调整图例与子图之间的距离，减少空白
    plt.subplots_adjust(top=0.95)

    # 调整布局
    plt.tight_layout(rect=[0, 0.05, 1, 0.92], h_pad=2)

    # 保存图表为高分辨率PDF
    output_pdf = os.path.join(folder_path, 'metrics_plot.pdf')
    plt.savefig(output_pdf, format='pdf', bbox_inches='tight', dpi=300)

    # 显示图表
    plt.show()

if __name__ == "__main__":
    # 获取命令行输入的文件夹路径
    if len(sys.argv) != 2:
        print("Usage: python script_name.py <folder_path>")
        sys.exit(1)

    folder_path = sys.argv[1]

    if not os.path.isdir(folder_path):
        print(f"Error: {folder_path} is not a valid directory.")
        sys.exit(1)

    # 加载数据并绘图
    data = load_data_from_folder(folder_path)
    plot_performance(data, folder_path)
