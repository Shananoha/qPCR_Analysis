import pandas as pd
import numpy as np


def load_data(file_path):
    """读取qPCR数据文件"""
    import os
    if not os.path.exists(file_path):
        raise FileNotFoundError(f'数据文件未找到：{file_path}\n请检查：\n1. 文件路径是否正确\n2. 当前工作目录是否为项目根目录')
    return pd.read_csv(file_path)


def preprocess(data, ref_gene, remove_outliers=True):
    """从数据预处理到最终筛选的全部逻辑"""
    # 数据预处理
    processed_rows = []

    # 分离参考基因和目标基因数据
    ref_data = data[data['Target'] == ref_gene][['Sample', 'Cq']]
    target_data = data[data['Target'] != ref_gene][['Sample', 'Target', 'Cq']]

    # 按目标基因分组处理
    for target, target_group in target_data.groupby('Target'):
        # 按样本分组
        for sample, sample_group in target_group.groupby('Sample'):
            # 获取当前样本的参考基因Ct值列表
            ref_ct = ref_data[ref_data['Sample'] == sample]['Cq'].tolist()
            target_ct = sample_group['Cq'].tolist()

            # 对齐数据长度并生成记录
            max_length = max(len(ref_ct), len(target_ct))
            for i in range(max_length):
                processed_rows.append({
                    'groupName': sample,
                    'targetName': target,
                    'Ctt': target_ct[i] if i < len(target_ct) else 'none',
                    'refName': ref_gene,
                    'Ctf': ref_ct[i] if i < len(ref_ct) else 'none',
                    'control_group': sample in ['NC', 'CON', 'CTRL']
                })

    # 创建最终DataFrame并清洗数据
    df_pre_process = pd.DataFrame(processed_rows)

    # 清理空值（将'none'替换为NaN后删除）
    df_pre_process.replace('none', pd.NA, inplace=True)
    df_pre_process.dropna(subset=['Ctt', 'Ctf'], how='any', inplace=True)

    # 计算ΔCT值
    df_process = df_pre_process.copy()
    df_process['deltaCT'] = df_process['Ctt'] - df_process['Ctf']

    # IQR 过滤异常值
    df_clear = df_process.copy()
    df_clear['control_group'] = df_clear['control_group'].astype('bool')

    # 对每个目标基因计算IQR范围并过滤异常值
    df_clear = df_clear.groupby('targetName', group_keys=False).apply(lambda x:
                                                    x[ (x.control_group & (
                                                                (x.deltaCT >= (x[x.control_group].deltaCT.quantile(0.25) - 1.5 * (
                                                        x[x.control_group].deltaCT.quantile(0.75) - x[x.control_group].deltaCT.quantile(0.25)))) &
                                                          (x.deltaCT <= (x[x.control_group].deltaCT.quantile(0.75) + 1.5 * (
                                                                  x[x.control_group].deltaCT.quantile(0.75) - x[x.control_group].deltaCT.quantile(0.25))))
                                                        ) | ~x.control_group ) ].reset_index(
        drop=True))

    # 计算相对表达量并标准化处理
    df_refExpr = df_clear.copy()
    df_refExpr['rel_expr'] = np.power(2, -df_refExpr['deltaCT'])
    control_mean = df_refExpr[df_refExpr.control_group] \
        .groupby('targetName')['rel_expr'] \
        .mean() \
        .rename('control_mean')
    df_refExpr = df_refExpr.merge(control_mean, on='targetName')
    df_refExpr['norm_expr'] = df_refExpr['rel_expr'] / df_refExpr['control_mean']

    print(f'执行非对照组IQR筛选，参数状态：remove_outliers={remove_outliers}')
    print(f'执行非对照组IQR筛选，参数状态：remove_outliers={remove_outliers}')
    print(f'筛选前数据量：{len(df_refExpr)}行')
    if remove_outliers:
        # 对非对照组数据进行IQR异常值筛选
        pre_count = len(df_refExpr)
        
        df_final = df_refExpr.groupby(['targetName', 'groupName'], group_keys=False).apply(lambda x:
            x[
                (x.control_group) | (
                    ~x.control_group & (
                        (x.norm_expr >= (
                            x[~x.control_group].norm_expr.quantile(0.25) -
                            1.5 * (
                                x[~x.control_group].norm_expr.quantile(0.75) -
                                x[~x.control_group].norm_expr.quantile(0.25)
                            )
                        )) & (
                            x.norm_expr <= (
                                x[~x.control_group].norm_expr.quantile(0.75) +
                                1.5 * (
                                    x[~x.control_group].norm_expr.quantile(0.75) -
                                    x[~x.control_group].norm_expr.quantile(0.25)
                                )
                            )
                        )
                    )
                )
            ]
        ).reset_index(drop=True)
    else:
        df_final = df_refExpr

    print(f'筛选后数据量：{len(df_final)}行 (减少{pre_count - len(df_final)}行)')
    
    return {
        'final_data': df_final,
        'statistics': df_final.describe(),
        'filter_applied': remove_outliers
    }


if __name__ == '__main__':
    # 数据导入
    df_row = load_data('testData/Biorad_test.csv')

    # 调用预处理函数
    ref_gene = 'Actin'  # 参考基因名称
    df_final = preprocess(df_row, ref_gene, remove_outliers=True)['final_data']

    print('最终数据列：', df_final.columns.tolist())
    print('数据样本：\n', df_final.head(2))

    # 输出最终结果
    df_final.to_csv('Biorad_final.csv', index=False, float_format='%.2f')