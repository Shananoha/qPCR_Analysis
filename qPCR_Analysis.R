rm(list = ls())

# 加载必要的包
library(dplyr)
library(car) # 用于Levene检验
library(broom) # 用于提取检验结果
library(rstatix) # 用于事后检验

# 声明全局变量
utils::globalVariables(c(
    "Target", "Sample", "ct", "line", "id", "Type", "Group",
    "sampleID", "targetName", "CT", "refCT", "deltaCt",
    "mean_deltaCt_control", "deltaDeltaCt", "relativeExpression",
    "Q1", "Q3", "IQR", "lower_bound", "upper_bound", "is_outlier",
    "."
))

# 1. 数据读取函数
read_data <- function(file_path) {
    if (!file.exists(file_path)) {
        stop("文件不存在，请检查文件路径。")
    }

    # 定义列名映射
    col_mapping <- list(
        csv = c(Target = "Target", Sample = "Sample", ct = "Cq"),
        txt = c(Target = "Target.Name", Sample = "Sample.Name", ct = "Cp")
    )

    # 根据文件类型读取数据
    file_type <- tools::file_ext(file_path)
    if (file_type == "csv") {
        data <- read.csv(file_path)
    } else if (file_type == "txt") {
        data <- read.delim(file_path, skip = 1)
    } else {
        stop("不支持的文件格式。请提供CSV或TXT文件。")
    }

    # 统一列名
    colnames(data) <- sapply(colnames(data), function(col) {
        if (col %in% col_mapping[[file_type]]) {
            names(col_mapping[[file_type]])[col_mapping[[file_type]] == col]
        } else {
            col
        }
    })

    # 检查数据完整性
    required_cols <- c("Target", "Sample", "ct")
    missing_cols <- setdiff(required_cols, colnames(data))
    if (length(missing_cols) > 0) {
        stop(paste("数据文件中缺少必要的列：", paste(missing_cols, collapse = ", ")))
    }

    if (!is.numeric(data$ct)) {
        tryCatch(
            {
                data$ct <- as.numeric(data$ct)
            },
            warning = function(w) {
                stop("ct列包含非数值数据，无法转换为数值型。")
            }
        )
    }

    # 添加文件信息
    attr(data, "file_type") <- file_type
    attr(data, "file_path") <- file_path

    return(data)
}

# 辅助函数：保存结果到文件
save_results <- function(data, file_name, suffix) {
    output_file <- paste0(file_name, "_", suffix, ".csv")
    write.csv(data, file = output_file, row.names = FALSE)
    message("结果已保存到：", output_file)
}

# 2. 数据预处理函数
preprocess_data <- function(data) {
    data <- data %>%
        mutate(line = row_number()) %>%
        arrange(Target, Sample) %>%
        select(Target, Sample, ct, line) %>%
        group_by(Target, Sample) %>%
        mutate(id = paste0(Sample, "_", row_number())) %>%
        ungroup() %>%
        mutate(
            Type = ifelse(tolower(Target) %in% c("actin", "graphd"), "Reference", "Test"),
            Group = ifelse(tolower(Sample) %in% c("nc", "con", "wt"), "control", "exper")
        ) %>%
        filter(!is.na(ct) & ct != 0)

    return(data)
}

# 3. 填充 refCT 和 refID 字段的函数
fill_ref_fields <- function(process_data, data) {
    process_data <- process_data %>%
        rowwise() %>%
        mutate(
            refCT = {
                ref_row <- data %>%
                    filter(Type == "Reference" & id == sampleID)
                if (nrow(ref_row) == 0) NA else ref_row$ct
            },
            refID = {
                ref_row <- data %>%
                    filter(Type == "Reference" & id == sampleID)
                if (nrow(ref_row) == 0) NA else ref_row$id
            }
        ) %>%
        ungroup()

    return(process_data)
}

# 4.1 从原始数据计算相对表达量的函数
calculate_relative_expression <- function(process_data) {
    results <- process_data %>%
        group_by(targetName) %>%
        mutate(
            deltaCt = CT - refCT,
            mean_deltaCt_control = mean(deltaCt[Group == "control"], na.rm = TRUE),
            deltaDeltaCt = deltaCt - mean_deltaCt_control,
            relativeExpression = 2^(-deltaDeltaCt)
        ) %>%
        ungroup()

    return(results)
}

# 4.2 从已有deltaCt计算相对表达量的函数
calc_rel_expr_from_deltact <- function(results) {
    results <- results %>%
        group_by(targetName) %>%
        mutate(
            mean_deltaCt_control = mean(deltaCt[Group == "control"], na.rm = TRUE),
            deltaDeltaCt = deltaCt - mean_deltaCt_control,
            relativeExpression = 2^(-deltaDeltaCt)
        ) %>%
        ungroup()

    return(results)
}

# 5. 使用 IQR 方法去除异常值的函数
remove_outliers_iqr <- function(results) {
    filtered_results <- results %>%
        group_by(targetName, Sample) %>%
        mutate(
            Q1 = quantile(deltaCt, 0.25, na.rm = TRUE), # 第一四分位数
            Q3 = quantile(deltaCt, 0.75, na.rm = TRUE), # 第三四分位数
            IQR = Q3 - Q1, # 四分位距
            lower_bound = Q1 - 1.5 * IQR, # 下界
            upper_bound = Q3 + 1.5 * IQR, # 上界
            is_outlier = deltaCt < lower_bound | deltaCt > upper_bound # 判断是否为异常值
        ) %>%
        filter(!is_outlier) %>% # 过滤掉异常值
        select(-Q1, -Q3, -IQR, -lower_bound, -upper_bound, -is_outlier) %>% # 删除临时列
        ungroup()

    return(filtered_results)
}

# 6. 显著性检验函数
perform_significance_tests <- function(results) {
    significance_results <- data.frame(
        Target = character(),
        Group1 = character(),
        Group2 = character(),
        p_value = numeric(),
        stringsAsFactors = FALSE
    )

    for (target in unique(results$targetName)) {
        target_data <- results %>% filter(targetName == target)
        sample_groups <- unique(target_data$Sample)

        if (length(sample_groups) == 2) {
            group1 <- target_data$relativeExpression[target_data$Sample == sample_groups[1]]
            group2 <- target_data$relativeExpression[target_data$Sample == sample_groups[2]]

            t_test_result <- t.test(group1, group2)
            p_value <- t_test_result$p.value

            significance_results <- rbind(significance_results, data.frame(
                Target = target,
                Group1 = sample_groups[1],
                Group2 = sample_groups[2],
                p_value = p_value
            ))
        }

        if (length(sample_groups) > 2) {
            anova_result <- aov(relativeExpression ~ Sample, data = target_data)
            tukey_result <- TukeyHSD(anova_result)
            tukey_p_values <- tukey_result$Sample[, "p adj"]

            for (i in 1:length(tukey_p_values)) {
                group_pair <- rownames(tukey_result$Sample)[i] %>%
                    strsplit("-") %>%
                    .[[1]]
                significance_results <- rbind(significance_results, data.frame(
                    Target = target,
                    Group1 = group_pair[1],
                    Group2 = group_pair[2],
                    p_value = tukey_p_values[i]
                ))
            }
        }
    }

    return(significance_results)
}

# 主函数
main <- function(fileName) {
    # 构建文件路径
    file_path <- file.path(getwd(), paste0(fileName, ".csv"))

    # 读取数据
    data <- read_data(file_path)

    # 数据预处理
    data <- preprocess_data(data)

    # 筛选 Type 为 Test 的行，提取 id, Target, ct, Type, Group 列
    process_data <- data %>%
        filter(Type == "Test") %>%
        select(Group, Sample, sampleID = id, targetName = Target, CT = ct)

    # 填充 refCT 和 refID 字段
    process_data <- fill_ref_fields(process_data, data)

    # 删除 ct 和 refCT 字段中的 NA 数据
    process_data <- process_data %>%
        filter(!is.na(CT) & !is.na(refCT))

    # 正常计算相对表达量
    results_normal <- calculate_relative_expression(process_data)

    # 保存结果
    save_results(results_normal, fileName, "result_normal")

    # 执行并保存显著性检验
    significance_results_normal <- perform_significance_tests(results_normal)
    save_results(significance_results_normal, fileName, "significance_normal")

    # 处理过滤数据
    results_filtered <- remove_outliers_iqr(results_normal) %>%
        calc_rel_expr_from_deltact()
    save_results(results_filtered, fileName, "result_filtered")

    significance_results_filtered <- perform_significance_tests(results_filtered)
    save_results(significance_results_filtered, fileName, "significance_filtered")

    # 打印结果
    print("正常计算结果：")
    print(results_normal)
    print("正常显著性检验结果：")
    print(significance_results_normal)
    print("过滤后计算结果：")
    print(results_filtered)
    print("过滤后显著性检验结果：")
    print(significance_results_filtered)
}

# 调用主函数
main("admin_test2")
