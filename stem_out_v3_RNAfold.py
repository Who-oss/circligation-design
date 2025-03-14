import streamlit as st
import random
from typing import List, Tuple, Dict
import re
import subprocess
import tempfile
import os
import matplotlib.pyplot as plt
import io
from PIL import Image
from pathlib import Path


# 设置页面配置
st.set_page_config(
    page_title="CircLigase设计工具",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# 移除 from Bio.SeqUtils import GC

# 添加自定义 GC 含量计算函数
def calculate_gc_content(sequence):
    """计算序列的 GC 含量百分比"""
    gc_count = sequence.count('G') + sequence.count('C')
    total = len(sequence)
    return (gc_count / total * 100) if total > 0 else 0

# 更新互补碱基字典，添加G-T配对
complement = {
    'A': 'T', 
    'T': 'A', 
    'C': 'G', 
    'G': 'C',
    # G-T wobble配对
    'GT': 'TG',
    'TG': 'GT'
}

def validate_dangling_length(five_length, three_length):
    """
    验证悬垂序列长度是否符合要求
    Args:
        five_length: 5'端悬垂长度
        three_length: 3'端悬垂长度
    Returns:
        errors: 错误信息列表
    """
    errors = []
    if not 3 <= five_length <= 6:
        errors.append("5'端悬垂长度需在3-6 nt之间")
    if not 3 <= three_length <= 10:
        errors.append("3'端悬垂长度需在3-10 nt之间")
    return errors

def generate_dangling_sequence(length, is_five_prime=True, at_rich=False, terminal_type="5'-G/3'-T"):
    """
    生成优化的悬垂序列
    Args:
        length: 悬垂序列长度
        is_five_prime: 是否为5'端悬垂
        at_rich: 是否AT富集
        terminal_type: 末端类型，可选 "5'-G/3'-T" 或 "5'-A/3'-G"
    """
    if length == 0:
        return ""
        
    # 设置末端碱基
    if terminal_type == "5'-G/3'-T":
        five_end = "G"
        three_end = "T"
    else:  # "5'-A/3'-G"
        five_end = "A"
        three_end = "G"
    
    def generate_inner_sequence(inner_length):
        if at_rich:
            return ''.join(random.choices(['A','T'], weights=[0.7, 0.3], k=inner_length))
        return ''.join(random.choices(['A','T','C','G'], k=inner_length))
    
    # 根据5'或3'端生成对应序列
    if is_five_prime:
        if length == 1:
            return five_end
        return five_end + generate_inner_sequence(length-1)
    else:
        if length == 1:
            return three_end
        return generate_inner_sequence(length-1) + three_end

def validate_stem_params(stem_length, gc_target):
    """参数验证（来源：知识库设计规则）"""
    errors = []
    if not 11 <= stem_length <= 26:
        errors.append("茎长度需在11-26 bp之间")
    if stem_length > 14 and gc_target < 50:
        errors.append("长茎(>14bp)建议GC含量≥50%")
    return errors

def generate_stem_pair(stem_length, gc_target):
    """生成严格互补的茎区对（5'和3'端）"""
    # 生成5'端茎序列
    stem_5 = []
    gc_count = 0
    for _ in range(stem_length):
        if random.random() < gc_target/100:
            nt = random.choice(['G','C'])
            gc_count += 1
        else:
            nt = random.choice(['A','T'])
        stem_5.append(nt)
    
    # 确保达到目标GC含量
    while (gc_count / stem_length) < (gc_target/100):
        pos = random.randint(0, stem_length-1)
        if stem_5[pos] in ['A','T']:
            stem_5[pos] = random.choice(['G','C'])
            gc_count += 1
    
    stem_5 = ''.join(stem_5)
    # 生成互补的3'端
    stem_3 = ''.join([complement[nt] for nt in stem_5[::-1]])
    return stem_5, stem_3

def is_base_pair(base1, base2):
    """检查两个碱基是否可以配对（包括G-T配对）"""
    if base1 + base2 in ['GT', 'TG']:  # G-T wobble配对
        return True
    return base1 == complement.get(base2)

def check_complementarity(seq1, seq2):
    """检查两个序列的互补程度（包括G-T配对）
    Returns:
        (matches, total): 匹配数和总碱基数
    """
    matches = sum(1 for a, b in zip(seq1, seq2) if is_base_pair(a, b))
    return matches, len(seq1)

def check_three_prime_end(sequence, length=5):
    """检查3'端序列质量（考虑G-T配对）"""
    three_prime = sequence[-length:]
    issues = []
    
    # 检查3'端是否以A结尾
    if three_prime[-1] == 'A':
        issues.append("3'端以A结尾可能导致错误引发")
    
    # 检查连续的G/C和G-T配对
    for i in range(len(three_prime)-2):
        window = three_prime[i:i+3]
        if 'GGG' in window or 'CCC' in window:
            issues.append("3'端含有连续的G/C序列")
        # 检查连续的G-T配对
        gt_count = sum(1 for j in range(len(window)-1) 
                      if window[j:j+2] in ['GT', 'TG'])
        if gt_count >= 2:
            issues.append("3'端含有过多的G-T配对")
    
    # 计算3'端GC含量
    gc_content = (three_prime.count('G') + three_prime.count('C')) / len(three_prime) * 100
    if gc_content < 40 or gc_content > 60:
        issues.append(f"3'端GC含量不适宜 ({gc_content:.1f}%)")
    
    return (len(issues) == 0, issues)

def check_hairpin(sequence, min_stem=4):
    """检查发夹结构（考虑G-T配对）"""
    for i in range(len(sequence)-min_stem):
        for j in range(i+min_stem+3, len(sequence)):  # +3为最小环大小
            stem_len = min(j-i, len(sequence)-j)
            seq1 = sequence[i:i+stem_len]
            seq2 = sequence[j:j+stem_len][::-1]  # 反向比对
            
            matches, total = check_complementarity(seq1, seq2)
            if matches >= total * 0.75:  # 允许25%错配
                return False, (f"发现潜在发夹结构: "
                             f"{seq1}...{seq2} "
                             f"({matches}/{total} 配对)")
    return True, ""

def check_primer_dimer(primer1, primer2, min_match=6):
    """检查引物二聚体（考虑G-T配对）"""
    def check_overlap(seq1, seq2, min_len=min_match):
        """检查两个序列的重叠区域是否形成二聚体"""
        for i in range(len(seq1)-min_len+1):
            for j in range(len(seq2)-min_len+1):
                subseq1 = seq1[i:i+min_len]
                subseq2 = seq2[j:j+min_len][::-1]  # 反向比对
                matches, total = check_complementarity(subseq1, subseq2)
                if matches >= total * 0.8:  # 允许20%错配
                    return False, (f"发现潜在二聚体: "
                                 f"{subseq1}...{subseq2} "
                                 f"({matches}/{total} 配对)")
        return True, ""
    
    # 检查3'端重叠（更严格）
    end_len = 8
    p1_end = primer1[-end_len:]
    p2_end = primer2[-end_len:]
    matches, total = check_complementarity(p1_end, p2_end[::-1])
    if matches >= 6:  # 3'端要求更严格
        return False, f"3'端可能形成二聚体 ({matches}/{total} 配对)"
    
    # 检查全长重叠
    ok, msg = check_overlap(primer1, primer2)
    if not ok:
        return False, msg
    
    # 检查自身重叠
    ok1, msg1 = check_overlap(primer1, primer1)
    if not ok1:
        return False, f"正向引物: {msg1}"
    ok2, msg2 = check_overlap(primer2, primer2)
    if not ok2:
        return False, f"反向引物: {msg2}"
    
    return True, ""

def check_junction_safety(junction, stem_seq):
    """
    严格检查连接序列的安全性
    Args:
        junction: 连接序列
        stem_seq: 茎序列
    Returns:
        bool: 是否安全
    """
    # 获取茎区前3bp及其所有可能的互补序列
    stem_3bp = stem_seq[:3]
    
    # 检查直接互补
    comp_3bp = ''.join(complement[nt] for nt in stem_3bp)
    
    # 检查所有可能的G-T配对组合
    unsafe_patterns = []
    # 添加标准互补序列
    unsafe_patterns.append(comp_3bp)
    
    # 生成所有可能的G-T配对变体
    for i in range(3):
        if stem_3bp[i] in 'GC':
            variant = list(comp_3bp)
            if stem_3bp[i] == 'G':
                variant[i] = 'T'
            elif stem_3bp[i] == 'C':
                variant[i] = 'G'
            unsafe_patterns.append(''.join(variant))
    
    # 检查连接序列是否包含任何不安全模式
    for pattern in unsafe_patterns:
        # 检查完全匹配
        if pattern in junction:
            return False
        # 检查部分匹配（2/3碱基互补）
        for i in range(len(junction)-2):
            subseq = junction[i:i+3]
            matches = sum(1 for a, b in zip(subseq, pattern) 
                        if is_base_pair(a, b))
            if matches >= 2:  # 如果有2个或更多碱基互补
                return False
    
    return True

def generate_junction(stem_5, stem_3, stem_gc):
    """生成安全的茎前连接序列"""
    max_attempts = 200  # 增加尝试次数
    
    for _ in range(max_attempts):
        # 高GC茎区强制AT富集
        if stem_gc > 60:
            junction = generate_dangling_sequence(3, at_rich=True)
        else:
            junction = generate_dangling_sequence(3)
        
        # 严格检查连接序列安全性
        if check_junction_safety(junction, stem_5) and \
           check_junction_safety(junction, stem_3):
            return junction
    
    # 如果无法生成安全序列，使用保守的AT富集序列
    return 'TAT'  # 更改默认安全序列为TAT

def is_complementary(seq1, seq2):
    """检查两个序列是否互补（考虑G-T配对）"""
    if len(seq1) != len(seq2):
        return False
    return all(is_base_pair(a, b) for a, b in zip(seq1, seq2))

def reverse_complement(sequence):
    """获取序列的反向互补序列，忽略非核苷酸字符"""
    complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    # 只处理核苷酸字符
    return ''.join(complement_dict[base] if base in complement_dict else base 
                  for base in reversed(sequence))

def calculate_tm(sequence):
    """计算引物的Tm值
    使用简化公式：4(G+C) + 2(A+T)
    """
    gc_count = sequence.count('G') + sequence.count('C')
    at_count = sequence.count('A') + sequence.count('T')
    return 4 * gc_count + 2 * at_count - 5

def analyze_primers(forward_primer, reverse_primer):
    """综合分析引物对的质量（更新版本）"""
    results = []
    
    # 基本参数检查
    if not (18 <= len(forward_primer) <= 30 and 18 <= len(reverse_primer) <= 30):
        results.append("引物长度应在18-30bp之间")
    
    # Tm值分析
    tm_f = calculate_tm(forward_primer)
    tm_r = calculate_tm(reverse_primer)
    tm_diff = abs(tm_f - tm_r)
    if tm_diff > 5:
        results.append(f"引物Tm值差异过大 ({tm_diff}°C)")
    if not (40 <= tm_f <= 75 and 40 <= tm_r <= 75):
        results.append("引物Tm值应在40-75°C范围内")
    
    # 3'端分析（包括G-T配对）
    f_three_ok, f_issues = check_three_prime_end(forward_primer)
    r_three_ok, r_issues = check_three_prime_end(reverse_primer)
    if not f_three_ok:
        results.extend([f"正向引物3'端: {issue}" for issue in f_issues])
    if not r_three_ok:
        results.extend([f"反向引物3'端: {issue}" for issue in r_issues])
    
    # 二级结构分析
    hairpin_ok, hairpin_msg = check_hairpin(forward_primer)
    if not hairpin_ok:
        results.append(f"正向引物发夹: {hairpin_msg}")
    hairpin_ok, hairpin_msg = check_hairpin(reverse_primer)
    if not hairpin_ok:
        results.append(f"反向引物发夹: {hairpin_msg}")
    
    # 二聚体分析
    dimer_ok, dimer_msg = check_primer_dimer(forward_primer, reverse_primer)
    if not dimer_ok:
        results.append(f"引物二聚体: {dimer_msg}")
    
    return results

def plot_structure_simple(sequence, structure):
    """使用matplotlib绘制简化的DNA二级结构图"""
    import matplotlib.pyplot as plt
    import numpy as np
    import io
    
    # 创建图形
    plt.figure(figsize=(10, 6))
    
    # 绘制序列
    x = np.arange(len(sequence))
    plt.plot(x, np.zeros_like(x), 'k-', alpha=0.2)
    
    # 绘制碱基
    for i, (base, struct) in enumerate(zip(sequence, structure)):
        color = {'A': 'green', 'T': 'red', 'G': 'orange', 'C': 'blue'}.get(base, 'black')
        y = 0 if struct == '.' else (0.5 if struct == '(' else -0.5)
        plt.text(i, y, base, ha='center', va='center', color=color, fontsize=10)
    
    # 绘制配对连接
    stack = []
    for i, char in enumerate(structure):
        if char == '(':
            stack.append(i)
        elif char == ')' and stack:
            start = stack.pop()
            plt.plot([start, i], [0.3, 0.3], 'k:', alpha=0.5)
    
    # 设置图形属性
    plt.axis('off')
    plt.title('DNA Secondary Structure')
    
    # 保存图形
    buf = io.BytesIO()
    plt.savefig(buf, format='png', dpi=300, bbox_inches='tight')
    plt.close()
    buf.seek(0)
    return buf

def find_matching_pair(structure, i):
    """找到配对的括号位置"""
    count = 1
    for j in range(i + 1, len(structure)):
        if structure[j] == '(':
            count += 1
        elif structure[j] == ')':
            count -= 1
            if count == 0:
                return j
    return None

def analyze_sequence_features(sequence):
    """增强的序列特征分析"""
    features = {
        'length': len(sequence),
        'gc_content': calculate_gc_content(sequence),
        'repeats': find_repeats(sequence),
        'motifs': find_special_motifs(sequence),
        'stability': estimate_stability(sequence)
    }
    return features

def find_repeats(sequence, min_length=4):
    """查找重复序列"""
    repeats = []
    for i in range(len(sequence)-min_length+1):
        for j in range(i+min_length, len(sequence)+1):
            subseq = sequence[i:j]
            if sequence.count(subseq) > 1:
                repeats.append((subseq, sequence.count(subseq)))
    return sorted(set(repeats), key=lambda x: (-len(x[0]), -x[1]))

def find_special_motifs(sequence):
    """查找特殊序列模式"""
    motifs = {
        'GC_runs': len(re.findall('(?:G|C){4,}', sequence)),
        'AT_runs': len(re.findall('(?:A|T){4,}', sequence)),
        'palindromes': find_palindromes(sequence)
    }
    return motifs

def estimate_stability(sequence):
    """估计序列稳定性"""
    # 简化的稳定性评估
    stability = {
        'gc_stability': calculate_gc_content(sequence),
        'terminal_stability': analyze_terminal_stability(sequence)
    }
    return stability

class SequenceValidationError(Exception):
    """序列验证错误类"""
    pass

def validate_sequence_input(sequence, context=""):
    """增强的序列验证"""
    try:
        if not sequence:
            raise SequenceValidationError("序列不能为空")
        
        invalid_chars = set(sequence) - set('ATCG')
        if invalid_chars:
            raise SequenceValidationError(
                f"发现无效字符: {', '.join(invalid_chars)}"
            )
        
        if len(sequence) < 10:
            raise SequenceValidationError(
                f"序列长度过短 ({len(sequence)} nt, 最小要求10 nt)"
            )
        
        # 添加更多验证规则...
        
        return True, None
        
    except SequenceValidationError as e:
        return False, f"{context}: {str(e)}"

# 主界面布局
st.title("🧬 辅助CircLigase设计工具 ")
st.caption("基于Analytical Chemistry (2021) DOI: 10.1021/acs.analchem.0c04668 的验证方案")

# 输入侧边栏
with st.sidebar:
    st.header("设计参数")
    terminal_type = st.radio("末端碱基对", 
                           ["5'-G/3'-T", "5'-A/3'-G"],
                           help="知识库推荐的高效末端组合")
    stem_length = st.slider("杂交茎长度 (bp)", 11, 26, 14,
                          help="短DNA：11-14 bp，长DNA：15-26 bp")
    stem_gc = st.slider("茎区GC含量 (%)", 30, 80, 55,
                      help="ΔG范围：-14.73至-19.95 kcal/mol @37°C")
    
    # 添加悬垂长度控制
    st.subheader("悬垂长度设置")
    dang_5_len = st.slider("5'端悬垂长度", 3, 6, 4, 
                          help="推荐范围：3-6 nt")
    dang_3_len = st.slider("3'端悬垂长度", 3, 10, 6, 
                          help="推荐范围：3-10 nt")

# 序列设计部分
st.header("1️⃣ 序列设计")

if st.button("生成优化设计"):
    # 参数验证
    errors = validate_stem_params(stem_length, stem_gc)
    dangling_errors = validate_dangling_length(dang_5_len, dang_3_len)
    errors.extend(dangling_errors)
    
    if errors:
        st.error("❗参数错误：\n" + "\n".join(errors))
        st.stop()
    
    # 生成茎区对
    stem_5, stem_3 = generate_stem_pair(stem_length, stem_gc)
    
    # 生成连接序列
    junction_5 = generate_junction(stem_5, stem_3, stem_gc)
    junction_3 = generate_junction(stem_3, stem_5, stem_gc)
    
    # 验证连接序列安全性并提供详细信息
    safe_5 = check_junction_safety(junction_5, stem_5) and \
             check_junction_safety(junction_5, stem_3)
    safe_3 = check_junction_safety(junction_3, stem_3) and \
             check_junction_safety(junction_3, stem_5)
    
    if not (safe_5 and safe_3):
        st.error("警告：检测到潜在的不安全连接序列")
        if not safe_5:
            st.write("5'端连接序列可能形成伪结构")
        if not safe_3:
            st.write("3'端连接序列可能形成伪结构")
        st.write("建议重新生成序列")
        st.stop()
    
    # 生成悬垂
    dang_5 = generate_dangling_sequence(dang_5_len, is_five_prime=True, 
                                      at_rich=stem_gc>60, terminal_type=terminal_type)
    dang_3 = generate_dangling_sequence(dang_3_len, is_five_prime=False, 
                                      at_rich=stem_gc>60, terminal_type=terminal_type)
    
    # 组装序列（不使用方括号）
    full_5 = f"{dang_5}{stem_5}{junction_5}"
    full_3 = f"{junction_3}{stem_3}{dang_3}"
    
    # 生成引物序列
    forward_primer = full_5
    reverse_primer = reverse_complement(full_3)
    
    # 为了显示目的，可以单独创建带标记的序列
    display_5 = f"{dang_5}[{stem_5}]{junction_5}"
    display_3 = f"{junction_3}[{stem_3}]{dang_3}"
    
    # 显示序列结构
    st.subheader("序列结构")
    st.code(f"5'端序列: {display_5}")
    st.code(f"3'端序列: {display_3}")
    
    # 存储生成的序列到session_state以供引物优化使用
    st.session_state.forward_primer = full_5
    st.session_state.reverse_primer = reverse_complement(full_3)
    st.session_state.display_5 = display_5
    st.session_state.display_3 = display_3

# 引物优化部分
st.header("2️⃣ 引物优化")

# 显示原始序列结构（如果存在）
if 'display_5' in st.session_state:
    st.subheader("原始序列结构")
    st.code(f"5'端序列: {st.session_state.display_5}")
    st.code(f"3'端序列: {st.session_state.display_3}")

    # 引物手动修改区域
    st.subheader("引物序列修改")
    col1, col2 = st.columns(2)
    
    with col1:
        modified_forward = st.text_input("修改正向引物",
                                       value=st.session_state.forward_primer,
                                       help="只允许输入ATCG碱基")
        
    with col2:
        modified_reverse = st.text_input("修改反向引物",
                                       value=st.session_state.reverse_primer,
                                       help="只允许输入ATCG碱基")
    
    # 验证输入的序列
    invalid_chars = set(modified_forward + modified_reverse) - set('ATCG')
    if invalid_chars:
        st.error(f"发现无效字符: {', '.join(invalid_chars)}")
    else:
        # 引物分析按钮
        if st.button("分析引物质量"):
            st.divider()
            st.subheader("引物分析结果")
            
            col1, col2 = st.columns(2)
            with col1:
                st.write("**正向引物**")
                st.code(f"5'-{modified_forward}-3'")
                st.write(f"长度: {len(modified_forward)}bp")
                st.write(f"Tm值: {calculate_tm(modified_forward)}°C")
                st.write(f"GC含量: {calculate_gc_content(modified_forward):.1f}%")
                
            with col2:
                st.write("**反向引物**")
                st.code(f"5'-{modified_reverse}-3'")
                st.write(f"长度: {len(modified_reverse)}bp")
                st.write(f"Tm值: {calculate_tm(modified_reverse)}°C")
                st.write(f"GC含量: {calculate_gc_content(modified_reverse):.1f}%")
            
            # 详细分析结果
            st.subheader("质量分析")
            
            # 3'端分析
            st.write("**3'端分析**")
            f_three_ok, f_issues = check_three_prime_end(modified_forward)
            r_three_ok, r_issues = check_three_prime_end(modified_reverse)
            
            if f_issues:
                st.warning("正向引物3'端问题：\n" + "\n".join(f"- {issue}" for issue in f_issues))
            if r_issues:
                st.warning("反向引物3'端问题：\n" + "\n".join(f"- {issue}" for issue in r_issues))
            if not f_issues and not r_issues:
                st.success("3'端序列符合要求")
            
            # 二级结构分析
            st.write("**二级结构分析**")
            f_hairpin_ok, f_hairpin_msg = check_hairpin(modified_forward)
            r_hairpin_ok, r_hairpin_msg = check_hairpin(modified_reverse)
            dimer_ok, dimer_msg = check_primer_dimer(modified_forward, modified_reverse)
            
            if not f_hairpin_ok:
                st.warning(f"正向引物: {f_hairpin_msg}")
            if not r_hairpin_ok:
                st.warning(f"反向引物: {r_hairpin_msg}")
            if not dimer_ok:
                st.warning(dimer_msg)
            if f_hairpin_ok and r_hairpin_ok and dimer_ok:
                st.success("未检测到明显的二级结构")
            
            # 综合建议
            st.subheader("优化建议")
            primer_issues = analyze_primers(modified_forward, modified_reverse)
            if primer_issues:
                st.warning("**建议修改：**\n" + "\n".join(f"- {issue}" for issue in primer_issues))
            else:
                st.success("✅ 引物设计符合所有要求")
else:
    st.info("请先点击'生成优化设计'按钮生成基础序列")

# 二级结构预测模块
st.header("3️⃣ DNA二级结构预测")

# 显示当前设计的序列
if 'display_5' in st.session_state and 'display_3' in st.session_state:
    st.subheader("当前序列结构")
    col1, col2 = st.columns(2)
    with col1:
        st.code(f"5'端: {st.session_state.display_5}")
    with col2:
        st.code(f"3'端: {st.session_state.display_3}")

# 中间序列输入
middle_seq = st.text_area("输入中间序列（ATCG字符）", 
                        height=100,
                        help="请输入连接5'和3'端的中间序列（10-100 nt）",
                        max_chars=100).upper()

# 添加预测参数
with st.expander("高级预测参数"):
    col1, col2 = st.columns(2)
    with col1:
        temperature = st.slider("温度 (°C)", 20, 60, 37)
    with col2:
        na_type = st.radio("核酸类型", ["DNA", "RNA"])

def run_rnafold(sequence, temp=37, dtype="DNA"):
    """使用RNAfold预测二级结构"""
    with tempfile.TemporaryDirectory() as tmpdir:
        # 准备输入文件
        input_file = os.path.join(tmpdir, "input.fa")
        with open(input_file, "w") as f:
            f.write(f">sequence\n{sequence}\n")
        
        # 构建命令参数
        cmd = ["RNAfold", "--noPS", "-T", str(temp)]
        if dtype == "DNA":
            cmd += ["--noconv", "-d2"]  # DNA模式参数
        
        try:
            # 运行RNAfold
            result = subprocess.run(
                cmd, 
                cwd=tmpdir,
                input=f"{sequence}\n",
                capture_output=True,
                text=True,
                timeout=30
            )
            
            # 解析输出
            if result.returncode == 0:
                output = result.stdout.split('\n')
                structure = output[1].split()[0]
                energy = output[1].split()[-1].strip('()')
                return structure, energy, None
            else:
                return None, None, f"RNAfold错误: {result.stderr}"
                
        except Exception as e:
            return None, None, str(e)

def plot_structure_vienna(sequence, structure):
    """使用ViennaRNA包绘制DNA结构图"""
    try:
        # 首先尝试使用matplotlib绘制
        return {'png': plot_structure_simple(sequence, structure).getvalue()}
    except Exception as e:
        st.error(f"结构图生成失败: {str(e)}")
        return None

if st.button("预测二级结构"):
    # 验证中间序列
    if not re.match("^[ATCG]+$", middle_seq):
        st.error("中间序列包含无效字符，只允许ATCG")
        st.stop()
    
    if len(middle_seq) < 10:
        st.error("中间序列过短（至少10 nt）")
        st.stop()
    
    # 组装完整序列
    full_5 = st.session_state.forward_primer
    full_3 = st.session_state.reverse_primer
    complete_seq = f"{full_5}{middle_seq}{reverse_complement(full_3)}"
    
    # 运行预测
    with st.spinner("正在运行RNAfold预测..."):
        structure, energy, error = run_rnafold(
            complete_seq,
            temp=temperature,
            dtype=na_type
        )
    
    if error:
        st.error(f"预测失败: {error}")
    else:
        st.subheader("预测结果")
        col1, col2 = st.columns([1, 2])
        
        with col1:
            st.markdown("**基本参数**")
            st.write(f"总长度: {len(complete_seq)} nt")
            st.write(f"自由能: {energy} kcal/mol")
            st.write(f"温度: {temperature}°C")
            
        with col2:
            st.markdown("**二级结构表示**")
            st.code(f"5' {structure} 3'")
            
            # 使用ViennaRNA生成结构图
            result = plot_structure_vienna(complete_seq, structure)
            if result:
                # 显示PNG图像
                st.image(result['png'], caption="DNA二级结构示意图", use_column_width=True)
                
                # 如果有PS文件，提供下载选项
                if result['ps']:
                    st.download_button(
                        "下载高清结构图(PS格式)",
                        result['ps'],
                        file_name="structure.ps",
                        mime="application/postscript"
                    )
            else:
                # 如果生成失败，使用matplotlib方案作为备选
                st.warning("结构图生成失败，使用简化图形显示")
                img_data = plot_structure_simple(complete_seq, structure)
                st.image(img_data, caption="DNA二级结构示意图", use_column_width=True)

# 帮助信息
with st.expander("设计原理说明"):
    st.markdown("""
    **核心设计规则（来自 Analytical Chemistry 2021）**：
    1. **末端互补强制配对**  
       通过生成完全互补的茎区序列，确保CircLigase有效连接（图2C）
    2. **茎前序列去互补化**  
       禁止连接区出现与茎区前3bp互补的序列，防止伪结形成（图S4）
    3. **动态AT富集策略**  
       高GC茎区自动生成AT富集的连接序列，增强单链灵活性（ΔG优化）
    4. **悬垂长度控制**  
       5'端3-6nt，3'端3-10nt悬垂，确保酶识别位点（图1B）
    """)

with st.expander("二级结构预测说明"):
    st.markdown("""
    **预测方法说明（基于RNAfold算法）**：
    1. **能量模型**：使用Nearest Neighbor热力学参数
    2. **DNA参数**：采用SantaLucia 2004 DNA热力学参数集
    3. **预测算法**：Zuker最小自由能算法
    4. **输出解释**：
       - 括号表示碱基配对（( ) = 标准配对，< > = G-T配对）
       - 自由能单位：kcal/mol（负值越大越稳定）
    """)

def create_sequence_input_section():
    """优化的序列输入界面"""
    st.subheader("序列输入")
    
    input_method = st.radio(
        "选择输入方式",
        ["手动输入", "文件上传", "示例序列"]
    )
    
    if input_method == "手动输入":
        sequence = st.text_area(
            "输入序列",
            help="请输入DNA序列（ATCG）",
            max_chars=1000
        )
    elif input_method == "文件上传":
        uploaded_file = st.file_uploader(
            "上传FASTA文件",
            type=['fa', 'fasta']
        )
        if uploaded_file:
            sequence = process_fasta_file(uploaded_file)
    else:
        sequence = load_example_sequence()
    
    return sequence.upper() if sequence else ""

def add_advanced_settings():
    """高级设置面板"""
    with st.expander("高级设置"):
        col1, col2 = st.columns(2)
        with col1:
            temperature = st.slider(
                "温度 (°C)",
                min_value=20,
                max_value=60,
                value=37,
                help="预测温度设置"
            )
        with col2:
            salt_conc = st.slider(
                "盐浓度 (mM)",
                min_value=0,
                max_value=1000,
                value=150,
                help="Na+浓度设置"
            )
    return temperature, salt_conc

def export_results(sequence, structure, analysis):
    """结果导出功能"""
    import json
    import csv
    from datetime import datetime
    
    # 创建导出数据
    export_data = {
        'timestamp': datetime.now().isoformat(),
        'sequence': sequence,
        'structure': structure,
        'analysis': analysis
    }
    
    # 提供多种导出格式
    export_format = st.selectbox(
        "选择导出格式",
        ["CSV", "JSON", "PDF"]
    )
    
    if st.button("导出结果"):
        if export_format == "CSV":
            return export_as_csv(export_data)
        elif export_format == "JSON":
            return export_as_json(export_data)
        else:
            return export_as_pdf(export_data)