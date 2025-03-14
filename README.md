# CircLigase设计工具

## 项目简介
这是一个基于Streamlit的CircLigase设计工具，可以帮助研究人员设计和优化CircLigase反应的引物序列。

## 在线访问
访问我们的[在线应用](https://share.streamlit.io)（部署后更新链接）

## 项目结构
├── stem_out_v3_RNAfold.py    # 主程序（含二级结构预测）
├── requirements.txt         # Python依赖
├── packages.txt            # 系统依赖
└── README.md                # 项目文档

## 新增功能模块
### DNA二级结构预测模块
- 输入：用户自定义中间序列
- 处理：结合设计的5'/3'端形成完整序列
- 算法：ViennaRNA RNAfold (最小自由能算法)
- 输出：
  - 二级结构点括号表示
  - 热力学自由能
  - 序列结构可视化

## 技术栈更新
- 新增依赖：ViennaRNA 2.5.0+
- 系统调用：subprocess管理RNAfold进程
- 临时文件处理：tempfile模块

## 功能特点
- 序列设计：自动生成优化的茎区和连接序列
- 引物优化：全面的引物质量分析
- 结构预测：使用ViennaRNA进行DNA二级结构预测
- 可视化：结构预测结果的图形化展示

## 本地部署
1. 克隆仓库
```bash
git clone https://github.com/Who-oss/circligation-design.git
cd circligation-design
```

2. 安装依赖
```bash
# 安装系统依赖
sudo apt-get update
sudo apt-get install viennarna imagemagick

# 安装Python依赖
pip install -r requirements.txt
```

3. 运行应用
```bash
streamlit run stem_out_v3_RNAfold.py
```

## 技术栈
- Streamlit 1.31.1
- ViennaRNA 2.5.0+
- ImageMagick
- Python 3.8+ 