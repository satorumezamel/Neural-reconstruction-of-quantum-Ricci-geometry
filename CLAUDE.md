# IDPC_Reproduction — Claude Instructions

## プロジェクト概要
EEG由来のNeural Ricci曲率と量子プロセッサ由来のQuantum Ricci曲率の
構造的対応を解析する研究リポジトリ。

論文: Intersection-Defined Phase Coordinates Reveal Localized Selection and a Non-Closed Observational Structure
著者: Satoru Watanabe (SIEL)
DOI: https://doi.org/10.5281/zenodo.19628769

## 環境
- Python 3.x
- 依存関係: requirements.txt参照
- メインノートブック: IDPC_Reproduction.ipynb
- デモノートブック: IDPC_Repro_Demo.ipynb

## データ構造
- IDPC_Reproduction/ : 被験者P1-P26のEEG・量子時系列CSV
  - PX_eeg_timeseries.csv : EEG時系列
  - PX_quantum_timeseries.csv : 量子時系列
  - PX_co_recon_features_W30.csv : co-reconstruction特徴量
  - PX_co_recon_features_W30_points.csv : 同、ポイント版
  - 集計CSVファイル各種
- IDPC_Reproduction_ricci/ : Ricci曲率時系列CSV（P1-P26）

## 解析を実行するとき
1. pip install -r requirements.txt で依存関係をインストール
2. IDPC_Reproduction.ipynb をJupyterで実行
3. 結果はIDPC_Reproduction/以下のCSVに出力される

## 主要な変数・指標
- ERicci(t): Neural Ricci曲率（EEGから算出）
- QRicci(t): Quantum Ricci曲率（量子測定から算出）
- φ(x,t): intersection variable（内部自由度）
- switch_gain: 対応構造の強度指標（観測値≈0.810, null≈0.565, p≈0.005）
- AUC ≈ 0.79, Accuracy ≈ 0.70: 状態分類性能
- mean abs Pearson ≈ 0.42: セッション間再構成一致度

## レポート出力ルール
- Markdown形式で出力
- 論文のセクション構造に従う（Abstract, Results, Discussion）
- 数値は小数点2桁、統計量はp値・効果量を必ず記載
- 図はmatplotlibで生成してレポートに埋め込む
- 再現解析の場合は元論文の数値と比較して差異を明記

## 追加実験・再現解析の指示例
- 「セッションP1-P10で再現解析を実行してレポートを出力して」
- 「switch_gainの閾値を0.7に変えて感度分析して」
- 「新しいセッションデータをZenodoからダウンロードして解析して」
