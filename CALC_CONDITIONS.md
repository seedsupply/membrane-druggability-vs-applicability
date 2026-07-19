# 確定計算条件（2026-07-13 検証）

## GPCR (20標的)
- 入力: vs_analysis/boltz2_outputs/{T}/... または boltz_outputs_{T}/...（1200化合物VS実行）
- **binder100_batch_out/ は使用禁止**（binderのみ別実行、値が違う）
- CB1R のみ: vs_analysis/results/CB1R_vs_results.csv
- pred_pKd = -affinity_pred_value
- 出力: gpcr_tail_FINAL_20.csv, gpcr_VS_FINAL_20.csv

## SLC (23標的 VS / 82標的 r)
- VS指標: slc_vs_all/*_vs_results.csv（列名 pred_pKd / predicted_pKd 混在）
- **r は別ソース**: vs_analysis/data/slc_merged_results.csv（binderのみ1391行/82標的）
- 出力: slc_VS_FINAL_23.csv

## 統一定義
- Sep = (binder_mean - nb_mean) / (binder_std + nb_std)   ※ std は ddof=1
- binder_top_margin = binder_max - nb_95pct
- multi-binder は除外しない
- VS成功 = AUROC >= 0.7

## 確定結果
- VS成功: GPCR 7/20, SLC 2/23
- Sep -> AUROC: GPCR 0.972, SLC 0.991
- Sep -> EF@1%: GPCR 0.922, SLC 0.573
- binder_top_margin -> EF@1%: GPCR 0.798, SLC 0.625
- nb_5pct -> EF@1%: 両クラスとも n.s.（下側裾仮説は棄却）
- 「SLCでSep-EFが弱い」= Sepレンジ圧縮による（GPCRをSep<0.55に絞ると0.655に低下）
- r（確定値, slc_merged_results.csv 由来）: SLC2A7 +0.729, SLC18A2 +0.709, SLC6A4 -0.216

## 過去の誤り
- 「OX1R は AUROC 0.088 で失敗」→ 誤り。正しくは 0.912 / EF 41.7 で成功。
  0.088 は論文 Table 2 の DRD3 の AP 値。
- CONFIRMED 辞書の Sep（ADRB1 1.501 等）は binder100_batch_out 混入版。
  正しくは gpcr_tail_FINAL_20.csv（ADRB1 1.466）。ただし順位・結論は不変。
