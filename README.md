# MicroNorm

**[Open MicroNorm](https://zbyloveyj.github.io/MicroNorm/)** · English / 中文 · Interface v0.2

A browser-local microbiome reference workbench with six pages:

- **Overview:** reference sample counts and analysis entry points.
- **Individual analysis:** complete profile input, M1–M3 coordinates, CMD and all 801 species percentiles, with search and exports.
- **Cohort comparison:** up to 250 samples in two groups, per-sample results, descriptive plots and eligible background-matched comparisons.
- **Reference atlas:** detection and positive-abundance curves from the deployed empirical model, age support and 473 frozen module memberships.
- **Methods:** equations, input requirements, matching rules, reference support and model scope.
- **Downloads:** complete synthetic examples, input templates, model parameters, membership tables, guides, offline web and Windows packages.

The first visit defaults to English. The language control remembers your selection and preserves calculated results while switching.

## Downloads

[Offline web package](downloads/MicroNorm_Web_v02.zip) · [Windows package](downloads/MicroNorm_Windows_v02.zip)

Unzip the web package and open `index.html` in a modern browser. For Windows, unzip the whole folder and run `MicroNorm.exe` (Microsoft Edge required), or open `app/index.html`. Application scoring works offline; links to GitHub and online downloads require an internet connection.

## Input

Single sample: a native MetaPhlAn species-level profile, or `species,abundance_percent` columns. Supply a complete detected species profile in percent and enter age (0–100) and recorded sex (female/male). Unlisted species are treated as zero.

Cohort: `sample_id,age,sex,group,study,species,abundance_percent`. One species per sample per row. `study` is optional for scoring and required for matching. At most two groups and 250 samples per batch.

The parser does not verify the upstream MetaPhlAn/database release. Format compatibility alone does not establish model compatibility.

## Model

Frozen model: `empirical_two_part_preview_v1`. It contains 11,105 broadly healthy first-visit reference samples, 801 species and three modules with 229, 138 and 106 members. The interface update does not retrain or replace the model.

This is an empirical age-bin research preview. Atlas shading denotes positive-abundance quantiles, not confidence intervals. Scores and CMD are reference positions, not clinical diagnoses or disease probabilities. See [model card](model-card.txt) and the Methods page for scope and validation limits.

## Privacy

Sample files are read and scored in the visitor's browser; the application does not submit sample files to a server. This repository contains the aggregated reference model and synthetic examples, not underlying participant profiles.

## v0.2 verification

Automated DOM tests covered all six pages in both languages, language switching with retained results, synthetic and file-input scoring, 801-row exports, cohort reports, reference-curve CSV, 473 memberships, invalid-input handling and deterministic 30-pair comparisons. Browser screenshot verification was unavailable. See [verification report](verification-v02.json).

## 中文说明

网站包含概览、单样本分析、队列比较、参考图谱、方法说明和下载中心。首次访问默认英文，右上角可切换中文。示例均为合成数据。离线网页包与 Windows 包使用同一参考模型和评分引擎；完整解压后即可使用。
