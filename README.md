# MFGRN: Multi-Method Fusion for Gene Regulatory Network Inference

## Overview

MFGRN (Multi-Method Fusion Gene Regulatory Network Inference) is a consensus-based framework for gene regulatory network (GRN) reconstruction that integrates multiple deep learning models to improve network inference accuracy, robustness, and biological reliability.

The framework combines three complementary GRN inference approaches:

- **3DCEMA** – 3D convolutional neural network based GRN inference
- **DeepSEM** – Variational autoencoder and structural equation model based GRN inference
- **DeepFGRN** – Directed graph embedding and GAN based GRN inference

By integrating the strengths of multiple methods and retaining only consensus-supported regulatory interactions, MFGRN generates high-confidence gene regulatory networks with improved robustness and reduced false positives.

This repository contains the complete MFGRN framework and its application to the study of transcriptional regulatory mechanisms in Sickle Cell Disease (SCD).

---

## Motivation

Gene regulatory network inference is a fundamental task in systems biology. However, different inference algorithms often generate inconsistent results due to methodological biases and dataset-specific characteristics.

Common challenges include:

- High false-positive rates
- Limited robustness across datasets
- Poor reproducibility
- Model-specific prediction bias

MFGRN addresses these limitations through a consensus fusion strategy that integrates multiple deep learning models and retains only high-confidence regulatory interactions supported by multiple methods.

---

## MFGRN Workflow

```text
RNA-seq Expression Matrix
           │
           ▼
 ┌─────────────────────┐
 │   3DCEMA Prediction │
 └─────────────────────┘
           │
           ▼
 ┌─────────────────────┐
 │ DeepSEM Prediction  │
 └─────────────────────┘
           │
           ▼
 ┌─────────────────────┐
 │ DeepFGRN Prediction │
 └─────────────────────┘
           │
           ▼
 Retain Top 10% Edges
           │
           ▼
 Consensus Voting
           │
           ▼
 High-confidence GRN
           │
           ├── Network topology analysis
           ├── Key TF identification
           ├── GO enrichment analysis
           ├── Module detection
           └── ChIP-seq validation
```

---

## Datasets

### Primary Dataset

**GSE244401**

Whole-blood bulk RNA-seq dataset containing:

- 23 pediatric SCD patients
- 17 healthy controls

Samples were collected before and after exercise.

### Experimental Groups

| Group | Description |
|---------|------------|
| CON_T1 | Healthy controls before exercise |
| CON_T2 | Healthy controls after exercise |
| SCD_T1 | SCD patients before exercise |
| SCD_T2 | SCD patients after exercise |

### External Validation Datasets

- GSE254951
- GSE117221

These datasets were used to evaluate the generalization capability of MFGRN.

---

## Methodology

### Step 1: Data Preprocessing

- Gene ID conversion
- Expression matrix normalization
- Removal of low-expression genes
- Human transcription factor annotation
- Construction of TF candidate sets

### Step 2: Independent GRN Inference

Three deep learning models independently infer TF-target regulatory relationships.

| Method | Description |
|----------|----------|
| 3DCEMA | CNN-based regulatory pattern learning |
| DeepSEM | Variational autoencoder + SEM |
| DeepFGRN | Directed graph embedding + GAN |

### Step 3: Consensus Fusion

For each model:

1. Rank regulatory interactions by confidence score
2. Retain top 10% predicted edges

MFGRN then:

- Integrates predictions from all models
- Performs consensus voting
- Retains shared regulatory interactions
- Constructs the final high-confidence GRN

### Step 4: Network Validation

Validation strategies include:

- ChIP-seq supported interactions
- AUROC evaluation
- Fisher's exact test
- Random network comparison
- External dataset validation

---

## Performance

### Internal Validation (GSE244401)

| Group | Best Single Model | MFGRN |
|---------|---------|---------|
| CON_T1 | 0.6598 | **0.7608** |
| CON_T2 | 0.5891 | **0.6361** |
| SCD_T1 | 0.7097 | **0.7109** |
| SCD_T2 | 0.7357 | **0.7851** |

MFGRN consistently achieved superior AUROC performance compared with individual methods.

### External Validation

| Dataset | AUROC |
|----------|----------|
| GSE254951 | 0.6211 |
| GSE117221 Control | 0.5492 |
| GSE117221 Disease | 0.5481 |

These results demonstrate improved robustness and transferability.

---

## Biological Findings

### Key Regulatory Factors

MFGRN identified several biologically relevant transcription factors associated with SCD, including:

- GATA1
- TAL1
- RUNX1
- SPI1
- MYB
- BCL11A

### Functional Characteristics

SCD-related regulatory networks were enriched in:

- Rhythmic process
- Respiratory system development
- Immune response
- Hematopoiesis
- Hypoxia-related regulation

### Network Topology

The inferred GRNs exhibited:

- Small-world characteristics
- High clustering coefficients
- Efficient information transfer
- Significant deviation from random networks

### Functional Modules

Identified modules were enriched in:

- Oxidative stress response
- Cytokine signaling
- Cellular senescence
- Metabolic regulation
- Inflammatory pathways

---

## Repository Structure

```text
MFGRN/
│
├── 3DCEMA/
│   └── 3DCEMA implementation
│
├── DeepSEM/
│   └── DeepSEM implementation
│
├── DeepFGRN/
│   └── DeepFGRN implementation
│
├── data/
│   └── Expression matrices and TF lists
│
├── results/
│   ├── Predicted networks
│   ├── Validation results
│   ├── GO enrichment
│   └── Module analysis
│
├── scripts/
│   ├── Data preprocessing
│   ├── Consensus fusion
│   ├── Network evaluation
│   └── Visualization
│
└── README.md
```

---

## Requirements

### Python

```bash
Python >= 3.8
```

### Main Dependencies

```bash
numpy
pandas
scikit-learn
torch
tensorflow
networkx
matplotlib
seaborn
scipy
```

Install dependencies:

```bash
pip install -r requirements.txt
```

---
## Applications

MFGRN can be applied to:

- Disease-specific GRN reconstruction
- Cancer systems biology
- Hematological disease studies
- Multi-omics integration
- Regulatory factor discovery
- Functional module identification

---

## Contact

**Zhenhao Zan**

GitHub: https://github.com/ZZH706

Research Interests:

- Gene Regulatory Networks
- Deep Learning for Bioinformatics
- Systems Biology
- Multi-omics Analysis
- Sickle Cell Disease

---

## License

This project is released under the MIT License.

---

⭐ If you find MFGRN useful, please consider giving this repository a star.
