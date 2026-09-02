# MFGRN: Multi-Algorithm Fusion for Constructing a Sickle Cell Disease-Specific Gene Regulatory Network
## Overview
MFGRN (Sickle Cell Disease–Multi-Algorithm Fusion Gene Regulatory Network) is a consensus-based framework designed for constructing high-confidence gene regulatory networks (GRNs) in Sickle Cell Disease (SCD).
The framework integrates three complementary deep learning-based GRN inference methods:
* **3DCEMA** – 3D convolutional neural network-based GRN inference
* **DeepSEM** – Variational autoencoder and structural equation model-based GRN inference
* **DeepFGRN** – Directed graph embedding and GAN-based GRN inference
By combining predictions from multiple algorithms and retaining consensus-supported regulatory interactions, SCD-MFGRN reduces false positives and improves the robustness and biological reliability of inferred regulatory networks.
The framework was developed using whole-blood RNA-seq datasets from pediatric SCD patients and healthy controls and was further validated using independent external datasets.
---

## Motivation

Sickle Cell Disease (SCD) is a hereditary hemoglobinopathy characterized by chronic hemolysis, inflammation, vascular dysfunction, and abnormal erythropoiesis.

Understanding transcriptional regulatory mechanisms underlying SCD is critical for identifying key regulators and potential therapeutic targets. However, different GRN inference algorithms often generate inconsistent regulatory networks due to methodological biases and dataset-specific characteristics.

To address these limitations, MFGRN integrates multiple deep learning-based GRN inference methods and constructs consensus-supported regulatory networks with improved reliability and biological interpretability.

---

## Research Workflow

<p align="center">
  <img src="Research process.svg" alt="SCD-MFGRN Workflow" width="1000">
</p>

---

## Datasets

### Primary Dataset

**GSE244401**

Whole-blood bulk RNA-seq dataset containing:

* 23 pediatric SCD patients
* 17 healthy controls

Samples were collected before and after exercise.

### Experimental Groups

| Group  | Description                      |
| ------ | -------------------------------- |
| CON_T1 | Healthy controls before exercise |
| CON_T2 | Healthy controls after exercise  |
| SCD_T1 | SCD patients before exercise     |
| SCD_T2 | SCD patients after exercise      |

### External Validation Datasets

* GSE254951
* GSE117221

These datasets were used to evaluate the generalization capability of SCD-MFGRN.

---

## Methodology

### Step 1: Data Preprocessing

* Gene ID conversion
* Expression matrix normalization
* Removal of low-expression genes
* Human transcription factor annotation
* Construction of TF candidate sets

### Step 2: Independent GRN Inference

Three deep learning models independently infer TF-target regulatory relationships.

| Method   | Description                                         |
| -------- | --------------------------------------------------- |
| 3DCEMA   | CNN-based regulatory pattern learning               |
| DeepSEM  | Variational autoencoder + Structural Equation Model |
| DeepFGRN | Directed graph embedding + GAN                      |

### Step 3: Consensus Fusion

For each model:

1. Rank regulatory interactions according to confidence scores
2. Retain the top 10% regulatory edges

SCD-MFGRN then:

* Integrates predictions from all methods
* Performs consensus voting
* Retains shared regulatory interactions
* Constructs a final high-confidence GRN

### Step 4: Network Validation

Validation strategies include:

* ChIP-seq supported interactions
* AUROC evaluation
* Fisher’s exact test
* Random network comparison
* External dataset validation

---

## Biological Insights from SCD-MFGRN

### Key Transcription Factors

SCD-MFGRN identified several transcription factors with established roles in hematopoiesis and erythroid differentiation:

* GATA1
* TAL1
* RUNX1
* SPI1
* MYB
* BCL11A

### Functional Enrichment

The inferred regulatory networks were significantly enriched in biological processes associated with:

* Hematopoiesis
* Immune regulation
* Hypoxia response
* Respiratory system development
* Circadian rhythm-related processes

### Network Characteristics

Network topology analysis demonstrated that SCD-MFGRN-generated networks exhibited:

* Small-world properties
* High clustering coefficients
* Efficient information transfer
* Significant deviation from random networks

### Functional Modules

Module analysis revealed subnetworks associated with:

* Oxidative stress response
* Cytokine signaling
* Cellular senescence
* Inflammatory regulation
* Metabolic processes

---

## Applications

SCD-MFGRN can be applied to:

* Disease-specific GRN reconstruction
* Hematological disease studies
* Cancer systems biology
* Regulatory factor discovery
* Multi-omics integration
* Functional module identification
* Systems biology research

---

## Contact

**Zhenhao Zan**

GitHub:

https://github.com/ZZH706

Research Interests:

* Gene Regulatory Networks
* Deep Learning for Bioinformatics
* Systems Biology
* Multi-omics Analysis
* Sickle Cell Disease

---

⭐ If you find SCD-MFGRN useful, please consider giving this repository a star.
