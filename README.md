# 🧬 qPCR Primer Designer  
A modern Streamlit-based web application for rapid and precise **qPCR primer design**, integrating NCBI sequence retrieval and Primer3-based primer generation.

<p align="center">
  <img src="assets/banner.png" width="80%">
</p>

---

## 🌟 Overview

The **qPCR Primer Designer** is a fully interactive web tool designed for molecular biology researchers who need reliable, publication-grade primers.  
This application automatically:

1. Fetches **mRNA sequences from NCBI (Entrez API)**
2. Designs **optimal primer pairs using Primer3**
3. Visualizes primer binding positions using **Plotly**
4. Allows users to **download primer tables and FASTA sequences**
5. Provides full control over primer parameters, GC%, Tm, amplicon size, and reaction conditions.

---

## 🔧 Features

### ✅ **Automated Workflow**
- Input a gene name (e.g., *TP53*), and the app automatically retrieves the correct mRNA sequence.
- Primer3 generates multiple primer pairs under your selected constraints.

### 🎚️ **Customizable Primer Design Parameters**
- Primer length (min/opt/max)
- Melting temperature (Tm)
- GC content range
- Amplicon size range
- Maximum Poly-X
- Max GC at 3' end
- Custom reaction buffer conditions

### 📊 **High-Quality Visualizations**
- Primer binding map along the mRNA sequence  
- Top primer pairs with Tm, GC%, and product size  
- Expandable details for each pair  

### 📥 **Downloads**
- Primer list (CSV)
- mRNA sequence (FASTA)
- Ready-to-upload sequences for NCBI Primer-BLAST

### 🖥️ **Clean & Modern UI**
- Fully themed with CSS for a professional look  
- Sidebar workflow indicator  
- Interactive tabs for parameters, advanced settings, and results  

---

## 📌 Installation

Clone the repository:

```bash
git clone https://github.com/your-repo/qPCR-Primer-Designer.git
cd qPCR-Primer-Designer
