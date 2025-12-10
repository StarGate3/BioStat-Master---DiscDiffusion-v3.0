# BioStat Master v3.0

**BioStat Master** is a specialized, modular Python application designed for the statistical analysis and visualization of **Disk Diffusion Assays** (Zone of Inhibition data). It automates the entire workflow from raw Excel data to publication-ready figures and reports, ensuring statistical rigor using standard biostatistical methods.

---

## 🚀 Key Features

### 📊 Statistical Analysis
*   **Automated Decision Tree**: Automatically selects between Parametric (ANOVA) and Non-Parametric (Kruskal-Wallis) tests based on Normality (Shapiro-Wilk) and Homogeneity of Variance (Levene’s test).
*   **Post-hoc Corrections**: Supports **Tukey HSD** (for ANOVA) and **Dunn’s Test** (for Kruskal-Wallis) with multiple correction methods:
    *   Holm-Bonferroni (Default)
    *   Benjamini-Hochberg (FDR)
    *   Bonferroni
*   **Effect Size**: Calculates **Cohen’s *d*** for all pairwise comparisons to determine the magnitude of differences.
*   **Outlier Detection**: Implements **Dixon’s Q Test** to identify and suggest removal of technical outliers in small sample sizes ($3 \le n \le 10$).

### 🎨 Scientific Visualization
Generates high-resolution, publication-quality figures using `Matplotlib` and `Seaborn`:
1.  **Main Comparison Plot**: Barplots, Boxplots, or Violinplots with significance asterisks.
2.  **Heatmaps**: Activity heatmaps and P-value significance matrices.
3.  **Dose-Response Trends**: Line plots with Spearman correlation coefficients and **MIC Estimation** (Minimal Inhibitory Concentration).
4.  **Effect Size Plot**: Lollipop charts visualizing the strength of differences (Cohen's d).
5.  **Cross-Species Comparison**: Summary view of activity across multiple bacterial strains.
6.  **PCA Analysis**: Principal Component Analysis to cluster bacterial strains based on their sensitivity profiles.

### 📝 Reporting
*   **PDF Reports**: detailed summary including methodology, statistical results, and embedded figures.
*   **Caption Generator**: Automatically generates scientific figure captions (e.g., "Figure 1. Antibacterial activity...") ready for copy-pasting into manuscripts.
*   **Excel Export**: Exports raw data, statistical summaries, and detailed post-hoc results.

---

## 🛠️ Installation & Requirements

Ensure you have **Python 3.8+** installed.

1.  **Clone/Download** the repository.
2.  **Install Dependencies**:
    ```bash
    pip install -r requirements.txt
    ```
    *Key libraries: `customtkinter`, `pandas`, `scipy`, `statsmodels`, `seaborn`, `scikit-posthocs`, `reportlab`.*

---

## 💻 Usage

1.  **Run the Application**:
    ```bash
    python main.py
    ```
2.  **Load Data**: Click "1. Wczytaj Excel". The Excel file should be formatted with columns for 'Bakterie' (Bacteria), 'Grupa' (Group/Substance), and 'Srednica_mm' (Zone Diameter).
3.  **Configure**:
    *   Select the bacterial strain to analyze.
    *   Choose a **Reference Group** (Negative Control) for comparisons.
    *   Select a Post-hoc correction method.
4.  **Run Analysis**: Click "URUCHOM ANALIZĘ".
5.  **Explore Results**: Switch between tabs to view different plots and the statistical log.
6.  **Export**: Save figures as high-res PNGs or generate a full PDF report.

---

## 🏗️ Architecture (v3.0 Modular)

The application follows a clean, modular Model-View-Controller (MVC) pattern for maintainability:

*   **`gui.py` (View/Controller)**: Handles the user interface using `customtkinter`. Orchestrates the application flow.
*   **`logic.py` (Model)**: Contains the `StatsEngine`. Pure Python class responsible for all statistical calculations (Shapiro, Levene, ANOVA/KW, Post-hoc). independent of the GUI.
*   **`plotting.py` (View)**: Contains the `Plotter` class. Encapsulates all `matplotlib` figure generation logic.
*   **`utils.py`**: Helper functions for outlier detection (Dixon), robust sorting, and string parsing.

---

## 🔬 Statistical Methods

The software adheres to standard biostatistical practices for small-sample biological data:
1.  **Normality Check**: Shapiro-Wilk test ($\alpha=0.05$).
2.  **Variance Check**: Levene’s test ($\alpha=0.05$).
3.  **Test Selection**:
    *   **Parametric**: If Normal AND Homogeneous Variances $\rightarrow$ **One-way ANOVA**.
    *   **Non-Parametric**: If Non-Normal OR Unequal Variances $\rightarrow$ **Kruskal-Wallis**.
4.  **Pairwise Comparison**:
    *   ANOVA $\rightarrow$ **Tukey HSD**.
    *   Kruskal-Wallis $\rightarrow$ **Dunn’s Test** (corrected).
