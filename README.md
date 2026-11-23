# Parodontite × Athérosclérose : gènes « crosstalk », ML explicable et causalité

Projet de Master Bioinformatique, binôme **Hammou** et **Maroune**.
Il valide trois modules : Big Data, Analyse intégrative multimodale, XAI et IA causale.

## 1. Article de référence

Dong W., Gong Y., Yang B., Li B. (2023). *Bioinformatics analysis and reveal potential crosstalk genetic and immune relationships between atherosclerosis and periodontitis.* Scientific Reports 13:10381.
DOI : 10.1038/s41598-023-37027-x (open access)
https://www.nature.com/articles/s41598-023-37027-x

**À lire en priorité** : *Materials and methods* et *Results* (10 minutes), puis les 7 figures.

### Ce que fait l'article (pipeline à reproduire)

| Étape | Outil | Paramètres / cibles de l'article |
|---|---|---|
| Données | GEO | PD : GSE16134 (découverte), GSE10334 (validation). AS : GSE100927 (découverte), GSE43292 (validation). IRGs : ImmPort |
| Infiltration immune | CIBERSORT (LM22) + Wilcoxon | 22 types cellulaires |
| DEGs | `limma` | \|log2FC\| > 1, p.adj < 0,05. Cibles : ~421 (AS), ~167 (PD) |
| WGCNA | `WGCNA` | `pickSoftThreshold`. β = 13 (AS), β = 15 (PD), R² ≈ 0,9. Gènes clés : \|MM\| > 0,8 et \|GS\| > 0,4. Cible : 14 gènes après intersection des 4 modules |
| ML | `caret` + `DALEX` | RF, SVM, GLM. GLM meilleur en AS, SVM en PD |
| Intersection | Venn | WGCNA ∩ DEGs ∩ IRGs = **BTK, ITGAL** |
| Validation | `rms`, `pROC`, `GSVA` | Nomogramme, calibration, DCA, ROC (AUC 0,81 à 0,94), ssGSEA, GO/KEGG |

Une reproduction exacte est peu probable (versions de packages, filtres non détaillés, permutations CIBERSORT). L'objectif est d'obtenir des valeurs **proches** et de **documenter les écarts**.

### Ce que nous ajoutons (notre contribution)

1. **Big Data** : ingestion GEO → Parquet (PySpark, MinIO), validation qualité, benchmark Spark vs pandas.
2. **Multimodal** : expression + fractions immunes par échantillon (MOFA+), par maladie.
3. **XAI** : validation croisée correcte, XGBoost, SHAP, calibration.
4. **Causal** : Mendelian randomization (`TwoSampleMR`) sur BTK et ITGAL.
5. **Généralisation** : refaire le pipeline sur d'autres jeux GEO de la même thématique.

## 2. Répartition

| | Hammou | Maroune |
|---|---|---|
| **Maladie** | Athérosclérose (AS) | Parodontite (PD) |
| **Big Data** | Docker (Spark + MinIO), ingestion AS, benchmark | Ingestion PD, `pandera`, `data_dictionary.yaml`, ETL GWAS/eQTL |
| **Multimodal** | CIBERSORT, `limma`, WGCNA AS, MOFA+ AS | CIBERSORT, `limma`, WGCNA PD, MOFA+ PD |
| **Demi-module final** | **Causal** : MR | **XAI** : ML, SHAP, nomogramme |
| **Commun** | Intersection finale, rapport, soutenance, revue croisée des PR | idem |

Hammou ajoute GO/KEGG et ssGSEA. Chacun ajuste `limma` sur l'âge et le sexe si les métadonnées le permettent.

## 3. Contrats de fichiers (même format pour les deux, pour travailler en parallèle)

```
data/silver/<GSE>/expr.parquet       # lignes = gènes (symbole HGNC), colonnes = échantillons
data/silver/<GSE>/pheno.parquet      # sample_id, group (control/disease), age, sex, platform
results/<AS|PD>/deg.csv              # gene, log2FC, pvalue, padj
results/<AS|PD>/wgcna_key_genes.csv  # gene, module, MM, GS
results/<AS|PD>/immune_fractions.csv # sample_id + 22 colonnes LM22
```

## 4. Structure du dépôt

```
perio-athero-crosstalk/
├── README.md
├── environment.yml            # Python
├── renv.lock                  # R
├── docker-compose.yml         # Spark master + 2 workers + MinIO
├── config/datasets.yaml       # accessions GEO, plateformes (GPL), URLs
├── data/                      # gitignoré
│   ├── bronze/  silver/  gold/
├── docs/
│   ├── data_dictionary.yaml
│   ├── article_reproduction.md   # tableau : valeur article vs notre valeur
│   └── meetings/
├── src/
│   ├── etl/        # geo_ingest.py, annotate.py, validate.py
│   ├── multimodal/ # cibersort.R, limma.R, wgcna.R, mofa.py
│   ├── xai/        # train.py, shap_analysis.py, nomogram.R
│   └── causal/     # mr_analysis.R
├── notebooks/
├── results/        # AS/  PD/  figures/  tables/
└── tests/
```

## 5. Installation

```bash
mamba create -n perio-athero python=3.11 r-base=4.3 -c conda-forge
mamba activate perio-athero
pip install GEOparse pandas pyarrow pandera pyspark scikit-learn xgboost shap \
            mofapy2 matplotlib seaborn
```

```r
install.packages(c("BiocManager", "remotes", "renv", "caret", "DALEX",
                   "rms", "pROC", "WGCNA", "data.table"))
BiocManager::install(c("GEOquery", "limma", "sva", "GSVA", "clusterProfiler",
                       "org.Hs.eg.db", "MOFA2"))
remotes::install_github("MRCIEU/TwoSampleMR")
```

**CIBERSORT** : le code source demande une inscription gratuite sur https://cibersortx.stanford.edu (non installable via CRAN). Alternative à évaluer : `immunedeconv`.

## 6. Démarrage rapide : ingestion d'un jeu GEO

```python
import GEOparse, pandas as pd

gse = GEOparse.get_GEO(geo="GSE16134", destdir="data/bronze")  # remplacer par ton GSE
expr = gse.pivot_samples("VALUE")                  # sondes × échantillons
print(gse.gpls)                                    # plateforme (GPL) : à noter dans config/datasets.yaml
pheno = gse.phenotype_data                         # métadonnées : group, age, sex...

expr.to_parquet("data/silver/GSE16134/expr_probes.parquet")
pheno.to_parquet("data/silver/GSE16134/pheno_raw.parquet")
```

Étapes suivantes : annotation sonde → gène via la table GPL (agrégation par médiane), harmonisation de `pheno` (colonne `group` = control/disease), validation `pandera`, écriture de `expr.parquet` et `pheno.parquet`.

## 7. Workflow Git

- Branches : `main` (protégée) ← `dev` ← `feat/<prénom>-<sujet>` (ex. `feat/maroune-pd-ingest`).
- Commits conventionnels : `feat:`, `fix:`, `docs:`.
- 1 Pull Request = 1 relecteur (l'autre membre du binôme).
- Tâches dans GitHub Projects (Kanban), 1 issue par tâche.
- Point de synchronisation **lundi et jeudi, 30 min**. Compte rendu dans `docs/meetings/AAAA-MM-JJ.md`.

## 8. Feuille de route (8 semaines)

| Sem. | Hammou | Maroune | Jalon |
|---|---|---|---|
| 1 | Repo, Docker, ingestion AS | Ingestion PD, `pandera` | Silver complet |
| 2 | CIBERSORT + `limma` AS | CIBERSORT + `limma` PD | ~421 / ~167 DEGs |
| 3 | WGCNA AS | WGCNA PD | 14 gènes retrouvés |
| 4 | Intersection, GO/KEGG, ssGSEA | RF/SVM/GLM + CV | BTK / ITGAL reproduits |
| 5 | MOFA+ AS, préparation MR | MOFA+ PD, XGBoost | Reproduction validée |
| 6 | MR (`TwoSampleMR`) | SHAP + nomogramme | Résultats XAI + causal |
| 7 | Nouveaux jeux GEO (AS) | Nouveaux jeux GEO (PD) | Test de généralisation |
| 8 | Rapport + soutenance | Rapport + soutenance | Livraison |

## 9. Checklist de la semaine 1

- [ ] Cloner le dépôt, créer sa branche `feat/<prénom>-...`
- [ ] Créer l'environnement (section 5)
- [ ] Lire l'article (Methods + Results)
- [ ] Télécharger ses 2 jeux GEO (AS pour Hammou : GSE100927, GSE43292 ; PD pour Maroune : GSE16134, GSE10334)
- [ ] Noter la plateforme (GPL) et les effectifs (contrôles / malades) dans `config/datasets.yaml`
- [ ] Produire `expr.parquet` et `pheno.parquet` selon les contrats
- [ ] Ouvrir une Pull Request et la faire relire par l'autre

## 10. Métriques à rapporter

- **Données** : lignes avant/après chaque étape, taux de complétude, taille Parquet vs CSV, temps Spark vs pandas.
- **WGCNA** : R² de l'indépendance d'échelle, corrélation module-trait.
- **ML** : AUROC, AUPRC, Brier, courbe de calibration, DCA, stabilité du top-k SHAP (Jaccard sur bootstrap).
- **MR** : OR (IVW), F > 10, test Q de Cochran, intercept d'Egger, médiane pondérée, Steiger.

## 11. Limites à discuter dans le rapport

- Confusion non contrôlée dans l'article (âge, sexe, traitements) ; nous ajustons quand les métadonnées le permettent.
- Les ROC « découverte » utilisent les données qui ont servi à choisir les gènes ; nous ajoutons une validation croisée et une validation externe.
- MR : instruments faibles probables pour la parodontite, eQTL de sang alors que l'analyse est tissulaire.
- MOFA+ est appliqué par maladie, car les échantillons AS et PD ne se recoupent pas.
