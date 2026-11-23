"""Ingestion GEO -> Silver (Parquet). Réutilisable pour AS et PD."""
import json
from pathlib import Path
import GEOparse
import pandas as pd


def download(gse_id, bronze="data/bronze"):
    return GEOparse.get_GEO(geo=gse_id, destdir=bronze)


# ---- Annotation sonde -> symbole (une fonction par type de plateforme) ----
def symbol_agilent(gpl):
    """GPL17077 (GSE100927)."""
    a = gpl[(gpl["CONTROL_TYPE"] == "FALSE") & gpl["GENE_SYMBOL"].notna()]
    a = a[a["GENE_SYMBOL"].str.strip() != ""].copy()
    a["ID"] = a["ID"].astype(str)
    return a.drop_duplicates("ID").set_index("ID")["GENE_SYMBOL"]


def _first_affy_st(s):
    if not isinstance(s, str) or s.strip() in ("", "---"):
        return None
    parts = s.split(" /// ")[0].split(" // ")
    if len(parts) > 1 and parts[1].strip() not in ("", "---"):
        return parts[1].strip()
    return None


def symbol_affy_st(gpl):
    """GPL6244 (GSE43292): symbole = 2e champ de gene_assignment."""
    a = gpl[gpl["category"] == "main"].copy()
    a["symbol"] = a["gene_assignment"].map(_first_affy_st)
    a = a.dropna(subset=["symbol"])
    a["ID"] = a["ID"].astype(str)
    return a.drop_duplicates("ID").set_index("ID")["symbol"]


def symbol_affy_u133(gpl, col="Gene Symbol"):
    """HG-U133 Plus 2 (à vérifier pour GSE16134 / GSE10334 : nom de la colonne)."""
    a = gpl[gpl[col].notna()].copy()
    a = a[~a[col].astype(str).str.strip().isin(["", "---"])]
    a["symbol"] = a[col].astype(str).str.split(" /// ").str[0].str.strip()
    a["ID"] = a["ID"].astype(str)
    return a.drop_duplicates("ID").set_index("ID")["symbol"]


# ---- Sondes -> gènes (médiane) ----
def collapse_to_genes(expr, mapping):
    expr = expr.copy()
    expr.index = expr.index.astype(str)
    e = expr.loc[expr.index.intersection(mapping.index)].copy()
    e["gene"] = e.index.map(mapping)
    g = e.groupby("gene").median()
    g.columns.name = None
    stats = {"probes_before": int(expr.shape[0]),
             "probes_annotated": int(e.shape[0]),
             "genes": int(g.shape[0])}
    return g, stats


# ---- Validation (contrat de fichiers) ----
def validate(expr_gene, pheno):
    need = {"sample_id", "group", "age", "sex", "platform"}
    assert need <= set(pheno.columns), f"colonnes manquantes: {need - set(pheno.columns)}"
    assert set(pheno["group"]) <= {"control", "disease"}, "group doit être control/disease"
    assert set(expr_gene.columns) == set(pheno["sample_id"]), "échantillons expr != pheno"
    assert expr_gene.index.is_unique, "gènes dupliqués"
    assert not expr_gene.isna().any().any(), "NaN dans expr"


def write_silver(gse_id, expr_gene, pheno, stats, silver="data/silver"):
    out = Path(silver) / gse_id
    out.mkdir(parents=True, exist_ok=True)
    expr_gene.reset_index().to_parquet(out / "expr.parquet")
    pheno.to_parquet(out / "pheno.parquet")
    lineage = {**stats, "gse": gse_id, "n_samples": int(pheno.shape[0]),
               "n_disease": int((pheno["group"] == "disease").sum()),
               "n_control": int((pheno["group"] == "control").sum())}
    (out / "lineage.json").write_text(json.dumps(lineage, indent=2))
    return out