#!/usr/bin/env python3
"""
Configuracion central del pipeline Track 1 (MVA Hackathon 2026).

Todas las rutas, umbrales y paneles de genes viven aqui para que el
pipeline sea reproducible y auditable. Ningun otro script debe llevar
constantes escritas a mano.
"""
from __future__ import annotations

import os
from pathlib import Path

# ---------------------------------------------------------------------------
# Rutas
# ---------------------------------------------------------------------------

PROJECT_ROOT = Path(os.environ.get("MVA_ROOT", Path.home() / "mva"))

DATA_DIR = PROJECT_ROOT / "data"
REF_DIR = PROJECT_ROOT / "ref"
DB_DIR = PROJECT_ROOT / "db"
WORK_DIR = PROJECT_ROOT / "work"
RESULTS_DIR = PROJECT_ROOT / "results"
LOGS_DIR = PROJECT_ROOT / "logs"
CONFIG_DIR = PROJECT_ROOT / "config"

# Entradas
RAW_VCF = DATA_DIR / "WGS_EX2312012_HGWCNDSX7.vcf.gz"
REFERENCE_FASTA = REF_DIR / "GRCh38.fa.gz"
GTF = REF_DIR / "Homo_sapiens.GRCh38.116.gtf.gz"
VEP_CACHE = DB_DIR / "vep"

# Intermedios
NORM_VCF = WORK_DIR / "proband.norm.vcf.gz"
VEP_TSV = WORK_DIR / "proband.vep.tsv"
VEP_STATS = WORK_DIR / "proband_vep_summary.html"
RARE_TSV = WORK_DIR / "rare_functional.tsv"

# Salidas
CANDIDATES_TSV = RESULTS_DIR / "recessive_candidates.tsv"
RANKED_TSV = RESULTS_DIR / "ranked_candidates.tsv"
SUBMISSION_CSV = RESULTS_DIR / "track1_submission.csv"

# ---------------------------------------------------------------------------
# Identificadores de la entrega
# ---------------------------------------------------------------------------

# IMPORTANTE: confirmar el identificador exacto en la plataforma del
# hackathon antes de enviar. El template usa el formato "PROBAND01".
PROBAND_ID = os.environ.get("MVA_PROBAND_ID", "PROBAND01")

TEAM_NAME = os.environ.get("MVA_TEAM_NAME", "TEAM_NAME_HERE")

# El VCF usa contigs sin prefijo (1, 2, ... X, Y) pero la plantilla de
# entrega pide "chr15". Este flag controla la conversion al escribir.
ADD_CHR_PREFIX_ON_SUBMIT = True

# ---------------------------------------------------------------------------
# Alcance del analisis
# ---------------------------------------------------------------------------

# Cromosomas primarios. Se excluyen mitocondrial, decoys (Un_JTFH01*),
# scaffolds no colocados y EBV: no aportan anotacion genica util y
# rompen la normalizacion contra la referencia de Ensembl.
PRIMARY_CHROMOSOMES = [str(i) for i in range(1, 23)] + ["X", "Y"]

# ---------------------------------------------------------------------------
# Umbrales de frecuencia poblacional
# ---------------------------------------------------------------------------

# Enfermedad ultra rara recesiva. 1% es deliberadamente permisivo:
# algunos alelos hipomorficos de BUB1B alcanzan frecuencias apreciables
# en poblaciones concretas, y un corte agresivo los eliminaria.
MAX_POPULATION_AF = 0.01

# Corte mas estricto usado solo para puntuar, no para filtrar.
ULTRA_RARE_AF = 0.0001

# ---------------------------------------------------------------------------
# Umbrales de calidad de genotipo
# ---------------------------------------------------------------------------

MIN_DEPTH = 10          # DP minima para confiar en el genotipo
MIN_GQ = 20             # calidad de genotipo minima
MIN_ALT_READS = 4       # lecturas que soportan el alelo alterno

# Balance alelico. Un heterocigoto real se situa cerca de 0.5.
# Este filtro es el que descarta el cluster artefactual de SERPINA1,
# donde 13 variantes consecutivas mostraban VAF de 0.15 a 0.30.
HET_VAF_MIN = 0.25
HET_VAF_MAX = 0.75
HOM_VAF_MIN = 0.85

# ---------------------------------------------------------------------------
# Deteccion de heterocigosis compuesta
# ---------------------------------------------------------------------------

# Un gen con demasiadas heterocigotas raras casi siempre indica un locus
# problematico (homopolimero, duplicacion segmental, pseudogen) y no
# biologia. Ademas evita la explosion combinatoria: 13 heterocigotas
# generan 78 pares.
MAX_HETS_PER_GENE = 5

# Dos variantes muy juntas suelen ser la misma alteracion representada de
# varias formas por el caller, no dos eventos independientes.
MIN_PAIR_DISTANCE_BP = 30

# ---------------------------------------------------------------------------
# Consecuencias consideradas funcionalmente relevantes
# ---------------------------------------------------------------------------

HIGH_IMPACT_CONSEQUENCES = {
    "transcript_ablation",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "stop_gained",
    "frameshift_variant",
    "stop_lost",
    "start_lost",
    "transcript_amplification",
}

MODERATE_IMPACT_CONSEQUENCES = {
    "inframe_insertion",
    "inframe_deletion",
    "missense_variant",
    "protein_altering_variant",
}

# Consecuencias de splicing que VEP clasifica como LOW pero que pueden
# ser relevantes. Se incluyen por subcadena.
SPLICE_KEYWORDS = ("splice",)

ACCEPTED_IMPACTS = {"HIGH", "MODERATE"}

# ---------------------------------------------------------------------------
# Paneles de genes (conocimiento previo, usado solo para puntuar)
# ---------------------------------------------------------------------------
#
# El nucleo del pipeline es agnostico: no filtra por gen. Estos paneles
# se aplican unicamente en la fase de ranking y su peso es configurable,
# de modo que se puede reportar el resultado con y sin conocimiento
# previo. Los jueces valoran poder distinguir ambas cosas.

MVA_CORE_GENES = {
    "BUB1B",    # MVA1, con predisposicion tumoral
    "CEP57",    # MVA2
    "TRIP13",   # MVA3
    "CENATAC",  # CCDC84
    "BUB1",
    "SLF2",
    "SMC5",
}

MITOTIC_CHECKPOINT_GENES = {
    "BUB3", "MAD1L1", "MAD2L1", "CDC20", "TTK", "ESPL1", "PLK4",
    "CENPE", "KNL1", "ZW10", "ZWILCH", "SPDL1", "NDC80", "AURKB",
    "INCENP", "SGO1", "PTTG1", "CEP192", "TUBB",
}

# Genes de predisposicion a rabdomiosarcoma y tumores embrionarios.
# Relevantes como hallazgos secundarios.
CANCER_PREDISPOSITION_GENES = {
    "TP53", "DICER1", "NF1", "HRAS", "KRAS", "NRAS", "PTPN11",
    "SOS1", "RB1", "WT1", "SMARCB1", "SMARCA4", "FH", "BRCA1", "BRCA2",
}

# ---------------------------------------------------------------------------
# Fenotipo del proband (HPO)
# ---------------------------------------------------------------------------
#
# Derivados del documento clinico del reto. Se usan para la priorizacion
# fenotipica (Exomiser) y no contienen informacion identificable.

HPO_TERMS = [
    "HP:0002859",  # Rabdomiosarcoma
    "HP:0000121",  # Nefrocalcinosis
    "HP:0004322",  # Talla baja
    "HP:0001508",  # Fallo de medro
    "HP:0003202",  # Atrofia muscular esqueletica
    "HP:0001622",  # Prematuridad
    "HP:0001518",  # Pequeno para edad gestacional
    "HP:0200067",  # Abortos espontaneos recurrentes (en los padres)
]

# ---------------------------------------------------------------------------
# Pesos del sistema de puntuacion
# ---------------------------------------------------------------------------

SCORE_WEIGHTS = {
    "impact_high": 30.0,
    "impact_moderate": 12.0,
    "splice_region": 8.0,
    "clinvar_pathogenic": 35.0,
    "clinvar_benign": -40.0,
    "absent_from_gnomad": 12.0,
    "ultra_rare": 8.0,
    "predictors_concordant_damaging": 10.0,
    "predictors_concordant_benign": -15.0,
    "phase_trans_confirmed": 25.0,
    "phase_cis_confirmed": -100.0,
    "model_compound_het": 10.0,
    "model_homozygous": 10.0,
    # Componente dependiente de conocimiento previo
    "panel_mva_core": 40.0,
    "panel_mitotic": 15.0,
    "panel_cancer": 5.0,
}

# Si es False, el ranking ignora los paneles y queda completamente
# agnostico. Util para generar la version comparativa del informe.
USE_GENE_PANELS = True

# ---------------------------------------------------------------------------
# Recursos de computo
# ---------------------------------------------------------------------------

THREADS = int(os.environ.get("MVA_THREADS", "10"))
VEP_BUFFER_SIZE = 20000


def ensure_directories() -> None:
    """Crea los directorios de trabajo si no existen."""
    for d in (WORK_DIR, RESULTS_DIR, LOGS_DIR, CONFIG_DIR):
        d.mkdir(parents=True, exist_ok=True)
