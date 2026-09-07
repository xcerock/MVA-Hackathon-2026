#!/usr/bin/env python3
"""
Paso 5: puntuacion, ranking y generacion del CSV de entrega.

El sistema de evaluacion del Track 1 concede 100 puntos si la variante
confirmada ocupa la primera posicion, credito parcial si aparece entre
las diez primeras, y medio credito si se acierta solo una de las dos
variantes de un par. El orden importa mas que la exhaustividad.

La puntuacion combina evidencia de varias fuentes independientes:

  Severidad de la consecuencia (HIGH sobre MODERATE).
  Clasificacion clinica previa en ClinVar.
  Rareza poblacional, con premio a las variantes ausentes de gnomAD.
  Concordancia entre SIFT y PolyPhen.
  Fase fisica: trans confirmada suma, cis confirmada elimina.
  Pertenencia a paneles de genes (conocimiento previo).

El componente de paneles se puede desactivar con --no-panels para
generar la version completamente agnostica del ranking. Reportar ambas
versiones demuestra que el hallazgo no depende de haber sabido de
antemano donde mirar.

Uso:
    python 05_rank_submission.py --proband-id PROBAND01
    python 05_rank_submission.py --no-panels --output-suffix agnostic
"""
from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

import config
from common import get_logger, parse_af, parse_prediction

W = config.SCORE_WEIGHTS


# ClinVar entrega una lista de clasificaciones separadas por comas, y
# algunas etiquetas contienen la palabra "pathogenic" sin afirmar
# patogenicidad. En particular "conflicting_classifications_of_
# pathogenicity" indica desacuerdo entre laboratorios, no evidencia a
# favor. La comparacion se hace por etiqueta completa, nunca por
# subcadena.

CLINVAR_PATHOGENIC_LABELS = {
    "pathogenic",
    "likely_pathogenic",
    "pathogenic/likely_pathogenic",
}

CLINVAR_BENIGN_LABELS = {
    "benign",
    "likely_benign",
    "benign/likely_benign",
}

CLINVAR_NEUTRAL_LABELS = {
    "uncertain_significance",
    "conflicting_classifications_of_pathogenicity",
    "conflicting_interpretations_of_pathogenicity",
    "not_provided",
    "other",
    "association",
    "risk_factor",
    "drug_response",
    "protective",
    "affects",
}


# ---------------------------------------------------------------------------
# Correccion de anotaciones de splicing espurias
# ---------------------------------------------------------------------------
#
# VEP marca como splice_donor_variant cualquier variante que solape la
# ventana de la region de splicing, incluidas deleciones largas que
# empiezan decenas de bases dentro del intron y no tocan el sitio donor
# canonico (posiciones +1 y +2, con relevancia funcional decreciente
# hasta +5). Sin esta correccion, deleciones intronicas profundas reciben
# impacto HIGH y desplazan a candidatos reales.

_HGVS_INTRON_OFFSET = re.compile(r"c\.-?\d+([+-])(\d+)")

# Distancia maxima al exon para considerar creible una anotacion de
# splicing. Donor: hasta +5. Aceptor: hasta -20, que cubre el punto de
# ramificacion y el tracto de polipirimidina.
MAX_DONOR_OFFSET = 5
MAX_ACCEPTOR_OFFSET = 20


def splice_annotation_is_credible(consequence: str, hgvsc: str) -> bool:
    """
    Comprueba si una anotacion de splicing es plausible por posicion.

    Devuelve True cuando la consecuencia no es de splicing, cuando no hay
    HGVSc para evaluar, o cuando el desplazamiento intronico cae dentro
    del rango funcionalmente relevante.
    """
    if "splice" not in consequence:
        return True
    if not hgvsc or hgvsc == "-":
        return True

    match = _HGVS_INTRON_OFFSET.search(hgvsc)
    if match is None:
        # Variante exonica anotada como splice_region: es creible.
        return True

    sign, offset = match.group(1), int(match.group(2))
    limit = MAX_DONOR_OFFSET if sign == "+" else MAX_ACCEPTOR_OFFSET
    return offset <= limit


def clinvar_component(clin_sig: str) -> tuple[float, list[str]]:
    """
    Puntua la clasificacion clinica previa.

    Solo cuentan las etiquetas que afirman patogenicidad o benignidad de
    forma inequivoca. Las clasificaciones conflictivas y las de
    significado incierto se tratan como ausencia de evidencia.
    """
    if not clin_sig or clin_sig == "-":
        return 0.0, []

    labels = {label.strip().lower()
              for label in clin_sig.split(",") if label.strip()}

    pathogenic = labels & CLINVAR_PATHOGENIC_LABELS
    benign = labels & CLINVAR_BENIGN_LABELS
    conflicting = labels & {
        "conflicting_classifications_of_pathogenicity",
        "conflicting_interpretations_of_pathogenicity",
    }

    # Clasificaciones contradictorias entre si: sin credito en ninguna
    # direccion, pero se deja constancia.
    if pathogenic and benign:
        return 0.0, ["clinvar_contradictorio"]
    if pathogenic:
        note = "clinvar_patogenico"
        if conflicting:
            # Hay evidencia de patogenicidad pero tambien desacuerdo:
            # se concede credito parcial.
            return W["clinvar_pathogenic"] * 0.4, [note + "_con_conflicto"]
        return W["clinvar_pathogenic"], [note]
    if benign:
        return W["clinvar_benign"], ["clinvar_benigno"]
    return 0.0, []


def frequency_component(max_af: str) -> tuple[float, list[str]]:
    af = parse_af(max_af)
    if af is None:
        return W["absent_from_gnomad"], ["ausente_gnomad"]
    if af < config.ULTRA_RARE_AF:
        return W["ultra_rare"], ["ultra_rara"]
    return 0.0, []


def impact_component(impact: str, consequence: str) -> tuple[float, list[str]]:
    if impact == "HIGH":
        return W["impact_high"], ["impacto_alto"]
    if impact == "MODERATE":
        score = W["impact_moderate"]
        notes = ["impacto_moderado"]
        if any(k in consequence for k in config.SPLICE_KEYWORDS):
            score += W["splice_region"]
            notes.append("region_splicing")
        return score, notes
    if any(k in consequence for k in config.SPLICE_KEYWORDS):
        return W["splice_region"], ["region_splicing"]
    return 0.0, []


def predictor_component(sift: str, polyphen: str) -> tuple[float, list[str]]:
    sift_label, _ = parse_prediction(sift)
    pp_label, _ = parse_prediction(polyphen)
    if sift_label is None or pp_label is None:
        return 0.0, []
    sift_bad = "deleterious" in sift_label
    pp_bad = "damaging" in pp_label
    if sift_bad and pp_bad:
        return W["predictors_concordant_damaging"], ["predictores_deletereos"]
    if not sift_bad and pp_label.startswith("benign"):
        return W["predictors_concordant_benign"], ["predictores_benignos"]
    return 0.0, []


def panel_component(gene: str) -> tuple[float, list[str]]:
    if gene in config.MVA_CORE_GENES:
        return W["panel_mva_core"], ["panel_mva"]
    if gene in config.MITOTIC_CHECKPOINT_GENES:
        return W["panel_mitotic"], ["panel_checkpoint_mitotico"]
    if gene in config.CANCER_PREDISPOSITION_GENES:
        return W["panel_cancer"], ["panel_predisposicion_tumoral"]
    return 0.0, []


def score_row(row: dict[str, str], use_panels: bool) -> tuple[float, list[str]]:
    score = 0.0
    notes: list[str] = []

    model = row.get("model", "")
    if model == "compound_het":
        score += W["model_compound_het"]
    elif model == "homozygous":
        score += W["model_homozygous"]

    phase = row.get("phase", "unknown")
    if phase == "trans":
        score += W["phase_trans_confirmed"]
        notes.append("fase_trans")
    elif phase == "cis":
        score += W["phase_cis_confirmed"]
        notes.append("fase_cis")

    # Ambos alelos contribuyen. El segundo aporta la mitad para no
    # duplicar el peso de un solo evento biologico.
    for suffix, weight in (("1", 1.0), ("2", 0.5)):
        if not row.get(f"chrom_{suffix}"):
            continue
        consequence = row.get(f"consequence_{suffix}", "")
        impact = row.get(f"impact_{suffix}", "")
        hgvsc = row.get(f"hgvsc_{suffix}", "")

        # Degrada las anotaciones de splicing que no son creibles por
        # posicion, para que no reciban credito de impacto alto.
        if not splice_annotation_is_credible(consequence, hgvsc):
            notes.append(f"splicing_intronico_profundo_{suffix}")
            impact = "MODIFIER"
            consequence = "intron_variant"

        for component in (
            impact_component(impact, consequence),
            clinvar_component(row.get(f"clin_sig_{suffix}", "")),
            frequency_component(row.get(f"max_af_{suffix}", "-")),
            predictor_component(row.get(f"sift_{suffix}", "-"),
                                row.get(f"polyphen_{suffix}", "-")),
        ):
            value, tags = component
            score += value * weight
            notes.extend(f"{t}_{suffix}" for t in tags)

    if use_panels:
        value, tags = panel_component(row.get("gene", ""))
        score += value
        notes.extend(tags)

    if "variantes_contiguas" in row.get("flags", ""):
        score -= 20.0
        notes.append("variantes_contiguas")

    return score, notes


def score_to_epcr(score: float, best: float) -> float:
    """
    Convierte la puntuacion en una probabilidad estimada dentro de (0, 1].

    El sistema de evaluacion solo usa el orden relativo, asi que basta
    una transformacion monotona. Se normaliza contra el mejor candidato
    y se acota para no declarar certeza absoluta.
    """
    if best <= 0:
        return 0.01
    value = max(score, 0.0) / best
    return round(min(0.95, max(0.01, value * 0.95)), 4)


def format_chrom(chrom: str) -> str:
    if not chrom:
        return ""
    if config.ADD_CHR_PREFIX_ON_SUBMIT and not chrom.startswith("chr"):
        return f"chr{chrom}"
    return chrom


def build_note(row: dict[str, str]) -> str:
    parts = [row.get("gene", "")]
    model = row.get("model", "")
    if model == "compound_het":
        parts.append("compound het")
        phase = row.get("phase", "unknown")
        if phase != "unknown":
            parts.append(f"phase {phase}")
    elif model == "homozygous":
        parts.append("homozygous")

    for suffix in ("1", "2"):
        hgvsp = row.get(f"hgvsp_{suffix}", "")
        consequence = row.get(f"consequence_{suffix}", "")
        if hgvsp and hgvsp != "-":
            protein = hgvsp.split(":")[-1]
            parts.append(f"{protein} ({consequence})")
        elif consequence:
            parts.append(consequence)

    clin = row.get("clin_sig_1", "")
    if clin and clin != "-":
        parts.append(f"ClinVar: {clin}")

    return "; ".join(p for p in parts if p)[:250]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", type=Path, default=config.CANDIDATES_TSV)
    parser.add_argument("--ranked-output", type=Path, default=config.RANKED_TSV)
    parser.add_argument("--submission-output", type=Path, default=config.SUBMISSION_CSV)
    parser.add_argument("--proband-id", default=config.PROBAND_ID)
    parser.add_argument("--top-n", type=int, default=10,
                        help="filas en el CSV de entrega (maximo 10)")
    parser.add_argument("--no-panels", action="store_true",
                        help="ranking sin conocimiento previo de genes")
    parser.add_argument("--output-suffix", default="",
                        help="sufijo para no sobrescribir salidas previas")
    args = parser.parse_args()

    config.ensure_directories()
    logger = get_logger("05_rank", config.LOGS_DIR / "05_rank.log")

    if args.proband_id == "PROBAND01":
        logger.warning("Usando el proband_id del template. Confirmar el "
                       "identificador real en la plataforma antes de enviar.")

    use_panels = config.USE_GENE_PANELS and not args.no_panels
    logger.info("Paneles de genes: %s", "activados" if use_panels else "desactivados")

    with open(args.input) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    logger.info("Candidatos leidos: %d", len(rows))

    scored = []
    for row in rows:
        score, notes = score_row(row, use_panels)
        record = dict(row)
        record["score"] = f"{score:.2f}"
        record["evidence"] = ";".join(notes)
        scored.append((score, record))

    scored.sort(key=lambda x: -x[0])

    ranked_path = args.ranked_output
    submission_path = args.submission_output
    if args.output_suffix:
        ranked_path = ranked_path.with_name(
            f"{ranked_path.stem}_{args.output_suffix}{ranked_path.suffix}")
        submission_path = submission_path.with_name(
            f"{submission_path.stem}_{args.output_suffix}{submission_path.suffix}")

    if scored:
        fieldnames = list(scored[0][1].keys())
        with open(ranked_path, "w", newline="") as out:
            writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            for _, record in scored:
                writer.writerow(record)
    logger.info("Ranking completo en %s", ranked_path)

    top = scored[:args.top_n]
    best_score = top[0][0] if top else 0.0

    with open(submission_path, "w", newline="") as out:
        writer = csv.writer(out)
        writer.writerow(["proband_id", "chrom_1", "pos_1", "ref_1", "alt_1",
                         "chrom_2", "pos_2", "ref_2", "alt_2",
                         "epcr", "finding_type", "notes"])
        for score, row in top:
            gene = row.get("gene", "")
            finding = "primary" if (
                gene in config.MVA_CORE_GENES
                or gene in config.MITOTIC_CHECKPOINT_GENES
            ) else "secondary"
            writer.writerow([
                args.proband_id,
                format_chrom(row.get("chrom_1", "")),
                row.get("pos_1", ""),
                row.get("ref_1", ""),
                row.get("alt_1", ""),
                format_chrom(row.get("chrom_2", "")),
                row.get("pos_2", ""),
                row.get("ref_2", ""),
                row.get("alt_2", ""),
                score_to_epcr(score, best_score),
                finding,
                build_note(row),
            ])

    logger.info("Entrega escrita en %s", submission_path)
    logger.info("Diez primeros candidatos:")
    for i, (score, row) in enumerate(top, 1):
        logger.info("  %2d. %-12s %-14s score=%7.2f  %s",
                    i, row.get("gene", ""), row.get("model", ""), score,
                    row.get("hgvsp_1", "")[:40])


if __name__ == "__main__":
    main()