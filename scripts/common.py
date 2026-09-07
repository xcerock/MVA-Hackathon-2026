#!/usr/bin/env python3
"""
Utilidades compartidas del pipeline.

Dos piezas requieren cuidado especial:

1. La salida tabular de VEP cambia de numero de columnas segun los flags
   usados (por ejemplo --hgvs anade HGVSc, HGVSp y HGVS_OFFSET), asi que
   todo acceso se hace por nombre de columna y nunca por indice fijo.

2. El campo Uploaded_variation contiene el ID del VCF cuando existe. Este
   VCF fue llamado con dbSNP como sitios conocidos, asi que la gran
   mayoria de las lineas traen un rsID en lugar del formato
   "chrom_pos_ref/alt". La identificacion de variantes usa por tanto la
   columna Location mas el alelo alterno, y recurre a Uploaded_variation
   solo como respaldo.
"""
from __future__ import annotations

import gzip
import logging
import subprocess
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------


def get_logger(name: str, log_file: Path | None = None) -> logging.Logger:
    logger = logging.getLogger(name)
    if logger.handlers:
        return logger
    logger.setLevel(logging.INFO)
    fmt = logging.Formatter("%(asctime)s [%(levelname)s] %(name)s: %(message)s",
                            datefmt="%Y-%m-%d %H:%M:%S")
    sh = logging.StreamHandler(sys.stderr)
    sh.setFormatter(fmt)
    logger.addHandler(sh)
    if log_file is not None:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        fh = logging.FileHandler(log_file)
        fh.setFormatter(fmt)
        logger.addHandler(fh)
    return logger


def run_command(cmd: list[str] | str, logger: logging.Logger,
                shell: bool = False) -> None:
    """Ejecuta un comando externo y aborta si falla."""
    printable = cmd if isinstance(cmd, str) else " ".join(map(str, cmd))
    logger.info("Ejecutando: %s", printable)
    result = subprocess.run(cmd, shell=shell)
    if result.returncode != 0:
        logger.error("El comando fallo con codigo %s", result.returncode)
        raise SystemExit(result.returncode)


def open_maybe_gzip(path: Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


# ---------------------------------------------------------------------------
# Lectura de la salida tabular de VEP
# ---------------------------------------------------------------------------


class VepTable:
    """Itera la salida --tab de VEP devolviendo diccionarios por nombre."""

    def __init__(self, path: Path):
        self.path = Path(path)
        self.columns: list[str] = []

    def __iter__(self) -> Iterator[dict[str, str]]:
        with open_maybe_gzip(self.path) as fh:
            for line in fh:
                if line.startswith("##"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if not self.columns:
                    self.columns = [c.lstrip("#") for c in fields]
                    continue
                yield dict(zip(self.columns, fields))

    def header(self) -> list[str]:
        if not self.columns:
            for _ in self:
                break
        return self.columns


def parse_location(value: str) -> tuple[str, int, int] | None:
    """
    Convierte el campo Location de VEP en (chrom, inicio, fin).

    Formatos posibles:
        "15:40209701"          sustitucion o insercion puntual
        "15:40209701-40209705" rango, para deleciones
        "15:40209702-40209701" rango invertido, para inserciones
    """
    if ":" not in value:
        return None
    chrom, _, rest = value.partition(":")
    if "-" in rest:
        start_s, _, end_s = rest.partition("-")
    else:
        start_s = end_s = rest
    try:
        return chrom, int(start_s), int(end_s)
    except ValueError:
        return None


def parse_uploaded_variation(value: str) -> tuple[str, str, str, str] | None:
    """
    Convierte Uploaded_variation en (chrom, pos, ref, alt).

    Solo funciona cuando VEP genero la clave a partir de coordenadas
    ("15_40209701_T/G"). Devuelve None cuando el campo contiene un rsID,
    que es el caso mayoritario en este dataset.
    """
    parts = value.split("_")
    if len(parts) < 3:
        return None
    pos = parts[-2]
    alleles = parts[-1]
    if "/" not in alleles or not pos.isdigit():
        return None
    chrom = "_".join(parts[:-2])
    ref, alt = alleles.split("/", 1)
    return chrom, pos, ref, alt


def variant_key(chrom: str, pos: str | int, ref: str, alt: str) -> str:
    return f"{chrom}:{pos}:{ref}:{alt}"


def position_key(chrom: str, pos: str | int) -> str:
    return f"{chrom}:{pos}"


def vep_allele_matches(vep_allele: str, ref: str, alt: str) -> bool:
    """
    Comprueba si el alelo que reporta VEP corresponde a un alelo del VCF.

    VEP normaliza a representacion minima: elimina la base ancla comun
    que el VCF conserva en indels. Una insercion "T -> TGGGG" en el VCF
    aparece en VEP como "GGGG", y una delecion "TGGG -> T" aparece como
    "-". Esta funcion cubre los tres casos.
    """
    if not vep_allele:
        return False
    if vep_allele == alt:
        return True

    # Insercion: VEP omite la base ancla inicial.
    if len(alt) > len(ref) and alt.startswith(ref):
        if vep_allele == alt[len(ref):]:
            return True

    # Delecion: VEP la representa con guion o con las bases eliminadas.
    if len(ref) > len(alt) and ref.startswith(alt):
        if vep_allele == "-":
            return True
        if vep_allele == ref[len(alt):]:
            return True

    return False


# ---------------------------------------------------------------------------
# Lectura de genotipos desde el VCF
# ---------------------------------------------------------------------------


@dataclass
class Genotype:
    chrom: str
    pos: str
    ref: str
    alt: str
    gt: str
    dp: int | None
    gq: int | None
    ad_ref: int | None
    ad_alt: int | None
    pid: str
    pgt: str

    @property
    def is_het(self) -> bool:
        return self.gt in ("0/1", "1/0")

    @property
    def is_hom_alt(self) -> bool:
        return self.gt == "1/1"

    @property
    def vaf(self) -> float | None:
        """Fraccion de lecturas que soportan el alelo alterno."""
        if self.ad_ref is None or self.ad_alt is None:
            return None
        total = self.ad_ref + self.ad_alt
        if total == 0:
            return None
        return self.ad_alt / total


def _to_int(value: str) -> int | None:
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


class GenotypeIndex:
    """
    Genotipos del VCF indexados por posicion y por clave completa.

    El indice por posicion permite localizar variantes cuando VEP solo
    aporta Location; en las posiciones con mas de una entrada (sitios
    multialelicos separados durante la normalizacion) se desambigua
    comparando el alelo alterno.
    """

    def __init__(self) -> None:
        self.by_key: dict[str, Genotype] = {}
        self.by_position: dict[str, list[Genotype]] = defaultdict(list)

    def add(self, genotype: Genotype) -> None:
        key = variant_key(genotype.chrom, genotype.pos,
                          genotype.ref, genotype.alt)
        self.by_key[key] = genotype
        self.by_position[position_key(genotype.chrom, genotype.pos)].append(genotype)

    def __len__(self) -> int:
        return len(self.by_key)

    def lookup_exact(self, chrom: str, pos: str, ref: str,
                     alt: str) -> Genotype | None:
        return self.by_key.get(variant_key(chrom, pos, ref, alt))

    def lookup_by_location(self, chrom: str, start: int, end: int,
                           vep_allele: str) -> Genotype | None:
        """
        Localiza la variante a partir del campo Location de VEP.

        VEP desplaza la coordenada de inicio en los indels porque elimina
        la base ancla, asi que se prueban la posicion informada y la
        inmediatamente anterior.
        """
        seen: set[int] = set()
        for candidate_pos in (start, start - 1, end, end - 1):
            if candidate_pos in seen or candidate_pos < 1:
                continue
            seen.add(candidate_pos)
            entries = self.by_position.get(position_key(chrom, candidate_pos))
            if not entries:
                continue
            for entry in entries:
                if vep_allele_matches(vep_allele, entry.ref, entry.alt):
                    return entry
            if len(entries) == 1 and vep_allele in ("", "-"):
                return entries[0]
        return None


def load_genotypes(vcf_path: Path, sample_index: int = 0) -> GenotypeIndex:
    """Carga los genotipos del VCF normalizado en un indice consultable."""
    index = GenotypeIndex()
    with open_maybe_gzip(vcf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 10:
                continue
            chrom, pos, ref, alt = p[0], p[1], p[3], p[4]
            fmt = p[8].split(":")
            values = p[9 + sample_index].split(":")
            d = dict(zip(fmt, values))

            ad_ref = ad_alt = None
            ad_raw = d.get("AD", "")
            if "," in ad_raw:
                pieces = ad_raw.split(",")
                ad_ref = _to_int(pieces[0])
                ad_alt = _to_int(pieces[1])

            index.add(Genotype(
                chrom=chrom, pos=pos, ref=ref, alt=alt,
                gt=d.get("GT", "./.").replace("|", "/"),
                dp=_to_int(d.get("DP", "")),
                gq=_to_int(d.get("GQ", "")),
                ad_ref=ad_ref, ad_alt=ad_alt,
                pid=d.get("PID", "."), pgt=d.get("PGT", "."),
            ))
    return index


# ---------------------------------------------------------------------------
# Fase fisica (PID / PGT del caller)
# ---------------------------------------------------------------------------


def resolve_phase(a_pid: str, a_pgt: str, b_pid: str, b_pgt: str) -> str:
    """
    Determina la relacion de fase entre dos variantes heterocigotas.

    El caller asigna el mismo PID a variantes que comparten lecturas y
    expresa el haplotipo en PGT. Cuando dos variantes comparten PID:
      - PGT identico  -> mismo cromosoma (cis)
      - PGT distinto  -> cromosomas opuestos (trans)

    Sin PID compartido la fase es desconocida, lo habitual cuando las
    variantes distan mas que la longitud del fragmento secuenciado.
    """
    if a_pid == "." or b_pid == "." or a_pid != b_pid:
        return "unknown"
    if a_pgt == "." or b_pgt == ".":
        return "unknown"
    return "cis" if a_pgt == b_pgt else "trans"


# ---------------------------------------------------------------------------
# Parseo de campos de anotacion
# ---------------------------------------------------------------------------


def parse_af(value: str) -> float | None:
    """MAX_AF vale '-' cuando la variante no aparece en gnomAD."""
    if value in ("-", "", None):
        return None
    try:
        return float(value)
    except ValueError:
        return None


def parse_prediction(value: str) -> tuple[str | None, float | None]:
    """
    SIFT y PolyPhen vienen como 'deleterious(0.01)' o 'benign(0.005)'.
    Devuelve (etiqueta, score).
    """
    if not value or value == "-":
        return None, None
    label = value
    score = None
    if "(" in value and value.endswith(")"):
        label, _, rest = value.partition("(")
        try:
            score = float(rest.rstrip(")"))
        except ValueError:
            score = None
    return label.strip(), score


def consequence_list(value: str) -> list[str]:
    return [c for c in value.split(",") if c] if value else []


def has_splice_consequence(value: str, keywords: tuple[str, ...]) -> bool:
    return any(k in value for k in keywords)