# MVA Hackathon 2026, Track 1: Variant Prediction

Pipeline reproducible para identificar la variante o variantes causales
a partir del genoma completo del proband y su descripcion fenotipica.

## Enfoque

El nucleo del pipeline es **agnostico respecto al gen**: no filtra por
panel ni por hipotesis diagnostica previa. Recorre las cerca de 4,7
millones de variantes del genoma, las anota funcionalmente, aplica el
modelo de herencia recesiva compatible con el MVA y ordena los
candidatos por evidencia acumulada.

El conocimiento previo sobre genes de MVA se usa unicamente en la fase
final de puntuacion, y es desactivable con `--no-panels`. Reportar el
ranking en ambas modalidades permite demostrar que el candidato principal
emerge del procedimiento y no de haber sabido de antemano donde mirar.

## Estructura

```
scripts/
  config.py               parametros, umbrales, paneles y rutas
  common.py               lectura de VEP y VCF, resolucion de fase
  01_normalize.py         filtro PASS, cromosomas primarios, normalizacion
  02_annotate.py          anotacion funcional con Ensembl VEP
  03_filter_variants.py   filtrado por frecuencia, consecuencia y calidad
  04_recessive_model.py   homocigotos y heterocigosis compuesta
  05_rank_submission.py   puntuacion, ranking y CSV de entrega
  run_pipeline.py         orquestador con registro de tiempos
```

Los datos del paciente, los archivos intermedios y las bases de datos de
anotacion no se versionan. El `.gitignore` los excluye explicitamente.

## Requisitos

Herramientas externas:

- `bcftools` 1.24 o superior
- Ensembl VEP 116 con cache indexado de GRCh38
- Python 3.11 (solo biblioteca estandar)

Datos de referencia:

- GRCh38 primary assembly de Ensembl, sin prefijo `chr` en los contigs,
  para que coincida con la nomenclatura del VCF de entrada
- Cache de VEP para homo_sapiens GRCh38, release 116

Instalacion de los entornos:

```bash
micromamba create -n mva -c conda-forge -c bioconda \
  bcftools samtools htslib tabix bedtools python=3.11
micromamba create -n vep -c conda-forge -c bioconda ensembl-vep=116.2

micromamba activate vep
mkdir -p ~/mva/db/vep
aria2c -x 16 -s 16 -d ~/mva/db/vep \
  https://ftp.ensembl.org/pub/release-116/variation/indexed_vep_cache/homo_sapiens_vep_116_GRCh38.tar.gz
tar xzf ~/mva/db/vep/homo_sapiens_vep_116_GRCh38.tar.gz -C ~/mva/db/vep
```

## Ejecucion

```bash
export MVA_ROOT=~/mva
export MVA_PROBAND_ID=<identificador del dataset>
export MVA_THREADS=10

micromamba activate mva
python scripts/run_pipeline.py
```

Pasos individuales:

```bash
python scripts/run_pipeline.py --from 3     # desde el filtrado
python scripts/run_pipeline.py --only 5     # solo el ranking
python scripts/05_rank_submission.py --no-panels --output-suffix agnostic
```

## Decisiones metodologicas

**Alcance genomico.** El analisis se restringe a los cromosomas
primarios (1-22, X, Y). Se excluyen el mitocondrial, los contigs decoy
`Un_JTFH01*` y los scaffolds no colocados: carecen de anotacion genica
util y su nomenclatura no coincide con la referencia de Ensembl.

**Umbral de frecuencia.** Se conservan variantes con MAX_AF por debajo
del 1%, un corte deliberadamente permisivo. Algunos alelos hipomorficos
de genes de MVA alcanzan frecuencias apreciables en poblaciones
concretas, y un filtro mas estricto los eliminaria.

**Control de balance alelico.** Los heterocigotos deben mostrar una
fraccion de lecturas alternas entre 0,25 y 0,75. Este criterio descarta
clusters artefactuales: regiones donde el caller produce muchas variantes
consecutivas con fraccion alelica sistematicamente baja, firma tipica de
lecturas mal mapeadas procedentes de duplicaciones segmentales o
pseudogenes.

**Limite de heterocigotos por gen.** Un gen con mas de cinco
heterocigotos raros se marca como locus problematico y se excluye del
analisis de pares. Ademas de eliminar falsos positivos, evita la
explosion combinatoria: trece heterocigotos generan setenta y ocho
combinaciones.

**Resolucion de fase sin progenitores.** GATK HaplotypeCaller asigna los
campos PID y PGT a variantes que comparten lecturas. Cuando dos
heterocigotos comparten PID, un PGT identico indica que estan en el mismo
cromosoma (cis, incompatible con heterocigosis compuesta) y un PGT
distinto indica cromosomas opuestos (trans, compuesto confirmado). Esta
distincion elimina pares espurios que un pipeline sin trio suele dar por
validos.

**Distancia minima entre pares.** Variantes separadas por menos de 30
pares de bases se penalizan, ya que suelen representar la misma
alteracion descrita de dos formas por el caller.

## Salidas

| Archivo | Contenido |
|---|---|
| `work/proband.norm.vcf.gz` | VCF normalizado |
| `work/proband.vep.tsv` | anotacion funcional completa |
| `work/rare_functional.tsv` | variantes raras que superan el control de calidad |
| `results/recessive_candidates.tsv` | candidatos bajo modelo recesivo |
| `results/ranked_candidates.tsv` | ranking con puntuacion y evidencia |
| `results/track1_submission.csv` | entrega en el formato del reto |

## Limitaciones

El analisis parte del VCF proporcionado, generado con GATK 4.2.4.0 y
filtros duros sin recalibracion VQSR. Las variantes que el caller no
detecto quedan fuera de alcance, en particular las variantes
estructurales y las deleciones exonicas, que requieren volver a los
FASTQ.

La fase de los pares separados por mas que la longitud del fragmento no
puede resolverse con los campos PID y PGT y queda como desconocida. Su
confirmacion requiere faseo por lecturas largas o analisis de haplotipos.

La ausencia de datos de los progenitores impide distinguir variantes de
novo de las heredadas y confirmar la segregacion.

## Licencia

Los resultados, el codigo y los informes se publican bajo CC BY 4.0,
conforme a las condiciones del hackathon.
