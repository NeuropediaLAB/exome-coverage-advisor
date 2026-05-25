"""
Exome Coverage Advisor – Backend v2
Datos: gnomAD r4 + MyGene.info (online). Kits: BED files aportados por el usuario.
"""

from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
from typing import Dict, List, Optional, Tuple
import bisect
import logging
import requests as http_requests
from functools import lru_cache

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

app = FastAPI(
    title="Exome Coverage Advisor API",
    description="Cobertura de exomas – gnomAD r4 + MyGene.info + BED files de kits",
    version="2.1.0",
)
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

GNOMAD_API = "https://gnomad.broadinstitute.org/api"
MYGENE_API = "https://mygene.info/v3"

# ── Kit BED registry (in-process memory) ─────────────────────────────────────

class KitBedIndex:
    """Efficient interval overlap index for a single exome kit BED file."""

    def __init__(self, name: str):
        self.name    = name
        self.regions: Dict[str, List[Tuple[int, int]]] = {}   # chrom → sorted [(start, stop)]
        self.n_regions = 0

    def load(self, bed_text: str):
        tmp: Dict[str, List[Tuple[int, int]]] = {}
        for line in bed_text.splitlines():
            if not line or line.startswith("#") or line.startswith("track") or line.startswith("browser"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            chrom = parts[0].lstrip("chr")
            try:
                start, stop = int(parts[1]), int(parts[2])
            except ValueError:
                continue
            tmp.setdefault(chrom, []).append((start, stop))
        for chrom, ivls in tmp.items():
            self.regions[chrom] = sorted(ivls)
        self.n_regions = sum(len(v) for v in self.regions.values())

    def covered_fraction(self, chrom: str, start: int, stop: int) -> float:
        """Fraction of [start, stop) covered by kit capture regions."""
        chrom = chrom.lstrip("chr")
        ivls  = self.regions.get(chrom, [])
        if not ivls:
            return 0.0
        total = 0
        # Binary search to skip irrelevant intervals
        lo = bisect.bisect_left(ivls, (start,))
        lo = max(0, lo - 1)
        for ivl_start, ivl_stop in ivls[lo:]:
            if ivl_start >= stop:
                break
            ov = min(stop, ivl_stop) - max(start, ivl_start)
            if ov > 0:
                total += ov
        return total / max(1, stop - start)

# Global kit registry
KIT_REGISTRY: Dict[str, KitBedIndex] = {}


# ── Models ────────────────────────────────────────────────────────────────────

class KitRegisterRequest(BaseModel):
    name:        str
    bed_content: str

class KitRegisterUrlRequest(BaseModel):
    name: str
    url:  str

class GeneInfo(BaseModel):
    gene:         str
    full_name:    str
    omim:         str
    chromosome:   str
    summary:      str
    pLI:          Optional[float] = None
    oe_lof_upper: Optional[float] = None
    lof_z:        Optional[float] = None

class ExonCoverage(BaseModel):
    exon:          int
    chrom:         str
    start:         int
    stop:          int
    length:        int
    mean_coverage: float
    pct_20x:       float
    pct_30x:       float
    status:        str    # excellent | good | low | poor | no_data
    kit_captured:  Optional[float] = None   # fraction captured by kit (None = no kit filter)

class TechnicalNote(BaseModel):
    category:       str
    severity:       str   # info | warning | critical
    message:        str
    recommendation: str

class CoverageReport(BaseModel):
    gene:           GeneInfo
    exons:          List[ExonCoverage]
    notes:          List[TechnicalNote]
    recommendation: str
    data_source:    str = "gnomAD r4 (~800k exomas)"
    kit:            Optional[str] = None
    n_exons_total:  int = 0
    n_exons_good:   int = 0
    n_exons_low:    int = 0
    n_exons_poor:   int = 0
    n_exons_uncaptured: int = 0

class GeneSummary(BaseModel):
    gene:                str
    full_name:           str
    omim:                str
    chromosome:          str
    pLI:                 Optional[float] = None
    oe_lof_upper:        Optional[float] = None
    constraint_category: str
    found:               bool = True

class BatchRequest(BaseModel):
    genes: List[str]

class BatchEntry(BaseModel):
    gene_name: str
    summary:   Optional[GeneSummary] = None
    found:     bool = True

class MatrixRequest(BaseModel):
    genes: List[str]
    kits:  List[str]   # must be registered kit names

class MatrixEntry(BaseModel):
    gene_name: str
    kit_name:  str
    report:    Optional[CoverageReport] = None
    found:     bool = True


# ── Online fetchers (cached) ──────────────────────────────────────────────────

@lru_cache(maxsize=5000)
def _fetch_mygene(symbol: str) -> dict:
    try:
        r = http_requests.get(
            f"{MYGENE_API}/query",
            params={"q": f"symbol:{symbol}", "species": "human",
                    "fields": "symbol,name,MIM,summary", "size": 1},
            timeout=10,
        )
        hits = r.json().get("hits", [])
        return hits[0] if hits else {}
    except Exception as e:
        logger.warning(f"MyGene {symbol}: {e}")
        return {}


@lru_cache(maxsize=5000)
def _fetch_gnomad_gene(symbol: str) -> dict:
    query = """
    query Gene($sym: String!) {
      gene(gene_symbol: $sym, reference_genome: GRCh38) {
        symbol omim_id chrom start stop
        gnomad_constraint { pLI lof_z mis_z syn_z oe_lof oe_lof_upper }
        exons { feature_type start stop }
      }
    }
    """
    try:
        r = http_requests.post(
            GNOMAD_API,
            json={"query": query, "variables": {"sym": symbol}},
            timeout=15,
        )
        data = r.json()
        if "errors" in data or not (data.get("data") or {}).get("gene"):
            return {}
        return data["data"]["gene"]
    except Exception as e:
        logger.warning(f"gnomAD gene {symbol}: {e}")
        return {}


@lru_cache(maxsize=20000)
def _fetch_gnomad_region(chrom: str, start: int, stop: int) -> Tuple[float, float, float]:
    """Returns (mean_coverage, pct_20x, pct_30x) for a genomic region."""
    query = """
    query Reg($chrom: String!, $start: Int!, $stop: Int!) {
      region(chrom: $chrom, start: $start, stop: $stop, reference_genome: GRCh38) {
        coverage(dataset: gnomad_r4) {
          exome { mean over_20 over_30 }
        }
      }
    }
    """
    try:
        r = http_requests.post(
            GNOMAD_API,
            json={"query": query, "variables": {"chrom": chrom, "start": start, "stop": stop}},
            timeout=20,
        )
        data = r.json()
        if "errors" in data:
            return (0.0, 0.0, 0.0)
        positions = (
            data.get("data", {}).get("region", {})
            .get("coverage", {}).get("exome", []) or []
        )
        if not positions:
            return (0.0, 0.0, 0.0)
        n = len(positions)
        return (
            round(sum(p["mean"]    for p in positions) / n, 1),
            round(sum(p["over_20"] for p in positions) / n * 100, 1),
            round(sum(p["over_30"] for p in positions) / n * 100, 1),
        )
    except Exception as e:
        logger.warning(f"gnomAD region {chrom}:{start}-{stop}: {e}")
        return (0.0, 0.0, 0.0)


def _constraint_category(pLI, oeup) -> str:
    if pLI is None and oeup is None:
        return "unknown"
    if (pLI or 0) >= 0.9 or (oeup or 1) <= 0.35:
        return "intolerant"
    if (pLI or 0) >= 0.5 or (oeup or 1) <= 0.6:
        return "moderate"
    return "tolerant"


def _gene_info_dict(symbol: str) -> Optional[dict]:
    mg = _fetch_mygene(symbol)
    gn = _fetch_gnomad_gene(symbol)
    if not mg.get("symbol") and not gn.get("symbol"):
        return None
    c = gn.get("gnomad_constraint") or {}
    ch = gn.get("chrom") or "?"
    return {
        "gene":         symbol,
        "full_name":    mg.get("name") or gn.get("symbol") or symbol,
        "omim":         str(mg.get("MIM") or gn.get("omim_id") or "N/A"),
        "chromosome":   f"chr{ch}" if ch != "?" else "?",
        "summary":      (mg.get("summary") or "")[:500],
        "pLI":          c.get("pLI"),
        "oe_lof_upper": c.get("oe_lof_upper"),
        "lof_z":        c.get("lof_z"),
    }


def _build_report(symbol: str, kit_name: Optional[str] = None) -> Optional[dict]:
    """
    Full per-exon coverage report.
    If kit_name is provided and registered, each exon also shows kit_captured fraction.
    """
    info = _gene_info_dict(symbol)
    if info is None:
        return None

    gnomad    = _fetch_gnomad_gene(symbol)
    chrom     = gnomad.get("chrom") or ""
    cds_exons = [e for e in (gnomad.get("exons") or []) if e.get("feature_type") == "CDS"]
    kit_idx   = KIT_REGISTRY.get(kit_name) if kit_name else None

    exons  = []
    n_good = n_low = n_poor = n_uncap = 0

    for i, exon in enumerate(cds_exons, 1):
        start, stop = exon["start"], exon["stop"]

        kit_frac = kit_idx.covered_fraction(chrom, start, stop) if kit_idx else None

        # If kit is specified and exon is not captured → mark as uncaptured, skip gnomAD call
        if kit_frac is not None and kit_frac == 0.0:
            n_uncap += 1
            exons.append({
                "exon": i, "chrom": chrom, "start": start, "stop": stop,
                "length": stop - start,
                "mean_coverage": 0.0, "pct_20x": 0.0, "pct_30x": 0.0,
                "status": "uncaptured", "kit_captured": 0.0,
            })
            continue

        mean_cov, pct_20x, pct_30x = _fetch_gnomad_region(chrom, start, stop)

        if pct_20x >= 95:
            status = "excellent"; n_good += 1
        elif pct_20x >= 80:
            status = "good";      n_good += 1
        elif pct_20x >= 50:
            status = "low";       n_low  += 1
        elif mean_cov > 0:
            status = "poor";      n_poor += 1
        else:
            status = "no_data"

        exons.append({
            "exon": i, "chrom": chrom, "start": start, "stop": stop,
            "length": stop - start,
            "mean_coverage": mean_cov, "pct_20x": pct_20x, "pct_30x": pct_30x,
            "status": status, "kit_captured": kit_frac,
        })

    # Notes
    notes = []
    pLI  = info.get("pLI") or 0
    oeup = info.get("oe_lof_upper") or 1
    if pLI > 0.9 or oeup < 0.35:
        notes.append({
            "category": "Constraint génico", "severity": "warning",
            "message":  f"pLI={pLI:.3f} · LOEUF={oeup:.3f} — Gen altamente intolerante a LoF.",
            "recommendation": "Asegurar cobertura completa. NGS ortogonal para exones con baja cobertura.",
        })
    if n_uncap > 0:
        notes.append({
            "category": "Captura del kit", "severity": "critical",
            "message":  f"{n_uncap} exón(es) CDS NO capturado(s) por el kit '{kit_name}'.",
            "recommendation": "Estos exones requieren un método ortogonal (panel dirigido, Sanger, WGS).",
        })
    if not cds_exons:
        notes.append({
            "category": "Estructura génica", "severity": "info",
            "message":  "No se encontraron exones CDS en gnomAD r4.",
            "recommendation": "Verificar símbolo del gen en HGNC/Ensembl.",
        })
    notes.append({
        "category": "Fuente", "severity": "info",
        "message":  "Cobertura: gnomAD r4 (~800k exomas, múltiples kits). Captura: BED del kit.",
        "recommendation": "Los valores de cobertura reflejan la media poblacional en regiones capturadas.",
    })

    n_total = len(exons)
    if n_total == 0:
        rec = "⚠️ SIN DATOS — Gen no encontrado en gnomAD r4."
    elif n_uncap == n_total:
        rec = f"❌ GEN NO CAPTURADO — Ningún exón CDS cubierto por el kit '{kit_name}'."
    elif n_uncap > 0:
        rec = f"❌ CAPTURA INCOMPLETA — {n_uncap}/{n_total} exones no capturados. Método ortogonal necesario."
    elif n_poor > 0:
        rec = f"❌ COBERTURA INSUFICIENTE — {n_poor}/{n_total} exón(es) <50% bases ≥20x."
    elif n_low > 0:
        rec = f"⚠️ COBERTURA PARCIAL — {n_low}/{n_total} exón(es) con cobertura subóptima."
    else:
        rec = f"✅ COBERTURA ADECUADA — {n_good}/{n_total} exones bien cubiertos."

    return {
        "gene": info, "exons": exons, "notes": notes,
        "recommendation": rec, "data_source": "gnomAD r4 (~800k exomas)",
        "kit": kit_name,
        "n_exons_total": n_total, "n_exons_good": n_good,
        "n_exons_low": n_low, "n_exons_poor": n_poor,
        "n_exons_uncaptured": n_uncap,
    }


# ── Endpoints ─────────────────────────────────────────────────────────────────

@app.get("/")
def root():
    return {"status": "ok", "service": "Exome Coverage Advisor API v2.1",
            "data_source": "gnomAD r4 + MyGene.info", "kits_registered": len(KIT_REGISTRY)}


# -- Kit management -----------------------------------------------------------

@app.post("/kits/register")
def register_kit(req: KitRegisterRequest):
    """Register a kit BED file (uploaded by user). Stored in process memory."""
    kit = KitBedIndex(req.name)
    kit.load(req.bed_content)
    if kit.n_regions == 0:
        raise HTTPException(status_code=400, detail="BED file vacío o sin regiones válidas.")
    KIT_REGISTRY[req.name] = kit
    logger.info(f"Kit registered: '{req.name}' with {kit.n_regions} regions")
    return {"name": req.name, "regions": kit.n_regions, "chroms": len(kit.regions)}


@app.post("/kits/register-url")
def register_kit_from_url(req: KitRegisterUrlRequest):
    """Download a BED file from a URL and register it as a kit."""
    MAX_BYTES = 200 * 1024 * 1024  # 200 MB
    try:
        r = http_requests.get(req.url, timeout=60, stream=True,
                              headers={"User-Agent": "ExomeCoverageAdvisor/2.2"})
        r.raise_for_status()
    except Exception as e:
        raise HTTPException(status_code=400, detail=f"Error al descargar el BED: {e}")

    chunks = []
    total = 0
    for chunk in r.iter_content(chunk_size=65536):
        total += len(chunk)
        if total > MAX_BYTES:
            raise HTTPException(status_code=413, detail="El archivo BED supera el límite de 200 MB.")
        chunks.append(chunk)

    raw = b"".join(chunks)

    # Decompress gzip if needed
    if req.url.endswith(".gz") or raw[:2] == b"\x1f\x8b":
        import gzip
        try:
            raw = gzip.decompress(raw)
        except Exception as e:
            raise HTTPException(status_code=400, detail=f"Error al descomprimir gzip: {e}")

    # Handle interval_list format (GATK / Picard) → convert to BED
    text = raw.decode("utf-8", errors="replace")
    if "@HD" in text[:500] or "@SQ" in text[:500]:
        lines = []
        for line in text.splitlines():
            if line.startswith("@"):
                continue
            parts = line.split("\t")
            if len(parts) >= 3:
                chrom, start, stop = parts[0], int(parts[1]) - 1, int(parts[2])
                lines.append(f"{chrom}\t{start}\t{stop}")
        text = "\n".join(lines)

    kit = KitBedIndex(req.name)
    kit.load(text)
    if kit.n_regions == 0:
        raise HTTPException(status_code=400, detail="El archivo descargado no contiene regiones BED válidas.")
    KIT_REGISTRY[req.name] = kit
    logger.info(f"Kit registered from URL: '{req.name}' ({kit.n_regions} regions) from {req.url}")
    return {"name": req.name, "regions": kit.n_regions, "chroms": len(kit.regions), "source_url": req.url}


@app.post("/kits/generate-gencode")
def generate_gencode_bed():
    """
    Genera un BED de regiones CDS de Gencode hg38 (básico) desde UCSC Table Browser.
    Útil como referencia pública cuando no se dispone del BED del kit.
    """
    UCSC_URL = (
        "https://api.genome.ucsc.edu/getData/track"
        "?genome=hg38&track=knownGene&fields=chrom,chromStart,chromEnd,name,score,strand,thickStart,thickEnd"
        "&maxItemsOutput=500000"
    )
    try:
        r = http_requests.get(UCSC_URL, timeout=60)
        r.raise_for_status()
        data = r.json()
    except Exception as e:
        raise HTTPException(status_code=502, detail=f"Error al contactar UCSC: {e}")

    rows = data.get("knownGene", [])
    bed_lines = []
    for row in rows:
        chrom = row.get("chrom", "")
        start = row.get("thickStart", row.get("chromStart", 0))
        stop  = row.get("thickEnd",  row.get("chromEnd",   0))
        if stop > start and chrom:
            bed_lines.append(f"{chrom}\t{start}\t{stop}")

    text = "\n".join(bed_lines)
    kit  = KitBedIndex("Gencode hg38 CDS (referencia)")
    kit.load(text)
    if kit.n_regions == 0:
        raise HTTPException(status_code=502, detail="UCSC no devolvió regiones válidas.")
    KIT_REGISTRY["Gencode hg38 CDS (referencia)"] = kit
    return {"name": "Gencode hg38 CDS (referencia)", "regions": kit.n_regions, "chroms": len(kit.regions)}


@app.delete("/kits/{kit_name}")
def delete_kit(kit_name: str):
    if kit_name not in KIT_REGISTRY:
        raise HTTPException(status_code=404, detail=f"Kit '{kit_name}' no registrado.")
    del KIT_REGISTRY[kit_name]
    return {"deleted": kit_name}


@app.get("/kits")
def list_kits():
    return [
        {"name": name, "regions": idx.n_regions, "chroms": len(idx.regions)}
        for name, idx in KIT_REGISTRY.items()
    ]


# -- Gene info ----------------------------------------------------------------

@app.get("/genes")
def get_genes():
    """No hay lista fija – todo el genoma disponible vía gnomAD/MyGene.info."""
    return []


@app.get("/gene/{gene_symbol}", response_model=GeneInfo)
def get_gene_info(gene_symbol: str):
    info = _gene_info_dict(gene_symbol.upper())
    if not info:
        raise HTTPException(status_code=404, detail=f"Gen {gene_symbol} no encontrado")
    return GeneInfo(**info)


# -- Coverage -----------------------------------------------------------------

@app.get("/coverage/{gene_symbol}", response_model=CoverageReport)
@app.get("/coverage/{gene_symbol}/{kit_name}", response_model=CoverageReport)
def get_coverage(gene_symbol: str, kit_name: Optional[str] = None):
    """
    Per-exon coverage from gnomAD r4.
    If kit_name is given, each exon also shows capture fraction from the registered BED.
    """
    if kit_name and kit_name not in KIT_REGISTRY:
        raise HTTPException(status_code=404, detail=f"Kit '{kit_name}' no registrado. Sube el BED primero.")
    report = _build_report(gene_symbol.upper(), kit_name)
    if not report:
        raise HTTPException(status_code=404, detail=f"Gen {gene_symbol} no encontrado")
    return CoverageReport(
        gene=GeneInfo(**report["gene"]),
        exons=[ExonCoverage(**e) for e in report["exons"]],
        notes=[TechnicalNote(**n) for n in report["notes"]],
        recommendation=report["recommendation"],
        data_source=report["data_source"],
        kit=report["kit"],
        n_exons_total=report["n_exons_total"],
        n_exons_good=report["n_exons_good"],
        n_exons_low=report["n_exons_low"],
        n_exons_poor=report["n_exons_poor"],
        n_exons_uncaptured=report["n_exons_uncaptured"],
    )


@app.post("/summary/batch", response_model=List[BatchEntry])
def get_summary_batch(request: BatchRequest):
    """Lightweight summary (no per-exon data) for N genes — fast."""
    results = []
    for gene in request.genes:
        sym = gene.strip().upper()
        if not sym:
            continue
        mg = _fetch_mygene(sym)
        gn = _fetch_gnomad_gene(sym)
        if not mg.get("symbol") and not gn.get("symbol"):
            results.append(BatchEntry(gene_name=sym, summary=None, found=False))
            continue
        c    = (gn.get("gnomad_constraint") or {})
        pLI  = c.get("pLI")
        oeup = c.get("oe_lof_upper")
        ch   = gn.get("chrom") or "?"
        results.append(BatchEntry(
            gene_name=sym, found=True,
            summary=GeneSummary(
                gene=sym,
                full_name=mg.get("name") or gn.get("symbol") or sym,
                omim=str(mg.get("MIM") or gn.get("omim_id") or "N/A"),
                chromosome=f"chr{ch}" if ch != "?" else "?",
                pLI=pLI, oe_lof_upper=oeup,
                constraint_category=_constraint_category(pLI, oeup),
            ),
        ))
    return results


@app.post("/coverage/matrix", response_model=List[MatrixEntry])
def get_coverage_matrix(request: MatrixRequest):
    """
    genes × kits matrix. Each kit must be registered via /kits/register.
    """
    results = []
    for gene in request.genes:
        sym = gene.strip().upper()
        for kit in request.kits:
            if kit not in KIT_REGISTRY:
                results.append(MatrixEntry(gene_name=sym, kit_name=kit, report=None, found=False))
                continue
            report = _build_report(sym, kit)
            if not report:
                results.append(MatrixEntry(gene_name=sym, kit_name=kit, report=None, found=False))
            else:
                results.append(MatrixEntry(
                    gene_name=sym, kit_name=kit, found=True,
                    report=CoverageReport(
                        gene=GeneInfo(**report["gene"]),
                        exons=[ExonCoverage(**e) for e in report["exons"]],
                        notes=[TechnicalNote(**n) for n in report["notes"]],
                        recommendation=report["recommendation"],
                        data_source=report["data_source"],
                        kit=report["kit"],
                        n_exons_total=report["n_exons_total"],
                        n_exons_good=report["n_exons_good"],
                        n_exons_low=report["n_exons_low"],
                        n_exons_poor=report["n_exons_poor"],
                        n_exons_uncaptured=report["n_exons_uncaptured"],
                    ),
                ))
    return results


if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)
