"""
Exome Coverage Advisor – Frontend v2.2
Streamlit UI: genes libres + PanelApp + kits BED propios + comparativa por kit.
Incluye tutorial interactivo de cuellos de botella NGS en la primera página.
"""

import streamlit as st
import requests
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from io import StringIO

API = "http://exome-advisor-api:8000"
PANELAPP = "https://panelapp.genomicsengland.co.uk/api/v1"

st.set_page_config(
    page_title="Exome Coverage Advisor",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

# ── helpers ────────────────────────────────────────────────────────────────────

def api_get(path, **params):
    try:
        r = requests.get(f"{API}{path}", params=params, timeout=120)
        r.raise_for_status()
        return r.json()
    except Exception as e:
        return None

def api_post(path, payload):
    try:
        r = requests.post(f"{API}{path}", json=payload, timeout=120)
        r.raise_for_status()
        return r.json()
    except Exception as e:
        return None

def api_coverage(gene: str, kit: str = None):
    path = f"/coverage/{gene}" if not kit else f"/coverage/{gene}/{kit}"
    return api_get(path)

def api_batch_summary(genes: list):
    return api_post("/summary/batch", {"genes": genes})

def api_coverage_matrix(genes: list, kits: list):
    return api_post("/coverage/matrix", {"genes": genes, "kits": kits})

def api_register_kit(name: str, bed_content: str):
    try:
        r = requests.post(f"{API}/kits/register", json={"name": name, "bed_content": bed_content}, timeout=30)
        r.raise_for_status()
        return r.json()
    except Exception as e:
        return None

def api_register_kit_url(name: str, url: str):
    try:
        r = requests.post(f"{API}/kits/register-url", json={"name": name, "url": url}, timeout=90)
        r.raise_for_status()
        return r.json()
    except Exception as e:
        return None

def api_generate_gencode():
    try:
        r = requests.post(f"{API}/kits/generate-gencode", timeout=90)
        r.raise_for_status()
        return r.json()
    except Exception as e:
        return None


    try:
        r = requests.delete(f"{API}/kits/{name}", timeout=10)
        r.raise_for_status()
        return True
    except:
        return False

def api_list_kits():
    return api_get("/kits") or []

# ── status colours ─────────────────────────────────────────────────────────────

STATUS_COLOR = {
    "excellent":  "#1a7c3a",
    "good":       "#4caf50",
    "low":        "#ff9800",
    "poor":       "#f44336",
    "no_data":    "#9e9e9e",
    "uncaptured": "#607d8b",
}
STATUS_LABEL = {
    "excellent":  "✅ Excelente (≥95% @20x)",
    "good":       "🟢 Buena (≥80% @20x)",
    "low":        "🟡 Baja (≥50% @20x)",
    "poor":       "🔴 Pobre (<50% @20x)",
    "no_data":    "⚫ Sin datos gnomAD",
    "uncaptured": "🔵 No capturado por el kit",
}
CONSTRAINT_COLOR = {
    "intolerant": "#f44336",
    "moderate":   "#ff9800",
    "tolerant":   "#4caf50",
    "unknown":    "#9e9e9e",
}

# ── PanelApp ───────────────────────────────────────────────────────────────────

@st.cache_data(ttl=3600)
def fetch_panels():
    try:
        r = requests.get(f"{PANELAPP}/panels/?format=json&page_size=500", timeout=20)
        data = r.json()
        panels = data.get("results", [])
        # paginate
        nxt = data.get("next")
        while nxt:
            r   = requests.get(nxt, timeout=20)
            d   = r.json()
            panels.extend(d.get("results", []))
            nxt = d.get("next")
        return [
            {"id": p["id"], "name": p["name"],
             "version": p.get("version", "?"), "n_genes": p.get("stats", {}).get("number_of_genes", "?")}
            for p in panels
        ]
    except Exception as e:
        return []

@st.cache_data(ttl=3600)
def fetch_panel_genes(panel_id: int):
    genes = []
    url = f"{PANELAPP}/panels/{panel_id}/genes/?format=json&page_size=200"
    while url:
        try:
            r    = requests.get(url, timeout=20)
            data = r.json()
            for g in data.get("results", []):
                sym = g.get("gene_data", {}).get("gene_symbol") or g.get("entity_name")
                if sym:
                    genes.append(sym)
            url = data.get("next")
        except:
            break
    return sorted(set(genes))


# ── coverage chart (per-exon bar chart) ───────────────────────────────────────

def exon_chart(report: dict, title: str = ""):
    exons = report.get("exons", [])
    if not exons:
        st.info("Sin datos de exones CDS.")
        return

    labels  = [f"E{e['exon']}" for e in exons]
    pct_20  = [e["pct_20x"]       for e in exons]
    pct_30  = [e["pct_30x"]       for e in exons]
    means   = [e["mean_coverage"] for e in exons]
    statuses= [e["status"]        for e in exons]
    colors  = [STATUS_COLOR.get(s, "#9e9e9e") for s in statuses]
    captured= [e.get("kit_captured") for e in exons]

    fig = go.Figure()

    # Bar: pct @20x
    fig.add_trace(go.Bar(
        name="%@20x", x=labels, y=pct_20,
        marker_color=colors,
        customdata=list(zip(means, pct_30, statuses,
                            [f"{c*100:.0f}%" if c is not None else "N/A" for c in captured])),
        hovertemplate=(
            "<b>%{x}</b><br>"
            "Mean cov: %{customdata[0]:.1f}x<br>"
            "%@20x: %{y:.1f}%<br>"
            "%@30x: %{customdata[1]:.1f}%<br>"
            "Estado: %{customdata[2]}<br>"
            "Captura kit: %{customdata[3]}<extra></extra>"
        ),
    ))

    # Line: mean coverage (secondary axis)
    fig.add_trace(go.Scatter(
        name="Media cobertura (x)", x=labels, y=means,
        mode="lines+markers", line=dict(color="#1976d2", width=2),
        yaxis="y2", hoverinfo="skip",
    ))

    # Reference lines
    for y_val, label, color in [(80, "80%", "#ff9800"), (95, "95%", "#4caf50")]:
        fig.add_hline(y=y_val, line_dash="dash", line_color=color, annotation_text=label,
                      annotation_position="right")

    fig.update_layout(
        title=title or f"{report['gene']['gene']} — Cobertura por exón CDS",
        xaxis_title="Exón",
        yaxis=dict(title="% bases ≥ 20x", range=[0, 105]),
        yaxis2=dict(title="Media cobertura (x)", overlaying="y", side="right"),
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
        height=380, plot_bgcolor="#ffffff",
        bargap=0.2,
    )
    st.plotly_chart(fig, use_container_width=True)


def exon_table(report: dict):
    exons = report.get("exons", [])
    if not exons:
        return
    rows = []
    for e in exons:
        cap = e.get("kit_captured")
        rows.append({
            "Exón":      f"E{e['exon']}",
            "Coords":    f"chr{e['chrom']}:{e['start']:,}-{e['stop']:,}",
            "Long (pb)": e["length"],
            "Media x":   e["mean_coverage"],
            "% ≥20x":    e["pct_20x"],
            "% ≥30x":    e["pct_30x"],
            "Captura %": f"{cap*100:.0f}%" if cap is not None else "—",
            "Estado":    STATUS_LABEL.get(e["status"], e["status"]),
        })
    df = pd.DataFrame(rows)
    st.dataframe(df, use_container_width=True, hide_index=True)


# ── single-gene view ───────────────────────────────────────────────────────────

def show_single(report: dict, kit: str = None):
    gene = report["gene"]
    rec  = report["recommendation"]

    # Header
    col1, col2, col3, col4 = st.columns(4)
    col1.metric("Gen", gene["gene"])
    col2.metric("Cromosoma", gene["chromosome"])
    col3.metric("OMIM", gene["omim"])
    pLI  = gene.get("pLI")
    oeup = gene.get("oe_lof_upper")
    col4.metric("pLI / LOEUF",
                f"{pLI:.3f} / {oeup:.3f}" if pLI is not None else "N/A")

    if gene.get("summary"):
        st.caption(gene["summary"])

    st.info(rec)

    # Counters
    tot  = report.get("n_exons_total", 0)
    good = report.get("n_exons_good",  0)
    low  = report.get("n_exons_low",   0)
    poor = report.get("n_exons_poor",  0)
    unc  = report.get("n_exons_uncaptured", 0)
    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Total exones CDS", tot)
    c2.metric("✅ Buenos", good)
    c3.metric("🟡 Bajos", low)
    c4.metric("🔴 Pobres / No capturados", poor + unc)

    exon_chart(report)
    exon_table(report)

    # Notes
    with st.expander("📋 Notas técnicas", expanded=False):
        for note in report.get("notes", []):
            icon = {"critical": "🔴", "warning": "🟡", "info": "ℹ️"}.get(note["severity"], "")
            st.markdown(f"**{icon} {note['category']}**: {note['message']}")
            st.caption(f"→ {note['recommendation']}")


# ── multi-gene summary view ────────────────────────────────────────────────────

def show_multi_summary(entries: list, kit: str = None):
    # Summary table
    rows = []
    for e in entries:
        if not e.get("found"):
            rows.append({"Gen": e["gene_name"], "Estado": "❌ No encontrado",
                         "Cromosoma": "—", "OMIM": "—", "pLI": None, "LOEUF": None,
                         "Constraint": "—", "Nombre completo": "—"})
            continue
        s = e.get("summary") or {}
        rows.append({
            "Gen":           s.get("gene",   e["gene_name"]),
            "Nombre completo": s.get("full_name", "—"),
            "Cromosoma":     s.get("chromosome", "—"),
            "OMIM":          s.get("omim", "—"),
            "pLI":           s.get("pLI"),
            "LOEUF":         s.get("oe_lof_upper"),
            "Constraint":    s.get("constraint_category", "unknown"),
            "Estado":        "✅ Encontrado",
        })
    df = pd.DataFrame(rows)
    st.dataframe(df.drop(columns=["Estado"] if "Estado" in df.columns else []),
                 use_container_width=True, hide_index=True)

    # pLI / LOEUF chart
    plot_df = df[df["pLI"].notna()].copy()
    if len(plot_df) > 0:
        fig = px.bar(
            plot_df.sort_values("pLI", ascending=False),
            x="Gen", y="pLI", color="Constraint",
            color_discrete_map={
                "intolerant": "#f44336", "moderate": "#ff9800",
                "tolerant":   "#4caf50", "unknown":  "#9e9e9e",
            },
            hover_data=["Nombre completo", "LOEUF"],
            title="pLI por gen (restricción a variantes pérdida de función)",
            labels={"pLI": "pLI", "Constraint": "Nivel"},
        )
        fig.add_hline(y=0.9, line_dash="dash", line_color="#f44336",
                      annotation_text="pLI=0.9", annotation_position="right")
        fig.update_layout(height=350, plot_bgcolor="#ffffff")
        st.plotly_chart(fig, use_container_width=True)

    # Per-gene detail (lazy load)
    st.subheader("Detalle por gen (cobertura completa)")
    for e in entries:
        if not e.get("found"):
            st.warning(f"❌ {e['gene_name']}: gen no encontrado en gnomAD/MyGene.info")
            continue
        sym  = e["gene_name"]
        key  = f"report_{sym}_{kit or 'none'}"
        btn_key = f"btn_{sym}_{kit or 'none'}"
        with st.container():
            col_a, col_b = st.columns([3, 1])
            col_a.markdown(f"### 🧬 {sym}")
            if key not in st.session_state:
                if col_b.button("Cargar cobertura", key=btn_key):
                    with st.spinner(f"Cargando {sym}…"):
                        st.session_state[key] = api_coverage(sym, kit)
            if key in st.session_state:
                r = st.session_state[key]
                if r:
                    show_single(r, kit)
                else:
                    st.error(f"{sym}: error al cargar datos de cobertura.")
            st.divider()


# ── kit comparison view ────────────────────────────────────────────────────────

def show_kit_comparison(matrix_entries: list, genes: list, kits: list):
    """
    Shows comparison heatmap (coverage % or capture %) for genes × kits.
    """
    # Build DataFrames: pct_20 and pct_uncaptured
    pct_data    = {kit: {} for kit in kits}
    uncap_data  = {kit: {} for kit in kits}

    for entry in matrix_entries:
        g = entry["gene_name"]
        k = entry["kit_name"]
        if not entry.get("found") or not entry.get("report"):
            pct_data[k][g]   = None
            uncap_data[k][g] = None
            continue
        rep  = entry["report"]
        tot  = rep.get("n_exons_total", 1) or 1
        good = rep.get("n_exons_good", 0)
        unc  = rep.get("n_exons_uncaptured", 0)
        exons = rep.get("exons", [])
        if exons:
            avg20 = sum(e["pct_20x"] for e in exons if e["status"] != "uncaptured") / max(1, len([e for e in exons if e["status"] != "uncaptured"]))
        else:
            avg20 = 0.0
        pct_data[k][g]   = round(avg20, 1)
        uncap_data[k][g] = round(unc / tot * 100, 1)

    df_cov  = pd.DataFrame(pct_data, index=genes).fillna(-1)
    df_unc  = pd.DataFrame(uncap_data, index=genes).fillna(-1)

    tab1, tab2 = st.tabs(["Cobertura @20x", "Exones no capturados"])

    with tab1:
        fig = go.Figure(data=go.Heatmap(
            z=df_cov.values.tolist(),
            x=list(df_cov.columns),
            y=list(df_cov.index),
            colorscale=[
                [0, "#f44336"], [0.5, "#ff9800"],
                [0.8, "#4caf50"], [1, "#1a7c3a"],
            ],
            zmin=0, zmax=100,
            text=[[f"{v:.1f}%" if v >= 0 else "N/A" for v in row] for row in df_cov.values.tolist()],
            texttemplate="%{text}",
            hovertemplate="Gen: %{y}<br>Kit: %{x}<br>%@20x: %{z:.1f}%<extra></extra>",
        ))
        fig.update_layout(
            title="Cobertura media @20x por gen y kit (%)",
            height=max(300, 80 + 40 * len(genes)),
            xaxis_title="Kit comercial",
            yaxis_title="Gen",
        )
        st.plotly_chart(fig, use_container_width=True)

    with tab2:
        fig2 = go.Figure(data=go.Heatmap(
            z=df_unc.values.tolist(),
            x=list(df_unc.columns),
            y=list(df_unc.index),
            colorscale=[
                [0, "#1a7c3a"], [0.2, "#4caf50"],
                [0.5, "#ff9800"], [1, "#f44336"],
            ],
            zmin=0, zmax=100,
            text=[[f"{v:.0f}%" if v >= 0 else "N/A" for v in row] for row in df_unc.values.tolist()],
            texttemplate="%{text}",
            hovertemplate="Gen: %{y}<br>Kit: %{x}<br>% no capturado: %{z:.0f}%<extra></extra>",
        ))
        fig2.update_layout(
            title="% exones CDS NO capturados por kit",
            height=max(300, 80 + 40 * len(genes)),
            xaxis_title="Kit comercial",
            yaxis_title="Gen",
        )
        st.plotly_chart(fig2, use_container_width=True)

    # Detailed per-kit per-gene (expandable)
    for kit in kits:
        with st.expander(f"📦 {kit} — detalle por gen"):
            for entry in [e for e in matrix_entries if e["kit_name"] == kit]:
                g = entry["gene_name"]
                if not entry.get("found") or not entry.get("report"):
                    st.warning(f"{g}: sin datos para este kit.")
                    continue
                key = f"cmp_{g}_{kit}"
                st.session_state.setdefault(key, entry["report"])
                r = st.session_state[key]
                st.markdown(f"**{g}** — {r.get('recommendation','')}")
                exon_chart(r, f"{g} · {kit}")
                st.divider()


# ── NGS bottleneck tutorial ────────────────────────────────────────────────────

# Each bottleneck: (id, icon, título, subtítulo, impacto, descripción, ejemplo, impacto_color)
BOTTLENECKS = [
    {
        "id": "kit",
        "icon": "🪤",
        "title": "Librería de captura (kit comercial)",
        "subtitle": "EL cuello de botella más crítico",
        "impact": "5 – 25% de exones no capturados",
        "color": "#d32f2f",
        "badge": "CRÍTICO",
        "description": (
            "Todo exoma comercial define un conjunto de sondas de hibridación que determinan "
            "**qué regiones del genoma son enriquecidas** antes de la secuenciación. "
            "Cada fabricante diseña su kit de forma diferente: Twist Exome, Agilent SureSelect, "
            "Illumina TruSeq, IDT xGen y KAPA HyperExome **no cubren exactamente las mismas regiones**. "
            "Un exón que está en el kit A puede estar completamente fuera del kit B."
        ),
        "example": (
            "**Caso real:** El gen *TTN* (titina, 364 exones) tiene una cobertura muy variable entre kits. "
            "Algunos kits excluyen exones de expresión no cardíaca por decisión de diseño. "
            "En epilepsia, *KCNQ2* exón 2 está mal capturado en ciertos kits de primera generación."
        ),
        "action": "Sube el archivo BED de tu kit en el panel lateral y verifica gen por gen.",
    },
    {
        "id": "gc",
        "icon": "🧪",
        "title": "Sesgo de contenido GC",
        "subtitle": "Limitación física de la PCR e hibridación",
        "impact": "2 – 15% pérdida de cobertura uniforme",
        "color": "#f57c00",
        "badge": "ALTO",
        "description": (
            "La PCR y la hibridación de sondas son menos eficientes en regiones con "
            "contenido GC extremo (< 30% o > 70%). Estas regiones forman estructuras secundarias "
            "que dificultan la desnaturalización y la extensión de la polimerasa. "
            "El resultado es una **cobertura sistemáticamente baja** en ciertos exones, "
            "independientemente de la profundidad media del experimento."
        ),
        "example": (
            "**Genes con GC alto problemático:** *GC*, *KCNQ1* (65–70% GC en algunos exones), "
            "*FMR1* (expansión CGG en 5'UTR, no detectable por WES). "
            "**Genes con GC bajo:** algunos exones de *DMD* y *NF1*."
        ),
        "action": "Revisa el gráfico de cobertura por exón: busca exones con media < 20x que coincidan con GC extremo.",
    },
    {
        "id": "depth",
        "icon": "📊",
        "title": "Profundidad de cobertura desigual",
        "subtitle": "La media oculta zonas frías",
        "impact": "Cobertura media 100x ≠ todos los exones a 100x",
        "color": "#1565c0",
        "badge": "IMPORTANTE",
        "description": (
            "Aunque el informe de secuenciación indique '100x de cobertura media', "
            "esta estadística **enmascara una distribución muy heterogénea**. "
            "La cobertura sigue una distribución aproximadamente log-normal: "
            "algunos exones pueden estar a 300x mientras otros están a 5x. "
            "El estándar clínico es evaluar el **% de bases a ≥ 20x** (o ≥ 30x para mosaicismos), "
            "no la media."
        ),
        "example": (
            "**Regla práctica:** Un exón con < 80% de bases a ≥ 20x tiene riesgo real de "
            "falsos negativos. Para variantes en mosaico (epilepsia, cáncer somático) "
            "se recomienda ≥ 100x con ≥ 95% bases a ≥ 50x."
        ),
        "action": "Este advisor evalúa % de bases a ≥ 20x y ≥ 30x por exón usando datos de gnomAD r4.",
    },
    {
        "id": "repeats",
        "icon": "🔁",
        "title": "Regiones repetidas y pseudogenes",
        "subtitle": "Reads que mapean en múltiples loci",
        "impact": "Hasta 100% de fallos en regiones altamente homólogas",
        "color": "#6a1b9a",
        "badge": "ALTO",
        "description": (
            "Las reads de Illumina (150 pb) no pueden resolverse unívocamente "
            "cuando mapean a más de una región del genoma. Los **pseudogenes** "
            "(copias no funcionales con > 95% de homología) y las **repeticiones en tándem** "
            "hacen que el alineador asigne calidad de mapeo Q=0, descartando esas reads. "
            "El resultado es cobertura artificialmente cero aunque el exón esté en el kit."
        ),
        "example": (
            "**Genes problemáticos:** *SMN1/SMN2* (atrofia muscular espinal: "
            "WES estándar NO puede diagnosticar esta enfermedad), "
            "*PMS2* (pseudogén *PMS2CL* con 98% homología), "
            "*STRC* (cojera hereditaria), *CYP21A2* (hiperplasia adrenal congénita)."
        ),
        "action": "El advisore señala genes con limitaciones conocidas de mapeo (pLI + notas técnicas).",
    },
    {
        "id": "svs",
        "icon": "✂️",
        "title": "Variantes estructurales y expansiones de repetición",
        "subtitle": "El punto ciego del exoma de lectura corta",
        "impact": "WES estándar detecta < 30% de las SVs clínicamente relevantes",
        "color": "#2e7d32",
        "badge": "IMPORTANTE",
        "description": (
            "El WES con reads de 150 pb detecta bien SNVs e indels pequeños (< 20 pb), "
            "pero tiene **limitaciones severas** para: "
            "(1) Grandes deleciones/duplicaciones (CNVs > 50 pb) — detección parcial via CNV calling, "
            "(2) Inversiones y translocaciones — indetectables con WES estándar, "
            "(3) Expansiones de repetición trinucleotídica (FMR1, HTT, DMPK, C9orf72) — "
            "requieren PCR específica o long-read sequencing (Oxford Nanopore, PacBio)."
        ),
        "example": (
            "**Diagnósticos imposibles con WES solo:** Síndrome X frágil (*FMR1*), "
            "Huntington (*HTT*), ELA familiar (*C9orf72*), Distrofia miotónica (*DMPK*), "
            "Ataxia de Friedreich (*FXN*). Para estos genes WES es **complementario**, "
            "no diagnóstico."
        ),
        "action": "Este advisor indica para cada gen si WES es la tecnología adecuada o se necesita otra estrategia.",
    },
    {
        "id": "pipeline",
        "icon": "⚙️",
        "title": "Pipeline bioinformático",
        "subtitle": "Alineador + variant caller + anotación",
        "impact": "2 – 8% de variantes reales no detectadas según pipeline",
        "color": "#37474f",
        "badge": "MODERADO",
        "description": (
            "Incluso con cobertura perfecta, la cadena bioinformática introduce sus propios "
            "sesgos: el **alineador** (BWA-MEM2, DRAGMAP) puede colocar reads incorrectamente, "
            "el **variant caller** (GATK HaplotypeCaller, DeepVariant) tiene sensibilidades "
            "distintas para SNVs vs indels, y las herramientas de **anotación** "
            "(VEP, Annovar) tienen bases de datos que se actualizan continuamente. "
            "Un pipeline desactualizado puede perder variantes patogénicas conocidas."
        ),
        "example": (
            "**Ejemplo:** GATK HaplotypeCaller tiene sensibilidad ~99% para SNVs "
            "pero solo ~95% para indels de 10–50 pb. DeepVariant mejora indels pero "
            "requiere más recursos computacionales. Para mosaicismos, Mutect2 es necesario."
        ),
        "action": "Pregunta a tu bioinformático qué versión de pipeline usó y si está validado clínicamente.",
    },
]


def bottleneck_sankey():
    """Sankey diagram: flujo de regiones exónicas a través del pipeline NGS."""

    # ── Nodes ──
    # 0–5: pipeline principal (izq → der)
    # 6–9: nodos de pérdida (abajo)
    labels = [
        "Exoma teórico<br>(60–80 Mb · ~2% genoma)",    # 0
        "Capturado<br>por kit",                          # 1
        "Buena cobertura<br>GC",                         # 2
        "Cobertura<br>≥ 20x",                            # 3
        "Mapeo<br>único",                                 # 4
        "✅ Variantes<br>detectables",                   # 5
        "⚠ No capturado<br>por kit<br>(5–15%)",         # 6
        "⚠ Sesgo GC<br>extremo<br>(2–10%)",             # 7
        "⚠ Profundidad<br>insuficiente<br>(2–8%)",      # 8
        "⚠ Mapeo ambiguo<br>(pseudogenes)<br>(1–5%)",   # 9
    ]

    # Colors: pipeline = azules, pérdidas = rojos/naranjas
    node_colors = [
        "#1565c0",  # 0 input
        "#1976d2",  # 1
        "#1e88e5",  # 2
        "#42a5f5",  # 3
        "#64b5f6",  # 4
        "#1a7c3a",  # 5 output (verde)
        "#d32f2f",  # 6 loss kit
        "#f57c00",  # 7 loss GC
        "#1565c0",  # 8 loss depth
        "#6a1b9a",  # 9 loss repeats
    ]

    # ── Links (values represent % of original 100 target bases) ──
    # Pipeline: 0→1→2→3→4→5
    # Losses:   0→6, 1→7, 2→8, 3→9
    sources = [0, 0, 1, 1, 2, 2, 3, 3, 4]
    targets = [1, 6, 2, 7, 3, 8, 4, 9, 5]
    values  = [88, 12, 80, 8, 74, 6, 71, 3, 71]
    link_colors = [
        "rgba(25,118,210,0.4)",   # 0→1 ok
        "rgba(211,47,47,0.5)",    # 0→6 kit loss
        "rgba(30,136,229,0.4)",   # 1→2 ok
        "rgba(245,124,0,0.5)",    # 1→7 GC loss
        "rgba(66,165,245,0.4)",   # 2→3 ok
        "rgba(21,101,192,0.4)",   # 2→8 depth loss
        "rgba(100,181,246,0.4)",  # 3→4 ok
        "rgba(106,27,154,0.5)",   # 3→9 repeat loss
        "rgba(26,124,58,0.5)",    # 4→5 output
    ]

    hover_labels = [
        "88% de las regiones diana<br>son capturadas por el kit",
        "12% NO capturado por el kit<br>(varía entre kits)",
        "80% tienen buen perfil GC<br>y buena hibridación",
        "8% con sesgo GC extremo<br>(cobertura reducida)",
        "74% alcanza cobertura ≥20x<br>en todas sus bases",
        "6% con profundidad insuficiente<br>pese a estar capturado",
        "71% de regiones mapean<br>unívocamente",
        "3% con mapeo ambiguo<br>(pseudogenes, repeticiones)",
        "71% de bases exónicas<br>disponibles para llamado de variantes",
    ]

    fig = go.Figure(go.Sankey(
        arrangement="snap",
        node=dict(
            pad=25,
            thickness=28,
            line=dict(color="white", width=1),
            label=labels,
            color=node_colors,
            hovertemplate="<b>%{label}</b><extra></extra>",
        ),
        link=dict(
            source=sources,
            target=targets,
            value=values,
            color=link_colors,
            customdata=hover_labels,
            hovertemplate="%{customdata}<extra></extra>",
        ),
    ))

    fig.update_layout(
        title=dict(
            text="Algoritmo de pérdida de información en WES: del exoma teórico a variantes detectables",
            font=dict(size=15),
        ),
        font_size=12,
        height=440,
        paper_bgcolor="#f8f9fa",
        margin=dict(l=10, r=10, t=50, b=10),
    )
    st.plotly_chart(fig, use_container_width=True)

    st.caption(
        "💡 Los porcentajes son valores medios representativos (varía por gen, kit y laboratorio). "
        "El cuello de botella del **kit de captura** es el único controlable en el momento de la solicitud clínica."
    )


def bottleneck_funnel():
    """Embudo de detección: muestra la reducción progresiva de regiones informativas."""
    stages = [
        ("Exoma teórico", 100, "#1565c0"),
        ("Capturado por kit", 88, "#1976d2"),
        ("Sin sesgo GC severo", 80, "#1e88e5"),
        ("Cobertura ≥ 20x", 74, "#42a5f5"),
        ("Mapeo unívoco", 71, "#64b5f6"),
        ("Variantes detectables", 68, "#1a7c3a"),
    ]
    labels = [s[0] for s in stages]
    values = [s[1] for s in stages]
    colors = [s[2] for s in stages]

    fig = go.Figure(go.Funnel(
        y=labels,
        x=values,
        textposition="inside",
        textinfo="value+percent initial",
        texttemplate="%{value}% · %{percentInitial}",
        marker=dict(color=colors, line=dict(width=1, color="white")),
        connector=dict(line=dict(color="#cfd8dc", dash="dot", width=2)),
        hovertemplate="<b>%{y}</b><br>%{x}% del exoma teórico<extra></extra>",
    ))

    fig.update_layout(
        title="Embudo de detección: ¿cuánto del exoma es realmente informativo?",
        height=360,
        paper_bgcolor="#f8f9fa",
        plot_bgcolor="#f8f9fa",
        margin=dict(l=5, r=5, t=45, b=5),
        font=dict(size=12),
    )
    st.plotly_chart(fig, use_container_width=True)


def kit_impact_chart():
    """Gráfico de barras comparando cobertura de genes clave según tipo de kit."""
    genes   = ["SCN1A", "KCNQ2", "PCDH19", "SMN1*", "PMS2*", "FMR1*"]
    kit_a   = [96,       92,       88,       10,       45,       0]   # Kit well-designed
    kit_b   = [78,       65,       90,       8,        30,       0]   # Kit less coverage
    kit_c   = [99,       95,       91,       12,       60,       0]   # Premium kit

    fig = go.Figure()
    for name, vals, color in [
        ("Kit A (estándar)", kit_a, "#1976d2"),
        ("Kit B (básico)", kit_b, "#ff9800"),
        ("Kit C (premium)", kit_c, "#1a7c3a"),
    ]:
        fig.add_trace(go.Bar(
            name=name, x=genes, y=vals,
            marker_color=color,
            hovertemplate="<b>%{x}</b><br>%{fullData.name}: %{y}%@20x<extra></extra>",
        ))

    fig.add_hline(y=80, line_dash="dash", line_color="#d32f2f",
                  annotation_text="Umbral clínico 80%", annotation_position="right")
    fig.update_layout(
        title="Impacto del kit: % bases a ≥20x por gen (valores ilustrativos)",
        barmode="group",
        yaxis=dict(title="% bases ≥ 20x", range=[0, 105]),
        xaxis_title="Gen",
        height=360,
        paper_bgcolor="#f8f9fa",
        plot_bgcolor="#ffffff",
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
        margin=dict(l=5, r=5, t=70, b=5),
    )
    st.plotly_chart(fig, use_container_width=True)
    st.caption("*SMN1, PMS2, FMR1: genes con limitaciones estructurales (pseudogenes/repeticiones) — ningún kit los resuelve bien con WES estándar.")


def show_intro_tutorial():
    """Página 1: tutorial interactivo sobre los cuellos de botella del WES."""

    st.markdown("""
    ## 🔬 Cuellos de botella de la secuenciación de exoma (WES)

    El **Whole Exome Sequencing** no es una ventana perfecta al genoma. Entre el ADN del paciente
    y el resultado clínico final existe una cadena de pasos, cada uno con sus propias limitaciones.
    Entender estos cuellos de botella es esencial para interpretar correctamente un informe genético.
    """)

    # ── Algorithm overview ──
    st.subheader("El algoritmo de pérdida: del genoma al diagnóstico")

    col_info, col_metric = st.columns([3, 1])
    with col_info:
        st.markdown("""
        El siguiente diagrama muestra cómo el **exoma teórico** (el conjunto de todas las
        regiones codificantes del genoma) pierde información en cada paso del proceso WES.
        Cada nodo naranja/rojo representa un cuello de botella donde se pierden variantes reales.
        """)
    with col_metric:
        st.metric(
            label="Eficiencia diagnóstica media WES",
            value="68 – 85%",
            delta="del exoma teórico disponible",
            delta_color="off",
        )

    tab_sankey, tab_funnel, tab_kit_impact = st.tabs([
        "🌊 Flujo de información (Sankey)",
        "🔽 Embudo de detección",
        "📦 Impacto del kit por gen",
    ])

    with tab_sankey:
        bottleneck_sankey()
    with tab_funnel:
        bottleneck_funnel()
    with tab_kit_impact:
        kit_impact_chart()

    st.divider()

    # ── Bottleneck cards ──
    st.subheader("⚠️ Los 6 cuellos de botella del WES")
    st.markdown("Despliega cada uno para entender su mecanismo, impacto clínico y cómo detectarlo.")

    for bn in BOTTLENECKS:
        badge_color = {
            "CRÍTICO": "#d32f2f", "ALTO": "#f57c00",
            "IMPORTANTE": "#1565c0", "MODERADO": "#37474f",
        }.get(bn["badge"], "#666")

        header_md = (
            f"{bn['icon']} **{bn['title']}** &nbsp;&nbsp;"
            f"<span style='background:{badge_color};color:white;padding:2px 8px;"
            f"border-radius:4px;font-size:0.78em;'>{bn['badge']}</span>&nbsp;&nbsp;"
            f"<span style='color:#888;font-size:0.85em;'>{bn['impact']}</span>"
        )

        with st.expander(f"{bn['icon']} {bn['title']} — {bn['badge']} | {bn['impact']}"):
            st.markdown(f"### {bn['icon']} {bn['title']}")
            col_desc, col_eg = st.columns([3, 2])
            with col_desc:
                st.markdown(f"**¿Qué ocurre?**")
                st.markdown(bn["description"])
            with col_eg:
                st.info(bn["example"])
            st.success(f"🛠 **¿Cómo actuar?** {bn['action']}")

    st.divider()

    # ── Decision flowchart ──
    st.subheader("🗺️ ¿Es el WES la tecnología adecuada para tu caso?")

    col1, col2, col3 = st.columns(3)
    with col1:
        st.markdown("""
        **✅ WES es adecuado cuando:**
        - Fenotipo complejo con muchos genes candidatos
        - Variantes puntuales (SNVs, indels pequeños)
        - Genes con buena cobertura en el kit
        - Diagnóstico diferencial amplio
        - Análisis de novo (trío padres–hijo)
        """)
    with col2:
        st.markdown("""
        **⚠️ WES con precaución cuando:**
        - Gen con pseudogenes conocidos (*PMS2*, *STRC*)
        - Regiones GC-extremas (*KCNQ1*, *GJB2*)
        - Mosaicismos (necesitas ≥ 100x específico)
        - Exones grandes con cobertura heterogénea
        - Se sospecha variante intrónica profunda
        """)
    with col3:
        st.markdown("""
        **❌ WES no es suficiente cuando:**
        - *SMN1/SMN2* (atrofia muscular espinal)
        - *FMR1* (síndrome X frágil)
        - *HTT*, *DMPK*, *C9orf72* (expansiones)
        - Sospecha de reordenamiento cromosómico grande
        - Sospecha de mosaicismo < 10%
        """)

    st.divider()

    # ── Kit library ──
    st.subheader("📦 Biblioteca de kits — descarga directa")

    AZ_BASE = "https://raw.githubusercontent.com/AstraZeneca-NGS/reference_data/master/hg38/bed"

    KIT_LIBRARY = [
        # ── DESCARGA DIRECTA (AstraZeneca-NGS/reference_data, público en GitHub) ──
        {
            "name": "Agilent SureSelect V6 (hg38)",
            "vendor": "Agilent", "grch": "GRCh38",
            "regions": "244.655", "access": "público",
            "bed_url": f"{AZ_BASE}/Exome-Agilent_V6.bed",
            "source_url": "https://github.com/AstraZeneca-NGS/reference_data",
            "notes": "SureSelect Human All Exon V6 (2016). El más usado en clínica en Europa.",
        },
        {
            "name": "Agilent SureSelect V5 Clinical (hg38)",
            "vendor": "Agilent", "grch": "GRCh38",
            "regions": "231.176", "access": "público",
            "bed_url": f"{AZ_BASE}/Exome-Agilent_V5.bed",
            "source_url": "https://github.com/AstraZeneca-NGS/reference_data",
            "notes": "SureSelect Clinical Research Exome V5. Incluye regiones intrónicas clínicas.",
        },
        {
            "name": "Agilent SureSelect V4 (hg38)",
            "vendor": "Agilent", "grch": "GRCh38",
            "regions": "186.396", "access": "público",
            "bed_url": f"{AZ_BASE}/Exome-Agilent_V4.bed",
            "source_url": "https://github.com/AstraZeneca-NGS/reference_data",
            "notes": "Versión legacy, común en laboratorios con muestras históricas.",
        },
        {
            "name": "Agilent SureSelect V2 (hg38)",
            "vendor": "Agilent", "grch": "GRCh38",
            "regions": "172.243", "access": "público",
            "bed_url": f"{AZ_BASE}/Exome-Agilent_V2.bed",
            "source_url": "https://github.com/AstraZeneca-NGS/reference_data",
            "notes": "Versión más antigua. Útil para reanalizar datos históricos.",
        },
        {
            "name": "Agilent SureSelect OneSeq (hg38)",
            "vendor": "Agilent", "grch": "GRCh38",
            "regions": "91.930", "access": "público",
            "bed_url": f"{AZ_BASE}/Exome-Agilent-OneSeq.bed",
            "source_url": "https://github.com/AstraZeneca-NGS/reference_data",
            "notes": "Kit focado en genes de expresión cardíaca y panel clínico reducido.",
        },
        {
            "name": "IDT xGen Exome Research Panel V1 (hg38)",
            "vendor": "IDT", "grch": "GRCh38",
            "regions": "205.846", "access": "público",
            "bed_url": f"{AZ_BASE}/Exome-IDT_V1.bed",
            "source_url": "https://github.com/AstraZeneca-NGS/reference_data",
            "notes": "xGen V1 (2016). Para V2 (2019, ~39 Mb) descarga desde el portal IDT.",
        },
        {
            "name": "Roche NimbleGen MedExome (hg38)",
            "vendor": "Roche/NimbleGen", "grch": "GRCh38",
            "regions": "196.355", "access": "público",
            "bed_url": f"{AZ_BASE}/Exome-MedExome.bed",
            "source_url": "https://github.com/AstraZeneca-NGS/reference_data",
            "notes": "SeqCap EZ MedExome. Diseñado para diagnóstico clínico (~47 Mb incluyendo UTRs).",
        },
        {
            "name": "Roche NimbleGen SeqCap EZ Exome V3 (hg38)",
            "vendor": "Roche/NimbleGen", "grch": "GRCh38",
            "regions": "184.706", "access": "público",
            "bed_url": f"{AZ_BASE}/Exome-NGv3.bed",
            "source_url": "https://github.com/AstraZeneca-NGS/reference_data",
            "notes": "SeqCap EZ Exome V3 (legacy). Base para el posterior MedExome.",
        },
        {
            "name": "Gencode hg38 CDS (referencia teórica)",
            "vendor": "UCSC / Gencode", "grch": "GRCh38",
            "regions": "~37 Mb CDS", "access": "público",
            "bed_url": None,   # generated via API
            "source_url": "https://api.genome.ucsc.edu",
            "notes": "No es un kit real: son todas las regiones CDS del genoma. Referencia para detectar qué está fuera de cualquier kit.",
        },
        # ── REQUIEREN REGISTRO ──
        {
            "name": "Agilent SureSelect V7 / V8 (hg38)",
            "vendor": "Agilent", "grch": "GRCh38",
            "regions": "~48–50 Mb", "access": "registro",
            "bed_url": None,
            "source_url": "https://earray.chem.agilent.com/suredesign/",
            "notes": "V7 (2019) y V8 (2021) son las versiones más actuales. Requieren cuenta SureDesign.",
        },
        {
            "name": "IDT xGen Exome Research Panel V2 (hg38)",
            "vendor": "IDT", "grch": "GRCh38",
            "regions": "~39 Mb", "access": "registro",
            "bed_url": None,
            "source_url": "https://www.idtdna.com/pages/tools/bed-files",
            "notes": "V2 (2019): añade regiones mitocondriales y no codificantes clínicas respecto a V1.",
        },
        {
            "name": "Twist Human Core Exome / Comprehensive (hg38)",
            "vendor": "Twist", "grch": "GRCh38",
            "regions": "~33–90 Mb", "access": "registro",
            "bed_url": None,
            "source_url": "https://www.twistbioscience.com/account/login",
            "notes": "Core Exome 2.0 (~33 Mb) para uso clínico. Comprehensive (TCE) añade UTRs y no codificante.",
        },
        {
            "name": "Illumina TruSeq Exome / Nextera RC (hg38)",
            "vendor": "Illumina", "grch": "GRCh38",
            "regions": "~37–45 Mb", "access": "registro",
            "bed_url": None,
            "source_url": "https://support.illumina.com/downloads/truseq-exome-product-files.html",
            "notes": "TruSeq Exome (2020, ~45 Mb). Nextera Rapid Capture Exome (legacy, ~37 Mb).",
        },
        {
            "name": "KAPA HyperExome (hg38)",
            "vendor": "Roche", "grch": "GRCh38",
            "regions": "~44 Mb", "access": "registro",
            "bed_url": None,
            "source_url": "https://sequencing.roche.com/en/applications/library-preparation/dna-seq/targeted-sequencing/hyperexome.html",
            "notes": "Alta uniformidad de cobertura. Requiere cuenta Roche. Contactar con distribuidor.",
        },
    ]

    public_kits  = [k for k in KIT_LIBRARY if k["access"] == "público"]
    private_kits = [k for k in KIT_LIBRARY if k["access"] == "registro"]

    st.markdown(
        f"**{len(public_kits)} kits con descarga directa** (sin registro) procedentes del repositorio público "
        f"[AstraZeneca-NGS/reference_data](https://github.com/AstraZeneca-NGS/reference_data) "
        f"+ {len(private_kits)} kits que requieren cuenta en el portal del fabricante."
    )

    tab_pub, tab_priv = st.tabs([
        f"🔓 Descarga directa ({len(public_kits)})",
        f"🔐 Requieren registro ({len(private_kits)})",
    ])

    with tab_pub:
        st.caption(
            "Fuente: [AstraZeneca-NGS/reference_data](https://github.com/AstraZeneca-NGS/reference_data) "
            "(GitHub público). Kits en formato hg38 (chr1, chr2…). "
            "Kits históricos pero válidos para comparación de cobertura."
        )
        for kit_info in public_kits:
            col_name, col_info, col_btn = st.columns([3, 2, 2])
            with col_name:
                st.markdown(f"**{kit_info['name']}**")
                st.caption(kit_info["notes"])
            with col_info:
                if kit_info["regions"] != "~37 Mb CDS":
                    st.metric("Regiones", kit_info["regions"], label_visibility="collapsed")
                st.caption(f"🧫 {kit_info['vendor']} · {kit_info['grch']}")
            with col_btn:
                btn_key = f"install_{kit_info['name']}"
                if kit_info["bed_url"]:
                    if st.button("⬇️ Instalar", key=btn_key, use_container_width=True):
                        with st.spinner(f"Descargando {kit_info['name']}…"):
                            res = api_register_kit_url(kit_info["name"], kit_info["bed_url"])
                        if res:
                            st.success(f"✅ {res['regions']:,} regiones instaladas.")
                            st.rerun()
                        else:
                            st.error("Error al descargar. Comprueba conectividad.")
                else:
                    if st.button("⬇️ Generar desde UCSC", key=btn_key, use_container_width=True):
                        with st.spinner("Contactando UCSC…"):
                            res = api_generate_gencode()
                        if res:
                            st.success(f"✅ {res['regions']:,} regiones generadas.")
                            st.rerun()
                        else:
                            st.error("Error al contactar UCSC.")
            st.divider()

    with tab_priv:
        st.caption("Estos kits requieren cuenta en el portal del fabricante. Una vez descargado el BED, regístralo en el **panel lateral → Añadir kit**.")
        for kit_info in private_kits:
            with st.expander(f"📦 {kit_info['name']} — {kit_info['vendor']} · {kit_info['regions']}"):
                col_a, col_b = st.columns([3, 1])
                with col_a:
                    st.markdown(f"**Versión/cobertura:** {kit_info['regions']} · {kit_info['grch']}")
                    st.caption(kit_info["notes"])
                    st.markdown("**Pasos:**  \n1. Accede al portal del fabricante  \n2. Descarga el archivo `.bed` o `.bed.gz`  \n3. Súbelo en el panel lateral → *Añadir kit → Subir archivo BED*")
                with col_b:
                    if kit_info["source_url"]:
                        st.link_button(
                            f"🌐 Portal {kit_info['vendor']}",
                            kit_info["source_url"],
                            use_container_width=True,
                        )

    st.divider()

    # ── CTA ──
    st.info(
        "👉 **¿Tienes genes específicos?** Ve a la pestaña **🔬 Análisis** para evaluar la cobertura "
        "real exón a exón usando datos de gnomAD r4, y compara entre los kits comerciales de tu laboratorio."
    )

    with st.expander("📚 Referencias y fuentes de datos"):
        st.markdown("""
        - **gnomAD r4** (Genome Aggregation Database) — cobertura empírica exón a exón en 730.947 exomas
        - **MyGene.info** — información funcional, OMIM, resumen génico
        - **PanelApp GEL** (Genomics England) — paneles clínicos curados con evidencia
        - **[AstraZeneca-NGS/reference_data](https://github.com/AstraZeneca-NGS/reference_data)** — BED files públicos de Agilent V2–V6, IDT V1, NimbleGen MedExome/V3 (hg38)
        - Warr et al. (2015) *Exome sequencing: Current and future perspectives* — Genomics Proteomics Bioinformatics
        - Biesecker & Green (2014) *Diagnostic clinical genome and exome sequencing* — NEJM
        - Meynert et al. (2014) *Variant detection sensitivity and biases in whole genome and exome sequencing* — BMC Bioinformatics
        """)


# ── sidebar ────────────────────────────────────────────────────────────────────

def sidebar():
    st.sidebar.title("🧬 Exome Coverage Advisor")
    st.sidebar.caption("gnomAD r4 + MyGene.info + BED kits propios")

    # ── Kit management ──
    st.sidebar.subheader("📦 Kits comerciales (BED)")
    kits = api_list_kits()
    kit_names = [k["name"] for k in kits]

    with st.sidebar.expander("Añadir kit", expanded=not bool(kits)):
        reg_method = st.radio(
            "Método de registro",
            ["📁 Subir archivo BED", "🔗 Desde URL directa"],
            key="kit_reg_method",
            horizontal=True,
        )

        kit_name = st.text_input("Nombre del kit", key="new_kit_name",
                                 placeholder="p.ej. Twist Exome 2.0")

        if reg_method == "📁 Subir archivo BED":
            bed_file = st.file_uploader("Archivo BED del kit", type=["bed", "txt", "gz"],
                                        key="new_kit_bed")
            if st.button("Registrar kit", key="btn_register_kit"):
                if not kit_name:
                    st.error("Introduce un nombre para el kit.")
                elif not bed_file:
                    st.error("Sube el archivo BED del kit.")
                else:
                    content = bed_file.read().decode("utf-8", errors="replace")
                    result  = api_register_kit(kit_name, content)
                    if result:
                        st.success(f"✅ Kit '{kit_name}' registrado: {result['regions']:,} regiones.")
                        st.rerun()
                    else:
                        st.error("Error al registrar el kit.")

        else:  # URL
            bed_url = st.text_input(
                "URL directa del BED",
                key="new_kit_url",
                placeholder="https://ejemplo.com/kit.bed o .bed.gz",
            )
            st.caption("Acepta BED, BED.gz e interval_list de GATK/Picard.")
            if st.button("Descargar y registrar", key="btn_register_url"):
                if not kit_name:
                    st.error("Introduce un nombre para el kit.")
                elif not bed_url or not bed_url.startswith("http"):
                    st.error("Introduce una URL válida (http/https).")
                else:
                    with st.spinner("Descargando BED…"):
                        result = api_register_kit_url(kit_name, bed_url)
                    if result:
                        st.success(f"✅ Kit '{kit_name}' registrado: {result['regions']:,} regiones.")
                        st.rerun()
                    else:
                        st.error("No se pudo descargar o parsear el BED. Comprueba que la URL es accesible y el archivo es un BED válido.")

    if kits:
        st.sidebar.markdown("**Kits registrados:**")
        for k in kits:
            col_a, col_b = st.sidebar.columns([3, 1])
            col_a.markdown(f"- **{k['name']}** ({k['regions']:,} reg.)")
            if col_b.button("🗑", key=f"del_{k['name']}"):
                if api_delete_kit(k["name"]):
                    st.rerun()

    st.sidebar.divider()

    # ── PanelApp ──
    st.sidebar.subheader("🔬 PanelApp (genes por indicación)")
    panels = fetch_panels()
    if panels:
        panel_map  = {f"{p['name']} (v{p['version']} · {p['n_genes']} genes)": p["id"] for p in panels}
        panel_opts = list(panel_map.keys())
        sel_panel  = st.sidebar.selectbox("Seleccionar panel clínico", ["— ninguno —"] + panel_opts,
                                          key="sel_panel")
        if sel_panel != "— ninguno —":
            panel_id = panel_map[sel_panel]
            key_p = f"panel_genes_{panel_id}"
            if key_p not in st.session_state:
                with st.spinner("Cargando genes del panel…"):
                    st.session_state[key_p] = fetch_panel_genes(panel_id)
            genes_panel = st.session_state[key_p]
            st.sidebar.caption(f"Panel activo: {sel_panel}\n{len(genes_panel)} genes cargados")
            if genes_panel and st.sidebar.button("Cargar genes en el área de texto"):
                st.session_state["gene_input"] = "\n".join(genes_panel)
                st.rerun()
    else:
        st.sidebar.warning("PanelApp no disponible.")

    st.sidebar.divider()
    st.sidebar.caption("Fuentes: gnomAD r4 · MyGene.info · PanelApp GEL")

    return kit_names


# ── main ───────────────────────────────────────────────────────────────────────

def main():
    kit_names = sidebar()

    st.title("🧬 Exome Coverage Advisor")

    tab_intro, tab_analysis = st.tabs([
        "🏠 Introducción · Cuellos de botella NGS",
        "🔬 Análisis de cobertura",
    ])

    with tab_intro:
        show_intro_tutorial()

    with tab_analysis:
        st.markdown(
            "Evalúa la cobertura de genes clínicamente relevantes en exomas comerciales. "
            "Introduce genes manualmente o carga un panel PanelApp. "
            "Sube los BED de tus kits comerciales para ver qué exones quedan fuera de captura."
        )

        # ── Gene input ──
        col_genes, col_kits = st.columns([2, 1])

        with col_genes:
            st.subheader("Genes de interés")
            gene_text = st.text_area(
                "Genes (uno por línea o separados por comas)",
                value=st.session_state.get("gene_input", ""),
                height=180,
                placeholder="SCN1A\nKCNQ2\nPCDH19\nMECP2",
                key="gene_input_area",
            )
            genes = [g.strip().upper() for g in gene_text.replace(",", "\n").splitlines() if g.strip()]
            st.caption(f"{len(genes)} gen(es) introducido(s)")

        with col_kits:
            st.subheader("Kits a evaluar")
            if kit_names:
                selected_kits = st.multiselect(
                    "Seleccionar kits registrados",
                    options=kit_names,
                    default=kit_names[:1] if kit_names else [],
                    key="sel_kits",
                )
                st.caption("Sin selección → cobertura gnomAD sin filtro de captura")
            else:
                st.info("Sube un BED de kit en el panel lateral para comparar kits.")
                selected_kits = []

        if not genes:
            st.warning("Introduce al menos un gen para comenzar.")
            return

        run = st.button("🔍 Evaluar cobertura", type="primary", use_container_width=True)

        if run or "last_run" in st.session_state:
            if run:
                st.session_state["last_run"] = {"genes": genes, "kits": selected_kits}
                # Clear old cached reports
                for k in list(st.session_state.keys()):
                    if k.startswith("report_") or k.startswith("btn_") or k.startswith("cmp_"):
                        del st.session_state[k]

            run_genes = st.session_state.get("last_run", {}).get("genes", genes)
            run_kits  = st.session_state.get("last_run", {}).get("kits", selected_kits)

            n_genes = len(run_genes)
            n_kits  = len(run_kits)

            # ── Scenario A: single gene, single kit (or no kit) ──
            if n_genes == 1 and n_kits <= 1:
                gene = run_genes[0]
                kit  = run_kits[0] if run_kits else None
                key  = f"report_{gene}_{kit or 'none'}"
                if key not in st.session_state:
                    with st.spinner(f"Cargando cobertura de {gene}…"):
                        st.session_state[key] = api_coverage(gene, kit)
                r = st.session_state.get(key)
                if r:
                    st.subheader(f"{'📦 ' + run_kits[0] + ' · ' if run_kits else ''}🧬 {gene}")
                    show_single(r, kit)
                else:
                    st.error(f"Gen '{gene}' no encontrado en gnomAD / MyGene.info.")

            # ── Scenario B: single gene, multiple kits → per-exon kit comparison ──
            elif n_genes == 1 and n_kits > 1:
                gene = run_genes[0]
                st.subheader(f"🧬 {gene} — Comparativa de kits")
                key_mtx = f"matrix_{gene}_{','.join(run_kits)}"
                if key_mtx not in st.session_state:
                    with st.spinner(f"Calculando cobertura de {gene} para {n_kits} kits…"):
                        st.session_state[key_mtx] = api_coverage_matrix([gene], run_kits) or []
                entries = st.session_state.get(key_mtx, [])
                if entries:
                    show_kit_comparison(entries, [gene], run_kits)
                else:
                    st.error("Error al obtener la matriz de comparación.")

            # ── Scenario C: multiple genes, no kit (or one kit) ──
            elif n_genes > 1 and n_kits <= 1:
                kit = run_kits[0] if run_kits else None
                key_sum = f"summary_{','.join(run_genes)}_{kit or 'none'}"
                if key_sum not in st.session_state:
                    with st.spinner(f"Cargando resumen de {n_genes} genes…"):
                        st.session_state[key_sum] = api_batch_summary(run_genes) or []
                entries = st.session_state.get(key_sum, [])
                if entries:
                    st.subheader(f"Resumen — {n_genes} genes" + (f" · {kit}" if kit else ""))
                    show_multi_summary(entries, kit)
                else:
                    st.error("Error al obtener el resumen de genes.")

            # ── Scenario D: multiple genes, multiple kits → heatmap ──
            elif n_genes > 1 and n_kits > 1:
                st.subheader(f"🗺️ Mapa de cobertura — {n_genes} genes × {n_kits} kits")
                key_mtx = f"matrix_{','.join(run_genes)}_{','.join(run_kits)}"
                if key_mtx not in st.session_state:
                    with st.spinner(f"Calculando matriz {n_genes}×{n_kits}… (puede tardar varios minutos)"):
                        st.session_state[key_mtx] = api_coverage_matrix(run_genes, run_kits) or []
                entries = st.session_state.get(key_mtx, [])
                if entries:
                    show_kit_comparison(entries, run_genes, run_kits)
                else:
                    st.error("Error al obtener la matriz de comparación.")



if __name__ == "__main__":
    main()
