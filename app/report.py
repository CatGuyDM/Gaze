"""Gaze PDF report builder."""
import io
from datetime import datetime, timezone

from reportlab.lib.pagesizes import A4
from reportlab.lib.units import mm
from reportlab.lib import colors
from reportlab.lib.styles import ParagraphStyle
from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, HRFlowable,
)

VERSION = "5.4.0"

BLUE       = colors.HexColor("#1B4F8A")
BLUE_LIGHT = colors.HexColor("#EBF2FB")
RED        = colors.HexColor("#C0392B")
ORANGE     = colors.HexColor("#D35400")
GREEN      = colors.HexColor("#1E8449")
GREY       = colors.HexColor("#555555")
BORDER     = colors.HexColor("#CCCCCC")
ROW_ALT    = colors.HexColor("#F8F8F8")


def _styles():
    return {
        "title":   ParagraphStyle("title",   fontSize=16, textColor=BLUE, spaceAfter=4),
        "sub":     ParagraphStyle("sub",     fontSize=9,  textColor=GREY, spaceAfter=2),
        "section": ParagraphStyle("section", fontSize=10, textColor=BLUE,
                                  fontName="Helvetica-Bold", spaceBefore=8, spaceAfter=3),
        "body":    ParagraphStyle("body",    fontSize=8),
        "small":   ParagraphStyle("small",   fontSize=7,  textColor=GREY),
        "cell":    ParagraphStyle("cell",    fontSize=7,  wordWrap="LTR"),
        "cellb":   ParagraphStyle("cellb",   fontSize=7,  fontName="Helvetica-Bold"),
        "good":    ParagraphStyle("good",    fontSize=10, textColor=GREEN,
                                  fontName="Helvetica-Bold"),
        "bad":     ParagraphStyle("bad",     fontSize=10, textColor=RED,
                                  fontName="Helvetica-Bold"),
    }


def _risk_color(risk):
    return {"HIGH": RED, "MEDIUM": ORANGE, "LOW": GREEN}.get(risk, colors.black)


def _wrap(text, n=20):
    if not text or len(text) <= n:
        return text or ""
    return "\n".join(text[i:i+n] for i in range(0, len(text), n))


def build_pdf(result):
    buf    = io.BytesIO()
    styles = _styles()
    ts     = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")

    sequence  = result.get("sequence", "")
    matches   = result.get("matches", [])
    analyst   = result.get("analyst", "") or "—"
    fname     = result.get("fname", "") or "—"

    doc = SimpleDocTemplate(buf, pagesize=A4,
                            leftMargin=15*mm, rightMargin=15*mm,
                            topMargin=15*mm, bottomMargin=15*mm,
                            title=f"Gaze — {sequence[:24]}",
                            author=f"Gaze v{VERSION}")
    story = []


    story.append(Paragraph("Gaze SPPS impurity analysis", styles["title"]))
    story.append(Paragraph(
        f"Generated {ts}  ·  Analyst: {analyst}  ·  Version {VERSION}",
        styles["sub"]))
    story.append(HRFlowable(width="100%", thickness=1.2, color=BLUE, spaceAfter=6))


    story.append(Paragraph("Sample", styles["section"]))
    seq_disp = "\n".join(sequence[i:i+40] for i in range(0, len(sequence), 40)) if sequence else "—"
    info_rows = [
        ["Sequence",            seq_disp],
        ["Length",              f"{len(sequence)} residues"],
        ["C-terminal amide",    "yes" if result.get("c_amide") else "no"],
        ["Linker / fixed mod",  f"+{result.get('fixed_mod_da', 0.0):.4f} Da"],
        ["Predicted parent",    f"{result.get('parent_mono', 0):.4f} Da mono "
                                f"({result.get('parent_avg', 0):.2f} Da avg)"],
        ["Coupling reagent",    result.get("coupling", "—")],
        ["mzML file",           fname],
    ]
    info_t = Table(info_rows, colWidths=[42*mm, 130*mm])
    info_t.setStyle(TableStyle([
        ("FONTSIZE",     (0,0),(-1,-1), 8),
        ("FONTNAME",     (0,0),( 0,-1), "Helvetica-Bold"),
        ("TEXTCOLOR",    (0,0),( 0,-1), BLUE),
        ("VALIGN",       (0,0),(-1,-1), "TOP"),
        ("BOTTOMPADDING",(0,0),(-1,-1), 2),
        ("TOPPADDING",   (0,0),(-1,-1), 2),
    ]))
    story.append(info_t)
    story.append(Spacer(1, 4*mm))


    story.append(Paragraph("Apex scan", styles["section"]))
    n_coadded = result.get("n_scans_coadded", 1)
    coadd_label = (f"{n_coadded} scans co-added" if n_coadded > 1
                   else "single apex scan (no co-add)")
    apex_rows = [
        ["Scan number",        str(result.get("scan_number", "—"))],
        ["Retention time",     f"{result.get('rt_min', 0):.3f} min"],
        ["TIC at apex",        f"{result.get('tic_at_apex', 0):.3e}"],
        ["MS1 scans analysed", str(result.get("n_analytical_scans", "—"))],
        ["Integration",        coadd_label],
        ["Peaks detected",     f"{result.get('n_peaks', 0)} (≥0.3% of base)"],
        ["Base peak",          f"m/z {result.get('base_peak_mz', 0):.4f}, "
                               f"intensity {result.get('base_peak_int', 0):.0f}"],
        ["Mass tolerance",     f"{result.get('tolerance_ppm', 0)} ppm"],
    ]
    apex_t = Table(apex_rows, colWidths=[42*mm, 130*mm])
    apex_t.setStyle(TableStyle([
        ("FONTSIZE",     (0,0),(-1,-1), 8),
        ("FONTNAME",     (0,0),( 0,-1), "Helvetica-Bold"),
        ("TEXTCOLOR",    (0,0),( 0,-1), BLUE),
        ("VALIGN",       (0,0),(-1,-1), "TOP"),
        ("BOTTOMPADDING",(0,0),(-1,-1), 2),
        ("TOPPADDING",   (0,0),(-1,-1), 2),
    ]))
    story.append(apex_t)
    story.append(Spacer(1, 4*mm))


    story.append(Paragraph("Spectrum coverage", styles["section"]))
    parent_pct      = result.get("parent_envelope_pct", 0)
    matched_pct     = result.get("matched_impurity_pct", 0)
    unexplained_pct = result.get("unexplained_pct", 0)
    n_unexp_above1  = result.get("n_unexplained_above_1pct", 0)
    coverage_rows = [
        ["Parent envelope",       f"{parent_pct:.1f}%   (M0..M+10 across z=1..max)"],
        ["Matched impurities",    f"{matched_pct:.2f}%"],
        ["Unexplained intensity", f"{unexplained_pct:.1f}%"],
        ["Peaks >1% of base unexplained", str(n_unexp_above1)],
    ]
    cov_t = Table(coverage_rows, colWidths=[60*mm, 112*mm])
    cov_t.setStyle(TableStyle([
        ("FONTSIZE",     (0,0),(-1,-1), 8),
        ("FONTNAME",     (0,0),( 0,-1), "Helvetica-Bold"),
        ("TEXTCOLOR",    (0,0),( 0,-1), BLUE),
        ("VALIGN",       (0,0),(-1,-1), "TOP"),
        ("BOTTOMPADDING",(0,0),(-1,-1), 2),
        ("TOPPADDING",   (0,0),(-1,-1), 2),
    ]))
    story.append(cov_t)
    if unexplained_pct > 20:
        story.append(Spacer(1, 2*mm))
        story.append(Paragraph(
            f"<b>Note:</b> {unexplained_pct:.0f}% of detected intensity is not assigned "
            f"to the parent envelope or any predicted impurity. See "
            f"\"Top unexplained peaks\" below for the largest unaccounted-for signals.",
            styles["body"]))
    story.append(Spacer(1, 4*mm))


    if result.get("warnings"):
        story.append(Paragraph("Parser notes", styles["section"]))
        for w in result["warnings"]:
            story.append(Paragraph(f"• {w}", styles["small"]))
        story.append(Spacer(1, 3*mm))


    n_high = sum(1 for m in matches if m.get("risk") == "HIGH")
    big_high = sum(1 for m in matches
                   if m.get("risk") == "HIGH" and m.get("rel_int_pct", 0) >= 5.0)

    story.append(Paragraph("Headline", styles["section"]))
    if not matches:
        story.append(Paragraph(
            "No predicted impurities matched any peak above the intensity threshold. "
            "Either the synthesis is very clean, or the predicted parent mass does "
            "not correspond to the dominant peak in the spectrum (check sequence + linker).",
            styles["body"]))
    elif big_high:
        story.append(Paragraph(
            f"{big_high} HIGH-risk impurit{'y' if big_high == 1 else 'ies'} detected at "
            f"≥5% of base peak. Synthesis quality likely compromised — see table below.",
            styles["bad"]))
    elif n_high:
        story.append(Paragraph(
            f"{n_high} HIGH-risk impurit{'y' if n_high == 1 else 'ies'} detected at "
            f"trace level (<5% of base peak). Spectrum looks clean overall.",
            styles["good"]))
    else:
        story.append(Paragraph(
            "No HIGH-risk impurities detected above the intensity threshold. "
            "Spectrum looks clean.", styles["good"]))
    n_iso_warn = sum(1 for m in matches if m.get("iso_warning"))
    if n_iso_warn > 0:
        story.append(Paragraph(
            f"{n_iso_warn} impurit{'y' if n_iso_warn == 1 else 'ies'} reported at "
            f"≥5% intensity without isotope-envelope confirmation. "
            f"{'This has' if n_iso_warn == 1 else 'These have'} been downgraded to "
            f"MEDIUM risk pending manual verification.",
            styles["body"]))
    story.append(Spacer(1, 3*mm))


    n_med = sum(1 for m in matches if m.get("risk") == "MEDIUM")
    n_low = sum(1 for m in matches if m.get("risk") == "LOW")
    counts = [
        ["Total matched", "HIGH risk", "MEDIUM risk", "LOW risk", "Predicted total"],
        [str(len(matches)), str(n_high), str(n_med), str(n_low),
         str(result.get("n_predicted", "—"))],
    ]
    ct = Table(counts, colWidths=[34*mm]*5)
    ct.setStyle(TableStyle([
        ("BACKGROUND", (0,0),(-1,0), BLUE),
        ("TEXTCOLOR",  (0,0),(-1,0), colors.white),
        ("FONTSIZE",   (0,0),(-1,-1), 8),
        ("FONTNAME",   (0,1),(-1,1), "Helvetica-Bold"),
        ("ALIGN",      (0,0),(-1,-1), "CENTER"),
        ("VALIGN",     (0,0),(-1,-1), "MIDDLE"),
        ("BACKGROUND", (0,1),(-1,1), BLUE_LIGHT),
        ("BOX",        (0,0),(-1,-1), 0.5, BORDER),
        ("INNERGRID",  (0,0),(-1,-1), 0.25, BORDER),
        ("TOPPADDING", (0,0),(-1,-1), 3),
        ("BOTTOMPADDING",(0,0),(-1,-1), 3),
    ]))
    story.append(ct)
    story.append(Spacer(1, 4*mm))


    story.append(Paragraph("Matched impurities", styles["section"]))
    if not matches:
        story.append(Paragraph("No matches above threshold.", styles["body"]))
    else:
        headers = ["Impurity", "Pos", "z", "Obs m/z", "ppm", "% base", "Env", "Iso",
                   "⚠Iso", "Pers", "Conf", "Risk", "Likely cause", "Recommended action"]
        col_w   = [22*mm, 14*mm, 5*mm, 16*mm, 8*mm, 10*mm, 7*mm, 6*mm,
                   5*mm, 8*mm, 8*mm, 9*mm, 30*mm, 30*mm]
        tdata   = [headers]
        rstyles = []


        msorted = sorted(matches, key=lambda m: (-m.get("confidence", 0),
                                                  -m.get("rel_int_pct", 0)))

        for i, m in enumerate(msorted, 1):
            risk = m.get("risk", "LOW")
            iso  = "✓" if m.get("iso_confirmed") else "—"
            why  = (m.get("why") or "")[:160]
            fix  = (m.get("fix") or "")[:160]
            env  = m.get("envelope_score", 0.0)
            conf = m.get("confidence", 0.0)
            pos  = m.get("position", "") or "—"
            pers = m.get("scan_persistence", 1.0)
            iso_warn = "⚠" if m.get("iso_warning") else ""
            tdata.append([
                Paragraph(_wrap(m.get("impurity", ""), 14), styles["cell"]),
                Paragraph(_wrap(str(pos), 8), styles["cell"]),
                str(m.get("z", "")),
                f"{m.get('obs_mz', 0):.4f}",
                f"{m.get('ppm', 0):.2f}",
                f"{m.get('rel_int_pct', 0):.2f}",
                f"{env:.2f}",
                iso,
                iso_warn,
                f"{pers:.2f}",
                f"{conf:.1f}",
                Paragraph(risk if risk != "MEDIUM" else "MED", styles["cellb"]),
                Paragraph(why, styles["cell"]),
                Paragraph(fix, styles["cell"]),
            ])
            rstyles.append(("TEXTCOLOR", (11, i), (11, i), _risk_color(risk)))
            rstyles.append(("FONTNAME",  (11, i), (11, i), "Helvetica-Bold"))
            if i % 2 == 0:
                rstyles.append(("BACKGROUND", (0, i), (-1, i), ROW_ALT))

        rt = Table(tdata, colWidths=col_w, repeatRows=1)
        rt.setStyle(TableStyle([
            ("BACKGROUND", (0,0),(-1,0), BLUE),
            ("TEXTCOLOR",  (0,0),(-1,0), colors.white),
            ("FONTSIZE",   (0,0),(-1,-1), 7),
            ("FONTNAME",   (0,0),(-1,0), "Helvetica-Bold"),
            ("ALIGN",      (0,0),(-1,-1), "LEFT"),
            ("ALIGN",      (2,0),(10,-1), "CENTER"),
            ("VALIGN",     (0,0),(-1,-1), "TOP"),
            ("TOPPADDING", (0,0),(-1,-1), 2),
            ("BOTTOMPADDING",(0,0),(-1,-1), 2),
            ("BOX",        (0,0),(-1,-1), 0.5, BORDER),
            ("INNERGRID",  (0,0),(-1,-1), 0.25, BORDER),
        ] + rstyles))
        story.append(rt)


    top_unexp = result.get("top_unexplained", [])
    if top_unexp:
        story.append(Spacer(1, 4*mm))
        story.append(Paragraph("Top unexplained peaks", styles["section"]))
        story.append(Paragraph(
            "Peaks not attributed to the parent envelope and not matched to any "
            "predicted impurity. Inspect manually for impurities outside Gaze's "
            "catalogue, parent multimers, or co-eluting species.",
            styles["small"]))
        story.append(Spacer(1, 1.5*mm))
        unexp_data = [["m/z", "Intensity", "% base"]]
        for mz, intensity, pct in top_unexp[:10]:
            unexp_data.append([f"{mz:.4f}", f"{intensity:.0f}", f"{pct:.2f}"])
        ut = Table(unexp_data, colWidths=[35*mm, 35*mm, 25*mm])
        ut.setStyle(TableStyle([
            ("BACKGROUND", (0,0),(-1,0), BLUE),
            ("TEXTCOLOR",  (0,0),(-1,0), colors.white),
            ("FONTSIZE",   (0,0),(-1,-1), 7),
            ("FONTNAME",   (0,0),(-1,0), "Helvetica-Bold"),
            ("ALIGN",      (0,0),(-1,-1), "LEFT"),
            ("VALIGN",     (0,0),(-1,-1), "TOP"),
            ("TOPPADDING", (0,0),(-1,-1), 2),
            ("BOTTOMPADDING",(0,0),(-1,-1), 2),
            ("BOX",        (0,0),(-1,-1), 0.5, BORDER),
            ("INNERGRID",  (0,0),(-1,-1), 0.25, BORDER),
        ]))
        story.append(ut)

    story.append(Spacer(1, 5*mm))
    story.append(HRFlowable(width="100%", thickness=0.4, color=BORDER))
    story.append(Spacer(1, 2*mm))
    story.append(Paragraph(
        "<b>Limitations.</b> Gaze matches predicted impurities to observed peaks "
        "by mass alone (with isotope envelope fitting). It does not analyse MS/MS "
        "fragmentation, retention time, UV chromatograms, or chromatographic peak shape. "
        "Mass-only matching cannot distinguish position-degenerate impurities (e.g. "
        "different deletion sites with the same composition), and cannot detect "
        "impurities outside the predefined catalogue. Predicted impurities require "
        "orthogonal confirmation (MS/MS for structure, chiral HPLC for racemization, "
        "retention time for isobars) before any analytical conclusion. Gaze is research "
        "software, not a regulatory tool.",
        styles["small"]))
    story.append(Spacer(1, 1*mm))
    story.append(Paragraph(
        f"Gaze v{VERSION}  ·  {ts}", styles["small"]))

    doc.build(story)
    return buf.getvalue()
