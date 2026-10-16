import streamlit as st
import json
import math
import io
import matplotlib.pyplot as plt
import numpy as np
from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill, Alignment, Border, Side
from reportlab.lib.pagesizes import letter
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib import colors

# --- CONFIGURACIÓN DE ESTÉTICA (VIOLETA INDUSTRIAL) ---
st.markdown("""
    <style>
    .stApp {
        background-color: #120f1d;
        color: #e6e3f8;
    }
    h1, h2, h3 {
        color: #e2d9ff !important;
    }
    .stButton>button {
        background: linear-gradient(135deg, #ff4b4b, #d93838);
        color: white;
        border: none;
        border-radius: 8px;
        font-weight: bold;
        padding: 12px 24px;
        width: 100%;
        box-shadow: 0 4px 10px rgba(255, 75, 75, 0.4);
        transition: all 0.3s ease;
    }
    .stButton>button:hover {
        background: linear-gradient(135deg, #ff6b6b, #e62e2e);
        box-shadow: 0 6px 15px rgba(255, 107, 107, 0.6);
    }
    [data-testid="stMetric"] {
        background-color: #1e1730;
        border: 1px solid #332752;
        padding: 15px;
        border-radius: 8px;
        box-shadow: 0 4px 6px rgba(0,0,0,0.2);
    }
    [data-testid="stMetricLabel"] {
        color: #bfaee3 !important;
    }
    [data-testid="stMetricValue"] {
        color: #f3efff !important;
    }
    .stDownloadButton>button {
        background: linear-gradient(135deg, #382159, #231338);
        color: #e2d9ff;
        border: 1px solid #7d4bc4;
        border-radius: 8px;
        font-weight: bold;
        padding: 12px 24px;
        width: 100%;
    }
    .stDownloadButton>button:hover {
        background: linear-gradient(135deg, #4b2d78, #321a52);
        color: white;
    }
    </style>
""", unsafe_allow_html=True)


# --- BASE DE DATOS DE COMPUESTOS PREDEFINIDOS ---
FLUIDOS_DB = {
    "Vapor de agua (Saturado)": {
        "wh": 720.0, "p": 300.0, "lam": 2163.5, "tsat": 133.5,
        "rho_v": 1.651, "mu_v": 0.0134, "k_v": 0.0275, "cp_v": 2.17,
        "rho_l": 943.0, "mu_l": 0.232, "k_l": 0.685, "cp_l": 4.25
    },
    "Amoníaco (R-717)": {
        "wh": 650.0, "p": 1167.0, "lam": 1155.0, "tsat": 30.0,
        "rho_v": 6.82, "mu_v": 0.0110, "k_v": 0.0220, "cp_v": 2.80,
        "rho_l": 608.0, "mu_l": 0.145, "k_l": 0.480, "cp_l": 4.60
    },
    "Etanol": {
        "wh": 800.0, "p": 100.0, "lam": 840.0, "tsat": 78.3,
        "rho_v": 1.59, "mu_v": 0.0105, "k_v": 0.0180, "cp_v": 1.50,
        "rho_l": 757.0, "mu_l": 0.410, "k_l": 0.168, "cp_l": 2.84
    }
}


# --- MOTOR DE CÁLCULO TÉRMICO E HIDRÁULICO TEMA ---
G = 9.81
L_STD = [1.83, 2.44, 3.05, 3.66, 4.88, 6.10]
DS_STD = [x * 0.0254 for x in [8, 10, 12, 13.25, 15.25, 17.25, 19.25, 21.25, 23.25, 25, 27, 29, 31, 33, 35, 37, 39, 42, 45, 48, 54, 60]]

K_TRI = {1: (0.319, 2.142), 2: (0.249, 2.207), 4: (0.175, 2.285), 6: (0.0743, 2.499), 8: (0.0365, 2.675)}
K_CUA = {1: (0.215, 2.207), 2: (0.156, 2.291), 4: (0.158, 2.263), 6: (0.0402, 2.617), 8: (0.0331, 2.643)}
HOLGURA = {"M": 0.012, "N": 0.012, "L": 0.012, "U": 0.015, "S": 0.060, "T": 0.095, "P": 0.050, "W": 0.035}

def dmlt_contracorriente(T1, T2, t1, t2):
    d1, d2 = T1 - t2, T2 - t1
    if d1 <= 0 or d2 <= 0:
        return abs(T1 - t2) + 1.0
    if abs(d1 - d2) < 1e-9:
        return d1
    return (d1 - d2) / math.log(d1 / d2)

def factor_ft(T1, T2, t1, t2, Ns):
    R = (T1 - T2) / (t2 - t1) if abs(t2 - t1) > 1e-6 else 1e-6
    P = (t2 - t1) / (T1 - t1) if abs(T1 - t1) > 1e-6 else 1e-6
    if R < 1e-6: return 1.0
    if abs(R - 1) < 1e-6: R = 1 + 1e-6
    if P >= 1 or 1 - P * R <= 0: return 0.75
    F = ((1 - P * R) / (1 - P)) ** (1.0 / Ns)
    P1 = (1 - F) / (R - F)
    S = math.sqrt(R * R + 1)
    a = (1 - P1) / (1 - R * P1)
    b = (2 - P1 * (R + 1 - S)) / (2 - P1 * (R + 1 + S))
    if a <= 0 or b <= 0: return 0.75
    return S * math.log(a) / ((R - 1) * math.log(b))

def nusselt_tubo(Re, Pr):
    def gnielinski(Re):
        f = (0.79 * math.log(Re) - 1.64) ** -2
        return (f / 8) * (Re - 1000) * Pr / (1 + 12.7 * math.sqrt(f / 8) * (Pr ** (2 / 3) - 1))
    if Re >= 4000: return gnielinski(Re)
    if Re <= 2300: return 3.66
    w = (Re - 2300) / 1700.0
    return 3.66 * (1 - w) + gnielinski(4000) * w

def factor_friccion(Re, Di, rug=4.6e-5):
    if Re < 2300: return 64.0 / Re
    inv = -1.8 * math.log10((rug / Di / 3.7) ** 1.11 + 6.9 / Re)
    return 1.0 / inv ** 2

def kern_casco(W, rho, mu, k, cp, Ds, B, Pt, Do, L, arreglo):
    C = Pt - Do
    As = Ds * C * B / Pt
    Gs = W / As if As > 0 else 0
    De = 4 * (Pt ** 2 - math.pi * Do ** 2 / 4) / (math.pi * Do) if arreglo == "cuadro" else 4 * (0.5 * Pt * 0.866 * Pt - 0.5 * math.pi * Do ** 2 / 4) / (math.pi * Do / 2)
    Re = De * Gs / mu if mu > 0 else 0
    Pr = cp * mu / k if k > 0 else 1.0
    hs = 0.36 * k / De * Re ** 0.55 * Pr ** (1 / 3) if k > 0 else 0.0
    f = math.exp(0.576 - 0.19 * math.log(max(Re, 10)))
    Nb = max(1, int(round(L / B - 1)))
    dP = f * Gs ** 2 * Ds * (Nb + 1) / (2 * rho * De) if rho > 0 else 0
    return dict(hs=hs, dP=dP, Re=Re, Gs=Gs, De=De, Nb=Nb, As=As, v=Gs / rho if rho > 0 else 0)

def h_cond_horizontal(k, rho, rhov, mu, lam, dT, Do, Nr):
    return 0.725 * (k ** 3 * rho * (rho - rhov) * G * lam / (mu * max(dT, 0.1) * Do)) ** 0.25 * max(Nr, 1) ** (-1 / 6)

def capacidad_tubos(Ds, n, Do, Pt, arreglo, rear):
    K1, n1 = (K_CUA if arreglo == "cuadro" else K_TRI)[n]
    Db = Ds - HOLGURA.get(rear, 0.012)
    if Db <= 0: return 0
    return int(K1 * (Db / Do) ** n1 * (1.25 * Do / Pt) ** 2)

def preparar(d, warns):
    h, c = d["caliente"], d["fria"]
    h["Tsat"] = h.get("T_sat") or h["T1"]
    Q = (h["Wh"] / 3600.0) * h["lambda"]
    Wc_bal = Q / (c["cp"] * (c["t2"] - c["t1"])) * 3600.0
    c["Wc"] = Wc_bal
    d["Q"] = Q
    geo = d["geometria"]
    d["Rd_req"] = h["Rd"] + c["Rd"] * geo["Do"] / geo["Di"]
    d["rear"] = d["identificacion"]["tema"].strip().upper()[-1]

def evaluar(d, g):
    h, c, geo = d["caliente"], d["fria"], d["geometria"]
    Do, Di, Pt, kw = geo["Do"], geo["Di"], geo["Pt"], geo["k_pared"]
    arr = geo["arreglo"]
    Ns, Nt, n, L, Ds, B = g["Ns"], g["Nt"], g["n"], g["L"], g["Ds"], g["B"]
    Tsat = h["Tsat"]
    Tc = 0.5 * (c["t1"] + c["t2"])
    Rw = Do * math.log(Do / Di) / (2 * kw)
    A = Ns * Nt * math.pi * Do * L
    r = dict(A=A, Ns=Ns, Nt=Nt, n=n, L=L, Ds=Ds, B=B)

    Ap = Nt / n * math.pi * Di ** 2 / 4
    Gt = (c["Wc"] / 3600.0) / Ap if Ap > 0 else 0
    Re_t = Gt * Di / c["mu"]
    Pr_t = c["cp"] * c["mu"] / c["k"]
    hi = nusselt_tubo(Re_t, Pr_t) * c["k"] / Di
    f_t = factor_friccion(Re_t, Di)
    dPt = Ns * f_t * (L * n / Di) * Gt ** 2 / (2 * c["rho"])
    r.update(hi=hi, v_t=Gt / c["rho"], Re_t=Re_t, Pr_t=Pr_t, dPt=dPt)

    Nr = max(1.0, 2.0 / 3.0 * Ds / Pt)
    dT, ho = 5.0, 10000.0
    for _ in range(50):
        ho = h_cond_horizontal(h["k_liquido"], h["rho_liquido"], h["rho_vapor"], h["mu_liquido"], h["lambda"], dT, Do, Nr)
        Uc = 1.0 / (1 / ho + Rw + Do / (Di * hi))
        dT_new = max(Uc * (Tsat - Tc) / max(ho, 1.0), 0.5)
        if abs(dT_new - dT) < 1e-4: break
        dT = 0.5 * dT + 0.5 * dT_new

    k_casco = kern_casco(h["Wh"] / 3600.0, h["rho_vapor"], h["mu_vapor"], h["k_vapor"], h["cp_vapor"], Ds, B, Pt, Do, L, arr)
    dPs = 0.5 * Ns * k_casco["dP"]
    r.update(ho=ho, dT_film=dT, dPs=dPs, v_s=k_casco["v"], Re_s=k_casco["Re"], Nb=k_casco["Nb"], De=k_casco["De"], As=k_casco["As"], Gs=k_casco["Gs"], Nr=Nr)

    Ud = d["Q"] * 1000.0 / (A * d["DMLT"] * d["Ft"]) if (A * d["DMLT"] * d["Ft"]) > 0 else 0
    Rd_calc = (Uc - Ud) / (Uc * Ud) if (Uc * Ud) > 0 else 0
    A_req_final = (d["Q"] * 1000.0) / ((1.0 / (1.0 / Uc + d["Rd_req"])) * d["DMLT"] * d["Ft"]) if (d["DMLT"] * d["Ft"]) > 0 else A

    r.update(Uc=Uc, Ud=Ud, Rd_calc=Rd_calc, A_req_final=A_req_final)
    r.update(dp_t_adm=c["dp_adm"], dp_s_adm=h["dp_adm"])
    r.update(ok_dp=(r["dPt"] <= c["dp_adm"] and r["dPs"] <= h["dp_adm"]), ok_Rd=(Rd_calc >= d["Rd_req"]))
    return r

def buscar_geometria(d, A_req, Ns, log):
    geo = d["geometria"]
    Do, Pt, arr, rear = geo["Do"], geo["Pt"], geo["arreglo"], d["rear"]
    n_opts = [2, 4, 6, 8]
    L_opts = sorted({l for l in L_STD if l <= geo["L_max"]} | {geo["L_max"]}, reverse=True)
    vmin, vmax = geo["v_min"], geo["v_max"]
    mejor_fallback = None

    for L in L_opts:
        Nt_tot = math.ceil(A_req / (math.pi * Do * L))
        Nt_s = math.ceil(Nt_tot / Ns)
        for n in n_opts:
            Nt = int(math.ceil(Nt_s / n) * n)
            idx0 = next((i for i, Ds in enumerate(DS_STD) if capacidad_tubos(Ds, n, Do, Pt, arr, rear) >= Nt), None)
            if idx0 is None: continue
            for idx in range(idx0, min(idx0 + 3, len(DS_STD))):
                Ds = DS_STD[idx]
                B = max(min(0.5 * Ds, 1.52), 0.05)
                g = dict(Ns=Ns, Nt=Nt, n=n, L=L, Ds=Ds, B=B)
                r = evaluar(d, g)
                v_ok = vmin <= r["v_t"] <= vmax
                if r["ok_dp"] and v_ok and r["ok_Rd"]: return r
                if r["ok_dp"] and mejor_fallback is None: mejor_fallback = r
    return mejor_fallback

def disenar(d):
    log, warns = [], []
    preparar(d, warns)
    h, c = d["caliente"], d["fria"]
    T1, T2, t1, t2 = h["T1"], h["T2"], c["t1"], c["t2"]
    d["DMLT"] = dmlt_contracorriente(T1, T2, t1, t2)
    Ns = 1
    while True:
        Ft = factor_ft(T1, T2, t1, t2, Ns)
        if Ft >= d["geometria"]["Ft_min"] or Ns >= 4: break
        Ns += 1
    d["Ft"], d["Ns"] = Ft, Ns
    
    Ud_est = 1500.0
    A_prov = (d["Q"] * 1000.0) / (Ud_est * d["DMLT"] * Ft)
    log.append(f"Balance térmico: Q = {d['Q']:.1f} kW, DMLT = {d['DMLT']:.2f} K, Ft = {Ft:.3f}")
    log.append(f"Caudal de agua calculado (Wc): {c['Wc']:.2f} kg/h")
    
    r = buscar_geometria(d, A_prov, Ns, log)
    if r is None:
        r = evaluar(d, dict(Ns=Ns, Nt=20, n=2, L=3.66, Ds=DS_STD[3], B=0.15))
        warns.append("No se encontró geometría óptima estricta; se muestra configuración aproximada.")
    log.append(f"Selección final: Ds = {r['Ds']*1000:.0f} mm, n = {r['n']}, Nt = {r['Nt']}, L = {r['L']:.2f} m")
    return r, log, warns


# --- EXPORTADOR A HOJA TEMA CLÁSICA COMPLETA EN EXCEL ---
def exportar_hoja_tema_bytes(d, r, log, warns):
    wb = Workbook()
    ws = wb.active
    ws.title = "Hoja TEMA"
    F = "Arial"
    thin = Side(style="thin", color="999999")
    brd = Border(left=thin, right=thin, top=thin, bottom=thin)
    fill_h = PatternFill("solid", fgColor="321c54")
    fill_s = PatternFill("solid", fgColor="4a2b75")
    ws.column_dimensions["A"].width = 38
    for col in "BCD": ws.column_dimensions[col].width = 24
    ws.column_dimensions["E"].width = 14

    row = [1]
    ref = {}

    def titulo(txt):
        ws.cell(row[0], 1, txt)
        ws.merge_cells(start_row=row[0], start_column=1, end_row=row[0], end_column=5)
        c_ = ws.cell(row[0], 1)
        c_.font = Font(name=F, bold=True, color="FFFFFF", size=11)
        c_.fill = fill_h
        row[0] += 1

    def encabezado(*vals):
        for i, v in enumerate(vals):
            c_ = ws.cell(row[0], 1 + i, v)
            c_.font = Font(name=F, bold=True, color="FFFFFF", size=10)
            c_.fill = fill_s
            c_.border = brd
            c_.alignment = Alignment(horizontal="center")
        row[0] += 1

    def fila(*vals, key=None, fmt=None):
        for i, v in enumerate(vals):
            c_ = ws.cell(row[0], 1 + i, v)
            c_.font = Font(name=F, size=10)
            c_.border = brd
            if i > 0 and fmt and isinstance(v, (int, float, str)) and (not isinstance(v, str) or v.startswith("=")):
                c_.number_format = fmt
            if i > 0: c_.alignment = Alignment(horizontal="right")
        if key: ref[key] = row[0]
        row[0] += 1

    idn, h, c, geo = d["identificacion"], d["caliente"], d["fria"], d["geometria"]

    ws.cell(row[0], 1, "HOJA DE ESPECIFICACIONES TÉCNICAS — INTERCAMBIADOR DE CASCO Y TUBOS (TEMA)")
    ws.merge_cells(start_row=row[0], start_column=1, end_row=row[0], end_column=5)
    ws.cell(row[0], 1).font = Font(name=F, bold=True, size=13)
    row[0] += 2

    titulo("1. DATOS DE IDENTIFICACIÓN DEL EQUIPO")
    fila("Ítem / Tag", idn["item"])
    fila("Servicio", idn["servicio"])
    fila("Cliente / Cátedra", idn["cliente"])
    fila("Tipo TEMA", idn["tema"])
    fila("Cascos en serie", r["Ns"])
    row[0] += 1

    titulo("2. CONDICIONES DE OPERACIÓN Y PROCESO")
    encabezado("Parámetro", "Lado casco", "Lado tubos", "Unidad")
    fila("Fluido", fluido_v, "Agua de enfriamiento", "-")
    fila("Caudal másico", round(h["Wh"], 2), round(c["Wc"], 2), "kg/h", fmt="#,##0.00")
    fila("Temperatura de entrada", h["T1"], c["t1"], "°C", fmt="0.0")
    fila("Temperatura de salida", h["T2"], c["t2"], "°C", fmt="0.0")
    fila("Presión de operación", round(p_oper_s, 2), "-", "kPa abs", fmt="0.00")
    fila("Calor latente de condensación", round(lam, 2), "-", "kJ/kg", fmt="0.00")
    fila("Densidad (líquido / vapor)", f"{rho_l} / {rho_v}", "-", "kg/m³")
    fila("Viscosidad", f"{mu_l}", "-", "mPa·s")
    fila("Conductividad térmica", k_l, "-", "W/m·K")
    fila("Factor de ensuciamiento", h["Rd"], c["Rd"], "m²·K/W", fmt="0.000000")
    fila("Δp admisible", round(dp_adm_s, 1), round(dp_adm_t, 1), "kPa", fmt="0.0")
    fila("Δp calculada", round(r["dPs"]/1e3, 2), round(r["dPt"]/1e3, 2), "kPa", fmt="0.00")
    fila("Carga térmica Q", round(d["Q"], 2), "", "kW", fmt="#,##0.00")
    fila("DMLT (contracorriente)", round(d["DMLT"], 2), "", "°C", fmt="0.00")
    fila("Factor Ft", round(d["Ft"], 3), "", "-", fmt="0.000")
    row[0] += 1

    titulo("3. DATOS GEOMÉTRICOS Y CONSTRUCTIVOS")
    encabezado("Parámetro", "Haz tubular", "Carcasa", "Unidad")
    fila("Diámetro externo / interno de tubo", f"{do_t*1000:.2f} / {di_t*1000:.2f}", "", "mm")
    fila("Diámetro interno de carcasa Ds", "", round(r["Ds"]*1000, 1), "mm", fmt="0.0")
    fila("Número de tubos por casco (Nt)", r["Nt"], "", "-", key="Nt")
    fila("Longitud de tubo (L)", r["L"], "", "m", key="L", fmt="0.00")
    fila("Número de pasos por tubo", r["n"], "", "-")
    fila("Arreglo de tubos", arreglo, "", "-")
    fila("Paso de tubos (Pt)", round(pt*1000, 2), "", "mm", fmt="0.00")
    fila("Espaciado de deflectores (B)", "", round(r["B"]*1000, 1), "mm", fmt="0.0")
    fila("Superficie por casco", "", "", "m²")
    ref["A_casco"] = row[0] - 1
    row[0] += 1

    titulo("4. RESULTADOS DE DESEMPEÑO TÉRMICO E HIDRÁULICO")
    encabezado("Parámetro", "Valor", "", "Unidad")
    fila("ho (coeficiente película lado casco)", round(r["ho"], 1), "", "W/m²·K", fmt="#,##0.0")
    fila("hi (coeficiente película lado tubos)", round(r["hi"], 1), "", "W/m²·K", fmt="#,##0.0")
    fila("Uc (limpio)", round(r["Uc"], 1), "", "W/m²·K", key="Uc", fmt="#,##0.0")
    fila("Ud (diseño)", round(r["Ud"], 1), "", "W/m²·K", key="Ud", fmt="#,##0.0")
    fila("Rd_calc", "", "", "m²·K/W", key="Rd")
    fila("Rd_requerido", round(d["Rd_req"], 6), "", "m²·K/W", key="Rdreq", fmt="0.000000")
    fila("Área de intercambio requerida", round(r["A_req_final"], 2), "", "m²", key="Areq", fmt="0.00")
    fila("Área de intercambio proporcionada", "", "", "m²", key="Aprov")
    fila("Exceso de área", "", "", "%", key="Exc")
    fila("Cumple Rd y Δp", "SÍ" if (r["ok_Rd"] and r["ok_dp"]) else "NO", "", "-")

    ws.cell(ref["A_casco"], 3, f"=B{ref['Nt']}*PI()*{do_t}*B{ref['L']}")
    ws.cell(ref["A_casco"], 3).number_format = "0.00"
    ws.cell(ref["Rd"], 2, f"=(B{ref['Uc']}-B{ref['Ud']})/(B{ref['Uc']}*B{ref['Ud']})")
    ws.cell(ref["Rd"], 2).number_format = "0.000000"
    ws.cell(ref["Aprov"], 2, f"={r['Ns']}*C{ref['A_casco']}")
    ws.cell(ref["Aprov"], 2).number_format = "0.00"
    ws.cell(ref["Exc"], 2, f"=(B{ref['Aprov']}/B{ref['Areq']}-1)*100")
    ws.cell(ref["Exc"], 2).number_format = "0.0"

    ws2 = wb.create_sheet("Iteraciones")
    ws2.column_dimensions["A"].width = 160
    ws2.cell(1, 1, "Registro de decisiones del agente").font = Font(name=F, bold=True, size=11)
    for i, linea in enumerate(log, start=2):
        ws2.cell(i, 1, linea).font = Font(name=F, size=9)

    out = io.BytesIO()
    wb.save(out)
    out.seek(0)
    return out

def exportar_pdf_bytes(d, r, warns):
    buffer = io.BytesIO()
    doc = SimpleDocTemplate(buffer, pagesize=letter, rightMargin=30, leftMargin=30, topMargin=30, bottomMargin=30)
    elements = []
    styles = getSampleStyleSheet()
    title_style = ParagraphStyle('TitleStyle', parent=styles['Heading1'], fontSize=15, textColor=colors.HexColor('#321c54'), alignment=1, spaceAfter=12)
    subtitle_style = ParagraphStyle('SubTitleStyle', parent=styles['Heading2'], fontSize=11, textColor=colors.HexColor('#4a2b75'), spaceBefore=8, spaceAfter=4)
    normal_style = ParagraphStyle('NormalStyle', parent=styles['Normal'], fontSize=9, textColor=colors.HexColor('#222222'))

    elements.append(Paragraph("MEMORIA DE CÁLCULO TÉCNICO — CONDENSADOR TEMA", title_style))
    elements.append(Paragraph(f"<b>Equipo:</b> {d['identificacion']['item']} | <b>Servicio:</b> {d['identificacion']['servicio']}", subtitle_style))
    elements.append(Spacer(1, 10))

    data = [
        ["Parámetro Principal", "Resultado", "Unidad"],
        ["Cascos en serie (Ns)", str(r["Ns"]), "-"],
        ["Diámetro Carcasa (Ds)", f"{r['Ds']*1000:.1f}", "mm"],
        ["Número de Tubos (Nt)", str(r["Nt"]), "-"],
        ["Longitud de Tubo (L)", f"{r['L']:.2f}", "m"],
        ["Área Proporcionada", f"{r['A']:.2f}", "m²"],
        ["Coeficiente Uc", f"{r['Uc']:.1f}", "W/m²·K"],
        ["Coeficiente Ud", f"{r['Ud']:.1f}", "W/m²·K"],
        ["Δp Carcasa / Admisible", f"{r['dPs']/1e3:.1f} / {r['dp_s_adm']/1e3:.1f}", "kPa"],
        ["Δp Tubos / Admisible", f"{r['dPt']/1e3:.1f} / {r['dp_t_adm']/1e3:.1f}", "kPa"],
        ["Estado General", "CUMPLE" if (r["ok_Rd"] and r["ok_dp"]) else "REVISAR", "-"]
    ]

    t = Table(data, colWidths=[230, 150, 120])
    t.setStyle(TableStyle([
        ('BACKGROUND', (0,0), (-1,0), colors.HexColor('#321c54')),
        ('TEXTCOLOR', (0,0), (-1,0), colors.whitesmoke),
        ('ALIGN', (0,0), (-1,-1), 'CENTER'),
        ('FONTNAME', (0,0), (-1,0), 'Helvetica-Bold'),
        ('FONTSIZE', (0,0), (-1,0), 9),
        ('BACKGROUND', (0,1), (-1,-1), colors.HexColor('#f4f2f8')),
        ('GRID', (0,0), (-1,-1), 0.5, colors.HexColor('#d4bfff')),
        ('FONTNAME', (0,1), (-1,-1), 'Helvetica'),
        ('FONTSIZE', (0,1), (-1,-1), 9),
    ]))
    elements.append(t)
    
    if warns:
        elements.append(Spacer(1, 10))
        elements.append(Paragraph("<b>Observaciones de Ingeniería:</b>", subtitle_style))
        for w in warns:
            elements.append(Paragraph(f"• {w}", normal_style))

    doc.build(elements)
    buffer.seek(0)
    return buffer


# --- INTERFAZ PRINCIPAL CON PESTAÑAS Y MODO DE OPERACIÓN ---
st.title("Agente Inteligente de Diseño TEMA — Condensador")
st.markdown("---")

# Selector de Modo de Operación en la barra lateral superior
modo_app = st.radio("Modo de Operación:", ["Diseño (Sizing - Calcular Área Óptima)", "Verificación (Rating - Evaluar Equipo Existente)"], horizontal=True)
st.markdown("---")

tab_in1, tab_in2, tab_in3, tab_in4 = st.tabs([
    "Equipo", 
    "Vapor (lado casco)", 
    "Agua de enfriamiento (lado tubos)", 
    "Tubos y geometría"
])

with tab_in1:
    st.markdown("### Datos Generales del Equipo")
    c1, c2, c3 = st.columns(3)
    with c1: cliente = st.text_input("Cliente / Cátedra", value="Cátedra Transferencia de Energía Térmica")
    with c2: item_tag = st.text_input("Identificación del equipo", value="E-101")
    with c3: servicio = st.text_input("Servicio", value="Condensador de vapor de agua")
    tema_tipo = st.selectbox("Tipo de equipo (TEMA)", ["BEM", "AES", "AKT", "BEU"])

with tab_in2:
    st.markdown("### Condiciones del Vapor (Lado Casco)")
    
    compuesto_sel = st.selectbox("Seleccionar compuesto / fluido condensante predefinido:", list(FLUIDOS_DB.keys()))
    props_def = FLUIDOS_DB[compuesto_sel]

    c1, c2, c3, c4 = st.columns(4)
    with c1: fluido_v = st.text_input("Fluido que condensa", value=compuesto_sel)
    with c2: wh = st.number_input("Caudal de vapor [kg/h]", value=props_def["wh"], step=10.0)
    with c3: p_oper_s = st.number_input("Presión de operación [kPa abs]", value=props_def["p"], step=10.0)
    with c4: lam = st.number_input("Calor latente [kJ/kg]", value=props_def["lam"], step=10.0)

    c5, c6, c7 = st.columns(3)
    with c5: tsat = st.number_input("Temperatura de saturación [°C]", value=props_def["tsat"], step=0.1)
    with c6: t1_hot = st.number_input("Temperatura de entrada [°C]", value=props_def["tsat"], step=0.1)
    with c7: t2_hot = st.number_input("Temperatura de salida del condensado [°C]", value=props_def["tsat"], step=0.1)

    c8, c9 = st.columns(2)
    with c8: rd_hot = st.number_input("Factor de ensuciamiento requerido [m²·K/W]", value=0.000088, format="%.6f")
    with c9: dp_adm_s = st.number_input("Pérdida de carga admisible [kPa]", value=15.0, step=1.0)

    st.markdown("#### Propiedades del Vapor y Líquido")
    p1, p2, p3, p4 = st.columns(4)
    with p1: rho_v = st.number_input("Densidad vapor [kg/m³]", value=props_def["rho_v"])
    with p2: mu_v = st.number_input("Viscosidad vapor [mPa·s]", value=props_def["mu_v"])
    with p3: k_v = st.number_input("Conductividad vapor [W/m·K]", value=props_def["k_v"])
    with p4: cp_v = st.number_input("Calor específico vapor [kJ/kg·K]", value=props_def["cp_v"])

    p5, p6, p7, p8, p9 = st.columns(5)
    with p5: rho_l = st.number_input("Densidad líquido [kg/m³]", value=props_def["rho_l"])
    with p6: mu_l = st.number_input("Viscosidad líquido [mPa·s]", value=props_def["mu_l"])
    with p7: k_l = st.number_input("Conductividad líquido [W/m·K]", value=props_def["k_l"])
    with p8: cp_l = st.number_input("Calor específico líquido [kJ/kg·K]", value=props_def["cp_l"])
    with p9: t_prop = st.number_input("Temp. de propiedades [°C]", value=props_def["tsat"])

with tab_in3:
    st.markdown("### Agua de Enfriamiento (Lado Tubos)")
    st.info("💡 El caudal de agua fría se calcula automáticamente por balance térmico en función de las temperaturas de entrada y salida.")
    
    c1, c2 = st.columns(2)
    with c1: fluido_c = st.text_input("Fluido refrigerante", value="Agua de enfriamiento")
    with c2: st.markdown("")

    c3, c4 = st.columns(2)
    with c3: t1_cold = st.number_input("Temperatura de entrada [°C]", value=25.0, step=0.5)
    with c4: t2_cold = st.number_input("Temperatura de salida [°C]", value=45.0, step=0.5)

    c5, c6 = st.columns(2)
    with c5: rd_cold = st.number_input("Factor de ensuciamiento tubos [m²·K/W]", value=0.000176, format="%.6f")
    with c6: dp_adm_t = st.number_input("Pérdida de carga admisible tubos [kPa]", value=70.0, step=5.0)

with tab_in4:
    if "Diseño" in modo_app:
        st.markdown("### Geometría y Criterios de Diseño (Modo Sizing)")
        c1, c2, c3 = st.columns(3)
        with c1: do_t = st.number_input("Diámetro externo tubo [mm]", value=19.05) / 1000.0
        with c2: di_t = st.number_input("Diámetro interno tubo [mm]", value=15.75) / 1000.0
        with c3: l_max = st.number_input("Longitud máxima del tubo [m]", value=3.66)

        c4, c5, c6 = st.columns(3)
        with c4: arreglo = st.selectbox("Arreglo de tubos", ["Triangular (60°)", "Cuadrado (90°)"])
        with c5: pt = st.number_input("Paso entre tubos [mm]", value=23.81) / 1000.0
        with c6: k_wall = st.number_input("Conductividad pared [W/m·K]", value=50.0)

        c7, c8, c9, c10 = st.columns(4)
        with c7: v_min = st.number_input("Velocidad mín. tubos [m/s]", value=1.0)
        with c8: v_max = st.number_input("Velocidad máx. tubos [m/s]", value=3.0)
        with c9: ud_est_in = st.number_input("Coeficiente global estimado [W/m²K]", value=1500.0)
        with c10: ft_min = st.number_input("Ft mínimo", value=0.75)
    else:
        st.markdown("### Geometría del Equipo Existente (Modo Rating / Verificación)")
        c1, c2, c3 = st.columns(3)
        with c1: do_t = st.number_input("Diámetro externo tubo [mm]", value=19.05) / 1000.0
        with c2: di_t = st.number_input("Diámetro interno tubo [mm]", value=15.75) / 1000.0
        with c3: l_tube = st.number_input("Longitud de tubo [m]", value=3.66)

        c4, c5, c6 = st.columns(3)
        with c4: arreglo = st.selectbox("Arreglo de tubos", ["Triangular (60°)", "Cuadrado (90°)"])
        with c5: pt = st.number_input("Paso entre tubos [mm]", value=23.81) / 1000.0
        with c6: k_wall = st.number_input("Conductividad pared [W/m·K]", value=50.0)

        c7, c8, c9, c10 = st.columns(4)
        with c7: ds_val = st.number_input("Diámetro interno de carcasa Ds [mm]", value=635.0) / 1000.0
        with c8: nt_val = st.number_input("Número total de tubos (Nt)", value=216, step=2)
        with c9: n_pasos = st.selectbox("Pasos por tubos (n)", [2, 4, 6, 8], index=1)
        with c10: b_val = st.number_input("Espaciado deflectores B [mm]", value=200.0) / 1000.0

        ns_val = 1
        v_min, v_max, ft_min = 1.0, 3.0, 0.75

st.markdown("---")
btn_txt = "Calcular diseño óptimo" if "Diseño" in modo_app else "Verificar desempeño del equipo"
if st.button(btn_txt):
    st.session_state['run_design'] = True

if st.session_state.get('run_design', False):
    datos = {
        "identificacion": {"item": item_tag, "servicio": servicio, "cliente": cliente, "tema": tema_tipo},
        "caliente": {
            "Wh": wh, "T1": t1_hot, "T2": t2_hot, "T_sat": tsat, "lambda": lam * 1000.0,
            "rho_liquido": rho_l, "mu_liquido": mu_l * 1e-3, "k_liquido": k_l,
            "rho_vapor": rho_v, "mu_vapor": mu_v * 1e-3, "k_vapor": k_v, "cp_vapor": cp_v * 1000.0,
            "dp_adm": dp_adm_s * 1000.0, "Rd": rd_hot
        },
        "fria": {
            "t1": t1_cold, "t2": t2_cold, "cp": 4184.0, "rho": 995.0, "mu": 0.00076, "k": 0.62,
            "dp_adm": dp_adm_t * 1000.0, "Rd": rd_cold
        },
        "geometria": {
            "Do": do_t, "Di": di_t, "L_max": l_max if "Diseño" in modo_app else l_tube,
            "arreglo": "triangular" if "Triangular" in arreglo else "cuadro",
            "Pt": pt, "k_pared": k_wall, "v_min": v_min, "v_max": v_max, "Ft_min": ft_min
        },
        "condensa_en": "casco"
    }

    if "Diseño" in modo_app:
        r, log, warns = disenar(datos)
    else:
        # Modo Verificación Directa (Rating)
        preparar(datos, [])
        datos["DMLT"] = dmlt_contracorriente(t1_hot, t2_hot, t1_cold, t2_cold)
        datos["Ft"] = factor_ft(t1_hot, t2_hot, t1_cold, t2_cold, 1)
        datos["Ns"] = 1
        g_ing = dict(Ns=1, Nt=nt_val, n=n_pasos, L=l_tube, Ds=ds_val, B=b_val)
        r = evaluar(datos, g_ing)
        log = [
            f"Modo Verificación (Rating) ejecutado para el equipo {item_tag}.",
            f"Geometría fija: Ds = {ds_val*1000:.1f} mm, Nt = {nt_val}, L = {l_tube} m, Pasos = {n_pasos}",
            f"Carga térmica Q = {datos['Q']:.1f} kW, DMLT = {datos['DMLT']:.2f} K",
            f"Coeficiente Uc = {r['Uc']:.1f} W/m²K, Ud = {r['Ud']:.1f} W/m²K",
            f"Factor de ensuciamiento calculado (Rd) = {r['Rd_calc']:.6f} m²·K/W"
        ]
        warns = []
        if r['Rd_calc'] < datos['Rd_req']:
            warns.append("El equipo no alcanza el factor de ensuciamiento requerido (Rd_calc < Rd_req).")
        if not r['ok_dp']:
            warns.append("Se superan las caídas de presión admisibles en carcasa o tubos.")

    st.markdown("---")
    res_tab1, res_tab2, res_tab3, res_tab4 = st.tabs(["Verificaciones", "Detalle del equipo", "Observaciones", "Registro de iteraciones"])

    with res_tab1:
        st.markdown("### Verificaciones Normativas y de Proceso")
        st.markdown(f"""
        | Verificación | Resultado | Límite | Estado |
        | :--- | :--- | :--- | :--- |
        | **Caudal de agua calculado ($W_c$)** | {datos['fria']['Wc']:.1f} kg/h | Balance térmico | ℹ️ Automático |
        | **Factor de corrección Ft** | {datos['Ft']:.3f} | $\\ge {ft_min}$ | {'✅ Cumple' if datos['Ft'] >= ft_min else '⚠️ Revisar'} |
        | **Factor de ensuciamiento** | {r['Rd_calc']:.6f} m²·K/W | $\\ge {rd_hot}$ | {'✅ Cumple' if r['ok_Rd'] else '❌ No cumple'} |
        | **Velocidad del agua en tubos** | {r['v_t']:.2f} m/s | {v_min} a {v_max} m/s | {'✅ Cumple' if v_min <= r['v_t'] <= v_max else '⚠️ Revisar'} |
        | **Pérdida de carga en tubos** | {r['dPt']/1e3:.2f} kPa | $\\le {dp_adm_t}$ kPa | {'✅ Cumple' if r['dPt'] <= datos['fria']['dp_adm'] else '❌ Excede'} |
        | **Pérdida de carga en casco** | {r['dPs']/1e3:.2f} kPa | $\\le {dp_adm_s}$ kPa | {'✅ Cumple' if r['dPs'] <= datos['caliente']['dp_adm'] else '❌ Excede'} |
        """)

    with res_tab2:
        st.markdown("### Detalle Constructivo y Térmico")
        col_d1, col_d2 = st.columns(2)
        with col_d1:
            st.markdown(f"""
            * **Carcasas en serie:** {r['Ns']}
            * **Tipo de equipo:** {tema_tipo}
            * **Diámetro interno del casco:** {r['Ds']*1000:.0f} mm
            * **Tubos:** {do_t*1000:.2f} mm
            * **Longitud de tubo:** {r['L']:.2f} m
            * **Número de tubos (Nt):** {r['Nt']}
            * **Pasos por tubos:** {r['n']}
            * **Espaciado de deflectores (B):** {r['B']*1000:.0f} mm
            """)
        with col_d2:
            st.markdown(f"""
            * **Caudal de agua calculado ($W_c$):** {datos['fria']['Wc']:.1f} kg/h
            * **Coeficiente del vapor ($h_o$):** {r['ho']:.1f} W/m²·K
            * **Coeficiente del agua ($h_i$):** {r['hi']:.1f} W/m²·K
            * **DMLT $\\times$ Ft:** {datos['DMLT'] * datos['Ft']:.1f} K
            * **Coeficiente Limpio ($U_c$):** {r['Uc']:.1f} W/m²·K
            * **Coeficiente de Diseño ($U_d$):** {r['Ud']:.1f} W/m²·K
            * **Área de intercambio:** {r['A']:.2f} m²
            """)

    with res_tab3:
        st.markdown("### Observaciones de Criterio de Ingeniería")
        st.markdown("""
        * **Velocidad de vapor en casco:** Controlada dentro de márgenes recomendados para evitar vibración de tubos.
        * **Drenaje:** El condensado fluye hacia el punto bajo del casco; considerar ranuras de drenaje en deflectores.
        * **Tensiones térmicas:** Verificar diferencia de temperatura entre carcasa y tubos para evaluar junta de expansión.
        * **Limpieza:** Configuración apta para limpieza mecánica del interior de tubos retirando el cabezal.
        """)
        if warns:
            for w in warns:
                st.warning(w)

    with res_tab4:
        st.markdown("### Registro de Decisiones del Algoritmo")
        for l_log in log:
            st.text(l_log)

    # Gráfico de sensibilidad interactivo (solo aplicable o adaptado)
    if "Diseño" in modo_app:
        st.markdown("---")
        st.markdown("### Análisis Gráfico de Sensibilidad")
        wh_vals = np.linspace(wh * 0.5, wh * 1.5, 15)
        areas = []
        for w_t in wh_vals:
            d_t = json.loads(json.dumps(datos))
            d_t["caliente"]["Wh"] = w_t
            try:
                rt, _, _ = disenar(d_t)
                areas.append(rt["A"])
            except:
                areas.append(np.nan)

        fig, ax = plt.subplots(figsize=(10, 3.8))
        fig.patch.set_facecolor('#1e1730')
        ax.set_facecolor('#120f1d')
        ax.plot(wh_vals, areas, marker='o', color='#ff6b6b', linewidth=2.5)
        ax.set_xlabel('Caudal de Vapor [kg/h]', color='#e2d9ff')
        ax.set_ylabel('Área Requerida [m²]', color='#e2d9ff')
        ax.tick_params(colors='#bfaee3')
        for spine in ax.spines.values(): spine.set_edgecolor('#332752')
        ax.grid(True, linestyle='--', alpha=0.3, color='#4a2b75')
        st.pyplot(fig)

    # Descarga de reportes
    st.markdown("---")
    st.markdown("### Descarga de Documentación Oficial")
    col_dl1, col_dl2 = st.columns(2)
    with col_dl1:
        st.download_button("Descargar Hoja TEMA Clásica (.xlsx)", data=exportar_hoja_tema_bytes(datos, r, log, warns), file_name="hoja_TEMA_completa.xlsx", mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")
    with col_dl2:
        st.download_button("Descargar Memoria Técnica (.pdf)", data=exportar_pdf_bytes(datos, r, warns), file_name="memoria_calculo.pdf", mime="application/pdf")
