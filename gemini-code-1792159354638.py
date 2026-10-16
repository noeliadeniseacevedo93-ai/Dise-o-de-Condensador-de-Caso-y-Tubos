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

# --- DISEÑO Y ESTÉTICA (TONOS VIOLETA OSCURO) ---
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
        background: linear-gradient(135deg, #6b3fa0, #4a2575);
        color: white;
        border: none;
        border-radius: 8px;
        font-weight: bold;
        padding: 12px 24px;
        width: 100%;
        box-shadow: 0 4px 10px rgba(107, 63, 160, 0.4);
        transition: all 0.3s ease;
    }
    .stButton>button:hover {
        background: linear-gradient(135deg, #7d4bc4, #5a2f8a);
        box-shadow: 0 6px 15px rgba(125, 75, 196, 0.6);
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


# --- FUNCIONES DE CÁLCULO ORIGINALES ---
G = 9.81  # m/s2
L_STD = [1.83, 2.44, 3.05, 3.66, 4.88, 6.10]
DS_STD = [x * 0.0254 for x in [8, 10, 12, 13.25, 15.25, 17.25, 19.25, 21.25, 23.25, 25, 27, 29, 31, 33, 35, 37, 39, 42, 45, 48, 54, 60]]

K_TRI = {1: (0.319, 2.142), 2: (0.249, 2.207), 4: (0.175, 2.285), 6: (0.0743, 2.499), 8: (0.0365, 2.675)}
K_CUA = {1: (0.215, 2.207), 2: (0.156, 2.291), 4: (0.158, 2.263), 6: (0.0402, 2.617), 8: (0.0331, 2.643)}
HOLGURA = {"M": 0.012, "N": 0.012, "L": 0.012, "U": 0.015, "S": 0.060, "T": 0.095, "P": 0.050, "W": 0.035}

def dmlt_contracorriente(T1, T2, t1, t2):
    d1, d2 = T1 - t2, T2 - t1
    if d1 <= 0 or d2 <= 0:
        raise ValueError("Cruce de temperaturas: T1 <= t2 o T2 <= t1.")
    if abs(d1 - d2) < 1e-9:
        return d1
    return (d1 - d2) / math.log(d1 / d2)

def factor_ft(T1, T2, t1, t2, Ns):
    R = (T1 - T2) / (t2 - t1)
    P = (t2 - t1) / (T1 - t1)
    if R < 1e-6: return 1.0
    if abs(R - 1) < 1e-6: R = 1 + 1e-6
    if P >= 1 or 1 - P * R <= 0: return 0.0
    F = ((1 - P * R) / (1 - P)) ** (1.0 / Ns)
    P1 = (1 - F) / (R - F)
    S = math.sqrt(R * R + 1)
    a = (1 - P1) / (1 - R * P1)
    b = (2 - P1 * (R + 1 - S)) / (2 - P1 * (R + 1 + S))
    if a <= 0 or b <= 0: return 0.0
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
    Gs = W / As
    De = 4 * (Pt ** 2 - math.pi * Do ** 2 / 4) / (math.pi * Do) if arreglo == "cuadro" else 4 * (0.5 * Pt * 0.866 * Pt - 0.5 * math.pi * Do ** 2 / 4) / (math.pi * Do / 2)
    Re = De * Gs / mu
    Pr = cp * mu / k if cp else 1.0
    hs = 0.36 * k / De * Re ** 0.55 * Pr ** (1 / 3) if k else 0.0
    f = math.exp(0.576 - 0.19 * math.log(Re))
    Nb = max(1, int(round(L / B - 1)))
    dP = f * Gs ** 2 * Ds * (Nb + 1) / (2 * rho * De)
    return dict(hs=hs, dP=dP, Re=Re, Gs=Gs, De=De, Nb=Nb, As=As, v=Gs / rho)

def h_cond_horizontal(k, rho, rhov, mu, lam, dT, Do, Nr):
    return 0.725 * (k ** 3 * rho * (rho - rhov) * G * lam / (mu * dT * Do)) ** 0.25 * Nr ** (-1 / 6)

def capacidad_tubos(Ds, n, Do, Pt, arreglo, rear):
    K1, n1 = (K_CUA if arreglo == "cuadro" else K_TRI)[n]
    Db = Ds - HOLGURA.get(rear, 0.012)
    if Db <= 0: return 0
    return int(K1 * (Db / Do) ** n1 * (1.25 * Do / Pt) ** 2)

def preparar(d, warns):
    h, c = d["caliente"], d["fria"]
    h["Tsat"] = h.get("T_sat") or (h["T1"] if abs(h["T1"] - h["T2"]) < 1e-9 else h["T1"])
    Q = h["Wh"] * (h["lambda"] + h.get("cp_vapor", 0) * max(h["T1"] - h["Tsat"], 0) + h.get("cp_liquido", 0) * max(h["Tsat"] - h["T2"], 0))
    Wc_bal = Q / (c["cp"] * (c["t2"] - c["t1"]))
    c["Wc"] = Wc_bal
    d["Q"] = Q
    geo = d["geometria"]
    d["Rd_req"] = (h["Rd"] + c["Rd"] * geo["Do"] / geo["Di"]) if d["condensa_en"] == "casco" else (c["Rd"] + h["Rd"] * geo["Do"] / geo["Di"])
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

    if d["condensa_en"] == "casco":
        Ap = Nt / n * math.pi * Di ** 2 / 4
        Gt = c["Wc"] / Ap
        Re = Gt * Di / c["mu"]
        Pr = c["cp"] * c["mu"] / c["k"]
        hi = nusselt_tubo(Re, Pr) * c["k"] / Di
        f = factor_friccion(Re, Di)
        dPt = Ns * f * (L * n / Di) * Gt ** 2 / (2 * c["rho"])
        r.update(hi=hi, v_t=Gt / c["rho"], Re_t=Re, Pr_t=Pr, dPt=dPt)
        Nr = max(1.0, 2.0 / 3.0 * Ds / Pt)
        dT, ho = 10.0, 5000.0
        for _ in range(80):
            ho = h_cond_horizontal(h["k_liquido"], h["rho_liquido"], h["rho_vapor"], h["mu_liquido"], h["lambda"], dT, Do, Nr)
            Uc = 1.0 / (1 / ho + Rw + Do / (Di * hi))
            dT_new = max(Uc * (Tsat - Tc) / ho, 0.5)
            if abs(dT_new - dT) < 1e-4: break
            dT = 0.5 * dT + 0.5 * dT_new
        k_ = kern_casco(h["Wh"], h["rho_vapor"], h["mu_vapor"], 1.0, 1.0, Ds, B, Pt, Do, L, arr)
        dPs = 0.5 * Ns * k_["dP"]
        r.update(ho=ho, dT_film=dT, dPs=dPs, v_s=k_["v"], Re_s=k_["Re"], Nb=k_["Nb"], De=k_["De"], As=k_["As"], Gs=k_["Gs"], Nr=Nr)
    
    Ud = d["Q"] / (A * d["DMLT"] * d["Ft"])
    r.update(Uc=Uc, Ud=Ud, Rd_calc=(Uc - Ud) / (Uc * Ud), A_req_final=d["Q"] / ((1 / (1 / Uc + d["Rd_req"])) * d["DMLT"] * d["Ft"]))
    dp_t_adm = c["dp_adm"]
    dp_s_adm = h["dp_adm"]
    r.update(dp_t_adm=dp_t_adm, dp_s_adm=dp_s_adm, ok_dp=(r["dPt"] <= dp_t_adm and r["dPs"] <= dp_s_adm), ok_Rd=(r["Rd_calc"] >= d["Rd_req"]))
    return r

def buscar_geometria(d, A_req, Ns, log):
    geo = d["geometria"]
    Do, Pt, arr, rear = geo["Do"], geo["Pt"], geo["arreglo"], d["rear"]
    h = d["caliente"]
    isoterma = abs(h["T1"] - h["T2"]) < 1e-9
    n_opts = [1, 2, 4, 6, 8] if (isoterma and rear != "U") else [2, 4, 6, 8]
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
            for idx in range(idx0, min(idx0 + 4, len(DS_STD))):
                Ds = DS_STD[idx]
                frac = 0.40
                while frac <= 1.0001:
                    B = max(min(frac * Ds, 1.52), 0.05, Ds / 5)
                    g = dict(Ns=Ns, Nt=Nt, n=n, L=L, Ds=Ds, B=B)
                    r = evaluar(d, g)
                    if r["dPs"] <= r["dp_s_adm"]: break
                    frac += 0.05
                v_ok = vmin <= r["v_t"] <= vmax if d["condensa_en"] == "casco" else True
                if r["ok_dp"] and v_ok: return r
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
        if Ft >= 0.75 or Ns >= 6: break
        Ns += 1
    d["Ft"], d["Ns"] = Ft, Ns
    Ud_est = d["geometria"]["Ud_est"]
    factibles, ultimo = [], None
    for _ in range(30):
        A_prov = d["Q"] / (Ud_est * d["DMLT"] * Ft)
        r = buscar_geometria(d, A_prov, Ns, log)
        if r is None:
            Ud_est *= 0.9
            continue
        ultimo = r
        exceso = r["A"] / r["A_req_final"] - 1
        Ud_tgt = 1 / (1 / r["Uc"] + d["Rd_req"])
        if r["ok_Rd"]:
            factibles.append(r)
            if exceso <= 0.25: break
            Ud_est = Ud_est * Ud_tgt / r["Ud"] * 0.97
        else:
            Ud_est = Ud_est * Ud_tgt / r["Ud"] * 0.98
    final = min(factibles, key=lambda x: x["A"]) if factibles else ultimo
    return final, log, warns

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
    for col in "BCD":
        ws.column_dimensions[col].width = 24
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
            if i > 0:
                c_.alignment = Alignment(horizontal="right")
        if key:
            ref[key] = row[0]
        row[0] += 1

    idn, h, c, geo = d["identificacion"], d["caliente"], d["fria"], d["geometria"]
    cas = d["condensa_en"] == "casco"
    S = (h if cas else c)
    T = (c if cas else h)

    ws.cell(row[0], 1, "HOJA DE ESPECIFICACIONES TÉCNICAS — INTERCAMBIADOR DE CASCO Y TUBOS (TEMA)")
    ws.merge_cells(start_row=row[0], start_column=1, end_row=row[0], end_column=5)
    ws.cell(row[0], 1).font = Font(name=F, bold=True, size=13)
    row[0] += 2

    titulo("1. DATOS DE IDENTIFICACIÓN DEL EQUIPO")
    fila("Ítem / Tag", idn["item"])
    fila("Servicio", idn["servicio"])
    fila("Cliente / Proyecto", idn["cliente"])
    fila("Tipo TEMA", idn["tema"])
    fila("Orientación", idn["orientacion"])
    fila("Cascos en serie", r["Ns"])
    row[0] += 1

    titulo("2. CONDICIONES DE OPERACIÓN Y PROCESO")
    encabezado("Parámetro", "Lado casco", "Lado tubos", "Unidad")
    fila("Fluido", S["fluido"], T["fluido"], "-")
    wS = h["Wh"] if cas else c["Wc"]
    wT = c["Wc"] if cas else h["Wh"]
    fila("Caudal másico", round(wS, 3), round(wT, 3), "kg/s", fmt="0.000")
    fila("Temperatura de entrada", S.get("T1", S.get("t1")), T.get("T1", T.get("t1")), "°C", fmt="0.0")
    fila("Temperatura de salida", S.get("T2", S.get("t2")), T.get("T2", T.get("t2")), "°C", fmt="0.0")
    fila("Presión de operación", round(S["P_oper"] / 1e3, 1), round(T["P_oper"] / 1e3, 1), "kPa abs", fmt="0.0")
    if cas:
        fila("Calor latente de condensación", round(h["lambda"] / 1e3, 1), "-", "kJ/kg", fmt="0.0")
        fila("Densidad (líquido / vapor)", f"{h['rho_liquido']} / {h['rho_vapor']}", c["rho"], "kg/m³")
        fila("Viscosidad", f"{h['mu_liquido'] * 1e3:.4f}", f"{c['mu'] * 1e3:.4f}", "mPa·s")
        fila("Conductividad térmica", h["k_liquido"], c["k"], "W/m·K")
        fila("Calor específico", "-", c["cp"], "J/kg·K")
    else:
        fila("Calor latente de condensación", "-", round(h["lambda"] / 1e3, 1), "kJ/kg", fmt="0.0")
        fila("Densidad", c["rho"], f"{h['rho_liquido']} / {h['rho_vapor']}", "kg/m³")
        fila("Viscosidad", f"{c['mu'] * 1e3:.4f}", f"{h['mu_liquido'] * 1e3:.4f}", "mPa·s")
        fila("Conductividad térmica", c["k"], h["k_liquido"], "W/m·K")
        fila("Calor específico", c["cp"], "-", "J/kg·K")
    fila("Velocidad", round(r["v_s"], 2), round(r["v_t"], 2), "m/s", fmt="0.00")
    fila("Factor de ensuciamiento", S["Rd"], T["Rd"], "m²·K/W", fmt="0.00000")
    fila("Δp admisible", round(r["dp_s_adm"] / 1e3, 1), round(r["dp_t_adm"] / 1e3, 1), "kPa", fmt="0.0")
    fila("Δp calculada", round(r["dPs"] / 1e3, 2), round(r["dPt"] / 1e3, 2), "kPa", fmt="0.00")
    fila("Carga térmica Q", round(d["Q"] / 1e3, 1), "", "kW", fmt="#,##0.0")
    fila("DMLT (contracorriente)", round(d["DMLT"], 2), "", "°C", fmt="0.00")
    fila("Factor Ft", round(d["Ft"], 3), "", "-", fmt="0.000")
    row[0] += 1

    titulo("3. DATOS GEOMÉTRICOS Y CONSTRUCTIVOS")
    encabezado("Parámetro", "Haz tubular", "Carcasa", "Unidad")
    fila("Diámetro externo / interno de tubo", f"{geo['Do'] * 1000:.2f} / {geo['Di'] * 1000:.2f}", "", "mm")
    fila("Diámetro interno de carcasa Ds", "", round(r["Ds"] * 1000, 1), "mm", fmt="0.0")
    fila("Número de tubos por casco (Nt)", r["Nt"], "", "-", key="Nt")
    fila("Longitud de tubo (L)", r["L"], "", "m", key="L", fmt="0.00")
    fila("Número de pasos por tubo", r["n"], "", "-")
    fila("Arreglo de tubos", "Triangular 60°" if geo["arreglo"] != "cuadro" else "Cuadro 90°", "", "-")
    fila("Paso de tubos (Pt)", round(geo["Pt"] * 1000, 2), "", "mm", fmt="0.00")
    fila("Material", idn["material_tubos"], idn["material_casco"], "-")
    fila("Espaciado de deflectores (B)", "", round(r["B"] * 1000, 1), "mm", fmt="0.0")
    fila("Número de deflectores", "", r["Nb"], "-")
    fila("Corte de deflector", "", "25 % segmentado", "-")
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
    fila("Δp calculada lado casco", round(r["dPs"] / 1e3, 2), "", "kPa", fmt="0.00")
    fila("Δp calculada lado tubos", round(r["dPt"] / 1e3, 2), "", "kPa", fmt="0.00")
    fila("Cumple Rd y Δp", "SÍ" if (r["ok_Rd"] and r["ok_dp"]) else "NO", "", "-")

    ws.cell(ref["A_casco"], 3, f"=B{ref['Nt']}*PI()*{geo['Do']}*B{ref['L']}")
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
    title_style = ParagraphStyle('TitleStyle', parent=styles['Heading1'], fontSize=16, textColor=colors.HexColor('#321c54'), alignment=1, spaceAfter=15)
    subtitle_style = ParagraphStyle('SubTitleStyle', parent=styles['Heading2'], fontSize=12, textColor=colors.HexColor('#4a2b75'), spaceBefore=10, spaceAfter=6)
    normal_style = ParagraphStyle('NormalStyle', parent=styles['Normal'], fontSize=9, textColor=colors.HexColor('#222222'))
    
    elements.append(Paragraph("REPORTE TÉCNICO DE DISEÑO TEMA", title_style))
    elements.append(Paragraph(f"<b>Equipo:</b> {d['identificacion']['item']} | <b>Servicio:</b> {d['identificacion']['servicio']}", subtitle_style))
    elements.append(Spacer(1, 10))
    
    data = [
        ["Parámetro de Desempeño", "Valor", "Unidad"],
        ["Cascos en serie (Ns)", str(r["Ns"]), "-"],
        ["Diámetro de Carcasa (Ds)", f"{r['Ds']*1000:.1f}", "mm"],
        ["Número de Tubos (Nt)", str(r["Nt"]), "-"],
        ["Longitud del Tubo (L)", f"{r['L']:.2f}", "m"],
        ["Coeficiente Limpio (Uc)", f"{r['Uc']:.1f}", "W/m²·K"],
        ["Coeficiente de Diseño (Ud)", f"{r['Ud']:.1f}", "W/m²·K"],
        ["Área de Intercambio Proporcionada", f"{r['A']:.2f}", "m²"],
        ["Δp Casco Calculada", f"{r['dPs']/1e3:.2f}", "kPa"],
        ["Δp Tubos Calculada", f"{r['dPt']/1e3:.2f}", "kPa"],
        ["Estado General", "CUMPLE" if (r["ok_Rd"] and r["ok_dp"]) else "NO CUMPLE", "-"]
    ]
    
    t = Table(data, colWidths=[230, 150, 120])
    t.setStyle(TableStyle([
        ('BACKGROUND', (0,0), (-1,0), colors.HexColor('#321c54')),
        ('TEXTCOLOR', (0,0), (-1,0), colors.whitesmoke),
        ('ALIGN', (0,0), (-1,-1), 'CENTER'),
        ('FONTNAME', (0,0), (-1,0), 'Helvetica-Bold'),
        ('FONTSIZE', (0,0), (-1,0), 9),
        ('BOTTOMPADDING', (0,0), (-1,0), 6),
        ('BACKGROUND', (0,1), (-1,-1), colors.HexColor('#f4f2f8')),
        ('GRID', (0,0), (-1,-1), 0.5, colors.HexColor('#d4bfff')),
        ('FONTNAME', (0,1), (-1,-1), 'Helvetica'),
        ('FONTSIZE', (0,1), (-1,-1), 9),
    ]))
    
    elements.append(t)
    
    if warns:
        elements.append(Spacer(1, 15))
        elements.append(Paragraph("<b>AVISOS DE DISEÑO:</b>", subtitle_style))
        for w in warns:
            elements.append(Paragraph(f"• {w}", normal_style))
            
    doc.build(elements)
    buffer.seek(0)
    return buffer


# --- INTERFAZ PRINCIPAL ---
st.title("Agente de Diseño TEMA — Condensador de Casco y Tubos")
st.markdown("---")

col_btn1, col_btn2 = st.columns([2, 1])
with col_btn1:
    ejecutar = st.button("🚀 Ejecutar Diseño Óptimo y Generar Hoja TEMA")

st.markdown("### 🎛️ Parámetros de Diseño y Proceso")

tab1, tab2, tab3, tab4 = st.tabs([
    "⚙️ Datos Operativos", 
    "🔥 Datos de Vapor (Caliente)", 
    "💧 Líquido Frío", 
    "📐 Datos Geométricos"
])

with tab1:
    st.markdown("#### Configuración Operativa General")
    c1, c2 = st.columns(2)
    with c1:
        item_tag = st.text_input("Ítem / Tag del equipo", value="E-101")
        servicio_desc = st.text_input("Servicio", value="Condensador de vapor de agua")
    with c2:
        tema_tipo = st.text_input("Tipo TEMA", value="BEM")

with tab2:
    st.markdown("#### Condiciones del Vapor y Propiedades")
    c1, c2 = st.columns(2)
    with c1:
        wh = st.number_input("Caudal másico vapor (Wh) [kg/s]", value=2.0, step=0.1)
        t1_hot = st.number_input("Temp. entrada vapor [°C]", value=120.2, step=0.1)
        t2_hot = st.number_input("Temp. salida vapor [°C]", value=120.2, step=0.1)
        tsat = st.number_input("Temperatura de condensación / sat. [°C]", value=120.2, step=0.1)
        lam = st.number_input("Calor latente (lambda) [J/kg]", value=2201900.0, step=100.0)
        rd_hot = st.number_input("Factor de ensuciamiento lado caliente (Rd) [m²·K/W]", value=0.00009, format="%.6f")
    with c2:
        rho_liq = st.number_input("Densidad líquido [kg/m³]", value=942.9)
        mu_liq = st.number_input("Viscosidad líquido [Pa·s]", value=0.000232, format="%.6f")
        k_liq = st.number_input("Conductividad líquido [W/m·K]", value=0.685)
        rho_vap = st.number_input("Densidad vapor [kg/m³]", value=1.129)
        mu_vap = st.number_input("Viscosidad vapor [Pa·s]", value=0.000013, format="%.6f")

with tab3:
    st.markdown("#### Condiciones del Líquido Frío y Propiedades")
    c1, c2 = st.columns(2)
    with c1:
        t1_cold = st.number_input("Temp. entrada agua fría (t1) [°C]", value=25.0, step=0.1)
        t2_cold = st.number_input("Temp. salida agua fría (t2) [°C]", value=40.0, step=0.1)
        rd_cold = st.number_input("Factor de ensuciamiento lado frío (Rd) [m²·K/W]", value=0.00035, format="%.6f")
    with c2:
        cp_c = st.number_input("Calor específico líquido frío [J/kg·K]", value=4178.0)
        rho_c = st.number_input("Densidad líquido frío [kg/m³]", value=995.0)
        mu_c = st.number_input("Viscosidad líquido frío [Pa·s]", value=0.00076, format="%.5f")
        k_c = st.number_input("Conductividad líquido frío [W/m·K]", value=0.62)

with tab4:
    st.markdown("#### Parámetros Geométricos y Constructivos")
    c1, c2 = st.columns(2)
    with c1:
        do_tube = st.number_input("Diámetro externo tubo (Do) [m]", value=0.01905, format="%.5f")
        di_tube = st.number_input("Diámetro interno tubo (Di) [m]", value=0.01483, format="%.5f")
        l_max = st.number_input("Longitud máxima de tubo (L_max) [m]", value=4.88, step=0.1)
    with c2:
        pt_pitch = st.number_input("Paso de tubos (Pt) [m]", value=0.02381, format="%.5f")
        k_wall = st.number_input("Conductividad pared tubo [W/m·K]", value=50.0)

if ejecutar:
    st.session_state['run_design'] = True

if st.session_state.get('run_design', False):
    datos = {
        "identificacion": {"item": item_tag, "servicio": servicio_desc, "cliente": "Práctica Integradora", "tema": tema_tipo, "orientacion": "Horizontal", "material_tubos": "Acero al carbono", "material_casco": "Acero al carbono"},
        "caliente": {"fluido": "Vapor de agua saturado", "Wh": wh, "T1": t1_hot, "T2": t2_hot, "T_sat": tsat, "P_oper": 200000.0, "lambda": lam, "cp_vapor": 0.0, "cp_liquido": 0.0, "rho_liquido": rho_liq, "mu_liquido": mu_liq, "k_liquido": k_liq, "rho_vapor": rho_vap, "mu_vapor": mu_vap, "dp_adm": 20000.0, "Rd": rd_hot},
        "fria": {"fluido": "Agua de enfriamiento", "Wc": None, "t1": t1_cold, "t2": t2_cold, "P_oper": 400000.0, "cp": cp_c, "mu": mu_c, "k": k_c, "rho": rho_c, "dp_adm": 70000.0, "Rd": rd_cold},
        "geometria": {"Do": do_tube, "Di": di_tube, "L_max": l_max, "arreglo": "triangular", "Pt": pt_pitch, "k_pared": k_wall, "Ud_est": 2000.0, "v_min": 1.0, "v_max": 3.0},
        "condensa_en": "casco"
    }
    
    r, log, warns = disenar(datos)
    
    st.markdown("---")
    st.header("📊 Resultados del Desempeño Térmico e Hidráulico")
    
    col1, col2, col3, col4 = st.columns(4)
    with col1:
        st.metric("Cascos en serie (Ns)", r["Ns"])
        st.metric("Diámetro carcasa (Ds)", f"{r['Ds']*1000:.1f} mm")
    with col2:
        st.metric("Número de tubos (Nt)", r["Nt"])
        st.metric("Longitud del tubo (L)", f"{r['L']:.2f} m")
    with col3:
        st.metric("Uc (Limpio)", f"{r['Uc']:.1f} W/m²K")
        st.metric("Ud (Diseño)", f"{r['Ud']:.1f} W/m²K")
    with col4:
        st.metric("Área de Intercambio", f"{r['A']:.2f} m²")
        st.metric("Espac. Deflectores (B)", f"{r['B']*1000:.1f} mm")
        
    if warns:
        for w in warns:
            st.warning(w)

    st.markdown("---")
    st.header("📈 Perfil de Temperaturas del Condensador")
    
    x_vals = np.linspace(0, r['L'], 100)
    T_hot = [tsat] * len(x_vals)
    t_cold = np.linspace(t1_cold, t2_cold, len(x_vals))
    
    fig, ax = plt.subplots(figsize=(10, 4.5))
    fig.patch.set_facecolor('#1e1730')
    ax.set_facecolor('#120f1d')
    
    ax.plot(x_vals, T_hot, label='Fluido Caliente (Vapor Sat.)', color='#d4bfff', linewidth=3)
    ax.plot(x_vals, t_cold, label='Fluido Frío (Agua)', color='#7d4bc4', linewidth=3)
    
    ax.set_xlabel('Posición axial en el intercambiador [m]', color='#e2d9ff', fontsize=10)
    ax.set_ylabel('Temperatura [°C]', color='#e2d9ff', fontsize=10)
    ax.set_title('Perfil Térmico Axial', color='#e2d9ff', fontsize=12, fontweight='bold')
    ax.legend(facecolor='#1a1528', edgecolor='#332752', labelcolor='#e2d9ff')
    ax.tick_params(colors='#bfaee3')
    for spine in ax.spines.values():
        spine.set_edgecolor('#332752')
    
    st.pyplot(fig)

    st.markdown("---")
    st.header("📉 Gráfico Interactivo de Sensibilidad")
    st.markdown("Analizá cómo varía el **Área de Intercambio requerida** al modificar el caudal de vapor másico ($Wh$).")
    
    wh_rangos = np.linspace(wh * 0.5, wh * 1.5, 20)
    areas_sensibilidad = []
    
    for w_test in wh_rangos:
        d_test = json.loads(json.dumps(datos))
        d_test["caliente"]["Wh"] = w_test
        try:
            r_test, _, _ = disenar(d_test)
            areas_sensibilidad.append(r_test["A"])
        except:
            areas_sensibilidad.append(np.nan)
            
    fig2, ax2 = plt.subplots(figsize=(10, 4.2))
    fig2.patch.set_facecolor('#1e1730')
    ax2.set_facecolor('#120f1d')
    
    ax2.plot(wh_rangos, areas_sensibilidad, marker='o', color='#bfaee3', linewidth=2.5, markerfacecolor='#6b3fa0')
    ax2.set_xlabel('Caudal Másico de Vapor ($Wh$) [kg/s]', color='#e2d9ff', fontsize=10)
    ax2.set_ylabel('Área de Intercambio ($A$) [m²]', color='#e2d9ff', fontsize=10)
    ax2.set_title('Sensibilidad del Área ante Cambios en el Caudal', color='#e2d9ff', fontsize=12, fontweight='bold')
    ax2.tick_params(colors='#bfaee3')
    for spine in ax2.spines.values():
        spine.set_edgecolor('#332752')
    ax2.grid(True, linestyle='--', alpha=0.3, color='#4a2b75')
    
    st.pyplot(fig2)

    st.markdown("---")
    st.header("📝 Registro de Iteraciones del Agente")
    with st.expander("Ver log paso a paso de la búsqueda geométrica", expanded=False):
        for linea in log:
            st.text(linea)

    st.markdown("---")
    st.header("📥 Descarga de Reportes y Especificaciones")
    
    col_dl1, col_dl2 = st.columns(2)
    
    excel_file = exportar_hoja_tema_bytes(datos, r, log, warns)
    pdf_file = exportar_pdf_bytes(datos, r, warns)
    
    with col_dl1:
        st.download_button(
            label="📥 Descargar Hoja TEMA en Excel (.xlsx)",
            data=excel_file,
            file_name="hoja_TEMA_actualizada.xlsx",
            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
        )
        
    with col_dl2:
        st.download_button(
            label="📄 Descargar Reporte Técnico en PDF",
            data=pdf_file,
            file_name="reporte_condensador.pdf",
            mime="application/pdf"
        )
