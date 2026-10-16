import streamlit as st
import json
import math
import io
from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill, Alignment, Border, Side

# --- DISEÑO Y ESTÉTICA (TONOS VIOLETA OSCURO) ---
st.markdown("""
    <style>
    /* Fondo general de la aplicación */
    .stApp {
        background-color: #120f1d;
        color: #e6e3f8;
    }
    
    /* Barra lateral en violeta oscuro profundo */
    [data-testid="stSidebar"] {
        background-color: #1a1528;
        border-right: 1px solid #2d2445;
    }
    
    /* Textos y títulos en la barra lateral y principal */
    [data-testid="stSidebar"] h1, [data-testid="stSidebar"] h2, [data-testid="stSidebar"] h3 {
        color: #d4bfff !important;
    }
    h1, h2, h3 {
        color: #e2d9ff !important;
    }
    
    /* Estilo de los botones principales */
    .stButton>button {
        background: linear-gradient(135deg, #6b3fa0, #4a2575);
        color: white;
        border: none;
        border-radius: 8px;
        font-weight: bold;
        padding: 10px 20px;
        width: 100%;
        box-shadow: 0 4px 10px rgba(107, 63, 160, 0.4);
        transition: all 0.3s ease;
    }
    .stButton>button:hover {
        background: linear-gradient(135deg, #7d4bc4, #5a2f8a);
        box-shadow: 0 6px 15px rgba(125, 75, 196, 0.6);
    }
    
    /* Cajas de métricas y resultados */
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
    
    /* Botón de descarga de Excel */
    .stDownloadButton>button {
        background: linear-gradient(135deg, #382159, #231338);
        color: #e2d9ff;
        border: 1px solid #7d4bc4;
        border-radius: 8px;
        font-weight: bold;
        width: 100%;
    }
    .stDownloadButton>button:hover {
        background: linear-gradient(135deg, #4b2d78, #321a52);
        color: white;
    }
    </style>
""", unsafe_allow_html=True)


# --- TUS FUNCIONES ORIGINALES DE CÁLCULO ---
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
    n_opts = [1, 2, 4, 6, 8] if (isoterma and rear
