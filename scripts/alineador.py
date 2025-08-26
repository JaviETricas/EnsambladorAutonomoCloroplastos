#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from pathlib import Path
import platform
import re
import subprocess
import sys


def fasta_to_oneline_py(src: Path, dst: Path) -> None:
    """Convierte cualquier FASTA a formato 'una línea por secuencia' (cabecera + secuencia sin saltos)."""
    with src.open("r", encoding="utf-8", errors="replace") as fi, \
         dst.open("w", encoding="utf-8") as fo:
        seq_parts = []
        header = None
        for line in fi:
            if line.startswith(">"):
                if header is not None:
                    fo.write("".join(seq_parts) + "\n")
                header = line.rstrip("\n")
                fo.write(header + "\n")
                seq_parts = []
            else:
                seq_parts.append(line.strip())
        if header is not None:
            fo.write("".join(seq_parts) + "\n")


def fasta_to_oneline_awk(src: Path, dst: Path) -> None:
    """Versión AWK (Unix)."""
    cmd = [
        "awk",
        r'/^>/ {if (seq) print seq; print; seq=""; next} {seq=seq $0} END {print seq}',
        str(src)
    ]
    with dst.open("w", encoding="utf-8") as fo:
        subprocess.run(cmd, stdout=fo, check=True)


def get_oneliner() -> callable:
    """Elige implementación oneline (Python en Windows, AWK en Unix si está disponible)."""
    if platform.system() == "Windows":
        return fasta_to_oneline_py
    # intenta AWK, cae a Python si falla
    try:
        subprocess.run(["awk", "--version"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
        return fasta_to_oneline_awk
    except Exception:
        return fasta_to_oneline_py


def normalize_basename(name: str) -> str:
    """
    Limpia sufijos comunes de nombres generados (Novowrap, .fastq.gz, etc.)
    para usarlo como 'hint' de cabecera.
    """
    n = name
    n = re.sub(r"Option_\d+_", "", n)
    n = re.sub(r"_Novowrap.*", "", n)
    n = re.sub(r"_pileup", "", n)
    n = re.sub(r"\.fastq(?:\.gz)?_?", "_", n)
    n = re.sub(r"\.(fasta|fa|fna|tsv|txt|fai)$", "", n, flags=re.IGNORECASE)
    return n

def _smith_waterman_endpoints(a: str, b: str, match: int = 2, mismatch: int = -1, gap: int = -2) -> tuple[int, int, int, int] | None:
    """
    Smith–Waterman local entre 'a' y 'b'.
    Devuelve (a_start, a_end, b_start, b_end) en coordenadas 1-based relativas a 'a' y 'b',
    o None si no hay alineamiento local con puntuación > 0.
    Pensado para ventanas cortas (O(len(a)*len(b))).
    """
    la, lb = len(a), len(b)
    if la == 0 or lb == 0:
        return None

    H = [[0] * (lb + 1) for _ in range(la + 1)]
    T = [[0] * (lb + 1) for _ in range(la + 1)]  # 0=stop,1=diag,2=up,3=left

    best = 0
    bi = bj = 0

    for i in range(1, la + 1):
        ai = a[i - 1]
        for j in range(1, lb + 1):
            bj = b[j - 1]
            s = match if ai == bj else mismatch
            diag = H[i - 1][j - 1] + s
            up   = H[i - 1][j] + gap
            left = H[i][j - 1] + gap
            val  = max(0, diag, up, left)
            H[i][j] = val
            if val == 0:
                T[i][j] = 0
            elif val == diag:
                T[i][j] = 1
            elif val == up:
                T[i][j] = 2
            else:
                T[i][j] = 3
            if val > best:
                best = val
                bi, bj = i, j

    if best <= 0:
        return None

    i, j = bi, bj
    while i > 0 and j > 0 and H[i][j] > 0:
        move = T[i][j]
        if move == 1:
            i, j = i - 1, j - 1
        elif move == 2:
            i -= 1
        elif move == 3:
            j -= 1
        else:
            break

    a_start = i + 1
    b_start = j + 1
    a_end   = bi
    b_end   = bj
    return (a_start, a_end, b_start, b_end)


def _micro_refine_with_sw(cons_ung: str,
                          new_ung: str,
                          q_st: int,
                          q_en: int,
                          c_st: int,
                          c_en: int,
                          flank: int = 250,
                          post: int = 60) -> tuple[int, int] | None:
    """
    Afina los cortes (q_st, q_en) con SW local alrededor de los bordes.
    - Izquierdo: mueve q_st al primer nt DESPUÉS del mejor match local cerca del borde.
    - Derecho : mueve q_en al nt ANTERIOR al arranque del mejor match local del borde derecho.
    Devuelve (q_st', q_en') o None si no hay mejora válida.
    """
    nQ = len(new_ung)
    nC = len(cons_ung)
    if not (1 <= q_st < q_en <= nQ) or not (1 <= c_st < c_en <= nC):
        return None

    # Izquierdo
    cL0 = max(1, c_st - flank)
    qL0 = max(1, q_st - flank)
    cL1 = min(nC, c_st - 1 + post)
    qL1 = min(nQ, q_st - 1 + post)
    left_cons = cons_ung[cL0 - 1:cL1]
    left_new  = new_ung[qL0 - 1:qL1]

    left_end_q = None
    if len(left_cons) >= 20 and len(left_new) >= 20:
        sw = _smith_waterman_endpoints(left_cons, left_new)
        if sw:
            _, _, b_start, b_end = sw
            left_end_q = (qL0 - 1) + b_end  # 1-based

    # Derecho
    cR0 = max(1, c_en + 1 - post)
    qR0 = max(1, q_en + 1 - post)
    cR1 = min(nC, c_en + flank)
    qR1 = min(nQ, q_en + flank)
    right_cons = cons_ung[cR0 - 1:cR1]
    right_new  = new_ung[qR0 - 1:qR1]

    right_start_q = None
    if len(right_cons) >= 20 and len(right_new) >= 20:
        sw = _smith_waterman_endpoints(right_cons, right_new)
        if sw:
            _, _, b_start, _ = sw
            right_start_q = (qR0 - 1) + b_start  # 1-based

    new_q_st = (left_end_q + 1) if left_end_q is not None else q_st
    new_q_en = (right_start_q - 1) if right_start_q is not None else q_en

    if new_q_st < new_q_en and (new_q_en - new_q_st + 1) >= 200:
        return (new_q_st, new_q_en)
    return None




def load_mutations(tsv_path: Path) -> dict[int, str]:
    """
    Lee el TSV con cabecera:
    CHROM, POS, ., ,, A, C, G, T, N, +, -, ^, $, other
    Devuelve un dict pos(1-based) -> base (A/C/G/T) a corregir.
    """
    muts: dict[int, str] = {}
    with tsv_path.open("r", encoding="utf-8", errors="replace") as tsv:
        for line in tsv:
            if not line.strip() or line.startswith("#") or line.startswith("CHROM"):
                continue
            cols = line.rstrip().split("\t")
            try:
                pos   = int(cols[1])  # POS
                dot   = int(cols[2])  # '.'
                comma = int(cols[3])  # ','
                A = int(cols[4]); C = int(cols[5]); G = int(cols[6]); T = int(cols[7])
            except (ValueError, IndexError):
                continue

            total = dot + comma + A + C + G + T
            if total == 0:
                continue

            # si la referencia no domina, elegir mayoría A/C/G/T
            if (dot + comma == 0) or ((dot + comma) / total < 0.5):
                bases  = ("A", "C", "G", "T")
                counts = (A, C, G, T)
                muts[pos] = bases[max(range(4), key=lambda i: counts[i])]
    return muts
def align_with_reference_and_number(oneline_fasta: Path,
                                    min_identity: float = 0.80,
                                    enable_ssc_fix: bool = True) -> Path | None:
    """
    Alinea 'oneline_fasta' vs Fusion_Cloroplastos_<max>.fasta (o referencia base).
    Guarda SOLO el resultado ONELINE (sin sufijo).
    Si enable_ssc_fix:
      - Detecta bloque con alta desidentidad (umbral adaptativo + fallback Kadane).
      - Refina límites con anclas MEM (±1500 nt; k_init=27→k_min=17) o, en su defecto, anclas simples (±1000; k=25).
      - Micro‑afina con Smith–Waterman en los bordes.
      - Invierte, re‑alinea con MAFFT y compara Jaccard (k=21) antes/después.
      - Elige la mejor variante.
    """
    import os
    ROOT_DIR = Path(__file__).resolve().parents[1]
    LIB_DIR  = ROOT_DIR / "libreris"
    RES_DIR  = ROOT_DIR / "Resultados"
    REF_DIR  = ROOT_DIR / "temporalDocs" / "cloroplasto_referencia"
    DISCARD_DIR = ROOT_DIR / "temporalDocs" / "descartesalineamiento"

    RES_DIR.mkdir(parents=True, exist_ok=True)

    # --- referencia y numeración ---
    pattern = re.compile(r"^Fusion_Cloroplastos_(\d+)\.fasta$", re.IGNORECASE)
    existentes: list[tuple[int, Path]] = []
    for p in RES_DIR.iterdir():
        if p.is_file():
            m = pattern.match(p.name)
            if m:
                try:
                    existentes.append((int(m.group(1)), p))
                except ValueError:
                    pass
    if existentes:
        existentes.sort(key=lambda x: x[0])
        max_n, ref_fasta = existentes[-1]
        next_n = max_n + 1
    else:
        ref_fasta = REF_DIR / "Fusion_Cloroplastos_1.fasta"
        if not ref_fasta.exists():
            print(f"[WARN] No hay previos en {RES_DIR} y no existe referencia: {ref_fasta}")
            return None
        next_n = 1

    out_path = RES_DIR / f"Fusion_Cloroplastos_{next_n}.fasta"

    for label, f in (("referencia", ref_fasta), ("query", oneline_fasta)):
        if not f.exists() or f.stat().st_size == 0:
            print(f"[WARN] FASTA de {label} inexistente o vacío: {f}")
            return None

    # --- rotación circular ---
    rotated_query = rotate_to_reference_anchor(ref_fasta, oneline_fasta, k_first=200, max_mismatches=10)
    cleanup_rotated = (rotated_query != oneline_fasta)

    # --- helpers MAFFT ---
    def resolve_mafft() -> list[str]:
        cand = LIB_DIR / "mafft"
        if cand.exists() and cand.is_file() and os.access(cand, os.X_OK):
            return [str(cand)]
        return ["mafft"]

    def run_mafft_input(input_fa: Path) -> tuple[int, str, str]:
        cmd = resolve_mafft() + ["--auto", str(input_fa)]
        print(f"[MAFFT] {' '.join(cmd)} -> {out_path}")
        res = subprocess.run(cmd, text=True, capture_output=True)
        return res.returncode, res.stdout, res.stderr

    # --- primer alineado ---
    tmp_in = RES_DIR / f".tmp_mafft_input_{next_n}.fasta"
    try:
        with tmp_in.open("w", encoding="utf-8") as fo:
            fo.write(ref_fasta.read_text(encoding="utf-8", errors="replace")); fo.write("\n")
            fo.write(rotated_query.read_text(encoding="utf-8", errors="replace")); fo.write("\n")
        rc, out, err = run_mafft_input(tmp_in)
        if rc != 0 and resolve_mafft()[0] != "mafft":
            cmd = ["mafft", "--auto", str(tmp_in)]
            print(f"[MAFFT] {' '.join(cmd)} -> {out_path}")
            p = subprocess.run(cmd, text=True, capture_output=True)
            rc, out, err = p.returncode, p.stdout, p.stderr
    finally:
        try: tmp_in.unlink(missing_ok=True)
        except Exception: pass

    if rc != 0:
        print(f"[WARN] MAFFT falló (rc={rc}). stderr:\n{err}")
        if cleanup_rotated:
            try: rotated_query.unlink(missing_ok=True)
            except Exception: pass
        return None

    # --- SSC-fix opcional ---
    flipped_query_path = None
    out_best = out
    if enable_ssc_fix:
        entries = _parse_fasta_string(out)
        if len(entries) >= 2:
            prev_aln = [seq for _h, seq in entries[:-1]]
            new_aln  = entries[-1][1]
            try:
                cons_aln = _consensus_from_entries(prev_aln)
            except ValueError:
                cons_aln = None

            if cons_aln and len(cons_aln) == len(new_aln):
                block = _detect_inversion_block(cons_aln, new_aln, window=400, thr=0.60, min_span=2000)
                print(f"[SSC] Bloque detectado (cols) = {block}")
                if block:
                    s_col, e_col = block

                    # Mapeos base columna->pos para disponer SIEMPRE de c_st/c_en y q_st/q_en
                    col2q = _map_alncol_to_pos(new_aln)
                    col2c = _map_alncol_to_pos(cons_aln)
                    q_st_base = next((col2q[i] for i in range(s_col, e_col + 1) if col2q[i] != -1), None)
                    q_en_base = next((col2q[i] for i in range(e_col, s_col - 1, -1) if col2q[i] != -1), None)
                    c_st      = next((col2c[i] for i in range(s_col, e_col + 1) if col2c[i] != -1), None)
                    c_en      = next((col2c[i] for i in range(e_col, s_col - 1, -1) if col2c[i] != -1), None)

                    # Iniciar cortes con los base (por seguridad)
                    q_st, q_en = q_st_base, q_en_base

                    # 1) Refinado por MEM anchors (preferido)
                    refined = _refine_cutpoints_by_mem_anchors(cons_aln, new_aln, s_col, e_col,
                                                               flank_nt=1500, k_init=27, k_min=17, search_span=400)
                    if refined:
                        q_st, q_en = refined
                        print(f"[SSC] Cortes tras MEM-anchors -> query 1-based: {q_st}-{q_en}")
                    else:
                        # 2) Fallback: anclas simples
                        refined_simple = _refine_cutpoints_by_anchors(cons_aln, new_aln, s_col, e_col,
                                                                      flank_nt=1000, k=25)
                        if refined_simple:
                            q_st, q_en = refined_simple
                            print(f"[SSC] Cortes tras anclas simples -> query 1-based: {q_st}-{q_en}")
                        else:
                            print(f"[SSC] Sin anclas; uso cortes base -> query 1-based: {q_st}-{q_en}")

                    # Validación de cortes antes de micro‑SW
                    if None not in (q_st, q_en, c_st, c_en) and q_st < q_en and c_st < c_en:
                        # 3) Micro‑afinado SW en bordes
                        cons_ung = _ungap(cons_aln)
                        new_ung  = _ungap(new_aln)
                        micro = _micro_refine_with_sw(cons_ung, new_ung, q_st, q_en, c_st, c_en,
                                                      flank=250, post=60)
                        if micro:
                            q_st, q_en = micro
                            print(f"[SSC] Micro‑afinado SW -> query 1-based: {q_st}-{q_en}")
                    else:
                        print("[SSC] Cortes inválidos; se omite micro‑SW.")

                    # Si tenemos cortes válidos, probar flip + re‑MAFFT
                    if q_st is not None and q_en is not None and q_st < q_en:
                        jacc_before = _kmer_jaccard(_ungap(cons_aln), _ungap(new_aln), k=21)
                        print(f"[SSC] Jaccard antes del flip: {jacc_before:.4f}")

                        try:
                            flipped_query_path = _flip_segment_in_oneline(rotated_query, q_st, q_en)
                        except Exception as e:
                            print(f"[WARN] No se pudo invertir SSC: {e}")
                            flipped_query_path = None

                        if flipped_query_path and flipped_query_path.exists():
                            tmp_in2 = RES_DIR / f".tmp_mafft_input_{next_n}_flip.fasta"
                            try:
                                with tmp_in2.open("w", encoding="utf-8") as fo:
                                    fo.write(ref_fasta.read_text(encoding="utf-8", errors="replace")); fo.write("\n")
                                    fo.write(flipped_query_path.read_text(encoding="utf-8", errors="replace")); fo.write("\n")
                                rc2, out2, err2 = run_mafft_input(tmp_in2)
                                if rc2 != 0:
                                    cmd = ["mafft", "--auto", str(tmp_in2)]
                                    print(f"[MAFFT] {' '.join(cmd)} -> {out_path}")
                                    p2 = subprocess.run(cmd, text=True, capture_output=True)
                                    rc2, out2, err2 = p2.returncode, p2.stdout, p2.stderr
                            finally:
                                try: tmp_in2.unlink(missing_ok=True)
                                except Exception: pass

                            if rc2 == 0:
                                entries2 = _parse_fasta_string(out2)
                                if len(entries2) >= 2:
                                    prev2 = [seq for _h, seq in entries2[:-1]]
                                    new2  = entries2[-1][1]
                                    try:
                                        cons2 = _consensus_from_entries(prev2)
                                    except ValueError:
                                        cons2 = None
                                    if cons2 and len(cons2) == len(new2):
                                        jacc_after = _kmer_jaccard(_ungap(cons2), _ungap(new2), k=21)
                                        print(f"[SSC] Jaccard tras flip: {jacc_after:.4f}")
                                        out_best = out2 if jacc_after >= jacc_before else out
                                    else:
                                        print("[SSC] No se pudo construir consenso tras flip.")

    # limpieza
    if flipped_query_path:
        try: flipped_query_path.unlink(missing_ok=True)
        except Exception: pass
    if cleanup_rotated:
        try: rotated_query.unlink(missing_ok=True)
        except Exception: pass

    # --- filtro global y guardado oneline ---
    ok, identity, considered, L = _assess_alignment_quality(out_best, min_identity=min_identity)

    tmp_aln = RES_DIR / f".tmp_mafft_out_{next_n}.fasta"
    try:
        with tmp_aln.open("w", encoding="utf-8") as fo:
            fo.write(out_best)

        if ok:
            write_oneline_copy(tmp_aln, out_path)
            print(f"[OK] Alineado (oneline): {out_path}")
            return out_path
        else:
            DISCARD_DIR.mkdir(parents=True, exist_ok=True)
            discard_path = DISCARD_DIR / out_path.name
            write_oneline_copy(tmp_aln, discard_path)
            print(f"[REJECT] Identidad {identity:.3f} (< {min_identity:.2f}) en {considered}/{L} posiciones comparables.")
            print(f"[INFO] Alineado (oneline) guardado en descartes: {discard_path}")
            return None
    finally:
        try: tmp_aln.unlink(missing_ok=True)
        except Exception: pass


def _refine_cutpoints_by_mem_anchors(consensus_aln: str,
                                     new_aln: str,
                                     s_col: int,
                                     e_col: int,
                                     flank_nt: int = 1500,
                                     k_init: int = 27,
                                     k_min: int = 17,
                                     search_span: int = 400) -> tuple[int, int] | None:
    """
    Refina los cortes (q_st, q_en) buscando ANCLAS MEM (maximal exact matches) en los flancos.
    - Flanco IZQUIERDO: preferimos un MEM que ACABE pegado al borde (q_st / c_st).
    - Flanco DERECHO : preferimos un MEM que EMPIECE pegado al borde (q_en / c_en).

    Parámetros:
      flank_nt    : tamaño del flanco a inspeccionar en consenso y query.
      k_init      : tamaño inicial de semilla k-mer.
      k_min       : tamaño mínimo de semilla si no hay matches con k_init.
      search_span : en cada flanco, sólo buscamos semillas cerca del borde (últimos/primeros 'search_span' nt).

    Devuelve (q_start_1based, q_end_1based) refinados sobre la QUERY sin gaps,
    o None si no se pudo refinar.
    """
    # 1) Validaciones básicas y mapeos columna->pos
    if len(consensus_aln) != len(new_aln) or len(new_aln) == 0:
        return None

    def _map_alncol_to_pos(seq_aln: str) -> list[int]:
        col2pos = [-1] * len(seq_aln)
        pos = 0
        for i, ch in enumerate(seq_aln):
            u = ch.upper()
            if u in "ACGTN":
                pos += 1
                col2pos[i] = pos
        return col2pos

    col2q = _map_alncol_to_pos(new_aln)
    col2c = _map_alncol_to_pos(consensus_aln)

    q_st = next((col2q[i] for i in range(s_col, e_col + 1) if col2q[i] != -1), None)
    q_en = next((col2q[i] for i in range(e_col, s_col - 1, -1) if col2q[i] != -1), None)
    c_st = next((col2c[i] for i in range(s_col, e_col + 1) if col2c[i] != -1), None)
    c_en = next((col2c[i] for i in range(e_col, s_col - 1, -1) if col2c[i] != -1), None)
    if None in (q_st, q_en, c_st, c_en) or not (q_st < q_en and c_st < c_en):
        return None

    cons_ung = _ungap(consensus_aln)
    new_ung  = _ungap(new_aln)
    nQ = len(new_ung); nC = len(cons_ung)

    # 2) Construcción de flancos (slices 0-based)
    #    Izquierdo: [.. c_st], [.. q_st] con un pequeño post hacia dentro para dar margen
    post = 60
    cL0 = max(1, c_st - flank_nt)
    qL0 = max(1, q_st - flank_nt)
    cL1 = min(nC, c_st - 1 + post)
    qL1 = min(nQ, q_st - 1 + post)
    left_cons = cons_ung[cL0-1:cL1]
    left_new  = new_ung[qL0-1:qL1]

    #    Derecho: [c_en ..], [q_en ..] también con post hacia dentro
    cR0 = max(1, c_en + 1 - post)
    qR0 = max(1, q_en + 1 - post)
    cR1 = min(nC, c_en + flank_nt)
    qR1 = min(nQ, q_en + flank_nt)
    right_cons = cons_ung[cR0-1:cR1]
    right_new  = new_ung[qR0-1:qR1]

    def _build_index(s: str, k: int) -> dict[str, list[int]]:
        """Índice simple k-mer -> lista de posiciones (0-based) en 's'."""
        idx: dict[str, list[int]] = {}
        if len(s) < k:
            return idx
        for i in range(len(s) - k + 1):
            kmer = s[i:i+k]
            if "N" in kmer:
                continue
            if not all(c in "ACGTacgt" for c in kmer):
                continue
            idx.setdefault(kmer.upper(), []).append(i)
        return idx

    def _extend_mem(a: str, ia: int, b: str, ib: int, k: int) -> tuple[int, int, int]:
        """
        Extiende una semilla exacta a[ia:ia+k] == b[ib:ib+k] en ambas direcciones.
        Devuelve (start_a, start_b, total_len) en 0-based.
        """
        # Izquierda
        L = 0
        while ia-1-L >= 0 and ib-1-L >= 0:
            ca = a[ia-1-L].upper(); cb = b[ib-1-L].upper()
            if ca != cb or ca not in "ACGT":
                break
            L += 1
        # Derecha
        R = 0
        La = len(a); Lb = len(b)
        while ia+k+R < La and ib+k+R < Lb:
            ca = a[ia+k+R].upper(); cb = b[ib+k+R].upper()
            if ca != cb or ca not in "ACGT":
                break
            R += 1
        start_a = ia - L
        start_b = ib - L
        total = k + L + R
        return start_a, start_b, total

    # 3) Buscar MEM en flanco IZQUIERDO (preferimos que ACABE cerca del borde)
    left_anchor_q_end = None  # coordenada 1-based en la query global
    best_len_left = -1
    if len(left_cons) >= k_min and len(left_new) >= k_min:
        for k in range(k_init, k_min-1, -2):  # probar tamaños decrecientes
            idx = _build_index(left_cons, k)
            if not idx:
                continue
            # exploramos la parte final del flanco de la query (cerca del borde)
            start_i = max(0, len(left_new) - search_span - k)
            end_i   = max(0, len(left_new) - k)
            for i in range(start_i, end_i+1):
                kmer = left_new[i:i+k].upper()
                if "N" in kmer:
                    continue
                pos_list = idx.get(kmer)
                if not pos_list:
                    continue
                for j in pos_list:
                    sa, sb, tot = _extend_mem(left_cons, j, left_new, i, k)
                    # fin en la query local = sb + tot - 1
                    q_end_local = sb + tot - 1
                    # preferimos tot más largo; empate -> fin más cercano al borde (mayor q_end_local)
                    if tot > best_len_left or (tot == best_len_left and q_end_local > (0 if left_anchor_q_end is None else (left_anchor_q_end - (qL0-1)))):
                        best_len_left = tot
                        left_anchor_q_end = (qL0 - 1) + q_end_local
            if left_anchor_q_end is not None:
                break  # ya tenemos un buen ancla, no bajar más k

    # 4) Buscar MEM en flanco DERECHO (preferimos que EMPIECE cerca del borde)
    right_anchor_q_start = None
    best_len_right = -1
    if len(right_cons) >= k_min and len(right_new) >= k_min:
        for k in range(k_init, k_min-1, -2):
            idx = _build_index(right_cons, k)
            if not idx:
                continue
            # exploramos la parte inicial del flanco de la query (cerca del borde)
            start_i = 0
            end_i   = min(len(right_new) - k, search_span)
            for i in range(start_i, end_i+1):
                kmer = right_new[i:i+k].upper()
                if "N" in kmer:
                    continue
                pos_list = idx.get(kmer)
                if not pos_list:
                    continue
                for j in pos_list:
                    sa, sb, tot = _extend_mem(right_cons, j, right_new, i, k)
                    # inicio en la query local = sb
                    q_start_local = sb
                    # preferimos tot más largo; empate -> inicio más cercano al borde (menor q_start_local)
                    if tot > best_len_right or (tot == best_len_right and q_start_local < (len(right_new) if right_anchor_q_start is None else (right_anchor_q_start - (qR0-1)))):
                        best_len_right = tot
                        right_anchor_q_start = (qR0 - 1) + q_start_local
            if right_anchor_q_start is not None:
                break

    # 5) Ajustar cortes (si hay anclas)
    new_q_st = (left_anchor_q_end + 1) if left_anchor_q_end is not None else q_st
    new_q_en = (right_anchor_q_start - 1) if right_anchor_q_start is not None else q_en

    # Validaciones
    if not (1 <= new_q_st < new_q_en <= nQ):
        return None
    # Evitar recortes minúsculos
    if (new_q_en - new_q_st + 1) < max(500, k_min * 5):
        return None

    return new_q_st, new_q_en



def rotate_to_reference_anchor(ref_fasta: Path, query_fasta: Path,
                               k_first: int = 200, max_mismatches: int = 10) -> Path:
    """
    Rota la secuencia de 'query_fasta' para que empiece alineada con la referencia ('ref_fasta'),
    útil para genomas circulares.

    Estrategia:
      - Prueba orientación forward y reverse-complement, elige la que mejor encaje.
      - Usa anclas largas en el inicio de ref: first_k = ref[:k_first] y circ = ref[-k_first:]+ref[:k_first].
      - Tolera mismatches hasta 'max_mismatches' y trata 'N' como comodín.
      - Si no encuentra coincidencias claras, devuelve la query ORIGINAL.

    Devuelve:
      - Ruta a un archivo FASTA temporal rotado (si hubo rotación/orientación).
      - Si no hay cambio necesario, devuelve la ruta original.
    """
    def revcomp(seq: str) -> str:
        comp = str.maketrans("ACGTacgtnN", "TGCAtgcanN")
        return seq.translate(comp)[::-1]

    def read_first_record(path: Path) -> tuple[str, str]:
        hdr = None; parts = []
        with path.open("r", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if line.startswith(">"):
                    if hdr is None:
                        hdr = line.rstrip("\n")
                    else:
                        break
                else:
                    if hdr is not None:
                        parts.append(line.strip())
        if hdr is None:
            raise ValueError(f"FASTA sin cabecera: {path}")
        return hdr, "".join(parts)

    def match_with_mismatches(text: str, pattern: str, max_mm: int) -> int:
        """
        Busca 'pattern' en 'text' permitiendo hasta max_mm mismatches.
        'N' en cualquiera de las secuencias se considera comodín (no suma mismatch).
        Devuelve índice o -1 si no encontrado.
        """
        m = len(pattern); n = len(text)
        if m > n:
            return -1
        for i in range(0, n - m + 1):
            mm = 0
            for a, b in zip(pattern, text[i:i+m]):
                if a == 'N' or b == 'N':
                    continue
                if a != b:
                    mm += 1
                    if mm > max_mm:
                        break
            if mm <= max_mm:
                return i
        return -1

    # Leer referencia y query
    _, ref = read_first_record(ref_fasta)
    qh, qseq_orig = read_first_record(query_fasta)
    if len(ref) < k_first or len(qseq_orig) == 0:
        return query_fasta

    # Construir anclas
    ref = ref.upper()
    qseq_orig = qseq_orig.upper()
    first = ref[:k_first]
    last  = ref[-k_first:]
    circ  = last + first

    candidates = [("fwd", qseq_orig), ("rev", revcomp(qseq_orig))]
    best = None  # (orientation, cut_index, mm_used, rule)

    for orient, qseq in candidates:
        doubled = qseq + qseq

        # A) ¿ya empieza por first?
        mm_head = sum(0 if (a == b or a == 'N' or b == 'N') else 1
                      for a, b in zip(first, qseq[:k_first]))
        if mm_head <= max_mismatches:
            if best is None or mm_head < best[2]:
                best = (orient, 0, mm_head, "HEAD_FIRST")
            if mm_head == 0:
                break

        # B) patrón circular last+first
        pos = match_with_mismatches(doubled, circ, max_mismatches)
        if pos != -1:
            cut = (pos + k_first) % len(qseq)
            if best is None or max_mismatches < best[2]:
                best = (orient, cut, max_mismatches, "CIRC_PATTERN")

        # C) buscar solo first en circular
        pos2 = match_with_mismatches(doubled, first, max_mismatches)
        if pos2 != -1 and pos2 < len(qseq):
            cut = pos2 % len(qseq)
            if best is None or max_mismatches < best[2]:
                best = (orient, cut, max_mismatches, "FIRST_ONLY")

    if best is None:
        return query_fasta

    orient, cut, _, _ = best
    qseq = qseq_orig if orient == "fwd" else revcomp(qseq_orig)
    if cut == 0 and orient == "fwd":
        return query_fasta

    rotated = qseq[cut:] + qseq[:cut]
    out = query_fasta.with_name(query_fasta.stem + f"_rotated_k{k_first}.fasta")
    with out.open("w", encoding="utf-8") as fo:
        fo.write(qh + "\n")
        fo.write(rotated + "\n")
    return out


def write_oneline_copy(src_fasta: Path, dst_oneline: Path) -> None:
    oneliner = get_oneliner()
    dst_oneline.parent.mkdir(parents=True, exist_ok=True)
    oneliner(src_fasta, dst_oneline)


def _parse_fasta_string(s: str) -> list[tuple[str, str]]:
    """
    Parsea un texto FASTA (ej. stdout de MAFFT) y devuelve [(header, seq_alineada), ...].
    - Mantiene el orden de aparición (MAFFT por defecto conserva el orden de entrada).
    - 'header' incluye el símbolo '>' inicial.
    - 'seq_alineada' conserva gaps ('-' o '.').
    """
    entries = []
    header = None
    seq_parts = []
    for line in s.splitlines():
        if line.startswith(">"):
            if header is not None:
                entries.append((header, "".join(seq_parts)))
            header = line.rstrip("\n")
            seq_parts = []
        else:
            if header is not None:
                seq_parts.append(line.strip())
    if header is not None:
        entries.append((header, "".join(seq_parts)))
    return entries

def _assess_alignment_quality(mafft_output: str, min_identity: float = 0.80) -> tuple[bool, float, int, int]:
    """
    Evalúa la calidad del alineamiento de MAFFT midiendo la identidad de la ÚLTIMA secuencia
    (la nueva) contra el CONSENSO de todas las anteriores.

    Retorna:
      (ok, identidad, posiciones_consideradas, longitud_alineada)

    Criterios:
      - Consenso por columna considerando solo A/C/G/T de las secuencias previas (ignora gaps y 'N').
      - Para la nueva secuencia, solo se consideran posiciones donde ella tiene A/C/G/T (ignora gaps y 'N').
      - Identidad = matches / posiciones_consideradas.
      - ok = identidad >= min_identity. Si no hay posiciones_consideradas, ok=False.
    """
    entries = _parse_fasta_string(mafft_output)
    if len(entries) < 2:
        # Sin suficiente material para evaluar; por seguridad, marcar como no ok.
        return (False, 0.0, 0, 0)

    # Separar secuencias previas y la última (nueva)
    prev = [seq for _h, seq in entries[:-1]]
    new  = entries[-1][1]

    # Comprobar longitudes coherentes
    L = len(new)
    if L == 0 or any(len(s) != L for s in prev):
        return (False, 0.0, 0, L)

    # Consenso por columna entre previas
    def consensus_base(chars: list[str]) -> str | None:
        counts = {'A':0, 'C':0, 'G':0, 'T':0}
        for ch in chars:
            u = ch.upper()
            if u in ('A','C','G','T'):
                counts[u] += 1
        # elegir la base con mayor recuento
        best_base = None
        best_count = 0
        for b in ('A','C','G','T'):
            if counts[b] > best_count:
                best_base, best_count = b, counts[b]
        return best_base if best_count > 0 else None

    matches = 0
    considered = 0
    for i in range(L):
        cons = consensus_base([s[i] for s in prev])
        if cons is None:
            continue  # columna sin consenso (todo gaps/N)
        q = new[i].upper()
        if q in ('-','.', 'N'):
            continue  # no considerar gaps ni N en la query
        considered += 1
        if q == cons:
            matches += 1

    if considered == 0:
        return (False, 0.0, 0, L)

    identity = matches / considered
    ok = identity >= min_identity
    return (ok, identity, considered, L)


def apply_mutations_single_contig(fasta_in: Path, fasta_out: Path, muts: dict[int, str]) -> int:
    """
    Aplica mutaciones sobre el PRIMER contig del FASTA (típico cloroplasto circular 1-contig).
    Devuelve la longitud final de la secuencia.
    """
    with fasta_in.open("r", encoding="utf-8", errors="replace") as fi, \
         fasta_out.open("w", encoding="utf-8") as fo:
        header = None
        seq_parts = []
        for line in fi:
            if line.startswith(">"):
                header = line.rstrip("\n")
                fo.write(header + "\n")
                break
        for line in fi:
            if line.startswith(">"):
                # ignoramos contigs adicionales en esta versión
                break
            seq_parts.append(line.strip())

        seq = list("".join(seq_parts))
        for pos, base in muts.items():
            idx = pos - 1
            if 0 <= idx < len(seq):
                seq[idx] = base

        fo.write("".join(seq) + "\n")
        return len(seq)


def rename_header_inplace(fasta_path: Path, name_hint: str, species: str = "Hordeum vulgare") -> None:
    """
    Renombra el encabezado del PRIMER contig como:
        >{name_hint}, {species}, {len}
    (len = nº de nucleótidos de la secuencia del primer contig)
    """
    tmp = fasta_path.with_suffix(fasta_path.suffix + ".tmp")
    with fasta_path.open("r", encoding="utf-8", errors="replace") as fi, \
         tmp.open("w", encoding="utf-8") as fo:
        header = None
        seq_parts = []
        for line in fi:
            if line.startswith(">"):
                header = line.rstrip("\n")
                # Dejamos la construcción del nuevo header más adelante (cuando sepamos len)
                break
        for line in fi:
            if line.startswith(">"):
                # sólo primer contig
                break
            seq_parts.append(line.strip())

        seq = "".join(seq_parts)
        new_header = f">{name_hint}, {species}, {len(seq)}"
        fo.write(new_header + "\n")
        fo.write(seq + "\n")

    tmp.replace(fasta_path)


def _consensus_from_entries(aligned_seqs: list[str]) -> str:
    """
    Construye el consenso columna-a-columna a partir de una lista de secuencias ALINEADAS (con gaps).
    - Considera sólo A/C/G/T para el consenso (ignora gaps y N).
    - Si no hay base dominante en una columna, pone 'N'.
    - Todas las secuencias deben tener la MISMA longitud (longitud del alineamiento).
    """
    if not aligned_seqs:
        return ""
    L = len(aligned_seqs[0])
    for s in aligned_seqs:
        if len(s) != L:
            raise ValueError("Secuencias alineadas con longitudes distintas al construir el consenso.")
    out = []
    for i in range(L):
        counts = {'A':0,'C':0,'G':0,'T':0}
        for s in aligned_seqs:
            ch = s[i].upper()
            if ch in counts:
                counts[ch] += 1
        best, cnt = None, 0
        for b in ('A','C','G','T'):
            if counts[b] > cnt:
                best, cnt = b, counts[b]
        out.append(best if cnt > 0 else 'N')
    return "".join(out)

def _detect_inversion_block(consensus_aln: str,
                            new_aln: str,
                            window: int = 400,
                            thr: float = 0.60,
                            min_span: int = 2000) -> tuple[int, int] | None:
    """
    Detecta un bloque grande con ALTA desidentidad (posible SSC invertida).
    Mejora:
      - Umbral adaptativo según el mismatch global.
      - Fallback tipo 'Kadane' (subarray de máximo exceso de mismatches) si no hay ventanas >= umbral.
    Devuelve: (start_col, end_col) índices 0-based del ALINEAMIENTO, o None.
    """
    L = len(new_aln)
    if len(consensus_aln) != L or L == 0:
        return None

    comp = [0]*L   # comparable (A/C/G/T vs A/C/G/T)
    mis  = [0]*L   # mismatch si comparable
    for i in range(L):
        a = consensus_aln[i].upper()
        b = new_aln[i].upper()
        if a in "ACGT" and b in "ACGT":
            comp[i] = 1
            mis[i]  = 1 if a != b else 0

    # Estadística global de mismatches en posiciones comparables
    total_comp = sum(comp)
    total_mis  = sum(mis)
    global_ratio = (total_mis / total_comp) if total_comp > 0 else 0.0

    # Umbral adaptativo: baja si el alineamiento global es “bueno”
    # (p.ej., identidad 0.92 => mismatch 0.08 => local_thr ≈ 0.35)
    local_thr = max(0.35, min(thr, global_ratio + 0.25))

    from itertools import accumulate
    pref_comp = [0] + list(accumulate(comp))
    pref_mis  = [0] + list(accumulate(mis))

    bad_windows: list[tuple[int,int]] = []
    step = max(100, window//4)
    for s in range(0, L - window + 1, step):
        e = s + window
        comp_w = pref_comp[e] - pref_comp[s]
        if comp_w < window*0.5:
            continue
        mis_w  = pref_mis[e] - pref_mis[s]
        ratio  = mis_w/comp_w if comp_w > 0 else 0.0
        if ratio >= local_thr:
            bad_windows.append((s, e-1))

    if bad_windows:
        # unir ventanas y escoger el bloque con más posiciones comparables
        bad_windows.sort()
        cur_s, cur_e = bad_windows[0]
        merged: list[tuple[int,int]] = []
        for s, e in bad_windows[1:]:
            if s <= cur_e + step:
                cur_e = max(cur_e, e)
            else:
                merged.append((cur_s, cur_e))
                cur_s, cur_e = s, e
        merged.append((cur_s, cur_e))

        best = None
        best_comp = -1
        for s, e in merged:
            comp_blk = pref_comp[e+1] - pref_comp[s]
            if comp_blk >= min_span and comp_blk > best_comp:
                best_comp = comp_blk
                best = (s, e)
        if best:
            return best

    # --- Fallback: Kadane sobre exceso de mismatches (mismatch=+1, match=-1) ---
    score = [0]*L
    for i in range(L):
        if comp[i]:
            score[i] = 1 if mis[i] else -1
        else:
            score[i] = 0  # columnas no comparables no aportan

    max_sum = cur_sum = 0
    best_s = best_e = -1
    temp_s = 0
    for i, val in enumerate(score):
        cur_sum += val
        if cur_sum <= 0:
            cur_sum = 0
            temp_s = i + 1
        elif cur_sum > max_sum:
            max_sum = cur_sum
            best_s, best_e = temp_s, i

    if max_sum > 0 and best_s >= 0:
        comp_blk = pref_comp[best_e+1] - pref_comp[best_s]
        if comp_blk >= max(min_span, 8000):  # algo conservador para SSC
            return (best_s, best_e)

    return None



def _map_alncol_to_pos(seq_aln: str) -> list[int]:
    """
    Dada una secuencia ALINEADA (con gaps), devuelve un vector 'col2pos' tal que:
      - col2pos[i] = posición 1-based en la secuencia SIN GAPS correspondiente a la columna i,
      - o -1 si en esa columna hay gap (o carácter no nucleotídico) en esta secuencia.

    Nota: sustituye a _map_alncol_to_qpos (mismo propósito, más genérico).
    """
    col2pos = [-1] * len(seq_aln)
    pos = 0
    for i, ch in enumerate(seq_aln):
        u = ch.upper()
        if u in "ACGTN":
            pos += 1
            col2pos[i] = pos
    return col2pos

def _refine_cutpoints_by_anchors(consensus_aln: str,
                                 new_aln: str,
                                 s_col: int,
                                 e_col: int,
                                 flank_nt: int = 1000,
                                 k: int = 25) -> tuple[int, int] | None:
    """
    Refina los límites [s_col, e_col] (índices de COLUMNA del alineamiento) buscando anclas de k-mer
    EXACTAS en los 1.000 nt (por defecto) ANTES y DESPUÉS del bloque candidato.
    Devuelve (q_start_1based, q_end_1based) REFINADOS sobre la QUERY (SIN GAPS),
    o None si no se puede refinar.

    Estrategia:
      1) Convertir columnas de alineamiento a posiciones 1-based en consensus y query (sin gaps).
      2) Extraer los flancos: [q_st-flank, q_st) y (q_en, q_en+flank] (y lo mismo en consensus).
      3) Buscar un k-mer idéntico:
         - En el flanco IZQUIERDO: desde el borde hacia atrás (más cercano al borde).
         - En el flanco DERECHO: desde el borde hacia delante.
      4) Si se hallan anclas, ajustar:
         q_start' = (fin de ancla izquierda) + 1
         q_end'   = (inicio de ancla derecha) - 1
         (así no tocamos las zonas idénticas).
      5) Comprobar que q_start' < q_end' y que el tramo no colapsa. Si algo falla, None.

    Requisitos:
      - consensus_aln y new_aln deben tener la MISMA longitud (alineamiento).
    """
    if len(consensus_aln) != len(new_aln) or len(new_aln) == 0:
        return None

    # Mapas columna->posición sin gaps
    col2q = _map_alncol_to_pos(new_aln)
    col2c = _map_alncol_to_pos(consensus_aln)

    # Posiciones iniciales 1-based (query y consensus) derivadas de columnas
    q_st = next((col2q[i] for i in range(s_col, e_col + 1) if col2q[i] != -1), None)
    q_en = next((col2q[i] for i in range(e_col, s_col - 1, -1) if col2q[i] != -1), None)
    c_st = next((col2c[i] for i in range(s_col, e_col + 1) if col2c[i] != -1), None)
    c_en = next((col2c[i] for i in range(e_col, s_col - 1, -1) if col2c[i] != -1), None)
    if None in (q_st, q_en, c_st, c_en) or not (q_st < q_en and c_st < c_en):
        return None

    cons_ung = _ungap(consensus_aln)
    new_ung  = _ungap(new_aln)
    nQ = len(new_ung); nC = len(cons_ung)

    # Flancos (índices 1-based -> slices 0-based)
    q_left_start = max(1, q_st - flank_nt)
    c_left_start = max(1, c_st - flank_nt)
    left_new  = new_ung[q_left_start - 1: q_st - 1]
    left_cons = cons_ung[c_left_start - 1: c_st - 1]

    q_right_start = q_en + 1
    c_right_start = c_en + 1
    q_right_end = min(nQ, q_en + flank_nt)
    c_right_end = min(nC, c_en + flank_nt)
    right_new  = new_ung[q_right_start - 1: q_right_end]
    right_cons = cons_ung[c_right_start - 1: c_right_end]

    # Búsqueda de ancla izquierda (desde el borde hacia atrás)
    left_anchor_q_end = None  # pos 1-based en query donde TERMINA el k-mer
    if len(left_new) >= k and len(left_cons) >= k:
        for i in range(len(left_new) - k, -1, -1):  # cerca del borde primero
            kmer = left_new[i:i + k]
            if "N" in kmer:
                continue
            if left_cons.find(kmer) != -1:
                left_anchor_q_end = (q_left_start - 1) + i + k  # 1-based end
                break

    # Búsqueda de ancla derecha (desde el borde hacia delante)
    right_anchor_q_start = None  # pos 1-based en query donde EMPIEZA el k-mer
    if len(right_new) >= k and len(right_cons) >= k:
        for i in range(0, len(right_new) - k + 1):
            kmer = right_new[i:i + k]
            if "N" in kmer:
                continue
            if right_cons.find(kmer) != -1:
                right_anchor_q_start = (q_right_start - 1) + i + 1  # 1-based start
                break

    # Ajuste de límites
    new_q_st = (left_anchor_q_end + 1) if left_anchor_q_end is not None else q_st
    new_q_en = (right_anchor_q_start - 1) if right_anchor_q_start is not None else q_en

    # Validaciones mínimas
    if not (1 <= new_q_st < new_q_en <= nQ):
        return None
    # Evitar recortes patológicos (p.ej. región demasiado pequeña)
    if (new_q_en - new_q_st + 1) < max(500, k * 5):  # heurística conservadora
        return None

    return new_q_st, new_q_en


def _flip_segment_in_oneline(oneline_fa: Path, start_1based: int, end_1based: int) -> Path:
    """
    Invierte (reverse-complement) el segmento [start_1based, end_1based] (INCLUSIVO)
    sobre un FASTA oneline (1 sola secuencia). Devuelve la ruta al nuevo FASTA creado.
    """
    def revcomp(seq: str) -> str:
        tbl = str.maketrans("ACGTacgtnN", "TGCAtgcanN")
        return seq.translate(tbl)[::-1]

    header = None
    seq = []
    with oneline_fa.open() as fh:
        for line in fh:
            if line.startswith(">") and header is None:
                header = line.rstrip("\n")
            elif not line.startswith(">"):
                seq.append(line.strip())
    s = "".join(seq)
    a = s[:start_1based-1]
    b = s[start_1based-1:end_1based]
    c = s[end_1based:]
    flipped = a + revcomp(b) + c
    out = oneline_fa.with_name(oneline_fa.stem + f"_sscflip_{start_1based}_{end_1based}.fasta")
    with out.open("w") as fo:
        fo.write(header + "\n")
        fo.write(flipped + "\n")
    return out


def _ungap(s: str) -> str:
    """Elimina gaps de una secuencia alineada."""
    return s.replace("-", "").replace(".", "")


def _kmer_jaccard(s1: str, s2: str, k: int = 21) -> float:
    """
    Jaccard de conjuntos de k-mers (en mayúsculas, sin 'N'); devuelve 0 si union vacía.
    Se espera que s1 y s2 estén SIN GAPS previamente.
    """
    s1 = s1.upper(); s2 = s2.upper()
    if len(s1) < k or len(s2) < k:
        return 0.0
    def kmerset(s: str) -> set[str]:
        out = set()
        for i in range(len(s)-k+1):
            kmer = s[i:i+k]
            if "N" in kmer:
                continue
            if all(c in "ACGT" for c in kmer):
                out.add(kmer)
        return out
    A = kmerset(s1); B = kmerset(s2)
    if not A and not B:
        return 0.0
    inter = len(A & B)
    uni   = len(A | B)
    return inter/uni if uni else 0.0


def main():
    """
    Flujo:
      1) Cargar TSV y decidir mutaciones por mayoría.
      2) Aplicar mutaciones sobre primer contig y escribir FASTA corregido.
      3) Renombrar cabecera del primer contig.
      4) Generar copia oneline del FASTA corregido.
      5) Alinear con referencia y numerar en Resultados/, usando el umbral --confianza.
         * Opcionalmente, intentar corregir una SSC invertida (desactivable con --no-ssc-fix).
    """
    ap = argparse.ArgumentParser(
        description="Corrige un FASTA (TSV mpileup), genera oneline y alinea numerando contra referencia.",
        allow_abbrev=False
    )
    ap.add_argument("--fasta", required=True, type=Path, help="FASTA a corregir (1 contig).")
    ap.add_argument("--tsv", required=True, type=Path, help="TSV de recuentos generado por BAMtsv.py.")
    ap.add_argument("--out", type=Path, default=None, help="FASTA corregido (por defecto: *_corrected.fasta).")
    ap.add_argument("--species", type=str, default="Hordeum vulgare", help="Nombre de especie para el header.")
    ap.add_argument("--oneline", action="store_true",
                    help="(Compat) Aceptado, pero el oneline se genera siempre en un fichero aparte.")
    ap.add_argument("--confianza", dest="confianza", type=float, default=0.95,
                    help="Identidad mínima para aceptar el alineamiento (por defecto 0.95).")
    ap.add_argument("--no-ssc-fix", action="store_true",
                    help="Desactiva la detección/corrección automática de una posible SSC invertida.")

    args = ap.parse_args()

    fasta_in = args.fasta.resolve()
    tsv = args.tsv.resolve()
    fasta_out = (args.out.resolve() if args.out else fasta_in.with_name(fasta_in.stem + "_corrected.fasta"))

    # 1) Cargar mutaciones
    muts = load_mutations(tsv)

    # 2) Aplicar mutaciones (primer contig)
    length = apply_mutations_single_contig(fasta_in, fasta_out, muts)

    # 3) Renombrar cabecera
    name_hint = normalize_basename(fasta_in.stem)
    rename_header_inplace(fasta_out, name_hint=name_hint, species=args.species)
    print(f"[OK] FASTA corregido: {fasta_out} (len={length})")

    # 4) Copia oneline del corregido
    oneline_out = fasta_out.with_name(fasta_out.stem + "_oneline.fasta")
    write_oneline_copy(fasta_out, oneline_out)
    print(f"[OK] FASTA oneline: {oneline_out}")

    # 5) Alinear usando el umbral de confianza solicitado + SSC-fix opcional
    align_with_reference_and_number(
        oneline_out,
        min_identity=args.confianza,
        enable_ssc_fix=not args.no_ssc_fix
    )

if __name__ == "__main__":
    main()

