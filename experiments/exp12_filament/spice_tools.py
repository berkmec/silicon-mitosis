"""ngspice (WSL, Debian ngspice 44.2) helpers for exp12 with the official SKY130 nfet_01v8 tt card.
4 K is NOT modelled by BSIM4; it is emulated at 27 C with explicit, documented knobs:
  fii   : multiplies every impact-ionisation alpha term (alpha0/alpha1 and their l/w/p variants) -> I_sub x fii
  dvt   : threshold increase at 4 K, applied as a gate-voltage offset V_g,eff = V_g - dvt
  body  : floating (frozen) body = large R_sub to the substrate tap; body-source junction replaced by an explicit diode
          (BSIM4 junction areas set to zero) whose turn-on voltage V_on is a knob (~0.6 V at 300 K, ~1.0-1.1 V at 4 K)."""
import re
import subprocess
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
SP = HERE / "spice"
SP.mkdir(exist_ok=True)
CARDS = Path(__file__).resolve().parents[2] / "data" / "sky130_models"


def wsl(p: Path) -> str:
    s = str(p.resolve()).replace("\\", "/")
    return "/mnt/" + s[0].lower() + s[2:]


def make_card():
    """Copy of the official card with every alpha* term multiplied by {fii} and the junction saturation currents
    (jss, jsws; drain-side and gate-edge values default to these in BSIM4) multiplied by {fjun}. fjun = 1e-30 removes
    BSIM4's own 300 K body diodes so that the explicit 4 K diode (V_on knob) is the only body discharge path."""
    out = SP / "nfet_01v8_tt_fii.spice"
    if out.exists():
        return out
    txt = (CARDS / "nfet_01v8_tt.pm3.spice").read_text(encoding="utf-8")
    pat = re.compile(r"\b([lwp]?alpha[01])\s*=\s*([-+0-9.eE]+)")
    new, n = pat.subn(lambda m: f"{m.group(1)} = {{{m.group(2)}*fii}}", txt)
    assert n > 100, n
    pj = re.compile(r"\b(jss|jsws|jswgs|jsd|jswd|jswgd)\s*=\s*([-+0-9.eE]+)")
    new, nj = pj.subn(lambda m: f"{m.group(1)} = {{{m.group(2)}*fjun}}", new)
    assert nj > 100, nj
    out.write_text(MODIFIED_NOTICE + ".param fii = 1\n.param fjun = 1\n" + new, encoding="utf-8", newline="\n")
    return out


# Apache-2.0 section 4(b): a modified copy must say that it was changed. The original header below it is kept.
MODIFIED_NOTICE = (
    "* NOTICE: modified copy written by experiments/exp12_filament/spice_tools.py (Silicon Mitosis project, counfield,\n"
    "* https://counfield.com/, 2026). Changes: every alpha* term is multiplied by {fii}; the junction saturation currents\n"
    "* (jss, jsws, jswgs, jsd, jswd, jswgd) are multiplied by {fjun}. Original: sky130_fd_pr__nfet_01v8__tt.pm3.spice\n"
    "* (Copyright 2020 The SkyWater PDK Authors, Apache-2.0; see data/sky130_models/LICENSE-Apache-2.0.txt).\n")


def header(fii=1.0, fjun=1.0):
    card = make_card()
    return (f".option scale=1e-6\n.temp 27\n.include {wsl(CARDS / 'nfet_mismatch.spice')}\n"
            f".include {wsl(card)}\n.param fii = {fii}\n.param fjun = {fjun}\n")


def diode_model(von, n=1.0, i_on=1e-6, w_um=1.0):
    """Diode that carries i_on (A per um of width) at V = von, ideality n at 27 C: IS = i_on * W * exp(-von/(n kT/q))."""
    vt = 0.025852 * n
    return f".model dbody{von:g}_{w_um:g} d is={i_on * w_um * np.exp(-von / vt):.4e} n={n}\n"


def run(name: str, text: str, timeout=1200):
    cir = SP / f"{name}.cir"
    cir.write_text(text, encoding="utf-8", newline="\n")
    r = subprocess.run(["wsl", "-d", "CB-Local", "--", "sh", "-c", f"cd {wsl(SP)} && ngspice -b {name}.cir"],
                       capture_output=True, text=True, timeout=timeout)
    log = (r.stdout or "") + (r.stderr or "")
    (SP / f"{name}.log").write_text(log, encoding="utf-8")
    bad = [ln for ln in log.splitlines() if re.search(r"(?i)error|singular|timestep too small|no convergence", ln)]
    return log, bad


def load(name: str):
    return np.loadtxt(SP / f"{name}.dat")
