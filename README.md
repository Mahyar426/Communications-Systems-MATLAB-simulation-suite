# 📡 Communications Systems — MATLAB Simulation Suite

[![MATLAB](https://img.shields.io/badge/MATLAB-R2023b-orange?style=flat-square&logo=mathworks)](https://www.mathworks.com/)
[![Topics](https://img.shields.io/badge/topics-channel_coding%20%7C%20CDMA%20%7C%20OFDM%20%7C%20satellite-blue?style=flat-square)]()
[![Institution](https://img.shields.io/badge/Politecnico_di_Torino-M.Sc._Communications_Engineering-darkred?style=flat-square)]()

A hands-on simulation library covering the full digital communications stack — from binary channel coding and error correction to spread-spectrum CDMA systems and live LEO satellite link budgets. Every result is derived analytically **and** validated by Monte Carlo simulation.

---

## What's Inside

### 🔐 Module 1 — Channel Coding & Error Detection
**Files:** `A1/EX1.m`, `A1/EX2.m`, `A1/EX3.m`

Built from scratch: linear block codes characterised by their generator matrix **G** and parity-check matrix **H**. Derived the full distance spectrum, computed undetected-error probability analytically, then cross-validated against BSC Monte Carlo runs across six decade-spaced BEP values.

Key concepts implemented:
- Systematic (8,4) linear block codes — two competing code designs compared head-to-head
- Hamming weight distributions and minimum distance extraction
- CRC generator polynomial design and detection-capability evaluation
- Syndrome-based error detection over a Binary Symmetric Channel

---

### ⚡ Module 2 — Hard vs. Soft Decoding & Synchronisation
**Files:** `A2/EX1.m`, `A2/EX2.m`, `A2/EX3.m`

Side-by-side analytical and simulated comparison of **hard-decision** and **soft-decision** decoding. Derived the union bound and asymptotic bound on codeword error rate for soft decoding and showed the 2–3 dB gain over hard decoding empirically.

Key concepts implemented:
- Maximum Likelihood hard decoding with syndrome lookup table
- Soft-decision decoding: union bound and asymptotic approximation vs. simulation
- **Linear Feedback Shift Register (LFSR)** design — maximal-length m-sequences generation
- Synchronisation marker insertion and detection
- Autocorrelation and power spectral density via FFT

---

### 🌐 Module 3 — Spread Spectrum, CDMA & OFDM
**Files:** `A3/EX1.m`, `A3/EX2.m`, `A3/EX3.m`

Explored the world of spread-spectrum signalling — from sequence design theory to a full CDMA system simulation and adaptive OFDM bit allocation.

Key concepts implemented:
- **Gold sequence** generation (degree-7 polynomials, N=127 chips) — autocorrelation and cross-correlation characterised and plotted
- Welch and Sidelnikov bounds on sequence sets, evaluated across family sizes
- **CDMA** multi-user system: interference analysis and capacity limits
- **OFDM** adaptive bit loading via the **Hughes-Hartogs algorithm** — per-subcarrier SNR-aware allocation to maximise throughput over a frequency-selective channel

---

### 🛰️ Module 4 — LEO Satellite Link Budget
**File:** `A4/A4.m`

End-to-end satellite communication system modelled in MATLAB's Satellite Communications Toolbox. Simulated a sun-synchronous LEO orbit (500 km, 97.4° inclination) over a 7-day scenario with a 60-second propagation timestep.

Key concepts implemented:
- Orbital mechanics: satellite scenario creation with Keplerian parameters
- Four global ground stations modelled (Inuvik, Svalbard, Awarua, Troll)
- Gaussian dish antenna sizing (0.18 m satellite, 3.7 m ground) with aperture efficiency
- **Link budget**: EIRP, G/T, free-space path loss, Doppler, Eb/N₀ margin
- Access interval and visibility window analysis across the station network
- Channel degradation: atmospheric loss, multipath, and interference effects

---

## Tech Stack

| Tool | Role |
|---|---|
| **MATLAB** | Core simulation & visualisation engine |
| **Communications Toolbox** | BSC, AWGN channels, Gold sequence generator |
| **Satellite Communications Toolbox** | Orbital mechanics, link budget, ground station modelling |
| **Signal Processing Toolbox** | FFT-based correlation, PSD analysis |

---

## Skills Demonstrated

`Channel Coding` · `Error Correction` · `Monte Carlo Simulation` · `Spread Spectrum` · `CDMA` · `OFDM` · `LFSR / m-sequences` · `Gold Sequences` · `Satellite Link Budget` · `Orbital Mechanics` · `Analytical-Simulation Cross-Validation`

---

> Part of my M.Sc. in Communications Engineering at Politecnico di Torino.
> Connect: [LinkedIn](https://linkedin.com/in/mahyaronsori) · [mahyaronsori99@gmail.com](mailto:mahyaronsori99@gmail.com)
