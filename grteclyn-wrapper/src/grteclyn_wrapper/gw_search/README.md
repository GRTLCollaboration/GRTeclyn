# LIGO matched-filter search methodology (wormhole / exotic compact object)

This note describes how to search public LIGO data for a short, low-frequency burst consistent with an intermediate-mass Ellis–Bronnikov–type collapse, using numerical-relativity (NR) output as a template. It assumes use of data from the Gravitational Wave Open Science Center (GWOSC) and standard community tools (PyCBC, GWpy).

---

## 1. Scientific context

### 1.1 Mass–frequency scaling and detector sensitivity

Gravitational-wave frequencies scale inversely with source mass,

\[
f \propto \frac{1}{M}.
\]

For a quasi-normal ringdown frequency of order \(f_{\mathrm{code}} \approx 0.3\,M^{-1}\) in code units (geometrized time), the mapping to physical frequency fixes a narrow mass band if the signal is to fall in LIGO’s sensitive band.

| Approximate mass | Ringdown frequency (order of magnitude) | LIGO context |
|------------------|----------------------------------------|----------------|
| \(30\,M_\odot\)  | \(\sim 2000\,\mathrm{Hz}\)             | High frequency; often faint and dominated by shot noise. |
| \(3000\,M_\odot\) | \(\sim 20\,\mathrm{Hz}\)              | Near seismic wall; poor sensitivity. |
| \(500\)–\(2000\,M_\odot\) | \(\sim 30\)–\(120\,\mathrm{Hz}\) | Overlaps the bottom of LIGO’s most sensitive bucket. |

Thus a **targeted** search over intermediate masses is physically motivated and computationally lighter than a full binary black-hole template bank.

### 1.2 Why a targeted search is tractable

1. **Template-bank size**  
   Standard CBC searches use large template banks over mass and spin. Here the hypothesis is a **single** intermediate-mass object (e.g. \(\sim 500\)–\(2000\,M_\odot\)), so the bank can shrink from millions of waveforms to on the order of tens of templates.

2. **Connection to short, heavy events**  
   LIGO has observed short, low-frequency, high-mass–like transients (e.g. GW190521), which are debated in the literature (merger vs. exotic or non-standard scenarios). A collapse template that predicts a short burst in this band is a natural hypothesis to test **explicitly** against public data, including around such events.

3. **Instrumental discrimination**  
   Short broadband glitches (“blip glitches”) can resemble high-mass bursts in simple matched-filter outputs. Time–frequency representations (e.g. CWT spectrograms) help: blips often appear as vertical broadband features, whereas a QNM-like ringdown tends to show a more structured, frequency-resolved tail. Spectrogram-level comparison supports separating astrophysical candidates from noise.

### 1.3 Scope of this document

The remainder specifies **units**, **data handling**, **noise estimation**, and **matched filtering** with PyCBC, plus a minimal Python reference script. Extensions (template grids, multiple events, comparison to CBC SNR) are outlined at the end.

---

## 2. Methodology

### 2.1 Units and sampling

- NR time is in **geometrized** units (\(t/M\)). Physical time is obtained via  
  \[
  t_{\mathrm{phys}} = t/M \times (M/M_\odot) \times 4.925\times 10^{-6}\,\mathrm{s},
  \]
  using the standard solar-mass conversion to seconds.

- The waveform must be **interpolated** onto a uniform grid with sample rate matching the analysis data (commonly **4096 Hz** for many GWOSC releases).

### 2.2 Conditioning the template

- Apply a **taper** (e.g. Tukey / start–end taper) so the template goes smoothly to zero at the edges, reducing FFT spectral leakage.

### 2.3 Strain data and noise model

- Obtain calibrated strain \(h(t)\) for the chosen detector and segment (e.g. H1 or L1 from the GWOSC catalog API or files).

- Estimate the **power spectral density** \(\mathrm{PSD}(f)\) of the noise (e.g. from neighboring data), with interpolation and low-frequency handling consistent with the analysis (e.g. cutoff below \(\sim 15\)–\(20\,\mathrm{Hz}\) where ground motion dominates).

### 2.4 Matched filtering

- Construct the **complex SNR time series** by matched filtering the template against \(h(t)\), weighting by the inverse noise spectrum (whitening).

- Report **\(|\mathrm{SNR}(t)|\)**. A conventional single-detector threshold for claiming a “loud” peak is often quoted around \(\mathrm{SNR} \gtrsim 8\); rigorous inference requires coherent multi-detector statistics, background estimation, and FAR—this script is a **starting point**, not a full search pipeline.

### 2.5 Reference event

For a short, low-frequency burst, **GW190521** is a standard test case (heaviest, short, and extensively discussed). The methodology below downloads H1 strain for that event via PyCBC’s catalog interface as an illustration.

---

## 3. Software dependencies

```bash
pip install pycbc gwpy scipy numpy matplotlib
```

---

## 4. Reference implementation (Python)

The following script: (i) loads or mocks NR strain, (ii) scales to physical units and resamples to 4096 Hz, (iii) tapers, (iv) fetches GW190521 H1 data, (v) estimates a PSD, (vi) runs `matched_filter`, and (vii) plots \(|\mathrm{SNR}(t)|\) near the peak.

```python
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from pycbc.types import TimeSeries
from pycbc.filter import matched_filter
from pycbc.psd import interpolate, inverse_spectrum_truncation
from pycbc.catalog import Catalog
from pycbc.waveform.utils import taper_timeseries

# --- 1. Load and scale NR data (replace dummy arrays with your file) ---
nr_time_M = np.linspace(-25, 25, 5000)
nr_strain = np.sin(0.3 * nr_time_M) * np.exp(-(nr_time_M) ** 2 / 20)  # mock burst

M_wormhole = 1000.0  # solar masses
M_to_sec = 4.925e-6
time_sec = nr_time_M * M_wormhole * M_to_sec

sample_rate = 4096
dt_ligo = 1.0 / sample_rate
time_interp = np.arange(time_sec[0], time_sec[-1], dt_ligo)
interpolator = interp1d(
    time_sec, nr_strain, kind="cubic", fill_value=0.0, bounds_error=False
)
strain_interp = interpolator(time_interp)

template = TimeSeries(strain_interp, delta_t=dt_ligo)
template = taper_timeseries(template, tapermethod="startend")

# --- 2. LIGO strain (example: GW190521, H1) ---
print("Downloading LIGO data for GW190521...")
catalog_event = Catalog("GWTC-2").data["GW190521"]
data = catalog_event["H1"]

template.resize(len(data))
template = template.cyclic_time_shift(template.start_time)

# --- 3. PSD ---
print("Calculating detector noise (PSD)...")
psd = data.psd(4)
psd = interpolate(psd, data.delta_f)
psd = inverse_spectrum_truncation(
    psd, int(4 * data.sample_rate), low_frequency_cutoff=15.0
)

# --- 4. Matched filter ---
print("Running matched filter...")
snr = matched_filter(template, data, psd=psd, low_frequency_cutoff=20.0)
snr_timeseries = snr.abs()

peak_snr = snr_timeseries.max()
peak_time = snr_timeseries.sample_times[snr_timeseries.argmax()]
print(f"Maximum SNR: {peak_snr:.2f} at GPS time {peak_time:.2f}")

# --- 5. Plot ---
plt.figure(figsize=(10, 5))
plt.plot(snr_timeseries.sample_times, snr_timeseries, color="purple")
plt.axhline(8.0, color="red", linestyle="--", label="Reference threshold (SNR = 8)")
plt.ylabel("Signal-to-noise ratio (|SNR|)")
plt.xlabel("GPS time (s)")
plt.title(
    rf"Search for $1000\,M_\odot$ template in GW190521 (H1)\nPeak SNR = {peak_snr:.2f}"
)
plt.xlim(peak_time - 1.0, peak_time + 1.0)
plt.legend()
plt.grid(True, alpha=0.3)
plt.show()
```

---

## 5. Extensions toward a publication-quality search

- **Template bank:** Loop over a mass grid (e.g. \(M \in \{500, 600, \ldots, 2000\}\,M_\odot\)) and repeat Section 4 for each template.

- **Multiple events:** Apply the same pipeline to other short-duration candidates (e.g. GW190426, etc.), with consistent segment selection and PSD estimation.

- **Model comparison:** Compare peak \(|\mathrm{SNR}|\) (or coherent network SNR) for the exotic template to that of the best-fit CBC waveform where both are computed fairly (same data, PSD, and frequency cuts).

---

## 6. From \(\Psi_4\) to strain

LIGO templates use **strain** \(h(t)\). If the NR extraction gives the Weyl scalar \(\Psi_4\), then in the time domain \(\Psi_4 = \ddot{h}\) (up to radiative gauge conventions used in your code). In the frequency domain,

\[
\tilde{h}(f) = -\frac{\tilde{\Psi}_4(f)}{(2\pi f)^2}.
\]

Apply this in Fourier space, inverse-FFT to \(h(t)\), then use that as `nr_strain` after unit scaling and resampling as in Section 2.1.

---

## 7. Summary

Intermediate masses \(\sim 500\)–\(2000\,M_\odot\) map ringdown frequencies into LIGO’s sensitive band, enabling a **small** template family and desktop-scale matched filtering on open data. The main caveats are correct \(\Psi_4 \to h\) conversion, glitch rejection (e.g. spectrogram discrimination), and full multi-detector statistical validation beyond a single SNR plot.
