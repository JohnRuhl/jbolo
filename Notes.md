# Some notes about CMB-S4 models to fix
 - Found that CHLAT and SPLAT yamls in pbdr_v2 have different window thicknesses.
 Left them that way for the comparisons, but will correct it when done doing those checks.
 (one was 1/8", the other 3/16")
 - Found that SAT_LFMF optics default temp was 273K, where SPLAT is 250K.
 This affects window in SAT, should be 250.  Leaving it for purposes of comparison.
 - SAT windows have been at 30mm thick, should now be 20mm per PBDR.


# Checks to do
- Check that bolocalc photon noise formula is correct, ie factor of 2 in bose term?


# Features to add in funcs and physics:
 - Johnson noise
 - Loop gain (requires an alpha input)


Test done:
- Checked NETs and Popts for all CMB-S4 SAT + LAT bands.  agreement to within 2%.
- Gets same Poptical as bolo-calc to within 1% or so for SPT3g, and for SPLAT MF.
- Gets same photon NEP, both correlated and uncorrelated, as bolo-calc for SPLAT MF, to within a few percent.
  - this means the basis for the correlation factors is very close too, need to add in detector noise to get correlation factor (of total NET)
- Gets a reasonable NET_photon, but haven't compared yet.

- Phonon noise is right on bolo-calc's for SPLAT MF_1 and MF_2
- Got ~2% or better agreement for SPLAT MF_1 and MF_2 on:
 - NEP_phonon, NEP_photon (corr and uncorr)
 - NET_corr, NET_uncorr
 - P_optical
 - NET wafers

 - Make yamls and compare (Popt, NETs, NEPs, corr_factors) with bolo-calc for:
  - SPLAT [Done, agrees to better than a ~1 percent.]
  - CHLAT [Done, agrees to better than a ~1 percent.]
  - SAT_LFMF [Done, agrees to better than 1%]
  - SAT_HF [Done.]

 - [Done] Check that psat_factor (rather than fixed psats) works

*Values in tables:  jbolo/bolo-calc*

| SAT band |  Poptical  | NET_NC  | corr_factor |
|----------|------------|---------|-------------|
| LF_1   | 0.568 / 0.573|176.6 / 176.0 | 1.040 / 1.038 |
| LF_2   | 2.58 / 2.57 | 217.1 / 217.6 | 1.007 / 1.007 |
| MF_1_1 | 3.08 / 3.06 | 311.8 / 312.7 | 1.041 / 1.040 |
| MF_2_1 | 3.28 / 3.23 | 276.0 / 277.1 | 1.023 / 1.023 |
| MF_1_2 | 5.44 / 5.50 | 331.4 / 330.1 | 1.007 / 1.007 |
| MV_2_2 | 5.95 / 5.90 | 355.3 / 353.9 | 1.003 / 1.003 |
| HF_1   | 13.3 / 13.2 | 734 / 727     | 1.004 / 1.004 |
| HF_1   | 16.4 / 16.2 | 1770 / 1747   | 1.002 / 1.002 |

| SPLAT band |  Poptical  | NET_NC  | corr_factor |
|----------|------------|---------|-------------|
| ULF_1| 0.131 / 0.134 | 338.0 / 332.7 | 1.213 / 1.215 |
| LF_1 | 0.223 / 0.227 | 290.1 / 286.8 | 1.221 / 1.240 |
| LF_2 | 1.409 / 1.413 | 269.1 / 269.0 | 1.046 / 1.041 |
| MF_1 | 1.501 / 1.507 | 288.1 / 287.0 | 1.147 / 1.146 |
| MF_2 | 3.626 / 3.610 | 270.8 / 268.8 | 1.012 / 1.012 |
| HF_1 | 8.711 / 8.615 | 555.6 / 549.3 | 1.006 / 1.006 |
| HF_2 | 11.42 / 11.268| 1316. / 1296. | 1.005 / 1.005 |

| CHLAT band |  Poptical  | NET_NC  | corr_factor |
|----------|------------|---------|-------------|
| LF_1 | 0.267 / 0.271 | 317.5 / 314.2 | 1.262 / 1.281 |
| LF_2 | 1.368 / 1.372 | 255.6 / 255.6 | 1.045 / 1.039 |
| MF_1 | 1.549 / 1.559 | 292.1 / 291.5 | 1.151 / 1.149 |
| MF_2 | 4.701 / 4.731 | 327.3 / 327.5 | 1.015 / 1.015 |
| HF_1 | 12.45 / 12.517| 739.1 / 740.3 | 1.008 / 1.008 |
| HF_2 | 16.91 / 17.016| 1846. / 1850. | 1.007 / 1.007 |
