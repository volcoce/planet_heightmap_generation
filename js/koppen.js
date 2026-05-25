// Köppen climate classification using the "worldbuilding pasta" band-based
// methodology.  Two-season (summer/winter) data is used as a proxy for
// warmest/coldest month values.
//
// Approach:
//   Step 1 – Temperature bands  (tropical → temperate → continental → tundra → ice cap)
//   Step 2 – Arid zones (B)     dry in both seasons → desert core + steppe fringe
//   Step 3 – Precipitation subtypes within each band (A / C / D details)
//
// IMPORTANT: The simulation labels "summer" and "winter" are NH-centric
// (NH summer = June-Aug, NH winter = Dec-Feb).  For each cell we determine
// the LOCAL warm/cold season from temperature and use that to assign the
// correct precipitation pattern (s/w/f).  Without this, Mediterranean (Cs)
// and monsoon (Cw/Dw) climates are hemisphere-flipped.
//
// When windResult is supplied, three additional improvements activate:
//   1. Direct latitude from r_lat (no hemisphere-guessing from temperature)
//   2. windwardScore: dot(windDir, towardOcean) → relaxed 's' threshold on
//      westerly windward coasts so Mediterranean (Cs) types appear correctly
//   3. Continentality-aware aridity check: continental interiors receive a
//      small precipitation boost before the Pthresh test, mirroring the
//      FMG MIN_MOISTURE fix that unlocks Dfa/Dsa generation

import { smoothstep } from './wind.js';
import { elevToHeightKm } from './color-map.js';

/**
 * Köppen class definitions: ID → { code, name, color [r,g,b] 0-1 }.
 */
export const KOPPEN_CLASSES = [
    { code: 'Ocean',  name: 'Ocean',                              color: [0.29, 0.44, 0.65] },  // #4a6fa5
    { code: 'Af',     name: 'Tropical rainforest',                color: [0.00, 0.00, 1.00] },  // #0000FF
    { code: 'Am',     name: 'Tropical monsoon',                   color: [0.00, 0.47, 1.00] },  // #0077FF
    { code: 'Aw',     name: 'Tropical savanna',                   color: [0.27, 0.67, 0.98] },  // #46AAFA
    { code: 'BWh',    name: 'Hot desert',                         color: [1.00, 0.00, 0.00] },  // #FF0000
    { code: 'BWk',    name: 'Cold desert',                        color: [1.00, 0.59, 0.59] },  // #FF9696
    { code: 'BSh',    name: 'Hot steppe',                         color: [0.96, 0.65, 0.00] },  // #F5A500
    { code: 'BSk',    name: 'Cold steppe',                        color: [1.00, 0.86, 0.39] },  // #FFDB63
    { code: 'Cfa',    name: 'Humid subtropical',                  color: [0.78, 1.00, 0.31] },  // #C8FF50
    { code: 'Cfb',    name: 'Oceanic',                            color: [0.39, 1.00, 0.31] },  // #64FF50
    { code: 'Cfc',    name: 'Subpolar oceanic',                   color: [0.20, 0.78, 0.00] },  // #32C800
    { code: 'Csa',    name: 'Hot-summer Mediterranean',           color: [1.00, 1.00, 0.00] },  // #FFFF00
    { code: 'Csb',    name: 'Warm-summer Mediterranean',          color: [0.78, 0.78, 0.00] },  // #C8C800
    { code: 'Csc',    name: 'Cold-summer Mediterranean',          color: [0.59, 0.59, 0.00] },  // #969600
    { code: 'Cwa',    name: 'Humid subtropical (monsoon)',         color: [0.59, 1.00, 0.59] },  // #96FF96
    { code: 'Cwb',    name: 'Subtropical highland',               color: [0.39, 0.78, 0.39] },  // #63C764
    { code: 'Cwc',    name: 'Cold subtropical highland',          color: [0.20, 0.59, 0.20] },  // #329633
    { code: 'Dfa',    name: 'Hot-summer continental',             color: [0.00, 1.00, 1.00] },  // #00FFFF
    { code: 'Dfb',    name: 'Warm-summer continental',            color: [0.22, 0.78, 1.00] },  // #37C8FF
    { code: 'Dfc',    name: 'Subarctic',                          color: [0.00, 0.49, 0.49] },  // #007D7D
    { code: 'Dfd',    name: 'Extremely cold subarctic',           color: [0.00, 0.27, 0.37] },  // #00465F
    { code: 'Dsa',    name: 'Hot-summer continental (dry summer)', color: [0.90, 0.50, 1.00] },  // #E680FF
    { code: 'Dsb',    name: 'Warm-summer continental (dry summer)', color: [0.70, 0.35, 0.85] },  // #B359D9
    { code: 'Dsc',    name: 'Subarctic (dry summer)',              color: [0.50, 0.20, 0.65] },  // #8033A6
    { code: 'Dsd',    name: 'Extremely cold subarctic (dry summer)', color: [0.35, 0.10, 0.45] },  // #591A73
    { code: 'Dwa',    name: 'Hot-summer continental (monsoon)',    color: [0.67, 0.69, 1.00] },  // #ABB1FF
    { code: 'Dwb',    name: 'Warm-summer continental (monsoon)',   color: [0.43, 0.47, 0.78] },  // #6E77C8
    { code: 'Dwc',    name: 'Subarctic (monsoon)',                color: [0.29, 0.31, 0.78] },  // #4A50C8
    { code: 'Dwd',    name: 'Extremely cold subarctic (monsoon)', color: [0.20, 0.00, 0.53] },  // #320087
    { code: 'ET',     name: 'Tundra',                             color: [0.70, 0.70, 0.70] },  // #B2B2B2
    { code: 'EF',     name: 'Ice cap',                            color: [0.41, 0.41, 0.41] },  // #686868
];

// Lookup table: KOPPEN_CLASSES code → ID (built once at import time)
const CODE_TO_ID = {};
KOPPEN_CLASSES.forEach((c, i) => { CODE_TO_ID[c.code] = i; });

/**
 * Classify each region into a Köppen climate type.
 *
 * @param {object}       mesh         - SphereMesh
 * @param {Float32Array}  r_elevation  - per-region elevation (<=0 = ocean)
 * @param {object}        tempResult   - { r_temperature_summer, r_temperature_winter } (0-1 → -45..+45 C)
 * @param {object}        precipResult - { r_precip_summer, r_precip_winter } (0-1 p95-normalized)
 * @param {object|null}   windResult   - optional; unlocks windwardScore, direct latitude, continentality
 * @returns {Uint8Array}  r_koppen     - per-region class ID (index into KOPPEN_CLASSES)
 */
export function classifyKoppen(mesh, r_elevation, tempResult, precipResult, windResult = null) {
    const n = mesh.numRegions;
    const r_koppen = new Uint8Array(n);

    const tSummer = tempResult.r_temperature_summer;
    const tWinter = tempResult.r_temperature_winter;
    const pSummer = precipResult.r_precip_summer;
    const pWinter = precipResult.r_precip_winter;

    // Optional geography/wind arrays from windResult
    const r_lat           = windResult?.r_lat           ?? null;
    const r_lon           = windResult?.r_lon           ?? null;
    const r_continentality = windResult?.r_continentality ?? null;
    const r_coastDist     = windResult?.r_coastDistLand  ?? null;
    const wE_s            = windResult?.r_wind_east_summer  ?? null;
    const wN_s            = windResult?.r_wind_north_summer ?? null;
    const wE_w            = windResult?.r_wind_east_winter  ?? null;
    const wN_w            = windResult?.r_wind_north_winter ?? null;
    const { adjOffset, adjList } = mesh;

    // ── Precompute per-cell toward-ocean direction in lat/lon tangent space ──
    // Direction = toward the neighbor with the lowest r_coastDistLand (or toward
    // ocean neighbors directly for coastal cells).  Used for windwardScore.
    const r_coastDirE = new Float32Array(n);
    const r_coastDirN = new Float32Array(n);
    if (r_lat && r_lon && r_coastDist) {
        for (let r = 0; r < n; r++) {
            if (r_elevation[r] <= 0) continue;
            const lat0 = r_lat[r], lon0 = r_lon[r];
            const cosLat = Math.cos(lat0);
            const myDist = r_coastDist[r];
            let bE = 0, bN = 0, count = 0;
            const end = adjOffset[r + 1];
            for (let ni = adjOffset[r]; ni < end; ni++) {
                const nb = adjList[ni];
                const nbDist = r_elevation[nb] <= 0 ? -1 : r_coastDist[nb];
                // Include ocean neighbors (dist=-1) and any land neighbor closer to coast
                if (nbDist < myDist || (r_elevation[nb] <= 0)) {
                    let dLon = r_lon[nb] - lon0;
                    if (dLon >  Math.PI) dLon -= 2 * Math.PI;
                    if (dLon < -Math.PI) dLon += 2 * Math.PI;
                    const dLat = r_lat[nb] - lat0;
                    bE += dLon * cosLat;
                    bN += dLat;
                    count++;
                }
            }
            const len = Math.sqrt(bE * bE + bN * bN);
            if (len > 1e-10) { r_coastDirE[r] = bE / len; r_coastDirN[r] = bN / len; }
        }
    }

    // ── Per-cell classification ──
    for (let r = 0; r < n; r++) {
        // ── Ocean ──
        if (r_elevation[r] <= 0) {
            r_koppen[r] = 0;
            continue;
        }

        // ── Convert normalised values to physical units ──
        const Ts = -45 + Math.max(0, Math.min(1, tSummer[r])) * 90;
        const Tw = -45 + Math.max(0, Math.min(1, tWinter[r])) * 90;
        const Thot  = Math.max(Ts, Tw);
        const Tcold = Math.min(Ts, Tw);
        const Tann  = (Ts + Tw) / 2;

        // ── Highland swing reduction (port of FMG highland modifier) ──
        // temperature.js already applies a lapse-rate reduction to the mean
        // temperature, but does NOT reduce the seasonal swing.  In reality,
        // thin-atmosphere highlands experience proportionally less seasonal
        // swing: cold seasons aren't as deep and warm seasons cap out.
        // Apply a modest narrowing above ~1.5 km so that highland tropics
        // classify as Cwb/Cwc (subtropical highland) rather than Cwa.
        const elevKm = elevToHeightKm(r_elevation[r]);
        const highlandFactor = smoothstep(1.5, 4.5, elevKm);  // 0 at 1.5 km, 1 at 4.5 km
        const swingNarrow = highlandFactor * (Thot - Tcold) * 0.22;
        // Warm season cools slightly more; cold season warms slightly (thin atmosphere)
        const ThhotEff  = Thot  - swingNarrow * 0.4;
        const TcoldEff  = Tcold + swingNarrow * 0.6;
        const TannEff   = (ThhotEff + TcoldEff) / 2;
        const TshoulderEff = ThhotEff - (ThhotEff - TcoldEff) * (1.2 / 6);

        // Shoulder-month temperature proxy (2 months before peak summer)
        const Tshoulder = TshoulderEff;

        // ── Hemisphere-aware local seasons ──
        const localSummerIsSim = Ts >= Tw;

        // Latitude: use r_lat directly when available (more accurate than
        // inferring hemisphere from the temperature seasonal swing alone)
        const latRad = r_lat ? r_lat[r] : (localSummerIsSim ? 1 : -1) * Math.abs(Math.asin(Math.max(-1, Math.min(1, TannEff / 28))));
        const absLat = Math.abs(latRad) * (180 / Math.PI);

        // ── Precipitation in mm (p95-calibrated) ──
        const Ps = Math.max(0, pSummer[r]) * 1000;
        const Pw = Math.max(0, pWinter[r]) * 1000;
        const Pann = Ps + Pw;

        const PsummerLocal = localSummerIsSim ? Ps : Pw;
        const PwinterLocal = localSummerIsSim ? Pw : Ps;
        const PsMonthLocal = PsummerLocal / 6;
        const PwMonthLocal = PwinterLocal / 6;

        // Estimate driest individual month from 6-month averages.
        // At equal seasons (ratio=1) → driest ≈ 0.70× avg.
        // At strong monsoon (ratio≥4) → driest ≈ 0.35× avg.
        const seasonRatio = Math.max(PsMonthLocal, PwMonthLocal) / (Math.min(PsMonthLocal, PwMonthLocal) || 1);
        const driestFraction = 0.60 - 0.35 * smoothstep(1, 4, seasonRatio);
        const Pdry = Math.min(PsMonthLocal, PwMonthLocal) * driestFraction;

        // ── windwardScore (positive = ocean is upwind = windward coast) ──
        // Uses annual-mean wind direction dotted against the toward-ocean vector,
        // negated so that onshore wind gives a positive score (FMG convention).
        let windwardScore = 0;
        if (wE_s && wN_s && (r_coastDirE[r] !== 0 || r_coastDirN[r] !== 0)) {
            const wE = (wE_s[r] + (wE_w ? wE_w[r] : wE_s[r])) * 0.5;
            const wN = (wN_s[r] + (wN_w ? wN_w[r] : wN_s[r])) * 0.5;
            const wLen = Math.sqrt(wE * wE + wN * wN) || 1;
            // Negative of (windDir · towardOcean): positive when wind blows FROM ocean
            windwardScore = Math.max(-1, Math.min(1,
                -((wE / wLen) * r_coastDirE[r] + (wN / wLen) * r_coastDirN[r])
            ));
        }

        // Continentality 0=coast, 1=deep interior
        const contVal = r_continentality ? r_continentality[r] : 0;

        // ================================================================
        //  STEP 1 – TEMPERATURE BANDS
        // ================================================================

        let band;
        let tempSubBand = '';

        if (ThhotEff < 0) {
            band = 'EF';
        } else if (ThhotEff < 10) {
            band = 'ET';
        } else if (TcoldEff >= 18) {
            band = 'A';
        } else if (TcoldEff >= 0) {
            band = 'C';
            tempSubBand = ThhotEff >= 22 ? 'hotSummer' : 'coolSummer';
        } else {
            band = 'D';
            tempSubBand = Tshoulder >= 10 ? 'humidCont' : 'subarctic';
        }

        if (band === 'EF') { r_koppen[r] = CODE_TO_ID['EF']; continue; }
        if (band === 'ET') { r_koppen[r] = CODE_TO_ID['ET']; continue; }

        // ================================================================
        //  STEP 2 – ARID ZONES (B)
        // ================================================================

        const summerFrac = Pann > 0 ? PsummerLocal / Pann : 0.5;
        let Pthresh;
        if (summerFrac >= 0.7) {
            Pthresh = 20 * TannEff + 280;
        } else if (summerFrac <= 0.3) {
            Pthresh = 20 * TannEff;
        } else {
            Pthresh = 20 * TannEff + 140;
        }
        Pthresh = Math.max(0, Pthresh);

        // Continental interior moisture boost — analog of FMG MIN_MOISTURE fix.
        // Deep interiors tend to receive summer convective precipitation that the
        // BFS cap in precipitation.js underestimates; a small continental boost
        // prevents borderline-Dfa/Dsa cells from falling into BSk/BSh.
        const contBoost = 1.0 + contVal * 0.28;
        const PannAdj = Pann * contBoost;

        if (PannAdj < Pthresh) {
            const isHot = TannEff >= 18;
            if (PannAdj < Pthresh * 0.5) {
                r_koppen[r] = isHot ? CODE_TO_ID['BWh'] : CODE_TO_ID['BWk'];
            } else {
                r_koppen[r] = isHot ? CODE_TO_ID['BSh'] : CODE_TO_ID['BSk'];
            }
            continue;
        }

        // ================================================================
        //  STEP 3 – PRECIPITATION SUBTYPES WITHIN EACH BAND
        // ================================================================

        // ── Determine s / w / f precipitation pattern ──
        //
        // For 's' (dry local summer): relax the monthly threshold on windward
        // westerly coasts (28–52° latitude, contVal < 0.35, windwardScore > 0.3).
        // The subtropical high genuinely suppresses summer precipitation on these
        // coasts; the physics-based simulation encodes this but the threshold can
        // clip real Mediterranean patterns near the boundary.
        const localSummerDrier = PsummerLocal < PwinterLocal;
        let sSummerThresh = 50; // mm/month (standard Köppen ≈ 40; relaxed for 6-month averages)
        if (windwardScore > 0.3 && absLat > 28 && absLat < 52 && contVal < 0.35) {
            // Up to 80 mm/month for strongly windward coasts at peak Mediterranean latitudes
            sSummerThresh = 50 + windwardScore * 30 * smoothstep(28, 38, absLat) * smoothstep(52, 42, absLat);
        }

        let precipPattern;
        if (localSummerDrier && PsMonthLocal < sSummerThresh && PsMonthLocal < PwMonthLocal / 2) {
            precipPattern = 's';
        } else if (!localSummerDrier && PwMonthLocal < PsMonthLocal / 3) {
            precipPattern = 'w';
        } else {
            precipPattern = 'f';
        }

        // ── Temperature sub-letter (a / b / c / d) ──
        let tempLetter;
        if (ThhotEff >= 22) {
            tempLetter = 'a';
        } else if (Tshoulder >= 10) {
            tempLetter = 'b';
        } else if (TcoldEff >= -38) {
            tempLetter = 'c';
        } else {
            tempLetter = 'd';
        }

        // ── Band A: Tropical ──
        if (band === 'A') {
            if (Pdry >= 60) {
                r_koppen[r] = CODE_TO_ID['Af'];
            } else if (Pann >= 25 * (100 - Pdry)) {
                r_koppen[r] = CODE_TO_ID['Am'];
            } else {
                r_koppen[r] = CODE_TO_ID['Aw'];
            }
            continue;
        }

        // ── Band C: Temperate ──
        if (band === 'C') {
            const code = 'C' + precipPattern + tempLetter;
            const id = CODE_TO_ID[code];
            if (id !== undefined) {
                r_koppen[r] = id;
            } else {
                r_koppen[r] = CODE_TO_ID['Cfb'];
            }
            continue;
        }

        // ── Band D: Continental ──
        if (band === 'D') {
            const code = 'D' + precipPattern + tempLetter;
            const id = CODE_TO_ID[code];
            if (id !== undefined) {
                r_koppen[r] = id;
            } else {
                const fallback = 'Df' + tempLetter;
                r_koppen[r] = CODE_TO_ID[fallback] || CODE_TO_ID['Dfc'];
            }
            continue;
        }
    }

    return r_koppen;
}
