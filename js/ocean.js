// Ocean current simulation: wind-stress driven with Ekman-layer Coriolis deflection.
// Surface drift = actual wind stress rotated 45° rightward (NH) / leftward (SH),
// fading to wind-following at the equator where f → 0.
// Western boundary intensification is approximated by redirecting westward-blocked
// Ekman flow poleward along continental margins (Gulf Stream / Kuroshio mechanism).
// Warmth is derived from the poleward component of the resulting current rather
// than from hardcoded coast-type rules.

console.log('[ocean.js] Module loaded');
import { percentile } from './climate-util.js';

const DEG = Math.PI / 180;

// ── Coast distance & classification via BFS ─────────────────────────────────

function computeCoastFields(mesh, r_xyz, r_isOcean,
    r_eastX, r_eastY, r_eastZ) {
    const { adjOffset, adjList, numRegions } = mesh;

    const westSeeds = [];
    const eastSeeds = [];
    const allCoastSeeds = [];

    for (let r = 0; r < numRegions; r++) {
        if (!r_isOcean[r]) continue;

        let landDirX = 0, landDirY = 0, landDirZ = 0;
        let hasLandNeighbor = false;

        const end = adjOffset[r + 1];
        for (let ni = adjOffset[r]; ni < end; ni++) {
            const nb = adjList[ni];
            if (!r_isOcean[nb]) {
                hasLandNeighbor = true;
                landDirX += r_xyz[3 * nb] - r_xyz[3 * r];
                landDirY += r_xyz[3 * nb + 1] - r_xyz[3 * r + 1];
                landDirZ += r_xyz[3 * nb + 2] - r_xyz[3 * r + 2];
            }
        }

        if (!hasLandNeighbor) continue;

        allCoastSeeds.push(r);

        // Project land direction into tangent frame east component
        const normalE = landDirX * r_eastX[r] + landDirY * r_eastY[r] + landDirZ * r_eastZ[r];

        // normalE < -0.2 → land is to the west → western coast seed
        // normalE > +0.2 → land is to the east → eastern coast seed
        if (normalE < -0.2) {
            westSeeds.push(r);
        } else if (normalE > 0.2) {
            eastSeeds.push(r);
        } else {
            if (normalE <= 0) westSeeds.push(r);
            else eastSeeds.push(r);
        }
    }

    const bfsQueue = new Int32Array(numRegions);

    function bfsDistance(seeds) {
        const dist = new Int32Array(numRegions);
        dist.fill(-1);
        let qLen = 0;
        for (const s of seeds) {
            dist[s] = 0;
            bfsQueue[qLen++] = s;
        }
        let head = 0;
        while (head < qLen) {
            const r = bfsQueue[head++];
            const d = dist[r] + 1;
            const end = adjOffset[r + 1];
            for (let ni = adjOffset[r]; ni < end; ni++) {
                const nb = adjList[ni];
                if (r_isOcean[nb] && dist[nb] === -1) {
                    dist[nb] = d;
                    bfsQueue[qLen++] = nb;
                }
            }
        }
        return dist;
    }

    const r_coastDist     = bfsDistance(allCoastSeeds);
    const r_westCoastDist = bfsDistance(westSeeds);
    const r_eastCoastDist = bfsDistance(eastSeeds);

    return { r_coastDist, r_westCoastDist, r_eastCoastDist };
}

// ── Circumpolar channel detection ───────────────────────────────────────────

function hasCircumpolarChannel(r_lat, r_lon, r_isOcean, numRegions, targetLat, bandWidth) {
    const NUM_BINS = 72;
    const binHasOcean = new Uint8Array(NUM_BINS);
    const latMin = targetLat - bandWidth;
    const latMax = targetLat + bandWidth;

    for (let r = 0; r < numRegions; r++) {
        if (!r_isOcean[r]) continue;
        const lat = r_lat[r];
        if (lat < latMin || lat > latMax) continue;

        let bin = Math.floor(((r_lon[r] + Math.PI) / (2 * Math.PI)) * NUM_BINS);
        bin = ((bin % NUM_BINS) + NUM_BINS) % NUM_BINS;
        binHasOcean[bin] = 1;
    }

    for (let i = 0; i < NUM_BINS; i++) {
        if (!binHasOcean[i]) return false;
    }
    return true;
}

// ── Laplacian smoothing (ocean only) ────────────────────────────────────────

function smoothOcean(mesh, field, r_isOcean, passes) {
    const { adjOffset, adjList, numRegions } = mesh;
    const tmp = new Float32Array(numRegions);

    for (let pass = 0; pass < passes; pass++) {
        for (let r = 0; r < numRegions; r++) {
            if (!r_isOcean[r]) { tmp[r] = field[r]; continue; }

            let sum = field[r], count = 1;
            const end = adjOffset[r + 1];
            for (let ni = adjOffset[r]; ni < end; ni++) {
                const nb = adjList[ni];
                if (r_isOcean[nb]) {
                    sum += field[nb];
                    count++;
                }
            }
            tmp[r] = sum / count;
        }
        field.set(tmp);
    }
}

// ── Main entry point ────────────────────────────────────────────────────────

/**
 * Compute ocean surface currents from wind stress with Ekman-layer Coriolis deflection.
 *
 * Physical model (simplified for concept art):
 *  1. Wind stress → Ekman surface drift: rotate wind vector 45° rightward (NH) /
 *     leftward (SH).  Coriolis fades to zero at the equator (|lat| < 5°).
 *  2. Western boundary intensification: westward Ekman drift blocked by a continental
 *     margin is redirected poleward, approximating the geostrophic western boundary
 *     current (Gulf Stream / Kuroshio / Brazil / Agulhas).
 *  3. Eastern boundary cold current: eastward drift blocked by eastern coasts deflects
 *     equatorward (California / Humboldt / Benguela / Canary).
 *  4. Warmth from flow direction: poleward current = warm (advecting equatorial water),
 *     equatorward = cold (advecting polar water).
 *
 * @param {SphereMesh} mesh
 * @param {Float32Array} r_xyz
 * @param {Float32Array} r_elevation
 * @param {object} windResult - output of computeWind()
 * @returns {object} r_ocean_current_east/north, r_ocean_speed, r_ocean_warmth (summer + winter)
 */
export function computeOceanCurrents(mesh, r_xyz, r_elevation, windResult) {
    console.log('[ocean.js] computeOceanCurrents called, numRegions:', mesh.numRegions);
    const numRegions = mesh.numRegions;
    const avgEdgeKm = (Math.PI * 6371) / Math.sqrt(numRegions);
    const timing = [];

    const { r_lat, r_sinLat, r_isLand,
        r_eastX, r_eastY, r_eastZ } = windResult;

    // Ocean mask
    const r_isOcean = new Uint8Array(numRegions);
    for (let r = 0; r < numRegions; r++) r_isOcean[r] = r_isLand[r] ? 0 : 1;

    // r_lon (may already be in windResult)
    let t0 = performance.now();
    let r_lon = windResult.r_lon;
    if (!r_lon) {
        r_lon = new Float32Array(numRegions);
        for (let r = 0; r < numRegions; r++) {
            r_lon[r] = Math.atan2(r_xyz[3 * r], r_xyz[3 * r + 2]);
        }
    }
    timing.push({ stage: 'Ocean: setup', ms: performance.now() - t0 });

    // Step 1: Coast BFS (shared between seasons)
    t0 = performance.now();
    const { r_coastDist, r_westCoastDist, r_eastCoastDist } =
        computeCoastFields(mesh, r_xyz, r_isOcean, r_eastX, r_eastY, r_eastZ);
    timing.push({ stage: 'Ocean: coast BFS', ms: performance.now() - t0 });

    // Step 2: Circumpolar channel detection
    t0 = performance.now();
    const circumpolarNH = hasCircumpolarChannel(r_lat, r_lon, r_isOcean, numRegions,  60 * DEG, 5 * DEG);
    const circumpolarSH = hasCircumpolarChannel(r_lat, r_lon, r_isOcean, numRegions, -60 * DEG, 5 * DEG);
    console.log(`[ocean.js] Circumpolar: NH=${circumpolarNH}, SH=${circumpolarSH}`);
    timing.push({ stage: 'Ocean: circumpolar', ms: performance.now() - t0 });

    // Boundary intensification radius (~750 km)
    const coastThreshold = Math.max(5, Math.round(750 / avgEdgeKm));

    // Per-cell Coriolis factor: 0 at equator, 1 at |lat| ≥ 5°.
    // Prevents the 1/f divergence at the equator.
    const sin5 = Math.sin(5 * DEG);
    const r_coriolisFactor = new Float32Array(numRegions);
    for (let r = 0; r < numRegions; r++) {
        r_coriolisFactor[r] = Math.min(1, Math.abs(r_sinLat[r]) / sin5);
    }

    const result = {};

    for (const season of ['summer', 'winter']) {
        const r_windE = windResult[`r_wind_east_${season}`];
        const r_windN = windResult[`r_wind_north_${season}`];

        t0 = performance.now();
        const currentE = new Float32Array(numRegions);
        const currentN = new Float32Array(numRegions);

        // Step 3: Ekman surface drift
        // Rotate wind vector by -(π/4)·sign(lat) (counterclockwise convention):
        //   NH: −45° (clockwise) — rightward of wind direction
        //   SH: +45° (counterclockwise) — leftward of wind direction
        // At the equator the Coriolis factor fades to 0 → current follows wind directly.
        //
        // Rotation of (wE, wN) by angle α:
        //   E' = wE·cos α − wN·sin α
        //   N' = wE·sin α + wN·cos α
        for (let r = 0; r < numRegions; r++) {
            if (!r_isOcean[r]) continue;

            const sinLat = r_sinLat[r];
            const cf     = r_coriolisFactor[r];
            const alpha  = -(Math.PI / 4) * (sinLat >= 0 ? 1 : -1) * cf;
            const cosA   = Math.cos(alpha);
            const sinA   = Math.sin(alpha);
            const wE     = r_windE[r];
            const wN     = r_windN[r];

            currentE[r] = wE * cosA - wN * sinA;
            currentN[r] = wE * sinA + wN * cosA;
        }

        // Step 4: Boundary intensification
        // Western boundary (Gulf Stream / Kuroshio mechanism):
        //   Trade winds drive westward Ekman drift → water piles against western coast →
        //   geostrophic pressure gradient forces poleward jet (×2 intensification).
        // Eastern boundary (California / Humboldt mechanism):
        //   Westerlies drive eastward Ekman drift → blocked by eastern coast →
        //   redirected equatorward as cold upwelling current.
        for (let r = 0; r < numRegions; r++) {
            if (!r_isOcean[r]) continue;

            const sinLat  = r_sinLat[r];
            const poleSign = sinLat >= 0 ? 1 : -1;

            const wDist = r_westCoastDist[r];
            if (wDist >= 0 && wDist < coastThreshold) {
                const prox = 1 - wDist / coastThreshold;
                const prox2 = prox * prox;
                // Westward Ekman component blocked by the coast
                const blockedWest = Math.max(0, -currentE[r]);
                currentN[r] += poleSign * blockedWest * prox2 * 2.0; // western intensification
                currentE[r] *= 1 - prox2 * 0.7;
            }

            const eDist = r_eastCoastDist[r];
            if (eDist >= 0 && eDist < coastThreshold) {
                const prox = 1 - eDist / coastThreshold;
                const prox2 = prox * prox;
                // Eastward Ekman component blocked by the coast
                const blockedEast = Math.max(0, currentE[r]);
                currentN[r] -= poleSign * blockedEast * prox2 * 0.8; // equatorward cold current
                currentE[r] *= 1 - prox2 * 0.5;
            }

            // Circumpolar boost: open Southern/Arctic Ocean channels sustain
            // a strong eastward circumpolar current driven by unimpeded westerlies.
            const isCircumpolar = (sinLat > 0 && circumpolarNH) || (sinLat < 0 && circumpolarSH);
            if (isCircumpolar) {
                const absLatDeg = Math.abs(r_lat[r]) / DEG;
                if (absLatDeg >= 55 && absLatDeg <= 75) {
                    const cStr = 1 - Math.abs(absLatDeg - 65) / 10;
                    currentE[r] = currentE[r] * (1 - cStr) + 1.5 * cStr;
                    currentN[r] *= 1 - cStr * 0.8;
                }
            }
        }

        timing.push({ stage: `Ocean: Ekman + boundary (${season})`, ms: performance.now() - t0 });

        // Step 5: Smooth ~125 km (scale-invariant)
        t0 = performance.now();
        const smoothPasses = Math.max(2, Math.round(125 / avgEdgeKm));
        smoothOcean(mesh, currentE, r_isOcean, smoothPasses);
        smoothOcean(mesh, currentN, r_isOcean, smoothPasses);
        for (let r = 0; r < numRegions; r++) {
            if (!r_isOcean[r]) { currentE[r] = 0; currentN[r] = 0; }
        }
        timing.push({ stage: `Ocean: smooth current (${season})`, ms: performance.now() - t0 });

        // Step 6: Warmth from poleward component of current
        // A current flowing toward the pole carries warm equatorial water (warm = +1).
        // A current flowing toward the equator carries cold polar water (cold = −1).
        t0 = performance.now();
        const r_warmth = new Float32Array(numRegions);
        for (let r = 0; r < numRegions; r++) {
            if (!r_isOcean[r]) continue;
            const sinLat   = r_sinLat[r];
            const poleSign = sinLat >= 0 ? 1 : -1;
            const poleward = currentN[r] * poleSign;
            const spd      = Math.sqrt(currentE[r] * currentE[r] + currentN[r] * currentN[r]);
            r_warmth[r]    = spd > 0.01 ? poleward / spd : 0;
        }
        // Smooth heavily — warmth signal should blur across basin scales (~900 km)
        const warmthPasses = Math.max(3, Math.round(900 / avgEdgeKm));
        smoothOcean(mesh, r_warmth, r_isOcean, warmthPasses);
        timing.push({ stage: `Ocean: warmth (${season})`, ms: performance.now() - t0 });

        // Step 7: Normalize speed to [0, 1] via 95th percentile
        t0 = performance.now();
        const r_speed    = new Float32Array(numRegions);
        const speedsSq   = new Float32Array(numRegions);
        let   oceanCount = 0;
        for (let r = 0; r < numRegions; r++) {
            const spdSq = currentE[r] * currentE[r] + currentN[r] * currentN[r];
            r_speed[r] = spdSq;
            if (r_isOcean[r] && spdSq > 0) speedsSq[oceanCount++] = spdSq;
        }
        const p95Sq   = percentile(speedsSq.subarray(0, oceanCount), 0.95);
        const invP95Sq = 1 / p95Sq;
        for (let r = 0; r < numRegions; r++) {
            r_speed[r] = Math.min(1, Math.sqrt(r_speed[r] * invP95Sq));
        }
        timing.push({ stage: `Ocean: normalize speed (${season})`, ms: performance.now() - t0 });

        console.log(`[Ocean ${season}] coastThreshold=${coastThreshold}, p95Sq=${p95Sq?.toExponential(3)}, oceanCells=${oceanCount}`);

        result[`r_ocean_current_east_${season}`]  = currentE;
        result[`r_ocean_current_north_${season}`] = currentN;
        result[`r_ocean_speed_${season}`]         = r_speed;
        result[`r_ocean_warmth_${season}`]        = r_warmth;
    }

    result._oceanTiming = timing;
    return result;
}
