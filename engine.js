/*
 * Meta-Analysis Audit Tool — pure statistical core.
 *
 * Extracted VERBATIM from the inline <script> of index.html so the statistical
 * core is a single source of truth, importable under Node for testing and
 * loaded by the browser as a plain <script> (functions become globals).
 *
 * Pure functions only — no DOM, no globals beyond these declarations. The
 * DOM/audit/render layer (collectStudyData, runAudit, displayResults, the SVG
 * plotters, CSV/report I/O) stays inline in index.html and consumes these.
 *
 * Method: log-scale pooling of effect measures (OR/RR/HR pooled as log; RD/SMD
 * on natural scale). Random-effects via DerSimonian-Laird, REML, Paule-Mandel,
 * and Empirical Bayes tau^2 estimators; optional Hartung-Knapp-Sidik-Jonkman
 * (HKSJ) CI; Q-profile tau^2 CI; t_{k-2} prediction interval; Egger's
 * regression test; leave-one-out. Faithful to the shipped app.
 *
 * Correctness fix applied during the 2026-06 revival (see applyHKSJ):
 *   HKSJ variance FLOOR. The original computed seAdj = sqrt(Qstar/((k-1)*sumW))
 *   with no floor, so when Q < k-1 the Hartung-Knapp interval NARROWED below the
 *   ordinary random-effects interval (anti-conservative — the documented HKSJ
 *   failure mode). Now floored at the classical RE variance via the standard
 *   multiplier max(1, Qstar/(k-1)). t_{k-1} critical value was already correct.
 */

// ---- Standard normal CDF (Abramowitz & Stegun 7.1.26 erf approximation) ----
// Phi(0) = 0.5 exactly; |error| < 7.5e-8.
function pnorm(x) {
    const a1 = 0.254829592, a2 = -0.284496736, a3 = 1.421413741;
    const a4 = -1.453152027, a5 = 1.061405429, p = 0.3275911;
    const sign = x < 0 ? -1 : 1;
    x = Math.abs(x) / Math.sqrt(2);
    const t = 1 / (1 + p * x);
    const y = 1 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * Math.exp(-x * x);
    return 0.5 * (1 + sign * y);
}

// ---- Log-gamma (Lanczos) ----
function lgamma(x) {
    const cof = [76.18009172947146, -86.50532032941677, 24.01409824083091,
                -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5];
    let y = x, tmp = x + 5.5;
    tmp -= (x + 0.5) * Math.log(tmp);
    let ser = 1.000000000190015;
    for (let j = 0; j < 6; j++) ser += cof[j] / ++y;
    return -tmp + Math.log(2.5066282746310005 * ser / x);
}

// ---- Regularized lower incomplete gamma ----
function gammainc(a, x) {
    if (x < 0 || a <= 0) return 0;
    if (x === 0) return 0;

    const ITMAX = 200, EPS = 3e-7;
    const gln = lgamma(a);

    if (x < a + 1) {
        let ap = a, sum = 1 / a, del = sum;
        for (let n = 1; n <= ITMAX; n++) {
            ap++;
            del *= x / ap;
            sum += del;
            if (Math.abs(del) < Math.abs(sum) * EPS)
                return sum * Math.exp(-x + a * Math.log(x) - gln);
        }
    } else {
        let b = x + 1 - a, c = 1e30, d = 1 / b, h = d;
        for (let i = 1; i <= ITMAX; i++) {
            const an = -i * (i - a);
            b += 2;
            d = an * d + b;
            if (Math.abs(d) < 1e-30) d = 1e-30;
            c = b + an / c;
            if (Math.abs(c) < 1e-30) c = 1e-30;
            d = 1 / d;
            const del = d * c;
            h *= del;
            if (Math.abs(del - 1) < EPS)
                return 1 - Math.exp(-x + a * Math.log(x) - gln) * h;
        }
    }
    return 1;
}

// ---- Chi-square CDF ----
function pchisq(x, df) {
    if (x <= 0) return 0;
    return gammainc(df / 2, x / 2);
}

// ---- Chi-square quantile (inverse CDF) via Newton-Raphson ----
function qchisq(p, df) {
    if (p <= 0) return 0;
    if (p >= 1) return Infinity;

    let x = df;
    for (let i = 0; i < 50; i++) {
        const fx = pchisq(x, df) - p;
        if (Math.abs(fx) < 1e-10) break;
        const pdf = Math.exp((df / 2 - 1) * Math.log(x) - x / 2 - (df / 2) * Math.log(2) - lgamma(df / 2));
        if (pdf < 1e-15) break;
        x = Math.max(0.001, x - fx / pdf);
    }
    return x;
}

// ---- Student-t CDF (incomplete beta) ----
function pt(t, df) {
    const x = df / (df + t * t);
    const a = df / 2, b = 0.5;
    const bt = Math.exp(lgamma(a + b) - lgamma(a) - lgamma(b) + a * Math.log(x) + b * Math.log(1 - x));
    const betainc = x < (a + 1) / (a + b + 2)
        ? bt * betacf(x, a, b) / a
        : 1 - bt * betacf(1 - x, b, a) / b;
    return t < 0 ? 0.5 * betainc : 1 - 0.5 * betainc;
}

function betacf(x, a, b) {
    const MAXIT = 100, EPS = 3e-7;
    const qab = a + b, qap = a + 1, qam = a - 1;
    let c = 1, d = 1 - qab * x / qap;
    if (Math.abs(d) < 1e-30) d = 1e-30;
    d = 1 / d;
    let h = d;
    for (let m = 1; m <= MAXIT; m++) {
        const m2 = 2 * m;
        let aa = m * (b - m) * x / ((qam + m2) * (a + m2));
        d = 1 + aa * d; if (Math.abs(d) < 1e-30) d = 1e-30;
        c = 1 + aa / c; if (Math.abs(c) < 1e-30) c = 1e-30;
        d = 1 / d; h *= d * c;
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
        d = 1 + aa * d; if (Math.abs(d) < 1e-30) d = 1e-30;
        c = 1 + aa / c; if (Math.abs(c) < 1e-30) c = 1e-30;
        d = 1 / d;
        const del = d * c;
        h *= del;
        if (Math.abs(del - 1) < EPS) break;
    }
    return h;
}

// ---- Student-t quantile via Newton-Raphson ----
function qt(p, df) {
    if (p <= 0) return -Infinity;
    if (p >= 1) return Infinity;
    if (p === 0.5) return 0;

    let t = p > 0.5 ? 1 : -1;
    for (let i = 0; i < 50; i++) {
        const cdf = pt(t, df);
        const err = cdf - p;
        if (Math.abs(err) < 1e-10) break;
        const pdf = Math.exp(lgamma((df + 1) / 2) - lgamma(df / 2)) / (Math.sqrt(df * Math.PI) * Math.pow(1 + t * t / df, (df + 1) / 2));
        t -= err / pdf;
    }
    return t;
}

// ============================================
// META-ANALYSIS ENGINE
// studies: [{ study, effect, se, ... }] with effect/se already on the analysis
// (log) scale. Pure — no DOM.
// ============================================

// Fixed effect meta-analysis
function metaFE(studies) {
    const k = studies.length;
    if (k === 0) return null;

    const w = studies.map(s => 1 / (s.se * s.se));
    const sumW = w.reduce((a, b) => a + b, 0);
    const effect = studies.reduce((sum, s, i) => sum + w[i] * s.effect, 0) / sumW;
    const se = Math.sqrt(1 / sumW);

    // Q statistic
    const Q = studies.reduce((sum, s, i) => sum + w[i] * Math.pow(s.effect - effect, 2), 0);
    const df = k - 1;
    const pHet = 1 - pchisq(Q, df);
    const I2 = Math.max(0, (Q - df) / Q * 100);

    const weights = w.map(wi => wi / sumW * 100);

    return {
        effect, se,
        ciLower: effect - 1.96 * se,
        ciUpper: effect + 1.96 * se,
        Q, df, pHet, I2, tau2: 0,
        weights, k,
        method: 'FE'
    };
}

// DerSimonian-Laird
function metaDL(studies) {
    const fe = metaFE(studies);
    if (!fe) return null;

    const k = studies.length;
    const w = studies.map(s => 1 / (s.se * s.se));
    const sumW = w.reduce((a, b) => a + b, 0);
    const C = sumW - w.reduce((sum, wi) => sum + wi * wi, 0) / sumW;

    const tau2 = Math.max(0, (fe.Q - fe.df) / C);

    const wRE = studies.map(s => 1 / (s.se * s.se + tau2));
    const sumWRE = wRE.reduce((a, b) => a + b, 0);
    const effect = studies.reduce((sum, s, i) => sum + wRE[i] * s.effect, 0) / sumWRE;
    const se = Math.sqrt(1 / sumWRE);

    return {
        effect, se,
        ciLower: effect - 1.96 * se,
        ciUpper: effect + 1.96 * se,
        Q: fe.Q, df: fe.df, pHet: fe.pHet,
        I2: Math.max(0, (fe.Q - fe.df) / fe.Q * 100),
        tau2,
        weights: wRE.map(wi => wi / sumWRE * 100),
        k,
        method: 'DL'
    };
}

// REML estimator (iterative)
function metaREML(studies, maxIter = 100, tol = 1e-8) {
    const fe = metaFE(studies);
    if (!fe) return null;

    const k = studies.length;
    const y = studies.map(s => s.effect);
    const v = studies.map(s => s.se * s.se);

    let tau2 = metaDL(studies).tau2;

    for (let iter = 0; iter < maxIter; iter++) {
        const w = v.map(vi => 1 / (vi + tau2));
        const sumW = w.reduce((a, b) => a + b, 0);
        const mu = y.reduce((sum, yi, i) => sum + w[i] * yi, 0) / sumW;

        const dl = -0.5 * w.reduce((sum, wi, i) => sum + wi - wi * wi * Math.pow(y[i] - mu, 2), 0);
        const d2l = 0.5 * w.reduce((sum, wi) => sum + wi * wi, 0);

        const delta = dl / d2l;
        const tau2New = Math.max(0, tau2 + delta);

        if (Math.abs(tau2New - tau2) < tol) break;
        tau2 = tau2New;
    }

    const w = v.map(vi => 1 / (vi + tau2));
    const sumW = w.reduce((a, b) => a + b, 0);
    const effect = y.reduce((sum, yi, i) => sum + w[i] * yi, 0) / sumW;
    const se = Math.sqrt(1 / sumW);

    return {
        effect, se,
        ciLower: effect - 1.96 * se,
        ciUpper: effect + 1.96 * se,
        Q: fe.Q, df: fe.df, pHet: fe.pHet,
        I2: tau2 > 0 ? Math.max(0, (fe.Q - fe.df) / fe.Q * 100) : 0,
        tau2,
        weights: w.map(wi => wi / sumW * 100),
        k,
        method: 'REML'
    };
}

// Paule-Mandel estimator
function metaPM(studies, maxIter = 100, tol = 1e-8) {
    const fe = metaFE(studies);
    if (!fe) return null;

    const k = studies.length;
    const y = studies.map(s => s.effect);
    const v = studies.map(s => s.se * s.se);
    const df = k - 1;

    let tau2 = metaDL(studies).tau2;

    for (let iter = 0; iter < maxIter; iter++) {
        const w = v.map(vi => 1 / (vi + tau2));
        const sumW = w.reduce((a, b) => a + b, 0);
        const mu = y.reduce((sum, yi, i) => sum + w[i] * yi, 0) / sumW;
        const Q = y.reduce((sum, yi, i) => sum + w[i] * Math.pow(yi - mu, 2), 0);

        if (Math.abs(Q - df) < tol) break;

        const dQ = -w.reduce((sum, wi, i) => sum + wi * wi * Math.pow(y[i] - mu, 2), 0);
        if (Math.abs(dQ) < 1e-15) break;

        const tau2New = Math.max(0, tau2 - (Q - df) / dQ);
        if (Math.abs(tau2New - tau2) < tol) break;
        tau2 = tau2New;
    }

    const w = v.map(vi => 1 / (vi + tau2));
    const sumW = w.reduce((a, b) => a + b, 0);
    const effect = y.reduce((sum, yi, i) => sum + w[i] * yi, 0) / sumW;
    const se = Math.sqrt(1 / sumW);

    return {
        effect, se,
        ciLower: effect - 1.96 * se,
        ciUpper: effect + 1.96 * se,
        Q: fe.Q, df: fe.df, pHet: fe.pHet,
        I2: tau2 > 0 ? Math.max(0, (fe.Q - fe.df) / fe.Q * 100) : 0,
        tau2,
        weights: w.map(wi => wi / sumW * 100),
        k,
        method: 'PM'
    };
}

// Empirical Bayes estimator
function metaEB(studies) {
    const fe = metaFE(studies);
    if (!fe) return null;

    const k = studies.length;
    const y = studies.map(s => s.effect);
    const v = studies.map(s => s.se * s.se);

    const w = v.map(vi => 1 / vi);
    const sumW = w.reduce((a, b) => a + b, 0);

    const tau2 = Math.max(0, (fe.Q - fe.df) / sumW);

    const wRE = v.map(vi => 1 / (vi + tau2));
    const sumWRE = wRE.reduce((a, b) => a + b, 0);
    const effect = y.reduce((sum, yi, i) => sum + wRE[i] * yi, 0) / sumWRE;
    const se = Math.sqrt(1 / sumWRE);

    return {
        effect, se,
        ciLower: effect - 1.96 * se,
        ciUpper: effect + 1.96 * se,
        Q: fe.Q, df: fe.df, pHet: fe.pHet,
        I2: tau2 > 0 ? Math.max(0, (fe.Q - fe.df) / fe.Q * 100) : 0,
        tau2,
        weights: wRE.map(wi => wi / sumWRE * 100),
        k,
        method: 'EB'
    };
}

// Apply HKSJ adjustment.
// FIX (2026-06 revival): floor the Hartung-Knapp multiplier at 1 so the
// adjusted SE never narrows below the classical RE SE when Q < k-1.
function applyHKSJ(meta, studies) {
    if (!meta || meta.k < 3) return meta;

    const k = meta.k;
    const y = studies.map(s => s.effect);
    const w = studies.map(s => 1 / (s.se * s.se + meta.tau2));
    const sumW = w.reduce((a, b) => a + b, 0);

    // Q* statistic for HKSJ
    const Qstar = y.reduce((sum, yi, i) => sum + w[i] * Math.pow(yi - meta.effect, 2), 0);
    // Hartung-Knapp variance = max(1, Q*/(k-1)) * (1/sumW); floor keeps it
    // from shrinking below the random-effects variance.
    const seAdj = Math.sqrt(Math.max(1, Qstar / (k - 1)) / sumW);

    const tCrit = qt(0.975, k - 1);

    return {
        ...meta,
        se: seAdj,
        ciLower: meta.effect - tCrit * seAdj,
        ciUpper: meta.effect + tCrit * seAdj,
        hksj: true
    };
}

// Prediction interval (t_{k-2} convention)
function calcPredictionInterval(meta, studies) {
    if (!meta || meta.k < 3) return { piLower: NaN, piUpper: NaN };

    const tCrit = qt(0.975, meta.k - 2);
    const piSE = Math.sqrt(meta.se * meta.se + meta.tau2);

    return {
        piLower: meta.effect - tCrit * piSE,
        piUpper: meta.effect + tCrit * piSE
    };
}

// Q-profile CI for tau^2
function tau2CI(studies, tau2Point) {
    const k = studies.length;
    if (k < 2) return { lower: NaN, upper: NaN };

    const df = k - 1;
    const y = studies.map(s => s.effect);
    const v = studies.map(s => s.se * s.se);

    function Qfunc(tau2) {
        const w = v.map(vi => 1 / (vi + tau2));
        const sumW = w.reduce((a, b) => a + b, 0);
        const mu = y.reduce((sum, yi, i) => sum + w[i] * yi, 0) / sumW;
        return y.reduce((sum, yi, i) => sum + w[i] * Math.pow(yi - mu, 2), 0);
    }

    const qLower = qchisq(0.975, df);
    const qUpper = qchisq(0.025, df);

    function solveTau2(targetQ) {
        let lo = 0, hi = 10;
        while (Qfunc(hi) > targetQ && hi < 1000) hi *= 2;

        for (let i = 0; i < 50; i++) {
            const mid = (lo + hi) / 2;
            const Q = Qfunc(mid);
            if (Math.abs(Q - targetQ) < 0.001) return mid;
            if (Q > targetQ) lo = mid;
            else hi = mid;
        }
        return (lo + hi) / 2;
    }

    const tau2Upper = solveTau2(qUpper);

    const Q0 = Qfunc(0);
    const tau2Lower = Q0 >= qLower ? solveTau2(qLower) : 0;

    return { lower: tau2Lower, upper: tau2Upper };
}

// Egger's regression test for funnel-plot asymmetry
function eggersTest(studies) {
    if (studies.length < 3) return null;

    const n = studies.length;
    const x = studies.map(s => 1 / s.se);          // precision
    const y = studies.map(s => s.effect / s.se);   // standardized effect

    const sumX = x.reduce((a, b) => a + b, 0);
    const sumY = y.reduce((a, b) => a + b, 0);
    const sumXY = x.reduce((sum, xi, i) => sum + xi * y[i], 0);
    const sumX2 = x.reduce((sum, xi) => sum + xi * xi, 0);

    const slope = (n * sumXY - sumX * sumY) / (n * sumX2 - sumX * sumX);
    const intercept = (sumY - slope * sumX) / n;

    const yHat = x.map(xi => intercept + slope * xi);
    const residSS = y.reduce((sum, yi, i) => sum + Math.pow(yi - yHat[i], 2), 0);
    const mse = residSS / (n - 2);
    const seIntercept = Math.sqrt(mse * sumX2 / (n * sumX2 - sumX * sumX));

    const t = intercept / seIntercept;
    const pValue = 2 * (1 - pt(Math.abs(t), n - 2));

    return { intercept, seIntercept, t, pValue };
}

// Leave-one-out analysis
function leaveOneOut(studies, estimator = 'DL') {
    const metaFunc = { DL: metaDL, REML: metaREML, PM: metaPM, EB: metaEB, FE: metaFE }[estimator] || metaDL;

    return studies.map((_, i) => {
        const subset = studies.filter((_, j) => j !== i);
        const result = metaFunc(subset);
        return {
            excluded: studies[i].study,
            ...result
        };
    });
}

if (typeof module !== 'undefined' && module.exports) {
    module.exports = {
        pnorm, lgamma, gammainc, pchisq, qchisq, pt, betacf, qt,
        metaFE, metaDL, metaREML, metaPM, metaEB,
        applyHKSJ, calcPredictionInterval, tau2CI, eggersTest, leaveOneOut
    };
}
