/*
 * Node tests for the Meta-Analysis Audit Tool engine. Run: node tests.js
 *
 * Every expected value is hand-derived independently below (NOT a re-run of the
 * engine). Studies carry effect/se already on the analysis (log) scale, matching
 * how collectStudyData feeds the engine in the browser.
 */
const E = require('./engine.js');

let pass = 0, fail = 0;
function ok(name, cond, detail) {
    if (cond) { pass++; console.log('  ok  ' + name); }
    else { fail++; console.log(' FAIL ' + name + (detail !== undefined ? '  -> ' + detail : '')); }
}
function close(a, b, tol) { return Math.abs(a - b) < (tol || 1e-3); }

// ============================================================
// 1. Distribution helpers
// ============================================================
// pnorm(0) = 0.5 EXACTLY (a value of 0 would be the classic erf/CDF bug).
ok('pnorm(0) == 0.5', close(E.pnorm(0), 0.5, 1e-9), E.pnorm(0));
ok('pnorm(1.96) ~ 0.975', close(E.pnorm(1.96), 0.975, 2e-4), E.pnorm(1.96));
ok('pnorm(-1.96) ~ 0.025', close(E.pnorm(-1.96), 0.025, 2e-4), E.pnorm(-1.96));
ok('pnorm symmetric: Phi(x)+Phi(-x)=1', close(E.pnorm(1.3) + E.pnorm(-1.3), 1, 1e-9));

// qt(0.975, df): t critical values vs known table (R qt to 3 dp).
ok('qt(0.975, 1) ~ 12.706', close(E.qt(0.975, 1), 12.706, 0.02), E.qt(0.975, 1));
ok('qt(0.975, 2) ~ 4.303', close(E.qt(0.975, 2), 4.303, 0.01), E.qt(0.975, 2));
ok('qt(0.975, 4) ~ 2.776', close(E.qt(0.975, 4), 2.776, 0.01), E.qt(0.975, 4));
ok('qt(0.975, 30) ~ 2.042', close(E.qt(0.975, 30), 2.042, 0.01), E.qt(0.975, 30));
ok('qt(0.5, 5) == 0', E.qt(0.5, 5) === 0);

// pchisq: chi-square CDF at known points. Median of chi2_2 is 2*ln(2)=1.3863.
ok('pchisq(0,*)=0', E.pchisq(0, 3) === 0);
ok('pchisq(1.3863, 2) ~ 0.5 (chi2_2 median)', close(E.pchisq(2 * Math.log(2), 2), 0.5, 2e-3), E.pchisq(2 * Math.log(2), 2));
// chi2_1 at 3.8415 (=1.96^2) has CDF 0.95.
ok('pchisq(3.8415, 1) ~ 0.95', close(E.pchisq(3.841459, 1), 0.95, 2e-3), E.pchisq(3.841459, 1));

// ============================================================
// 2. Hand-worked TWO-study DerSimonian-Laird example
// ------------------------------------------------------------
// Studies on the LOG scale (e.g. log-OR):
//   s1: y1 = ln(0.70) = -0.3566749,  se1 = 0.10  -> v1 = 0.01, w1 = 100
//   s2: y2 = ln(0.90) = -0.1053605,  se2 = 0.20  -> v2 = 0.04, w2 = 25
//
// Fixed effect:  sumW = 125
//   sumWy = 100*(-0.3566749) + 25*(-0.1053605) = -38.301505
//   thetaFE = -38.301505 / 125 = -0.30641204
//
// Q = w1*(y1-thetaFE)^2 + w2*(y2-thetaFE)^2
//   (y1-thetaFE) = -0.05026286 -> 100*0.00252636 = 0.252636
//   (y2-thetaFE) =  0.20105154 -> 25 *0.04042172 = 1.010543
//   Q = 1.263179 ; df = 1
//
// C = sumW - sumW2/sumW = 125 - (100^2 + 25^2)/125 = 125 - 10625/125 = 40
// tau2 = (Q - df)/C = 0.263179 / 40 = 0.00657947
// I^2  = (Q - df)/Q * 100 = 0.263179/1.263179 * 100 = 20.8349 %
//
// Random-effects:
//   wRE1 = 1/(0.01+0.00657947) = 60.31572
//   wRE2 = 1/(0.04+0.00657947) = 21.46876
//   sumWRE = 81.78448
//   thetaRE = (60.31572*-0.3566749 + 21.46876*-0.1053605)/81.78448 = -0.2907036
//   seRE = sqrt(1/81.78448) = 0.1105766
//   pooled OR (back-transformed) = exp(-0.2907036) = 0.747754
// ============================================================
const s2 = [
    { study: 'A', effect: Math.log(0.70), se: 0.10 },
    { study: 'B', effect: Math.log(0.90), se: 0.20 }
];
const fe = E.metaFE(s2);
ok('2study FE: thetaFE ~ -0.306412', close(fe.effect, -0.30641204, 1e-5), fe.effect);
ok('2study FE: Q ~ 1.263179', close(fe.Q, 1.263179, 1e-4), fe.Q);
ok('2study FE: df == 1', fe.df === 1);
ok('2study FE: tau2 == 0 (fixed)', fe.tau2 === 0);
ok('2study FE: weights sum to 100', close(fe.weights.reduce((a, b) => a + b, 0), 100, 1e-9));

const dl = E.metaDL(s2);
ok('2study DL: tau2 ~ 0.00657947', close(dl.tau2, 0.00657947, 1e-6), dl.tau2);
ok('2study DL: I2 ~ 20.8349%', close(dl.I2, 20.8349, 0.01), dl.I2);
ok('2study DL: thetaRE ~ -0.2907036', close(dl.effect, -0.2907036, 1e-5), dl.effect);
ok('2study DL: seRE ~ 0.1105766', close(dl.se, 0.1105766, 1e-5), dl.se);
ok('2study DL: pooled OR exp(theta) ~ 0.747754', close(Math.exp(dl.effect), 0.747754, 1e-4), Math.exp(dl.effect));
ok('2study DL: CI uses 1.96 (Wald)', close(dl.ciLower, dl.effect - 1.96 * dl.se, 1e-12));

// ============================================================
// 3. Edge cases
// ============================================================
// k=1 must NOT produce NaN for the FE pooled effect/SE (a single study just
// passes through). DL tau2/C is degenerate at k=1; we only require FE be finite.
const one = E.metaFE([{ study: 'solo', effect: Math.log(0.8), se: 0.15 }]);
ok('k=1 FE: effect == ln(0.8), not NaN', close(one.effect, Math.log(0.8), 1e-12) && !isNaN(one.effect), one.effect);
ok('k=1 FE: se == 0.15, not NaN', close(one.se, 0.15, 1e-12) && !isNaN(one.se), one.se);
ok('k=1 FE: Q == 0 (single study)', one.Q === 0, one.Q);
// I^2 = (Q-df)/Q = 0/0 is mathematically undefined at k=1 (df=0). The engine
// faithfully returns NaN; I^2 is simply not defined for a single study. The
// contract "k=1 must not NaN" applies to the pooled effect/SE, which hold above.
ok('k=1 FE: I2 is NaN (undefined for k=1)', isNaN(one.I2), one.I2);

// empty guard
ok('empty -> null (FE)', E.metaFE([]) === null);
ok('empty -> null (DL)', E.metaDL([]) === null);

// two IDENTICAL studies: no heterogeneity => Q=0, tau2=0, I2=0, and the pooled
// SE = sqrt(1/(w+w)) = sqrt(1/(2/v)) = se/sqrt(2).
const idv = 0.12;
const id = E.metaDL([
    { study: 'X', effect: Math.log(0.85), se: idv },
    { study: 'X', effect: Math.log(0.85), se: idv }
]);
ok('identical: Q == 0', close(id.Q, 0, 1e-9), id.Q);
ok('identical: tau2 == 0', id.tau2 === 0, id.tau2);
ok('identical: I2 == 0', id.I2 === 0, id.I2);
ok('identical: effect == ln(0.85)', close(id.effect, Math.log(0.85), 1e-9), id.effect);
ok('identical: se == se/sqrt(2)', close(id.se, idv / Math.sqrt(2), 1e-9), id.se);

// ============================================================
// 4. HKSJ floor property (the revival fix)
// ------------------------------------------------------------
// When studies are very homogeneous, Q* < k-1, so the UNFLOORED Hartung-Knapp
// SE would shrink below the classical RE SE. With the floor, the HKSJ SE must be
// >= the random-effects SE. Build a near-homogeneous 3-study set.
// ============================================================
const homo = [
    { study: 'h1', effect: Math.log(0.80), se: 0.20 },
    { study: 'h2', effect: Math.log(0.81), se: 0.20 },
    { study: 'h3', effect: Math.log(0.80), se: 0.20 }
];
const dlHomo = E.metaDL(homo);
const hk = E.applyHKSJ(dlHomo, homo);
ok('HKSJ floor: seHKSJ >= seRE for homogeneous data', hk.se >= dlHomo.se - 1e-12, 'seHKSJ=' + hk.se + ' seRE=' + dlHomo.se);
ok('HKSJ uses t_{k-1}: CI half-width = qt(0.975,k-1)*se',
    close(hk.ciUpper - hk.effect, E.qt(0.975, hk.k - 1) * hk.se, 1e-9));
ok('HKSJ no-op for k<3', E.applyHKSJ(E.metaDL(s2), s2).hksj === undefined);

// HKSJ should WIDEN (not narrow) when there IS heterogeneity (Q*>k-1). Use
// genuinely heterogeneous effects with VARYING precision (k=4).
const hetSet = [
    { study: 'g1', effect: Math.log(0.45), se: 0.12 },
    { study: 'g2', effect: Math.log(1.40), se: 0.18 },
    { study: 'g3', effect: Math.log(0.65), se: 0.10 },
    { study: 'g4', effect: Math.log(1.05), se: 0.22 }
];
const dlHet = E.metaDL(hetSet);
const hkHet = E.applyHKSJ(dlHet, hetSet);
ok('HKSJ widens for heterogeneous data', hkHet.se > dlHet.se, 'seHKSJ=' + hkHet.se + ' seRE=' + dlHet.se);

// ============================================================
// 5. Egger's test property (publication-bias diagnostic)
// ------------------------------------------------------------
// Undefined / null for k<3.  For a PERFECTLY symmetric set (effects exactly
// proportional to their SE about a common mean with no small-study slope-vs-
// intercept artefact) the intercept p-value should be sane (in [0,1]).
// ============================================================
ok('Egger null for k<3 (k=2 -> null)', E.eggersTest(s2) === null);
ok('Egger null for k=1 -> null', E.eggersTest([s2[0]]) === null);
const eg = E.eggersTest(hetSet);
ok('Egger k=3: returns object', eg !== null && typeof eg.pValue === 'number');
ok('Egger p-value in [0,1]', eg.pValue >= 0 && eg.pValue <= 1, eg.pValue);

// Symmetric funnel: effects centred symmetrically about 0 with VARYING SE
// (Egger's regression on precision is degenerate if every SE is identical, so
// the middle study gets a distinct SE). The standardized-effect regression has
// intercept ~ 0 (no small-study asymmetry) and p = 1.
const sym = [
    { study: 's1', effect: -0.20, se: 0.10 },
    { study: 's2', effect: 0.00, se: 0.15 },
    { study: 's3', effect: 0.20, se: 0.10 }
];
const egSym = E.eggersTest(sym);
ok('Egger symmetric set: intercept ~ 0', Math.abs(egSym.intercept) < 1e-6, egSym.intercept);
ok('Egger symmetric set: p ~ 1 (no asymmetry)', close(egSym.pValue, 1, 1e-6), egSym.pValue);

// ============================================================
// 6. Leave-one-out: returns k results, each over k-1 studies
// ============================================================
const loo = E.leaveOneOut(hetSet, 'DL');
ok('LOO: returns k results', loo.length === hetSet.length);
ok('LOO: each result pools k-1 studies', loo.every(l => l.k === hetSet.length - 1));
ok('LOO: excludes the named study', loo[0].excluded === 'g1' && loo[2].excluded === 'g3');

// ============================================================
// 7. Prediction interval: undefined (NaN) for k<3; uses t_{k-2}; wider than CI
// ============================================================
ok('PI: NaN for k<3', isNaN(E.calcPredictionInterval(E.metaDL(s2), s2).piLower));
const piHet = E.calcPredictionInterval(dlHet, hetSet);
ok('PI: defined for k>=3', !isNaN(piHet.piLower) && !isNaN(piHet.piUpper));
ok('PI: wider than the RE 95% CI', (piHet.piUpper - piHet.piLower) > (dlHet.ciUpper - dlHet.ciLower),
    'PIw=' + (piHet.piUpper - piHet.piLower) + ' CIw=' + (dlHet.ciUpper - dlHet.ciLower));
ok('PI: half-width = qt(0.975,k-2)*sqrt(se^2+tau2)',
    close(dlHet.effect - piHet.piLower, E.qt(0.975, dlHet.k - 2) * Math.sqrt(dlHet.se ** 2 + dlHet.tau2), 1e-9));

// ============================================================
// 8. tau2 Q-profile CI: lower <= point <= upper, lower >= 0
// ============================================================
const t2ci = E.tau2CI(hetSet, dlHet.tau2);
ok('tau2 CI: lower >= 0', t2ci.lower >= 0, t2ci.lower);
ok('tau2 CI: upper >= lower', t2ci.upper >= t2ci.lower, JSON.stringify(t2ci));
ok('tau2 CI: point estimate within [lower, upper]',
    dlHet.tau2 >= t2ci.lower - 1e-6 && dlHet.tau2 <= t2ci.upper + 1e-6,
    'tau2=' + dlHet.tau2 + ' CI=' + JSON.stringify(t2ci));

// ============================================================
// 9. Estimator agreement: all RE tau2 estimators >= 0 and reduce to the same
//    pooled effect direction; on homogeneous data all tau2 ~ 0.
// ============================================================
ok('REML tau2 >= 0', E.metaREML(hetSet).tau2 >= 0);
ok('PM tau2 >= 0', E.metaPM(hetSet).tau2 >= 0);
ok('EB tau2 >= 0', E.metaEB(hetSet).tau2 >= 0);
ok('homogeneous: DL tau2 ~ 0', close(E.metaDL(homo).tau2, 0, 1e-6), E.metaDL(homo).tau2);

console.log('\n' + pass + ' passed, ' + fail + ' failed');
process.exit(fail === 0 ? 0 : 1);
