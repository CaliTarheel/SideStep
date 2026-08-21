//! Spherical geometry helpers, Fibonacci sampling, smooth noise, spatial hash.
use rand::Rng;
use std::collections::HashMap;
use std::f64::consts::PI;

pub const R_KM: f64 = 6371.0;
pub type V3 = [f64; 3];

#[inline] pub fn dot(a: V3, b: V3) -> f64 { a[0] * b[0] + a[1] * b[1] + a[2] * b[2] }
#[inline] pub fn cross(a: V3, b: V3) -> V3 {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}
#[inline] pub fn add(a: V3, b: V3) -> V3 { [a[0] + b[0], a[1] + b[1], a[2] + b[2]] }
#[inline] pub fn sub(a: V3, b: V3) -> V3 { [a[0] - b[0], a[1] - b[1], a[2] - b[2]] }
#[inline] pub fn scale(a: V3, s: f64) -> V3 { [a[0] * s, a[1] * s, a[2] * s] }
#[inline] pub fn norm(a: V3) -> f64 { dot(a, a).sqrt() }
#[inline] pub fn normalize(a: V3) -> V3 {
    let n = norm(a);
    if n > 1e-300 { scale(a, 1.0 / n) } else { [0.0, 1.0, 0.0] }
}
/// Chord distance between unit vectors (≈ angular distance in radians for small angles).
#[inline] pub fn dist(a: V3, b: V3) -> f64 { norm(sub(a, b)) }
#[inline] pub fn angle(a: V3, b: V3) -> f64 { dot(a, b).clamp(-1.0, 1.0).acos() }

/// Rotate unit vector `p` about angular-velocity vector `omega` (rad/Myr) for `dt` Myr.
pub fn rotate(p: V3, omega: V3, dt: f64) -> V3 {
    let w = norm(omega);
    let th = w * dt;
    if th < 1e-15 { return p; }
    let k = scale(omega, 1.0 / w);
    let (s, c) = th.sin_cos();
    let t1 = scale(p, c);
    let t2 = scale(cross(k, p), s);
    let t3 = scale(k, dot(k, p) * (1.0 - c));
    normalize(add(add(t1, t2), t3))
}

/// Unit tangent vector at `p` pointing toward `q`.
pub fn tangent_toward(p: V3, q: V3) -> V3 {
    normalize(sub(q, scale(p, dot(p, q))))
}

/// Surface velocity (km/Myr) of a point `p` on a plate rotating with `omega` (rad/Myr).
#[inline] pub fn surface_velocity(omega: V3, p: V3) -> V3 { scale(cross(omega, p), R_KM) }

/// Some unit tangent vector at `p` (deterministic, for degenerate cases).
pub fn any_tangent(p: V3) -> V3 {
    let trial = if p[0].abs() < 0.9 { [1.0, 0.0, 0.0] } else { [0.0, 1.0, 0.0] };
    normalize(cross(p, trial))
}

pub fn fibonacci_sphere(n: usize) -> Vec<V3> {
    let golden = PI * (3.0 - 5f64.sqrt());
    (0..n)
        .map(|i| {
            let y = 1.0 - 2.0 * (i as f64 + 0.5) / n as f64;
            let r = (1.0 - y * y).max(0.0).sqrt();
            let th = golden * i as f64;
            [r * th.cos(), y, r * th.sin()]
        })
        .collect()
}

/// y is the polar axis. Returns (lat, lon) in radians.
#[inline] pub fn to_latlon(p: V3) -> (f64, f64) { (p[1].clamp(-1.0, 1.0).asin(), p[2].atan2(p[0])) }
#[inline] pub fn from_latlon(lat: f64, lon: f64) -> V3 {
    let c = lat.cos();
    [c * lon.cos(), lat.sin(), c * lon.sin()]
}

pub fn random_unit<R: Rng>(rng: &mut R) -> V3 {
    loop {
        let v = [rng.gen::<f64>() * 2.0 - 1.0, rng.gen::<f64>() * 2.0 - 1.0, rng.gen::<f64>() * 2.0 - 1.0];
        let n = norm(v);
        if n > 1e-3 && n <= 1.0 { return scale(v, 1.0 / n); }
    }
}

/// Smooth band-limited noise on the sphere: a sum of random plane waves. Output roughly in [-1, 1].
pub struct Noise { terms: Vec<(V3, f64, f64, f64)>, norm: f64 }
impl Noise {
    pub fn new<R: Rng>(rng: &mut R, n_terms: usize, k_min: f64, k_max: f64) -> Self {
        let mut terms = Vec::with_capacity(n_terms);
        let mut norm = 0.0;
        for _ in 0..n_terms {
            let u = random_unit(rng);
            let k = rng.gen_range(k_min..k_max);
            let ph = rng.gen_range(0.0..2.0 * PI);
            let amp = 1.0 / k.sqrt();
            norm += amp * amp;
            terms.push((u, k, ph, amp));
        }
        Noise { terms, norm: (norm / 2.0).sqrt().max(1e-9) }
    }
    pub fn eval(&self, p: V3) -> f64 {
        let mut s = 0.0;
        for &(u, k, ph, amp) in &self.terms { s += amp * (k * dot(u, p) + ph).sin(); }
        (s / self.norm).clamp(-1.5, 1.5)
    }
}

/// Uniform-grid spatial hash over the unit cube, for neighbour queries on the sphere.
pub struct SpatialHash { cell: f64, dim: i32, map: HashMap<i32, Vec<u32>> }
impl SpatialHash {
    pub fn new(cell: f64) -> Self {
        let dim = (2.0 / cell).ceil() as i32 + 2;
        SpatialHash { cell, dim, map: HashMap::new() }
    }
    #[inline] fn coord(&self, v: f64) -> i32 { ((v + 1.0) / self.cell).floor() as i32 + 1 }
    #[inline] fn key(&self, ix: i32, iy: i32, iz: i32) -> i32 { ix + self.dim * (iy + self.dim * iz) }
    pub fn build<I: Iterator<Item = (u32, V3)>>(&mut self, pts: I) {
        self.map.clear();
        for (i, p) in pts {
            let k = self.key(self.coord(p[0]), self.coord(p[1]), self.coord(p[2]));
            self.map.entry(k).or_default().push(i);
        }
    }
    /// Visit every stored index whose cell intersects the ball of radius `r` around `p`.
    /// Callers must check the distance themselves.
    #[inline]
    pub fn query<F: FnMut(u32)>(&self, p: V3, r: f64, mut f: F) {
        let lo = [self.coord(p[0] - r), self.coord(p[1] - r), self.coord(p[2] - r)];
        let hi = [self.coord(p[0] + r), self.coord(p[1] + r), self.coord(p[2] + r)];
        for iz in lo[2]..=hi[2] {
            for iy in lo[1]..=hi[1] {
                for ix in lo[0]..=hi[0] {
                    if let Some(v) = self.map.get(&self.key(ix, iy, iz)) {
                        for &i in v { f(i); }
                    }
                }
            }
        }
    }
}

pub type M3 = [[f64; 3]; 3];

/// Rotation matrix for angular velocity `omega` applied for `dt` (Rodrigues).
pub fn rot_matrix(omega: V3, dt: f64) -> M3 {
    let w = norm(omega);
    let th = w * dt;
    if th < 1e-15 { return [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]; }
    let k = scale(omega, 1.0 / w);
    let (s, c) = th.sin_cos();
    let t = 1.0 - c;
    [
        [c + k[0] * k[0] * t, k[0] * k[1] * t - k[2] * s, k[0] * k[2] * t + k[1] * s],
        [k[1] * k[0] * t + k[2] * s, c + k[1] * k[1] * t, k[1] * k[2] * t - k[0] * s],
        [k[2] * k[0] * t - k[1] * s, k[2] * k[1] * t + k[0] * s, c + k[2] * k[2] * t],
    ]
}
pub fn mat_mul(a: M3, b: M3) -> M3 {
    let mut o = [[0.0; 3]; 3];
    for i in 0..3 { for j in 0..3 { for k in 0..3 { o[i][j] += a[i][k] * b[k][j]; } } }
    o
}
pub fn mat_t(a: M3) -> M3 { [[a[0][0], a[1][0], a[2][0]], [a[0][1], a[1][1], a[2][1]], [a[0][2], a[1][2], a[2][2]]] }
pub fn mat_apply(a: M3, v: V3) -> V3 { [dot(a[0], v), dot(a[1], v), dot(a[2], v)] }

/// Axis-angle of a rotation matrix as (axis unit vector, angle in radians).
pub fn mat_to_axis_angle(m: M3) -> (V3, f64) {
    let tr = (m[0][0] + m[1][1] + m[2][2] - 1.0) / 2.0;
    let th = tr.clamp(-1.0, 1.0).acos();
    if th < 1e-9 { return ([0.0, 1.0, 0.0], 0.0); }
    let axis = [m[2][1] - m[1][2], m[0][2] - m[2][0], m[1][0] - m[0][1]];
    let n = norm(axis);
    if n < 1e-12 {
        // 180-degree case: take the largest diagonal
        let i = if m[0][0] >= m[1][1] && m[0][0] >= m[2][2] { 0 } else if m[1][1] >= m[2][2] { 1 } else { 2 };
        let mut a = [0.0; 3];
        a[i] = ((m[i][i] + 1.0) / 2.0).max(0.0).sqrt();
        for j in 0..3 { if j != i { a[j] = m[i][j] / (2.0 * a[i]); } }
        return (normalize(a), th);
    }
    (scale(axis, 1.0 / n), th)
}
