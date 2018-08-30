
use std::f64::consts::PI;

/// convert croodinate from wgs84 to mercator projection
/// - https://en.wikipedia.org/wiki/Mercator_projection
/// # Examples
///
/// ```
/// use mercator::lnglat_to_mercator;
/// let k0:f64 = 0.9999;
/// let dx:f64 = 250000.0;
/// let lng:f64 = 120.982025;
/// let lat:f64 = 23.973875;
/// let (x, y) = lnglat_to_mercator(lng, lat, 121.0, k0, dx);
/// ```
#[allow(non_snake_case)]
pub fn lnglat_to_mercator(lng:f64, lat:f64, center_lng:f64, k0:f64, dx:f64) ->(f64, f64) {
    let a:f64 = 6378137.0;
    let b:f64 = 6356752.314245;
    let lng0:f64 = center_lng * PI / 180.0;
    let lng = (lng/180.0) * PI;
    let lat = (lat/180.0) * PI;
    
    //---------------------------------------------------------
    let e:f64 = (1.0 - b.powf(2.0) / a.powf(2.0)).powf(0.5);
    let e2:f64 = e.powf(2.0)/(1.0 - e.powf(2.0)); 
    let n:f64 = ( a - b ) / ( a + b );
    let nu:f64 = a / ((1.0 - e.powf(2.0) * lat.sin().powf(2.0))).powf(0.5);
    let p:f64 = lng - lng0;
    let A:f64 = a * (1.0 - n + (5.0/4.0) * (n.powf(2.0) - n.powf(3.0)) + (81.0/64.0) * (n.powf(4.0)  - n.powf(5.0)));
    let B:f64 = (3.0 * a * n/2.0) * (1.0 - n + (7.0/8.0)*(n.powf(2.0) - n.powf(3.0)) + (55.0/64.0)*(n.powf(4.0) - n.powf(5.0)));
    let C:f64 = (15.0 * a * (n.powf(2.0))/16.0)*(1.0 - n + (3.0/4.0)*(n.powf(2.0) - n.powf(3.0)));
    let D:f64 = (35.0 * a * (n.powf(3.0))/48.0)*(1.0 - n + (11.0/16.0)*(n.powf(2.0) - n.powf(3.0)));
    let E:f64 = (315.0 * a * (n.powf(4.0))/51.0)*(1.0 - n);
    
    let S:f64 = A * lat - B * (2.0 * lat).sin() +C * (4.0 * lat).sin() - D * (6.0 * lat).sin() + E * (8.0 * lat).sin();

    // get x
    let K1 = S*k0;
    let K2 = k0*nu*(2.0*lat).sin()/4.0;
    let K3 = (k0*nu*lat.sin()*lat.cos().powf(3.0)/24.0) * (5.0 - lat.tan().powf(2.0) + 9.0 * e2 * lat.cos().powf(2.0) + 4.0*(e2.powf(2.0))*(lat.cos().powf(4.0)));        
    let y = K1 + K2*p.powf(2.0) + K3*p.powf(4.0);

    // get y
    let K4 = k0*nu*lat.cos();
    let K5 = (k0*nu*lat.cos().powf(3.0)/6.0) * (1.0 - lat.tan().powf(2.0) + e2*(lat.cos().powf(2.0)));
    let x = K4 * p + K5 * p.powf(3.0) + dx;

    (x, y)
}

pub fn wgs84_to_twd97(lng:f64, lat:f64) -> (f64, f64) {
    let k0:f64 = 0.9999;
    let dx:f64 = 250000.0;
    return lnglat_to_mercator(lng, lat, 121.0, k0, dx)
}

pub fn wgs84_to_2degree_zone(lng:f64, lat:f64, center_lng:f64) -> (f64, f64) {
    let k0:f64 = 0.9999;
    let dx:f64 = 250000.0;
    return lnglat_to_mercator(lng, lat, center_lng, k0, dx)
}

pub fn wgs84_to_3degree_zone(lng:f64, lat:f64, center_lng:f64) -> (f64, f64) {
    let k0:f64 = 1.0;
    let dx:f64 = 350000.0;
    return lnglat_to_mercator(lng, lat, center_lng, k0, dx)
}

pub fn wgs84_to_6degree_zone(lng:f64, lat:f64, center_lng:f64) -> (f64, f64) {
    let k0:f64 = 0.9996;
    let dx:f64 = 500000.0;
    return lnglat_to_mercator(lng, lat, center_lng, k0, dx)
}

/// convert croodinate from mercator projection to wgs84
/// - https://en.wikipedia.org/wiki/Mercator_projection
/// # Examples
///
/// ```
/// use mercator::lnglat_to_mercator;
/// let k0:f64 = 0.9999;
/// let dx:f64 = 250000.0;
/// let x:f64 = 248170.82572;
/// let y:f64 = 2652129.9773;
/// let (lng, lat) = lnglat_to_mercator(x, y, 121.0, k0, dx);
/// ```
pub fn mercator_to_lnglat(x:f64, y:f64, center_lng:f64, k0:f64, dx:f64) -> (f64, f64) {
    let a:f64 = 6378137.0;
    let b:f64 = 6356752.314245;
    let lng0:f64 = center_lng * PI / 180.0;

    let dy:f64 = 0.0;
    let e:f64 = (1.0 - b.powf(2.0)/a.powf(2.0)).powf(0.5);

    let x:f64 = x - dx;
    let y:f64 = y - dy;

    // calculate the meridional arc
    let m:f64 = y/k0;

    // calculate Footprint Latitude
    let mu:f64 = m/(a*(1.0 - e.powf(2.0)/4.0 - 3.0*e.powf(4.0)/64.0 - 5.0*e.powf(6.0)/256.0));
    let e1:f64 = (1.0 - (1.0 - e.powf(2.0)).powf(0.5)) / (1.0 + (1.0 - e.powf(2.0)).powf(0.5));

    let j1 = 3.0*e1/2.0 - 27.0*e1.powf(3.0)/32.0;
    let j2 = 21.0*e1.powf(2.0)/16.0 - 55.0*e1.powf(4.0)/32.0;
    let j3 = 151.0*e1.powf(3.0)/96.0;
    let j4 = 1097.0*e1.powf(4.0)/512.0;

    let fp = mu + j1*(2.0*mu).sin() + j2*(4.0*mu).sin() + j3*(6.0*mu).sin() + j4*(8.0*mu).sin();

    // calculate Latitude and Longitude

    let e2 = (e*a/b).powf(2.0);
    let c1 = e2*fp.cos().powf(2.0);
    let t1 = fp.tan().powf(2.0);
    let r1 = a*(1.0-e.powf(2.0))/(1.0-e.powf(2.0)*fp.sin().powf(2.0)).powf(3.0/2.0);
    let n1 = a/(1.0-e.powf(2.0)*fp.sin().powf(2.0)).powf(0.5);

    let d = x/(n1*k0);

    // get lat
    let q1 = n1*fp.tan()/r1;
    let q2 = d.powf(2.0)/2.0;
    let q3 = (5.0 + 3.0*t1 + 10.0*c1 - 4.0*c1.powf(2.0) - 9.0*e2)*d.powf(4.0)/24.0;
    let q4 = (61.0 + 90.0*t1 + 298.0*c1 + 45.0*t1.powf(2.0) - 3.0*c1.powf(2.0) - 252.0*e2)*d.powf(6.0)/720.0;
    let lat = fp - q1*(q2 - q3 + q4);

    // get lng
    let q5 = d;
    let q6 = (1.0 + 2.0*t1 + c1)*d.powf(3.0)/6.0;
    let q7 = (5.0 - 2.0*c1 + 28.0*t1 - 3.0*c1.powf(2.0) + 8.0*e2 + 24.0*t1.powf(2.0))*d.powf(5.0)/120.0;
    let lng = lng0 + (q5 - q6 + q7)/fp.cos();

    let lat = (lat * 180.0) / PI;
    let lng = (lng * 180.0) / PI;
    return (lng, lat)
}

pub fn f2degree_zone_to_wgs84(x:f64, y:f64, center_lng:f64) -> (f64, f64) {
    let k0:f64 = 0.9999;
    let dx:f64 = 250000.0;
    return lnglat_to_mercator(x, y, center_lng, k0, dx)
}

pub fn f3degree_zone_to_wgs84(x:f64, y:f64, center_lng:f64) -> (f64, f64) {
    let k0:f64 = 1.0;
    let dx:f64 = 350000.0;
    return lnglat_to_mercator(x, y, center_lng, k0, dx)
}

pub fn f6degree_zone_to_wgs84(x:f64, y:f64, center_lng:f64) -> (f64, f64) {
    let k0:f64 = 0.9996;
    let dx:f64 = 500000.0;
    return lnglat_to_mercator(x, y, center_lng, k0, dx)
}
