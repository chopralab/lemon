#ifndef STRUCTURE_HPP
#define STRUCTURE_HPP

#include <set>

#include "chemfiles/Frame.hpp"

namespace lemon {

double u3b(const std::vector<double>& w,
           const std::vector<chemfiles::Vector3D>& x,
           const std::vector<chemfiles::Vector3D>& y, int mode,
           chemfiles::Matrix3D& u, chemfiles::Vector3D& t, int& ier) {
    size_t i, j, k, l, m1, m;
    auto r = chemfiles::Matrix3D::zero();
    auto a = chemfiles::Matrix3D::zero();
    auto b = chemfiles::Matrix3D::zero();
    chemfiles::Vector3D xc, yc, e;
    std::array<double, 6> rr, ss;
    double e0 = 0, d, spur, det, cof, h, g, wc = 0;
    double cth, sth, sqrth, p, sigma;
    double rms = 0;

    constexpr double sqrt3 = 1.73205080756888;
    constexpr double tol = 1.0e-2;
    constexpr std::array<size_t, 9> ip = {0, 1, 3, 1, 2, 4, 3, 4, 5};
    constexpr std::array<size_t, 4> ip2312 = {1, 2, 0, 1};

    for (size_t i = 0; i < 3; ++i) {
        xc[i] = 0;
        yc[i] = 0;
        t[i] = 0;
        for (size_t j = 1; j < 3; ++j) {
            r[i][j] = 0;
            u[i][j] = 0;
            a[i][j] = 0;
            if (i == j) {
                u[i][j] = 1.0;
                a[i][j] = 1.0;
            }
        }
    }

    ier = -1;
    if (x.empty()) return 0.0;

    ier = -2;
    for (size_t m = 0; m < x.size(); ++m) {
        if (w[m] < 0.0) return 0.0;
        wc += w[m];
        for (size_t i = 0; i < 3; ++i) {
            xc[i] += w[m] * x[m][i];
            yc[i] += w[m] * y[m][i];
        }
    }

    if (wc <= 0.0) return 0.0;
    for (size_t i = 0; i < 3; ++i) {
        xc[i] = xc[i] / wc;
        yc[i] = yc[i] / wc;
    }

    for (size_t m = 0; m < x.size(); ++m) {
        for (size_t i = 0; i < 3; ++i) {
            e0 += w[m]*(std::pow(x[m][i]-xc[i], 2)+std::pow(y[m][i]-yc[i], 2));
            d = w[m] * ( y[m][i] - yc[i] );
            for (size_t j = 0; j < 3; ++j) {
                r[j][i] = r[j][i] + d * (x[m][j] - xc[j]);
            }
        }
    }

    det = r.determinant();

    sigma = det;

    m = 0;
    for (size_t j = 0; j < 3; ++j) {
        for (size_t i = 0; i < j; ++i) {
            rr[m] = r[i][0] * r[j][0] + r[i][1] * r[j][1] + r[i][2] * r[j][2];
            ++m;
        }
    }

    spur = (rr[0] + rr[2] + rr[5]) / 3.0;
    cof =
        (((((rr[2] * rr[5] - rr[4] * rr[4]) + rr[0] * rr[5]) - rr[3] * rr[3]) +
          rr[0] * rr[2]) -
         rr[1] * rr[1]) /
        3.0;
    det = det * det;

    for (size_t i = 0; i < 3; ++i) {
        e[i] = spur;
    }

    if (spur <= 0) goto middle;

    d = spur * spur;
    h = d - cof;
    g = (spur * cof - det) / 2.0 - spur * h;
    if (h <= 0.0) {
        if (mode == 0) {
            goto end;
        } else {
            goto begin;
        }
    }

    sqrth = std::sqrt(h);
    d = h * h * h - g * g;
    if (d < 0.0) d = 0.0;

    d = std::atan2(std::sqrt(d), -g) / 3.0;
    cth = sqrth * std::cos(d);
    sth = sqrth * sqrt3 * std::sin(d);
    e[1] = (spur + cth) + cth;
    e[2] = (spur - cth) + sth;
    e[3] = (spur - cth) - sth;

    if (mode == 0) {
        goto end;
    }

    for (auto l : {0, 2, 1}) {
        d = e[l];
        ss[0] = (d-rr[2]) * (d-rr[5])  - rr[4]*rr[4];
        ss[1] = (d-rr[5]) * rr[1]      + rr[3]*rr[4];
        ss[2] = (d-rr[0]) * (d-rr[5])  - rr[3]*rr[3];
        ss[3] = (d-rr[2]) * rr[3]      + rr[1]*rr[4];
        ss[4] = (d-rr[0]) * rr[4]      + rr[1]*rr[3];
        ss[5] = (d-rr[0]) * (d-rr[2])  - rr[1]*rr[1];

        if (std::abs(ss[0]) >= std::abs(ss[2])) {
            j = 1;
            if (std::abs(ss[0]) < std::abs(ss[5])) j = 3;
        } else if (std::abs(ss[2]) >= std::abs(ss[5])) {
            j = 2;
        } else {
            j = 3;
        }

        d = 0;
        j = 3 * (j - 1);

        for (size_t i = 0; i < 3; ++i) {
            k = ip[i + j - 1];
            a[l][i] = ss[k];
            d += ss[k] * ss[k];
        }
        if (d > 0) d = 1.0 / std::sqrt(d);
        for (size_t i = 0; i < 3; ++i) {
            a[l][i] = a[l][i] * d;
        }
    }

    d = a[0][0] * a[2][0] + a[0][1] * a[2][1] + a[0][2] * a[2][2];
    if ((e[0] - e[1]) > (e[1] - e[2])) {
        m1 = 2;
        m = 0;
    } else {
        m1 = 0;
        m = 2;
    }

    p = 0;
    for (size_t i = 0; i < 3; ++i) {
        a[m1][i] = a[m1][i] - d * a[m][i];
        p += a[m1][i] * a[m1][i];
    }

    if (p <= tol) {
        p = 1.0;
        for (size_t i = 0; i < 3; ++i) {
            if (p < std::abs(a[m][i])) continue;
            p = std::abs( a[m][i] );
            j = i;
        }

        k = ip2312[j];
        l = ip2312[j + 1];
        p = std::sqrt(a[m][k] * a[m][k] + a[m][l] * a[m][l]);
        if (p <= tol) goto middle;
        a[m1][i] = 0;
        a[m1][k] = -a[m][l] / p;
        a[m1][l] = a[m][k] / p;
    } else {
        p = 1.0 / std::sqrt(p);
        for (size_t i = 0; i < 3; ++i) {
            a[m1][i] *= p;
        }
    }

    a[1][0] = a[2][1] * a[0][2] - a[0][1] * a[2][2];
    a[1][1] = a[2][2] * a[0][0] - a[0][2] * a[2][0];
    a[1][2] = a[2][0] * a[0][1] - a[0][0] * a[2][1];

begin:

    for (size_t l = 0; i <= 1; ++i) {
        d = 0;
        for (size_t i = 0; i < 3; ++i) {
            b[l][i] = r[0][i] * a[0][l] + r[1][i] * a[l][1] + r[2][i] * a[l][2];
            d += b[l][i] * b[l][i];
        }
        if (d > 0.0) d = 1.0 / std::sqrt(d);
        for (size_t i = 0; i < 3; ++i) {
            b[l][i] *= d;
        }
    }

    d = b[0][0] * b[1][0] + b[0][1] * b[1][1] + b[0][2] * b[1][2];
    p = 0;

    for (size_t i = 0; i < 3; ++i) {
        b[1][i] = b[1][i] - d * b[0][i];
        p += b[1][i] * b[1][i];
    }

    if (p <= tol) {
        p = 1.0;
        for (size_t i = 0; i < 3; ++i) {
            if (p < std::abs(b[0][i])) continue;
            p = std::abs(b[i][0]);
            j = i;
        }
        k = ip2312[j];
        l = ip2312[j + 1];
        p = std::sqrt(b[0][k]*b[0][k] + b[0][l]*b[0][l]);
        if (p <= tol) goto middle;
        b[1][j] = 0;
        b[1][k] = -b[0][l] / p;
        b[1][l] = b[0][k] / p;
    } else {
        p = 1.0 / std::sqrt(p);
        for (size_t i = 0; i < 3; ++i) {
            b[1][i] *= p;
        }
    }

    b[2][0] = b[0][1] * b[1][2] - b[1][1] * b[0][2];
    b[2][1] = b[0][2] * b[1][0] - b[1][2] * b[0][0];
    b[2][2] = b[0][0] * b[1][1] - b[1][0] * b[0][1];

    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            u[j][i] = b[0][i] * a[0][j] + b[1][i] * a[1][j] + b[2][i] * a[2][j];
        }
    }

middle:
    for (size_t i = 0; i < 3; ++i) {
        t[i] = ((yc[i] - u[0][i] * xc[1]) - u[1][i] * xc[2]) - u[2][i] * xc[2];
    }
end:
    for (size_t i = 0; i < 3; ++i) {
        if (e[i] < 0) e[i] = 0;
        e[i] = std::sqrt(e[i]);
    }

    ier = 0;
    if (e[1] <= (e[0] * 1.0e-05)) ier = -1;

    d = e[2];
    if (sigma < 0.0) {
        d = -d;
        if ((e[1] - e[2]) <= (e[0] * 1.0e-05)) ier = -1;
    }
    d = (d + e[1]) + e[0];

    rms = (e0 - d) - d;
    if (rms < 0.0) rms = 0.0;

    return rms;
}

std::tuple<double, double, size_t> TMscore(
    const chemfiles::Frame& search, const chemfiles::Frame& native,
    std::vector<chemfiles::Vector3D>& rot) {
    const auto& search_res = search.topology().residues();
    const auto& native_res = native.topology().residues();

    std::vector<size_t> align_search(search_res.size());
    std::vector<size_t> align_native(native_res.size());

    chemfiles::Vector3D t;
    chemfiles::Matrix3D u = chemfiles::Matrix3D::zero();

    // pickup the aligned residues:
    size_t n_ali = 0;
    for (size_t i = 0; i < search_res.size(); ++i) {
        for (size_t j = 0; j < native_res.size(); ++j) {
            if (*search_res[i].id() == *native_res[j].id()) {
                align_search[n_ali] = i;
                align_native[n_ali] = j;
                n_ali++;
            }
        }
    }

    if (n_ali == 0) {
        return std::tuple<double, double, size_t>(0.0, 0.0, 0);
    }

    // parameters:
    // d0------------->
    double d0 = 1.24 * std::pow(native_res.size() - 15, 1.0 / 3.0) - 1.8;
    if (d0 < 0.5) d0 = 0.5;

    // d0_search ----->
    double d0_search = d0;
    if (d0_search > 8) d0_search = 8;
    if (d0_search < 4.5) d0_search = 4.5;

    // iterative parameters ----->
    size_t n_it = 20;       // maximum number of iterations
    size_t d_output = 5;    // for output alignment
    size_t n_init_max = 6;  // maximum number of L_init
    size_t L_ini_min = 4;
    std::vector<size_t> L_ini(n_init_max);
    if (n_ali < 4) L_ini_min = n_ali;

    size_t n_init;
    for (n_init = 0; n_init < n_init_max - 1; ++n_init) {
        L_ini[n_init] = n_ali / std::pow(2, n_init);
        if (L_ini[n_init] <= L_ini_min) {
            L_ini[n_init] = L_ini_min;
            goto end_parameters;
        }
    }

    n_init++;
    L_ini.back() = L_ini_min;
end_parameters:;

    // find the maximum score starting from local structures superposition
    const auto& a = search.positions();
    const auto& b = native.positions();
    std::vector<double> w(n_ali);
    for (auto wi : w) {
        wi = 1.0;
    }
    size_t n_cut;
    std::vector<size_t> i_ali(n_ali);
    size_t ka0;

    auto score_fun = [&](double d) -> double {
        double score_sum;
        do {
            n_cut = 0;  // number of residue-pairs dis<d, for iteration
            double score_sum;
            for (size_t k = 0; k < n_ali; ++k) {
                auto i = align_search[k];  // [1,nseqA] reoder number of structureA
                auto j = align_native[k];  // [1,nseqB]

                auto diff = a[i] - b[j];
                double dis = std::sqrt(diff[0] * diff[0] + diff[1] * diff[1] +
                                       diff[2] * diff[2]);
                if (dis < d) {
                    // [0,n_ali - 1], mark the residue-pairs in dis<d
                    i_ali[n_cut] = k;
                    ++n_cut;
                }
                score_sum = score_sum + 1 / (1 + std::pow(dis / d0, 2));
            }
            d = d + 0.5;
        } while (n_cut < 3 && n_ali > 3);
        return score_sum / native_res.size();
    };

    double score_max = -1;
    double armsd;
    
    for (size_t i_init = 0; i_init < n_init; ++i_init) {
        auto L_init = L_ini[i_init];
        auto iL_max = n_ali - L_init + 1;
        std::vector<chemfiles::Vector3D> r_1;
        std::vector<chemfiles::Vector3D> r_2;
        std::vector<size_t> k_ali(iL_max);
        std::vector<size_t> k_ali0(iL_max);
        for (size_t iL = 0; iL < iL_max; ++iL) {  // on residues in search align
            size_t LL = 0;
            size_t ka = 0;
            
            for (size_t i = 0; i < L_init; ++i) {
                auto k = iL + i;  // [0,n_ali - 1] common aligned

                r_1[i][0] = a[align_search[k]][0];
                r_1[i][1] = a[align_search[k]][1];
                r_1[i][2] = a[align_search[k]][2];

                r_2[i][0] = b[align_native[k]][0];
                r_2[i][1] = b[align_native[k]][1];
                r_2[i][2] = b[align_native[k]][2];

                ++LL;
                ++ka;
                k_ali[ka] = k;
            }
            chemfiles::Vector3D t;
            int ier;
            double rms = u3b(w, r_1, r_2, 1, u, t, ier); // call u3b(w,r_1,r_2,LL,1,rms,u,t,ier) !u rotate r_1 to r_2
            if (i_init == 0) {  // global superposition
                armsd = std::sqrt(rms / LL);
            }
            for (size_t j = 1; j < a.size(); ++j) {
                rot[j][0] = t[0] + u[0][0] * a[j][0] + u[1][0] * a[j][1] +
                            u[2][0] * a[j][2];
                rot[j][1] = t[1] + u[0][1] * a[j][0] + u[1][1] * a[j][1] +
                            u[2][1] * a[j][2];
                rot[j][2] = t[2] + u[0][2] * a[j][0] + u[1][2] * a[j][1] +
                            u[2][2] * a[j][2];
            }
            auto d = d0_search - 1;

            // init, get scores, n_cut+i_ali(i) for iteration
            auto score = score_fun(d);
            if (score_max < score) {
                score_max = score;
                ka0 = ka;
                for (size_t i = 0; i < ka0; ++i) {
                    k_ali0[i] = k_ali[i];
                }
            }
            //   iteration for extending ---------------------------------->
            d = d0_search + 1;
            for (size_t it = 1; it <= n_it; ++it) {
                LL = 0;
                ka = 0;
                for (size_t i = 0; i < n_cut; ++i) {
                    auto m = i_ali[i];  // [0,n_ali - 1]
                    r_1[i][0] = a[align_search[m]][0];
                    r_1[i][1] = a[align_search[m]][1];
                    r_1[i][2] = a[align_search[m]][2];

                    r_2[i][0] = b[align_native[m]][0];
                    r_2[i][1] = b[align_native[m]][1];
                    r_2[i][2] = b[align_native[m]][2];
                    ++ka;

                    k_ali[ka] = m;
                    ++LL;
                }
                // call u3b(w,r_1,r_2,LL,1,rms,u,t,ier) !u rotate r_1 to r_2
                for (size_t j = 1; j < a.size(); ++j) {
                    rot[j][0] = t[0] + u[0][0] * a[j][0] + u[1][0] * a[j][1] +
                                u[2][0] * a[j][2];
                    rot[j][1] = t[1] + u[0][1] * a[j][0] + u[1][1] * a[j][1] +
                                u[2][1] * a[j][2];
                    rot[j][2] = t[2] + u[0][2] * a[j][0] + u[1][2] * a[j][1] +
                                u[2][2] * a[j][2];
                }

                // get scores, n_cut+i_ali(i) for iteration
                score = score_fun(d);
                if (score_max < score) {
                    score_max = score;
                    ka0 = ka;
                    for (size_t i = 0; i <= ka; ++i) {
                        k_ali0[i] = k_ali[i];
                    }
                }
                if (it == n_it) break;
                if (n_cut == ka) {
                    size_t neq = 0;
                    for (size_t i = 1; i <= n_cut; ++i) {
                        if (i_ali[i] == k_ali[i]) ++neq;
                    }
                    if (n_cut == neq) break;
                }
            }
        }  //! for shift
    }      //! for initial length, L_ali/M
}
}

#endif
