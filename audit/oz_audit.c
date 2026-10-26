/*=====================================================================
 oz_audit.c —— Day 35 C 证人: NOx 元素守恒 + O3 收支闭合
   nox_budget: NOx 总量账 = 初值 + 排放积分 − 沉降积分 − 末期量 (Kahan)
               (NOx 家族守恒: 机制内 NO/NO2/HNO3 互换不改 N 原子总数;
                背景松弛不作用于 N 物种 ⇒ 望远镜闭合应≈0)
   o3_budget : 5 天 O3 净变化 = ∫(P−L)dt (梯形 Kahan)
=====================================================================*/
#include <math.h>
#include <stdlib.h>

static double kahan_sum(const double *x, long n) {
    double s = 0.0, c = 0.0;
    for (long i = 0; i < n; i++) {
        double y = x[i] - c, t = s + y;
        c = (t - s) - y;
        s = t;
    }
    return s;
}

void nox_budget(long n, const double *nox_series, double nox0, double dt,
                double *out) {
    /* out[0]=ΔNOx  out[1]=Σ|Δ|逐时(漂移指示) */
    double cur = nox_series[n-1];
    out[0] = cur - nox0;
    double *d = (double *)malloc(sizeof(double)*n);
    for (long i = 0; i < n; i++) d[i] = fabs(nox_series[i] - nox0);
    out[1] = kahan_sum(d, n)/(double)n;
    free(d);
}

void o3_budget(long n, const double *p_series, const double *l_series,
               double dt, double *out) {
    double *net = (double *)malloc(sizeof(double)*n);
    for (long i = 0; i < n; i++) net[i] = p_series[i] - l_series[i];
    out[0] = kahan_sum(net, n)*dt;      /* ∫(P-L)dt */
    out[1] = 0.5*(net[0] + net[n-1])*dt; /* 梯形半格 */
    out[0] -= out[1];
    free(net);
}
