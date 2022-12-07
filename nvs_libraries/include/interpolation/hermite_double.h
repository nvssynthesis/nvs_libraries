// original
static inline double hermite1(double x, double y0, double y1, double y2, double y3)
{
    // 4-point, 3rd-order Hermite (x-form)
    double c0 = y1;
    double c1 = 0.5f * (y2 - y0);
    double c2 = y0 - 2.5 * y1 + 2 * y2 - 0.5 * y3;
    double c3 = 1.5 * (y1 - y2) + 0.5 * (y3 - y0);

    return ((c3 * x + c2) * x + c1) * x + c0;
}

// james mccartney
static inline double hermite2(double x, double y0, double y1, double y2, double y3)
{
    // 4-point, 3rd-order Hermite (x-form)
    double c0 = y1;
    double c1 = 0.5 * (y2 - y0);
    double c3 = 1.5 * (y1 - y2) + 0.5 * (y3 - y0);
    double c2 = y0 - y1 + c1 - c3;

    return ((c3 * x + c2) * x + c1) * x + c0;
}

// james mccartney
static inline double hermite3(double x, double y0, double y1, double y2, double y3)
{
        // 4-point, 3rd-order Hermite (x-form)
        double c0 = y1;
        double c1 = 0.5 * (y2 - y0);
        double y0my1 = y0 - y1;
        double c3 = (y1 - y2) + 0.5 * (y3 - y0my1 - y2);
        double c2 = y0my1 + c1 - c3;

        return ((c3 * x + c2) * x + c1) * x + c0;
}

// laurent de soras
static inline double hermite4(double frac_pos, double xm1, double x0, double x1, double x2)
{
   const double    c     = (x1 - xm1) * 0.5;
   const double    v     = x0 - x1;
   const double    w     = c + v;
   const double    a     = w + v + (x2 - x0) * 0.5;
   const double    b_neg = w + a;

   return ((((a * frac_pos) - b_neg) * frac_pos + c) * frac_pos + x0);
}
