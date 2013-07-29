#ifndef INCLUDE_HALTON_H
#define INCLUDE_HALTON_H

class Halton: public Rng
{
public:

    Halton(long i, int base)
    {
        _base = base;
        double f = inv_base = 1.0/base;
        value = 0.0;
        while (i>0)
        {
            value += f* (double)(i % base);
            i /= base;
            f *= inv_base;
        }
    }

    void seed(long i = 1)
    {
        double f = inv_base;
        value = 0.0;
        while (i>0)
        {
            value += f* (double)(i % _base);
            i /= _base;
            f *= inv_base;
        }
    }

    double next()
    {
        double r = 1.0 - value - 0.0000000001;
        if (inv_base < r)
            value += inv_base;
        else
        {
            double h = inv_base, hh;
            do
            {
                hh = h;
                h *= inv_base;
            } while (h >= r);
            value += hh + h - 1.0;
        }

        return value;  //*0.8 + rng()*0.2;
    }

    operator double()
    {
        return value;
    }

private:
  int _base;
  double value;
  double inv_base;
};

#endif
