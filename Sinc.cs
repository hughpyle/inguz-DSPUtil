using System;
using System.Collections.Generic;
using System.Text;

using DSPUtil;

namespace DSPUtil
{
    public class Sinc : SoundObj
    {
        int _size;
        double _fc;

        public Sinc(int size, double fC_Fraction)
        {
            _size = size;
            _fc = fC_Fraction;
        }

        public override IEnumerator<ISample> Samples
        {
            get
            {
                ushort nc = NumChannels;
                int m2 = _size/2;
                double k = 1 / (2 * Math.PI * _fc);
                for (int i = 0; i < _size; i++)
                {
                    int n = i - m2;
                    double val;
                    if (n == 0)
                    {
                        val = 1;
                    }
                    else
                    {
                        val = Math.Sin(2 * Math.PI * _fc * n) * k / n;
                    }
                    ISample s = (nc == 2) ? new Sample2() : new Sample(nc) as ISample;
                    for (int c = 0; c < nc; c++)
                    {
                        s[c] = val;
                    }
                    yield return s;
                }
            }
        }
    }
}
