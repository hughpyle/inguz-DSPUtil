using System;
using System.Collections.Generic;
using System.Text;

// Copyright (c) 2006 by Hugh Pyle, inguzaudio.com

namespace DSPUtil
{
    public interface ISample
    {
        ushort NumChannels { get; }
        double this[int arg] { get; set; }
    }

    [Serializable]
    public struct Sample2 : ISample
    {
        double _0;
        double _1;
        public Sample2(double a, double b)
        {
            _0 = a;
            _1 = b;
        }
        public ushort NumChannels { get { return 2; } }
        public double this[int arg]
        {
            get
            {
                return (arg == 0) ? _0 : _1;
            }
            set
            {
                if (arg == 0) _0 = value; else _1 = value;
            }
        }
        public override bool Equals(object obj)
        {
            ISample o = obj as ISample;
            if (o == null)
            {
                return false;
            }
            return o.NumChannels == 2 && o[0] == _0 && o[1] == _1;
        }
        public override int GetHashCode()
        {
            return _0.GetHashCode() ^ _1.GetHashCode();
        }
    }


    [Serializable]
    public struct Sample : ISample
    {
        // A multi-channel sample value.  Each data point is a double-precision float.
        double[] _value;

        public Sample(ushort numChannels)
        {
            _value = new double[numChannels];
        }

        // Single-channel constructor
        public Sample(double val)
        {
            _value = new double[1];
            _value[0] = val;
        }

        // Copy constructor
        public Sample(ISample sample)
        {
            ushort nc = sample.NumChannels;
            _value = new double[nc];
            for (int n = 0; n < nc; n++)
            {
                _value[n] = sample[n];
            }
        }

        public Sample(ISample sample, double gain)
        {
            ushort nc = sample.NumChannels;
            _value = new double[nc];
            for (int n = 0; n < nc; n++)
            {
                _value[n] = sample[n] * gain;
            }
        }

        public ushort NumChannels { get { return (ushort)_value.Length; } }

        public double this[int arg]
        {
            get
            {
                return _value[arg];
            }
            set
            {
                _value[arg] = value;
            }
        }

        public override bool Equals(object obj)
        {
            ISample o = obj as ISample;
            if (o == null)
            {
                return false;
            }
            ushort nc = NumChannels;
            if (o.NumChannels != nc)
            {
                return false;
            }
            for (int n = 0; n < nc; n++)
            {
                if (o[n] != _value[n])
                {
                    return false;
                }
            }
            return true;
        }

        public override int GetHashCode()
        {
            int h = 0;
            ushort nc = NumChannels;
            for (int n = 0; n < nc; n++)
            {
                h = h ^ _value[n].GetHashCode();
            }
            return h;
        }

        public static Sample operator +(Sample c1, Sample c2)
        {
            if (c1.NumChannels != c2.NumChannels)
            {
                throw new ArgumentException("Samples must have the same number of channels");
            }
            Sample s = new Sample(c1.NumChannels);
            for (int n = 0; n < c1.NumChannels; n++)
            {
                s[n] = c1[n] + c2[n];
            }
            return s;
        }

        public static Sample operator -(Sample c1, Sample c2)
        {
            if (c1.NumChannels != c2.NumChannels)
            {
                throw new ArgumentException("Samples must have the same number of channels");
            }
            Sample s = new Sample(c1.NumChannels);
            for (int n = 0; n < c1.NumChannels; n++)
            {
                s[n] = c1[n] - c2[n];
            }
            return s;
        }

        // Multiplication operator (needed for convolution)
        public static Sample operator *(Sample c1, Sample c2)
        {
            ushort nc = c1.NumChannels;
            if (nc != c2.NumChannels)
            {
                throw new ArgumentException("Samples must have the same number of channels");
            }
            Sample s = new Sample(nc);
            for (int n = 0; n < nc; n++)
            {
                s[n] = c1[n] * c2[n];
            }
            return s;
        }

        // Multiply by scalar (gain)
        public static Sample operator *(double gain, Sample c)
        {
            ushort nc = c.NumChannels;
            Sample s = new Sample( nc );
            for (int n = 0; n < nc; n++)
            {
                s[n] = c[n] * gain;
            }
            return s;
        }
    }

    /*
    [Serializable]
    public struct ComplexSample : ISample<Complex>
    {
        // A multi-channel sample value.  Each data point is a double-precision complex number.
        Complex[] _value;

        public ComplexSample(ushort numChannels)
        {
            _value = new Complex[numChannels];
        }

        // Copy constructor
        public ComplexSample(Sample sample)
        {
            ushort nc = sample.NumChannels;
            _value = new Complex[nc];
            for (int n = 0; n < nc; n++)
            {
                _value[n] = new Complex(sample[n],0);
            }
        }

        public ushort NumChannels { get { return (ushort)_value.Length; } }
        public Complex[] Value
        {
            get { return _value; }
        }

        public Complex this[int arg]
        {
            get
            {
                return _value[arg];
            }
            set
            {
                _value[arg] = value;
            }
        }

        public static ComplexSample operator +(ComplexSample c1, ComplexSample c2)
        {
            if (c1.NumChannels != c2.NumChannels)
            {
                throw new ArgumentException("Samples must have the same number of channels");
            }
            ComplexSample s = new ComplexSample(c1.NumChannels);
            for (int n = 0; n < c1.NumChannels; n++)
            {
                s[n] = c1[n] + c2[n];
            }
            return s;
        }

        public static ComplexSample operator -(ComplexSample c1, ComplexSample c2)
        {
            if (c1.NumChannels != c2.NumChannels)
            {
                throw new ArgumentException("Samples must have the same number of channels");
            }
            ComplexSample s = new ComplexSample(c1.NumChannels);
            for (int n = 0; n < c1.NumChannels; n++)
            {
                s[n] = c1[n] - c2[n];
            }
            return s;
        }

        // Multiplication operator (needed for convolution)
        public static ComplexSample operator *(ComplexSample c1, ComplexSample c2)
        {
            if (c1.NumChannels != c2.NumChannels)
            {
                throw new ArgumentException("Samples must have the same number of channels");
            }
            ComplexSample s = new ComplexSample(c1.NumChannels);
            for (int n = 0; n < c1.NumChannels; n++)
            {
                s[n] = c1[n] * c2[n];
            }
            return s;
        }

        // Multiply by scalar (gain)
        public static ComplexSample operator *(double gain, ComplexSample c)
        {
            ComplexSample s = new ComplexSample(c.NumChannels);
            for (int n = 0; n < c.NumChannels; n++)
            {
                s[n] = c[n] * gain;
            }
            return s;
        }
    }
    */
}
