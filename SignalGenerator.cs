using System;
using System.Collections.Generic;
using System.Text;

// Copyright (c) 2006, 2007 by Hugh Pyle, inguzaudio.com

namespace DSPUtil
{
    /// <summary>
    /// 0dbfs signal generators
    /// </summary>
    public class SignalGenerator : SoundObj
    {
        protected const double twopi = 2 * Math.PI;
        protected double _freq;
        protected double _gain;

        public SignalGenerator(ushort numChannels, uint sampleRate, double freq, double gain)
        {
            NumChannels = numChannels;
            SampleRate = sampleRate;
            _freq = freq;
            _gain = gain;
        }

        /// <summary>
        /// Get an iterator for samples
        /// </summary>
        public override IEnumerator<ISample> Samples
        {
            get
            {
                yield break;
            }
        }

        /// <summary> Number of iterations expected to do the signal processing </summary>
        public override int Iterations
        {
            get { return int.MaxValue; }
        }
    }

    public class SineGenerator : SignalGenerator
    {
        public SineGenerator(ushort numChannels, uint sampleRate, double freq, double gain)
            : base(numChannels, sampleRate, freq, gain)
        {
        }

        /// <summary>
        /// Get an iterator for samples
        /// </summary>
        public override IEnumerator<ISample> Samples
        {
            get
            {
                ushort nc = NumChannels;
                uint sr = SampleRate;
                double mul = twopi * _freq / sr;
                int n = 0;
                Sample s = new Sample(nc);
                while (true)
                {
                    double v = _gain * Math.Sin(n * mul);
                    for (int c = 0; c < nc; c++)
                    {
                        s[c] = v;
                    }
                    n++;
                    yield return s;
                }
            }
        }
    }

    public class SineQuadGenerator : SignalGenerator
    {
        public SineQuadGenerator(ushort numChannels, uint sampleRate, double freq, double gain)
            : base(numChannels, sampleRate, freq, gain)
        {
        }

        /// <summary>
        /// Get an iterator for samples
        /// </summary>
        public override IEnumerator<ISample> Samples
        {
            get
            {
                ushort nc = NumChannels;
                if (nc < 2)
                {
                    throw new ArgumentOutOfRangeException("NumChannels", "Quadrature requires two channels.");
                }
                uint sr = SampleRate;
                double mul = twopi * _freq / sr;
                int n = 0;
                Sample s = new Sample(nc);
                while (true)
                {
                    s[0] = _gain * Math.Cos(n * mul);
                    s[1] = _gain * Math.Sin(n * mul);
                    n++;
                    yield return s;
                }
            }
        }
    }

    public class SquareGenerator : SignalGenerator
    {
        public SquareGenerator(ushort numChannels, uint sampleRate, double freq, double gain)
            : base(numChannels, sampleRate, freq, gain)
        {
        }

        /// <summary>
        /// Get an iterator for samples
        /// </summary>
        public override IEnumerator<ISample> Samples
        {
            get
            {
                ushort nc = NumChannels;
                uint sr = SampleRate;
                double mul = 2 * _freq / sr;
                int n = 0;
                Sample s = new Sample(nc);
                while (true)
                {
                    double saw = ((n * mul) % 2) - 1;
                    double v = saw >0 ? _gain : -_gain;
                    for (int c = 0; c < nc; c++)
                    {
                        s[c] = v;
                    }
                    n++;
                    yield return s;
                }
            }
        }
    }

    public class Harmonic
    {
        public double Gain;
        public double Phase;
        public Harmonic(double gain)
        {
            Gain = gain;
            Phase = 0;
        }
        public Harmonic(double gain, double phase)
        {
            Gain = gain;
            Phase = phase;
        }
    }

    public abstract class HarmonicGenerator : SignalGenerator
    {
        public HarmonicGenerator(ushort numChannels, uint sampleRate, double freq, double gain)
            : base(numChannels, sampleRate, freq, gain)
        {
        }

        public virtual List<Harmonic> Harmonics
        {
            get
            {
                List<Harmonic> harmonics = new List<Harmonic>();
                return harmonics;
            }
        }

        /// <summary>
        /// Get an iterator for samples
        /// </summary>
        public override IEnumerator<ISample> Samples
        {
            get
            {
                ushort nc = NumChannels;
                uint sr = SampleRate;
                double mul = twopi * _freq / sr;
                int n = 0;
                List<Harmonic> lH = Harmonics;
                List<Sample> samples = new List<Sample>((int)sr);
                while(n<sr)
                {
                    double v = 0;
                    for (int nH = 0; nH < lH.Count; nH++)
                    {
                        Harmonic h = lH[nH];
                        if (h.Gain > 0) v += h.Gain * Math.Sin((nH * n * mul) + h.Phase);
                    }
                    Sample s = new Sample(nc);
                    for (int c = 0; c < nc; c++)
                    {
                        s[c] = v;
                    }
                    n++;
                    samples.Add(s);
                }
                while (true)
                {
                    foreach (Sample s in samples)
                    {
                        yield return s;
                    }
                }
            }
        }
    }

    public class BandLimitedSquareGenerator : HarmonicGenerator
    {
        public BandLimitedSquareGenerator(ushort numChannels, uint sampleRate, double freq, double gain)
            : base(numChannels, sampleRate, freq, gain)
        {
        }

        public override List<Harmonic> Harmonics
        {
            get
            {
                List<Harmonic> harmonics = new List<Harmonic>();
                uint sr = SampleRate;
                int nH = 0;
                double g = _gain * 4 / Math.PI;
                while (true)
                {
                    double fN = _freq * (nH + 1);
                    if (fN > sr / 2)
                    {
                        // Above Nyquist; stop
                        break;
                    }
                    Harmonic h = new Harmonic(nH % 2 == 0 ? 0 : g / nH);
                    harmonics.Add(h);
                    nH++;
                }
                return harmonics;
            }
        }
    }

    public class TriangleGenerator : SignalGenerator
    {
        public TriangleGenerator(ushort numChannels, uint sampleRate, double freq, double gain)
            : base(numChannels, sampleRate, freq, gain)
        {
        }

        /// <summary>
        /// Get an iterator for samples
        /// </summary>
        public override IEnumerator<ISample> Samples
        {
            get
            {
                ushort nc = NumChannels;
                uint sr = SampleRate;
                double mul = 2 * _freq / sr;
                int n = 0;
                Sample s = new Sample(nc);
                while (true)
                {
                    double saw = ((n * mul) % 2);
                    double v = 2 * saw;
                    if (v > 1) v = 2 - v;
                    if (v < -1) v = -2 - v;
                    for (int c = 0; c < nc; c++)
                    {
                        s[c] = _gain * v;
                    }
                    n++;
                    yield return s;
                }
            }
        }
    }

    public class BandLimitedTriangleGenerator : HarmonicGenerator
    {
        public BandLimitedTriangleGenerator(ushort numChannels, uint sampleRate, double freq, double gain)
            : base(numChannels, sampleRate, freq, gain)
        {
        }

        public override List<Harmonic> Harmonics
        {
            get
            {
                List<Harmonic> harmonics = new List<Harmonic>();
                uint sr = SampleRate;
                int nH = 0;
                double g = _gain * 8 / (Math.PI * Math.PI);
                while (true)
                {
                    double fN = _freq * (nH + 1);
                    if (fN > sr / 2)
                    {
                        // Above Nyquist; stop
                        break;
                    }
                    Harmonic h = new Harmonic(nH % 2 == 0 ? 0 : g / (nH * nH), nH % 4 == 1 ? Math.PI : 0);
                    harmonics.Add(h);
                    nH++;
                }
                return harmonics;
            }
        }
    }

    public class SawtoothGenerator : SignalGenerator
    {
        public SawtoothGenerator(ushort numChannels, uint sampleRate, double freq, double gain)
            : base(numChannels, sampleRate, freq, gain)
        {
        }

        /// <summary>
        /// Get an iterator for samples
        /// </summary>
        public override IEnumerator<ISample> Samples
        {
            get
            {
                ushort nc = NumChannels;
                uint sr = SampleRate;
                double mul = 2 * _freq / sr;
                int n = 0;
                Sample s = new Sample(nc);
                while (true)
                {
                    double saw = ((n * mul) % 2) - 1;
                    double v = _gain * saw;
                    for (int c = 0; c < nc; c++)
                    {
                        s[c] = v;
                    }
                    n++;
                    yield return s;
                }
            }
        }
    }

    public class BandLimitedSawtoothGenerator : HarmonicGenerator
    {
        public BandLimitedSawtoothGenerator(ushort numChannels, uint sampleRate, double freq, double gain)
            : base(numChannels, sampleRate, freq, gain)
        {
        }

        public override List<Harmonic> Harmonics
        {
            get
            {
                List<Harmonic> harmonics = new List<Harmonic>();
                uint sr = SampleRate;
                int nH = 1;
                double g = _gain * 2 / Math.PI;
                harmonics.Add(new Harmonic(0,0));
                while (true)
                {
                    double fN = _freq * (nH + 1);
                    if (fN > sr / 2)
                    {
                        // Above Nyquist; stop
                        break;
                    }
                    Harmonic h = new Harmonic(g / nH, Math.PI);
                    harmonics.Add(h);
                    nH++;
                }
                return harmonics;
            }
        }
    }

}
