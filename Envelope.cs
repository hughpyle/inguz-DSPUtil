using System;
using System.Collections.Generic;
using System.Text;

// Copyright (c) 2006, 2007 by Hugh Pyle, inguzaudio.com

namespace DSPUtil
{
    [Serializable]
    public class LinearEnvelope : SoundObj
    {
        // Linear envelope:
        // apply a linear-time linear-gain change from startGain (dB) to endGain (dB)
        // over the course of nSamples

        private double _startGain;
        private double _endGain;
        private int _nSamples;

        public LinearEnvelope(double dBstartGain, double dBendGain, int nSamples)
        {
            _startGain = MathUtil.gain(dBstartGain);
            _endGain = MathUtil.gain(dBendGain);
            _nSamples = nSamples;
        }

        /// <summary>
        /// Get an iterator for samples
        /// </summary>
        public override IEnumerator<ISample> Samples
        {
            get
            {
                if (_input == null)
                {
                    yield break;
                }
                int n = 0;
                // Multiply each input sample
                // by the instantaneous gain we want
                foreach (ISample sample in Input)
                {
                    double frac = (double)n / (double)_nSamples;
                    double gain = _startGain + frac * (_endGain - _startGain);
                    Sample s = new Sample( sample, gain );
                    n++;
                    yield return s;
                }
            }
        }
    }

    [Serializable]
    public class LinearDbEnvelope : SoundObj
    {
        // Linear envelope:
        // apply a linear-time log-gain change from startGain (dB) to endGain (dB)
        // over the course of nSamples

        private double _startGain;
        private double _endGain;
        private int _nSamples;

        public LinearDbEnvelope(double dBstartGain, double dBendGain, int nSamples)
        {
            _startGain = dBstartGain;
            _endGain = dBendGain;
            _nSamples = nSamples;
        }

        /// <summary>
        /// Get an iterator for samples
        /// </summary>
        public override IEnumerator<ISample> Samples
        {
            get
            {
                if (_input == null)
                {
                    yield break;
                }
                int n = 0;
                // Multiply each input sample
                // by the instantaneous gain we want
                foreach (ISample sample in Input)
                {
                    double frac = (double)n / (double)_nSamples;
                    double gain = _startGain + frac * (_endGain - _startGain);
                    Sample s = new Sample(sample, MathUtil.gain(gain));;
                    n++;
                    yield return s;
                }
            }
        }
    }


    [Serializable]
    public class CosWindow : SoundObj
    {
        // A cosine-based window function
        // centered on sample # 'center'
        // returns the full input stream, muted outside the windowed region

        private int _center;
        private int _sidewidth;
        private int _halfplateauwidth;
        private double _c0;
        private double _c1;
        private double _c2;
        private double _c3;

        public CosWindow(int center, int sidewidth, int halfplateauwidth, double c0, double c1, double c2, double c3)
        {
            _center = center;
            _sidewidth = sidewidth;
            _halfplateauwidth = halfplateauwidth;
            _c0 = c0;
            _c1 = c1;
            _c2 = c2;
            _c3 = c3;
        }

        public IEnumerator<double> Gains
        {
            get
            {
                int n = 0;
                int nstart = _center - _halfplateauwidth - _sidewidth;
                int nend = _center + _halfplateauwidth + _sidewidth;
                int plateaustart = _center - _halfplateauwidth;
                int plateauend = _center + _halfplateauwidth;

                while (true)
                {
                    double gain = 0;
                    if (n < nstart)
                    {
                        //gain = 0;
                    }
                    else if (n > nend)
                    {
                        //gain = 0;
                    }
                    else if (n > plateaustart && n < plateauend)
                    {
                        gain = 1;
                    }
                    else
                    {
                        // Raised Cosine: 0.5* ( cos(phi) + 1 ), from phi=pi to 2pi
                        // 
                        long nn = (long)n - _center;
                        if (nn < 0) { nn += (long)_halfplateauwidth; } else { nn -= (long)_halfplateauwidth; }
                        long N = (long)_sidewidth;
                        if (Math.Abs(nn) < 2 * N)
                        {
                            double phi = Math.PI * ((double)nn / N);
                            gain = _c0 + (_c1 * Math.Cos(phi)) + (_c2 * Math.Cos(2 * phi)) + (_c3 * Math.Cos(3 * phi));
                        }
                    }
                    n++;
                    yield return gain;
                }
            }
        }

        /// <summary>
        /// Get an iterator for samples
        /// </summary>
        public override IEnumerator<ISample> Samples
        {
            get
            {
                if (_input == null)
                {
                    yield break;
                }

                IEnumerator<double> gains = Gains;
                foreach (ISample sample in _input)
                {
                    gains.MoveNext();
                    double gain = gains.Current;
                    Sample s = new Sample(sample, gain);
                    yield return s;
                }
            }
        }
    }

    // The usual types of cosine-based windows, centered, with no additional plateau
    // and acting as filters, i.e. only returning the portion of the stream which lies within the window
    [Serializable]
    public class Hann : CosWindow
    {
        public Hann(int center, int sidewidth)
            : base(center, sidewidth, 0, 0.5, 0.5, 0, 0)
        {
        }
    }

    [Serializable]
    public class Hamming : CosWindow
    {
        public Hamming(int center, int sidewidth)
            : base(center, sidewidth, 0, 0.54, 0.46, 0, 0)
        {
        }

        public Hamming(int center, int sidewidth, int halfplateauwidth)
            : base(center, sidewidth, halfplateauwidth, 0.54, 0.46, 0, 0)
        {
        }
    }

    [Serializable]
    public class Blackman : CosWindow
    {
        public Blackman(int center, int sidewidth)
            : base(center, sidewidth, 0, 0.42, 0.5, 0.08, 0)
        {
        }

        public Blackman(int center, int sidewidth, int halfplateauwidth)
            : base(center, sidewidth, halfplateauwidth, 0.42, 0.5, 0.08, 0)
        {
        }
    }

    [Serializable]
    public class BlackmanHarris : CosWindow
    {
        public BlackmanHarris(int center, int sidewidth)
            : base(center, sidewidth, 0, 0.35875, 0.48829, 0.14128, 0.01168)
        {
        }

        public BlackmanHarris(int center, int sidewidth, int halfplateauwidth)
            : base(center, sidewidth, halfplateauwidth, 0.35875, 0.48829, 0.14128, 0.01168)
        {
        }
    }

    [Serializable]
    public class NormalWindow : SoundObj
    {
        private int _center;
        private int _plateau;
        private double _sigma;

        /// <summary>
        /// A normal distribution around the center, with std deviation of 'sigma' samples
        /// (e.g. use center=N/2, sigma=N/6 to cover most of the curve).
        /// The window hiehgt is set such that the center sample is unchanged (per most windows).
        /// rather than normalizing area under the curve.
        /// </summary>
        public NormalWindow(int center, double sigma)
        {
            _center = center;
            _plateau = 0;
            _sigma = sigma;
        }

        public NormalWindow(int center, int plateau, double sigma)
        {
            _center = center;
            _plateau = plateau;
            _sigma = sigma;
        }

        /// <summary>
        /// Get an iterator for samples
        /// </summary>
        public override IEnumerator<ISample> Samples
        {
            get
            {
                if (_input == null)
                {
                    yield break;
                }
                int n = 0;
//              double mul = 1/(_sigma * Math.Sqrt(2*Math.PI));
                double sig2 = 2 * _sigma * _sigma;
                foreach (ISample sample in Input)
                {
                    long nn = (long)_center - (long)n;
                    double gain = 1;
                    // tbd: plateau
                    gain = Math.Exp(-(nn * nn / sig2)); // * mul; to normalize
                    Sample s = new Sample(sample, gain);
                    n++;
                    yield return s;
                }
            }
        }

    }

    [Serializable]
    public class ShelfEnvelope : SoundObj
    {
        // "Raised sine shelf" envelope:
        // given
        //  start gain: dB
        //  end gain: dB
        //  start sample #
        //  end sample #
        // apply a raised-sine envelope between start and end.

        private double _startGain;
        private double _endGain;
        private int _startSample;
        private int _endSample;

        public ShelfEnvelope(double dBstartGain, double dBendGain, int startSample, int endSample)
        {
            _startGain = dBstartGain;
            _endGain = dBendGain;
            _startSample = startSample;
            _endSample = endSample;
        }

        /// <summary>
        /// Get an iterator for samples
        /// </summary>
        public override IEnumerator<ISample> Samples
        {
            get
            {
                if (_input == null)
                {
                    yield break;
                }
                int n = 0;
                // Multiply each input sample
                // by the instantaneous gain we want
                foreach (ISample sample in Input)
                {
                    double gain;
                    if (n < _startSample)
                    {
                        gain = _startGain;
                    }
                    else if (n > _endSample)
                    {
                        gain = _endGain;
                    }
                    else
                    {
                        // Raised Cosine: 0.5* ( cos(phi) + 1 ), from phi=pi to 2pi
                        // 
                        double frac = (double)(n - _startSample) / (double)(_endSample - _startSample);
                        double phi = Math.PI * (1 + frac);
                        double rcos = ( 1 + Math.Cos(phi) )/2;
                        gain = _startGain + rcos * (_endGain - _startGain);
                    }
                    Sample s = new Sample( sample, MathUtil.gain(gain) );
                    n++;
                    yield return s;
                }
            }
        }
    }


    [Serializable]
    public class LogBasisShelfEnvelope : SoundObj
    {
        // "Raised sine shelf" envelope:
        // given
        //  start gain: dB
        //  end gain: dB
        //  start sample #
        //  end sample #
        // apply a raised-sine envelope between start and end, with log scale
        // (e.g. suitable for use on frequency-domain data)

        private double _startGain;
        private double _endGain;
        private int _startSample;
        private int _endSample;

        public LogBasisShelfEnvelope(double dBstartGain, double dBendGain, int startSample, int endSample)
        {
            _startGain = dBstartGain;
            _endGain = dBendGain;
            _startSample = startSample;
            _endSample = endSample;
        }

        /// <summary>
        /// Get an iterator for samples
        /// </summary>
        public override IEnumerator<ISample> Samples
        {
            get
            {
                if (_input == null)
                {
                    yield break;
                }
                int n = 0;

                // Add 1 so we don't blow up when starting from zero
                double logStart = Math.Log(1+_startSample);
                double logEnd = Math.Log(1+_endSample);

                // Multiply each input sample
                // by the instantaneous gain we want
                foreach (ISample sample in Input)
                {
                    double gain;
                    if (n < _startSample)
                    {
                        gain = _startGain;
                    }
                    else if (n > _endSample)
                    {
                        gain = _endGain;
                    }
                    else
                    {
                        // Raised Cosine: 0.5* ( cos(phi) + 1 ), from phi=pi to 2pi
                        // 
                        // fraction 0 to 1, linear basis
//                      double frac = (double)(n - _startSample) / (double)(_endSample - _startSample);
                        // fraction 0 to 1, log basis
                        double frac = (Math.Log(n + 1) - logStart) / (logEnd - logStart);
                        double phi = Math.PI * (1 + frac);
                        double rcos = (1 + Math.Cos(phi)) / 2;
                        gain = _startGain + rcos * (_endGain - _startGain);
                    }
                    Sample s = new Sample( sample, MathUtil.gain(gain) );
                    n++;
                    yield return s;
                }
            }
        }
    }

}
