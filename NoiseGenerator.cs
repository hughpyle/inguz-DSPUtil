using System;
using System.Collections.Generic;
using System.Text;

// Copyright (c) 2006 by Hugh Pyle, inguzaudio.com

namespace DSPUtil
{
    public enum NoiseType
    {
        SILENCE = 0,
        WHITE = 1,
        PINK = 2,
        EWEIGHTED = 3,
        WHITE_FLAT = 4,
        ARBITRARY = 5,
        DIRAC = 6
    }

    [Serializable]
    public class NoiseGenerator : SoundObj
    {
        const double twopi = 2 * Math.PI;

        NoiseType _type;
        int _lengthSamples;
        double _gain;
        bool _mono;

        public NoiseGenerator(NoiseType type, ushort numChannels, double lengthSecs, uint sampleRate, double gain, bool mono)
        {
            init(type, numChannels, (int)(sampleRate * lengthSecs), sampleRate, gain, mono);
        }

        public NoiseGenerator(NoiseType type, ushort numChannels, int lengthSamples, uint sampleRate, double gain, bool mono)
        {
            init(type, numChannels, lengthSamples, sampleRate, gain, mono);
        }

        private void init(NoiseType type, ushort numChannels, int lengthSamples, uint sampleRate, double gain, bool mono)
        {
            _type = type;
            _lengthSamples = lengthSamples;
            NumChannels = numChannels;
            SampleRate = sampleRate;
            _gain = gain;
            _mono = mono;
            if (_gain == 0) _gain = 1;
            _random = new Random();
        }

        // List of coefficients, only used for "arbitrary" noise
        private FilterProfile _coeffs;
        public FilterProfile Coefficients
        {
            get { return _coeffs; }
            set { _coeffs = value; }
        }


        // Cached random-number generator
        private static Random _random = new Random();
        // Random value from zero to 1
        private static double NextRandom()
        {
            lock (_random)
            {
                return _random.NextDouble();
            }
        }
        // Random value from -1 to +1
        private static double NextRandom2()
        {
            lock (_random)
            {
                return 2 * _random.NextDouble() - 1;
            }
        }

        private static IEnumerator<double> Silence
        {
            get
            {
                while (true)
                {
                    yield return 0;
                }
            }
        }

        private static IEnumerator<double> White
        {
            get
            {
                while (true)
                {
                    yield return NextRandom2();
                }
            }
        }

        private static IEnumerator<double> Pink
        {
            get
            {
                double[] b = new double[7];
                while (true)
                {
                    /*
                    Filter to make pink noise from white
                    ------------------------------------
                    This is an approximation to a -10dB/decade filter using a weighted sum
                    of first order filters. It is accurate to within +/-0.05dB above 9.2Hz 
                    (44100Hz sampling rate). Unity gain is at Nyquist, but can be adjusted
                    by scaling the numbers at the end of each line.

                    If 'white' consists of uniform random numbers, such as those generated
                    by the rand() function, 'pink' will have an almost gaussian level 
                    distribution.

                    An 'economy' version with accuracy of +/-0.5dB is also available.

                      b0 = 0.99765 * b0 + white * 0.0990460;
                      b1 = 0.96300 * b1 + white * 0.2965164;
                      b2 = 0.57000 * b2 + white * 1.0526913;
                      pink = b0 + b1 + b2 + white * 0.1848;
                    ---
                    paul.kellett@maxim.abel.co.uk
                    http://www.abel.co.uk/~maxim/
                    
                    (hfp note: this is bitrate-independent, although the low-freq limit will change)
                    */
                    double white = NextRandom2();
                    b[0] = 0.99886 * b[0] + white * 0.0555179;
                    b[1] = 0.99332 * b[1] + white * 0.0750759;
                    b[2] = 0.96900 * b[2] + white * 0.1538520;
                    b[3] = 0.86650 * b[3] + white * 0.3104856;
                    b[4] = 0.55000 * b[4] + white * 0.5329522;
                    b[5] = -0.7616 * b[5] - white * 0.0168980;
                    double pink = b[0] + b[1] + b[2] + b[3] + b[4] + b[5] + b[6] + white * 0.5362;
                    b[6] = white * 0.115926;

                    yield return pink / 5;
                }
            }
        }

        private static IEnumerator<double> EWeighted
        {
            get
            {
                // Also used in the dither filter, but here we return E-weighted noise instead of adding E-weight dither

                // These Coefficients from "Psychoacoustically Optimal Noise Shaping," by Robert A. Wannamaker
                // (From AES Journal Volume 40, No. 7/8, 1992 July/August)

                // !!!NOTE: This FIR is only good for 44.1kHz sampling rate!!!

                double[] _shaped = {  2.412f, -3.370f, 3.937f, -4.174f, 3.353f, -2.205f, 1.281f, -0.569f, 0.0847f};
                int ORDER = 9;
                int _HistPos = ORDER-1;
                double[] _EH = new double[2 * ORDER];
                while (true)
                {
                    double samp = NextRandom2();

                    // Convolve history with the FIR...
                    for (int x = 0; x < ORDER; x++)
                    {
                        samp += _shaped[x] * _EH[_HistPos + x];
                    }

                    _HistPos--;
                    if (_HistPos < 0) _HistPos += ORDER;

                    // Update buffer (both copies)
                    _EH[_HistPos + 9] = _EH[_HistPos] = samp;

                    yield return samp / 16;
                }
            }
        }

        private static IEnumerator<double> WhiteFlat(int length)
        {
            // Generate "truly white" noise (flat spectrum)
            // by generating random-phase flat-magnitude in the frequency domain, then IFFT.
            // This noise is periodic - length (next power of 2 from 'n')

            int n = MathUtil.NextPowerOfTwo(length);
            Complex[] data = new Complex[n];

            double logn = Math.Log(n);
            for (int j = 0; j < n; j++)
            {
                // Create a random phase value from 0 to 2pi
                double phi = NextRandom() * 2 * Math.PI;
                // Magnitude is 1, so just trig
                double re = Math.Cos(phi);
                double im = Math.Sin(phi);
                data[j] = new Complex(logn*re,logn*im);
            }

            // IFFT
            Fourier.IFFT((int)n, data);

            // Return the real component
            int k = 0;
            while (true)
            {
                yield return data[k].Re * logn;
                k++;
                if (k >= n)
                {
                    k = 0;
                }
            }
        }

        private static IEnumerator<double> Arbitrary(int length, FilterProfile coeffs, uint sampleRate)
        {
            // Generate random-phase with specified magnitudes in the frequency domain, then IFFT.
            // This noise is periodic - length (next power of 2 from 'n')

            int l = MathUtil.NextPowerOfTwo(length);
            Complex[] data = new Complex[l];

            int n = 0;
            double freq1 = coeffs[n].Freq * 2 * l / sampleRate;
            double freq2 = coeffs[n + 1].Freq * 2 * l / sampleRate;
            double gain1 = coeffs[n].Gain;
            double gain2 = coeffs[n + 1].Gain;

            double logn = Math.Log(l);
            for (int j = 0; j < l; j++)
            {
                double gainDb;
                double gain;
                if (j > freq2)
                {
                    // Move to the next coefficient
                    n++;
                    if (n < coeffs.Count-1)
                    {
                        freq1 = coeffs[n].Freq * 2 * l / sampleRate;
                        freq2 = coeffs[n + 1].Freq * 2 * l / sampleRate;
                        gain1 = coeffs[n].Gain;
                        gain2 = coeffs[n + 1].Gain;
                    }
                }
                if (j < freq1)
                {
                    gainDb = gain1;
                }
                else if (j > freq2)
                {
                    gainDb = gain2;
                }
                else
                {
                    // Raised Cosine: 0.5* ( cos(phi) + 1 ), from phi=pi to 2pi
                    // 
                    double frac = (double)(j - freq1) / (double)(freq2 - freq1);
                    double ph = Math.PI * (1 + frac);
                    double rcos = (1 + Math.Cos(ph)) / 2;
                    gainDb = gain1 + rcos * (gain2 - gain1);
                }
                gain = MathUtil.gain(gainDb);

                // Create a random phase value from 0 to 2pi
                double phi = NextRandom() * 2 * Math.PI;
                // Magnitude is 1, so just trig
                double re = Math.Cos(phi);
                double im = Math.Sin(phi);
                data[j] = new Complex(gain * logn * re, gain * logn * im);
            }

            // IFFT
            Fourier.IFFT((int)l, data);

            // Return the real component
            int k = 0;
            while (true)
            {
                yield return data[k].Re * logn;
                k++;
                if (k >= l)
                {
                    k = 0;
                }
            }
        }

        private static IEnumerator<double> Dirac
        {
            get
            {
                yield return 1;
                while (true)
                {
                    yield return 0;
                }
            }
        }

        private IEnumerator<double> Noise
        {
            get
            {
                switch (_type)
                {
                    case NoiseType.SILENCE:
                        return Silence;
                    case NoiseType.WHITE:
                        return White;
                    case NoiseType.PINK:
                        return Pink;
                    case NoiseType.EWEIGHTED:
                        return EWeighted;
                    case NoiseType.WHITE_FLAT:
                        return WhiteFlat(_lengthSamples);
                    case NoiseType.ARBITRARY:
                        return Arbitrary(_lengthSamples,_coeffs,SampleRate);
                    case NoiseType.DIRAC:
                        return Dirac;
                    default:
                        throw new ArgumentException("Unknown noise type");
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
                ushort nc = _mono ? (ushort)1 : NumChannels;
                IEnumerator<double>[] noise = new IEnumerator<double>[nc];
                for (int c = 0; c < nc; c++)
                {
                    noise[c] = Noise;
                    noise[c].MoveNext();
                }
                for (int j = 0; j < _lengthSamples; j++)
                {
                    Sample sample = new Sample(NumChannels);
                    if (_mono)
                    {
                        double v = _gain * noise[0].Current;
                        noise[0].MoveNext();
                        for (int c = 0; c < NumChannels; c++)
                        {
                            sample[c] = v;
                        }
                    }
                    else
                    {
                        for (int c = 0; c < NumChannels; c++)
                        {
                            sample[c] = _gain * noise[c].Current;
                            noise[c].MoveNext();
                        }
                    }
                    yield return sample;
                }
            }
        }

        /// <summary> Number of iterations expected to do the signal processing </summary>
        public override int Iterations
        {
            get { return (_lengthSamples); }
        }

    }
}
