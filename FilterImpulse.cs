using System;
using System.Collections.Generic;
using System.Text;

// Copyright (c) 2006, 2007 by Hugh Pyle, inguzaudio.com

namespace DSPUtil
{
    // Constructs an impulse for parametric EQ filter

    [Serializable]
    public struct FreqGain
    {
        private double _freq;
        /// <summary>
        /// Frequency, Hz
        /// </summary>
        public double Freq
        {
            get { return _freq; }
            set { _freq = value; }
        }

        private double _gain;
        /// <summary>
        /// Gain, dB
        /// </summary>
        public double Gain
        {
            get { return _gain; }
            set { _gain = value; }
        }

        /// <summary>
        /// Frequency/Gain pair
        /// </summary>
        /// <param name="freq">Frequency, Hz</param>
        /// <param name="gain">Gain, dB</param>
        public FreqGain(double freq, double gain)
        {
            _freq = freq;
            _gain = gain;
        }
    }

    public class FilterProfile : List<FreqGain>
    {
        public FilterProfile()
            : base()
        {
        }

        /// <summary>
        /// Copy constructor
        /// </summary>
        /// <param name="lfg"></param>
        public FilterProfile(FilterProfile lfg)
            : base(lfg)
        {
        }

        /// <summary>
        /// Compute the filter profile of an impulse.
        /// NOTE: very memory-intensive!
        /// </summary>
        /// <param name="impulse">The impulse file</param>
        /// <param name="fractionsOfERB">1.0 for ERB-spaced bands. Smaller for less smoothing</param>
        public FilterProfile(ISoundObj impulse, double fractionsOfERB)
            : base()
        {
            /*
            double[] muff = magbands(impulse, 300);

            // smooth the muff
            double[] smoo = ERB.smooth(muff, (int)(38 / fractionsOfERB));

            // sample frequency response at ERB band centers
            FilterProfile lfg = ERB.profile(smoo, impulse.SampleRate, fractionsOfERB);
            */
            FilterProfile lfg = Smoothing.Profile(impulse, SmoothingType.OCTAVE, fractionsOfERB);
            AddRange(lfg);
        }


        private static double[] magbands(ISoundObj impulse, double bins)
        {
            uint nSR = impulse.SampleRate;
            uint nSR2 = nSR / 2;

            int nn = (int)bins + 1;
            double[] muff = new double[nn];

            ushort nChannels = impulse.NumChannels;
            for (ushort c = 0; c < nChannels; c++)
            {
                // Read channel into a buffer
                SingleChannel channel = impulse.Channel(c);
                SoundBuffer buff = new SoundBuffer(channel);
                buff.ReadAll();

                // And then double in length to prevent wraparound
                buff.PadTo(buff.Count * 2);
                // Pad to next higher power of two
                buff.PadToPowerOfTwo();
                // Read out into array of complex
                Complex[][] data = buff.ToComplexArray();
                Complex[] cdata = data[0];

                // Then we're done with the buffer for this channel
                buff = null;
                GC.Collect();

                // FFT in place
                Fourier.FFT(cdata.Length, cdata);

                int n = cdata.Length / 2;

                // Drop the FFT magnitudes into the 'muff' array
                // according to an ERB-based scale (near-logarithmic).
                // Then smoothing is easy.
                double binw = (nSR2 / (double)n);
                int prevbin = 0;
                int nbin = 0;
                double v = 0;
                for (int j = 0; j < n; j++)
                {
                    double f = (double)j * binw;    // equiv freq, Hz
                    int bin = (int)ERB.f2bin(f, nSR, bins); // the bin we drop this sample in
                    v += cdata[j].Magnitude;
                    nbin++;

                    if ((bin > prevbin) || (j == n - 1))
                    {
                        muff[prevbin] += (v / nbin);
                        v = 0;
                        nbin = 0;
                        prevbin = bin;
                    }
                }
            }

            // Now muff is sum(all channels) of average-magnitude-per-bin.
            // Divide it all by the number of channels, so our gains are averaged...
            for (int j = 0; j < muff.Length; j++)
            {
                muff[j] = muff[j] / nChannels;
            }

            return muff;
        }


        /// <summary>
        /// Find the maximum gain.
        /// </summary>
        /// <returns></returns>
        public double MaxGain()
        {
            double g = double.MinValue;
            foreach (FreqGain fg in this)
            {
                if (fg.Gain > g)
                {
                    g = fg.Gain;
                }
            }
            return g;
        }

        /// <summary>
        /// Find the minimum gain.
        /// </summary>
        /// <returns></returns>
        public double MinGain()
        {
            double g = double.MaxValue;
            foreach (FreqGain fg in this)
            {
                if (fg.Gain < g)
                {
                    g = fg.Gain;
                }
            }
            return g;
        }
        
        /// <summary>
        /// Create a subset of the profile, within the given frequency range.
        /// </summary>
        /// <param name="fMin"></param>
        /// <param name="fMax"></param>
        /// <returns></returns>
        public FilterProfile FreqRange(double fMin, double fMax)
        {
            FilterProfile lfg = new FilterProfile();
            foreach (FreqGain fg in this)
            {
                if (fg.Freq >= fMin && fg.Freq <= fMax)
                {
                    lfg.Add(fg);
                }
            }
            return lfg;
        }


        /// <summary>
        /// Invert the profile.
        /// </summary>
        /// <returns></returns>
        public FilterProfile Inverse()
        {
            return Inverse(0);
        }

        public FilterProfile Inverse(double maxGain)
        {
            FilterProfile lfg = new FilterProfile();
            foreach (FreqGain fg in this)
            {
                double newGain = -fg.Gain;
                if (maxGain != 0)
                {
                    // Don't make more than +/-maxGain dB difference!
                    newGain = Math.Min(newGain, maxGain);
                    newGain = Math.Max(newGain, -maxGain);
                    lfg.Add(new FreqGain(fg.Freq, newGain));
                }
            }
            return lfg;
        }

        /// <summary>
        /// Multiply the profile by the given gain factor.
        /// </summary>
        /// <param name="scale"></param>
        /// <returns></returns>
        public FilterProfile Scale(double scale)
        {
            FilterProfile lfg = new FilterProfile();
            for (int j = 0; j < lfg.Count; j++)
            {
                FreqGain fg = lfg[j];
                lfg.Add(new FreqGain(fg.Freq, fg.Gain * scale));
            }
            return lfg;
        }

        public static FilterProfile operator *(FilterProfile c1, double num)
        {
            return c1.Scale(num);
        }

        public static FilterProfile operator *(double num, FilterProfile c1)
        {
            return c1.Scale(num);
        }

        /// <summary>
        /// Subtract one profile from another.
        /// Only valid if they share the same frequency set.
        /// </summary>
        /// <param name="c1"></param>
        /// <param name="c2"></param>
        /// <returns></returns>
        public static FilterProfile operator -(FilterProfile c1, FilterProfile c2)
        {
            if (c1.Count != c2.Count)
            {
                return null;
            }
            FilterProfile newProfile = new FilterProfile();
            for (int j = 0; j < c1.Count; j++)
            {
                if (c1[j].Freq != c2[j].Freq)
                {
                    return null;
                }
                double g1 = c1[j].Gain;
                double g2 = c2[j].Gain;
                newProfile.Add(new FreqGain(c1[j].Freq, g1 - g2));
            }
            return newProfile;
        }

        /// <summary>
        /// Convert to a JSON string
        /// </summary>
        /// <param name="name"></param>
        /// <param name="description"></param>
        /// <returns></returns>
        public string ToJSONString(string name, string description)
        {
            StringBuilder sb = new StringBuilder();
            sb.AppendLine("\"" + name + "\": \"" + description + "\",");
            sb.Append("\"" + name + "_loop\":[ ");
            bool first = true;
            foreach (FreqGain fg in this)
            {
                double gain = fg.Gain;
                if (!Double.IsNaN(gain) && !Double.IsInfinity(gain))
                {
                    if (first)
                    {
                        first = false;
                    }
                    else
                    {
                        sb.Append(", ");
                    }
                    sb.Append("{\"");
                    sb.AppendFormat("{0:f2}", fg.Freq);
                    sb.Append("\":");
                    sb.AppendFormat("{0:f2}", gain);
                    sb.Append("}");
                }
            }
            sb.AppendLine(" ]");
            return sb.ToString();
        }
    }

    // Comparer, compares frequency portion only
    public class FreqGainComparer : Comparer<FreqGain>
    {
        public override int Compare(FreqGain x, FreqGain y)
        {
            return x.Freq.CompareTo(y.Freq);
        }
    }

    public enum FilterInterpolation
    {
        COSINE = 0,
        SPLINE = 1
    }

    [Serializable]
    public class FilterImpulse : SoundObj
    {
        int _nSamples;
        FilterProfile _coeffs;
        bool _allZero;
        FilterInterpolation _int;

        /// <summary>
        /// Create a linear-phase filter impulse.
        /// </summary>
        /// <param name="nSize">Size: determines length of the final filter - suggest 4096 or 8192 (or 0 to auto-set)</param>
        /// <param name="coefficients">List of {frequency Hz, gain dB}</param>
        /// <param name="interpolation">type of interpolation (cosine)</param>
        public FilterImpulse(int nSize, FilterProfile coefficients, FilterInterpolation interpolation, uint sampleRate)
        {
            if (sampleRate == 0)
            {
                throw new Exception("Sample rate cannot be zero");
            }

            _nSamples = nSize;
            _coeffs = coefficients;
            _int = interpolation;
            SampleRate = sampleRate;

            // Check the coefficients are sorted and non-overlapping
            _coeffs.Sort(new FreqGainComparer());
            FilterProfile co = new FilterProfile();
            _allZero = true;
            double prevFreq = -1;
            double freqDiff = double.MaxValue;
            foreach (FreqGain fg in _coeffs)
            {
                if (fg.Freq > prevFreq && fg.Freq >= 0)
                {
                    co.Add(fg);
                    freqDiff = Math.Min(freqDiff, fg.Freq - prevFreq);
                    prevFreq = fg.Freq;
                }
                if (fg.Gain != 0)
                {
                    _allZero = false;
                }
            }
            _coeffs = co;

            if (_nSamples == 0)
            {
                // Make up a length from 4k to 32k, based on the closest-together frequencies
                _nSamples = _allZero ? 1024 : Math.Max(4096, (int)Math.Min(32768.0, 327680 / freqDiff));
            }
            _nSamples = MathUtil.NextPowerOfTwo(_nSamples);
        }

#region Overrides
        /// <summary> Number of iterations expected to do the signal processing </summary>
        public override int Iterations
        {
            get { return _nSamples; }
        }

        public override ushort NumChannels
        {
            get
            {
                ushort nc = base.NumChannels;
                if (nc == 0) return 1;
                return nc;
            }
            set
            {
                base.NumChannels = value;
            }
        }

#endregion

        public override IEnumerator<ISample> Samples
        {
            get
            {
                int j;

                // Make a Dirac impulse
                Complex[] data = new Complex[_nSamples];
                int centerpoint = (_nSamples / 2) -1;
                data[centerpoint] = new Complex(1, 0);

                Complex[][] ddata = new Complex[1][];
                ddata[0] = data;

                ISoundObj cbr;
                if (!_allZero)
                {
                    // FFT it
                    Fourier.FFT((int)_nSamples, data);

                    // Make a buffer from the FFT data
                    cbr = new ComplexBufferReader(ddata, 1, 0, _nSamples / 2, ComplexBufferFlags.Both);
                    ISoundObj shape = cbr;

                    cbr = new CallbackSource(1, SampleRate, delegate(long i){
                        if (i > _nSamples / 2)
                        {
                            return null;
                        }
                        return new Sample(1.0);
                    });

                    if (_int == FilterInterpolation.COSINE)
                    {
                        // Construct a series of shelf envelopes over the FFT data to match the coefficients,
                        // with each shelf "section" feeding from the previous.
                        // So for example
                        // (-12@60, -3*1k,+12@15k) gives
                        // shelf(0, 15) <- shelf(-12, -3)  <- cbr
                        LogBasisShelfEnvelope rcsePrev = null;
                        for (j = 1; j < _coeffs.Count; j++)
                        {
                            double freq1 = _coeffs[j - 1].Freq;
                            double freq2 = _coeffs[j].Freq;

                            // gains in dB
                            double gain1 = _coeffs[j - 1].Gain;
                            double gain2 = _coeffs[j].Gain;

                            // Since we're chaining these filters together,
                            // we actually want to adjust for pre-sweep == 0dB, then end-sweep=(difference)
                            double gainEnd = gain2 - gain1;

                            int startSample = (int)Math.Round(freq1 * _nSamples / SampleRate);
                            int endSample = (int)Math.Round(freq2 * _nSamples / SampleRate);
                            if (endSample <= startSample)
                            {
                                endSample = startSample + 1;
                            }

                            LogBasisShelfEnvelope rcse = new LogBasisShelfEnvelope((j == 1) ? gain1 : 0, (j == 1) ? gain2 : gainEnd, startSample, endSample);

                            // Input of this shelf is output of the previous shelf
                            // (or, for the first shelf, the sweep)
                            if (rcsePrev == null)
                            {
                                rcse.Input = cbr;
                            }
                            else
                            {
                                rcse.Input = rcsePrev;
                            }
                            rcsePrev = rcse;
                        }
                        shape = (rcsePrev == null) ? cbr : rcsePrev as ISoundObj;
                    }
                    else if (_int == FilterInterpolation.SPLINE)
                    {
                        throw new Exception("Not implemented");
                    }

                    // Apply the shelf filters to the data
                    j = 0;
                    foreach (ISample sample in shape)
                    {
                        if (j + j < _nSamples)
                        {
                            double val = sample[0];
                            data[j].mul(val); // = new Complex(sample[0], sample[1]);
                            data[_nSamples - j - 1].mul(val); // = new Complex(sample[0], sample[1]);
                        }
                        j++;
                    }

                    // IFFT to get the filter's impulse response
                    Fourier.IFFT((int)_nSamples, data);

                    // Scale to unity gain (ignoring the imaginary portion)
                    double mul = 2; // Math.Log(_nSamples / 4);
                    for (j = 0; j < _nSamples; j++)
                    {
                        data[j].mul(mul);
                    }
                }

                // Make a buffer of the real portion
                cbr = new ComplexBufferReader(ddata, 1, 0, _nSamples, ComplexBufferFlags.RealOnly);

                // Window to attenuate the extremities
                // SoundObj envelope = new NormalWindow(centerpoint, _nSamples / 14);
//                SoundObj envelope = new BlackmanHarris(centerpoint, centerpoint);
//                envelope.Input = cbr;

                ushort nc = NumChannels;
                if (nc <= 1)
                {
                    // Yield samples from the windowed impulse
                    foreach (ISample sample in cbr /* envelope */)
                    {
                        yield return sample;
                    }
                }
                else
                {
                    foreach (ISample sample in cbr /* envelope */)
                    {
                        ISample s = (nc == 2) ? new Sample2() : new Sample(nc) as ISample;
                        for (int c = 0; c < nc; c++)
                        {
                            s[c] = sample[0];
                        }
                        yield return s;
                    }
                }
            }
        }
    }
}
