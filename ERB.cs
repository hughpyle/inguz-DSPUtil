using System;
using System.Collections.Generic;
using System.Text;

// Copyright (c) 2006, 2007 by Hugh Pyle, inguzaudio.com

namespace DSPUtil
{
    public enum SmoothingType
    {
        OCTAVE = 1,
        ERB
    }

    public class Smoothing
    {
        public static double WindowLength(int bin, int bins, double sampleRate, SmoothingType type, double resolution, out int bin0, out int bin1)
        {
            double len = 0;
            bin0 = bin;
            bin1 = bin;

            // Frequency (Hz) of the center of the given bin:
            // Bins range from 0 through sampleRate-1,
            // so the width of each bin is (sampleRate/bins) Hz,
            // and the center of bin <N> is (<N> * (sampleRate/bins)) + half a bin
            double binw = sampleRate / bins;
            double freq = (0.5 + bin) * binw;

            // The window goes from f0=(freq/X) to f1=(freq*X)
            // where the width of the window is (resolution)*(octaves or ERB bands)
            double bw = 0;
            if (type == SmoothingType.OCTAVE)
            {
                // Fraction X of an octave is f(x) = 2^x 
                bw = freq * Math.Pow(2, resolution);
            }
            else if(type==SmoothingType.ERB)
            {
                bw = ERB.ERBWidth(freq)*resolution;
            }

            // f.x - f/x = bw
            // f.x.x - bw.x - f = 0
            // f.x.(x-(bw/f)) = f
            // x.(x-(bw/f)) = 1
            // x.x - bw.x - 1/

            if (type == SmoothingType.OCTAVE)
            {
                // One octave is frequency*2.  Two is freq*4.
                // Fraction X of an octave is f(x) = 2^x 
//                double f1 = freq * (1 - Math.Pow(2, r2));
//                double f2 = freq * (1 - Math.Pow(2, r2));
            }

            return len * resolution;
        }
        
        public static FilterProfile Profile(ISoundObj impulse, SmoothingType type, double resolution)
        {
            uint nSR = impulse.SampleRate;
            uint nSR2 = nSR / 2;

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

                // Now we have an array of complex, from 0Hz to Nyquist and back again.
                // We really only care about the first half of the cdata buffer, but
                // treat it as circular anyway (i.e. wrap around for negative values).
                //
                // We're only working with magnitudes from here on,
                // so we can save some space by computing mags right away and storing them in the
                // real part of the complex array; then we can use the imaginary portion for the
                // smoothed data.
                for (int j = 0; j < cdata.Length; j++)
                {
                    cdata[j].Re = cdata[j].Magnitude;
                    cdata[j].Im = 0;
                }

                // Take a rectangular window of width (resolution)*(octave or ERB band)
                // Add up all magnitudes falling within this window
                //
                // Move the window forward by one thingummajig
                //double wMid = 0;    // center of the window
                //double wLen = 0;
            }
            return new FilterProfile(); // temp
        }
    }

    public class ERB
    {
        /*
         * Misc stuff for dealing with ERB scales.
         */


        /// <summary>
        /// Compute the equivalent rectangular bandwidth (ERB) at frequency f (Hz)
        /// </summary>
        /// <param name="f">Frequency (Hz)</param>
        /// <returns>Bandwidth (Hz)</returns>
        public static double ERBWidth(double f)
        {
            // from http://www.ling.su.se/STAFF/hartmut/bark.htm
            double erb = (6.23e-6 * (f * f)) + (9.339e-2 * f) + 28.52;
            return erb;
        }

        /// <summary>
        /// ERB value for a given frequency
        /// </summary>
        /// <param name="f">Frequency (Hz)</param>
        /// <returns></returns>
        public static double ERBVal(double f)
        {
            // from http://en.wikipedia.org/wiki/Equivalent_rectangular_bandwidth
            double i = 1 + (46.06538 * f / (f + 14678.49));
            double v = 11.17268 * Math.Log(i);
            return v;
        }

        public static double invERBVal(double v)
        {
            // from http://en.wikipedia.org/wiki/Equivalent_rectangular_bandwidth
            double f = (676170.4 / (47.06538 - Math.Exp(0.08950404 * v))) - 14678.49;
            return f;
        }

        /// <summary>
        /// Inverse of ERB(f)
        /// </summary>
        /// <param name="b">bandwidth</param>
        /// <returns>f (Hz)</returns>
        public static double invERB(double b)
        {
            const double c0 = 6.23e-6;
            const double c1 = 9.339e-2 / (2.0 * c0);
            const double c2 = (c1 * c1) * c0;
            double f = Math.Sqrt((b + c2 - 28.52) / c0) - c1;
            return f;
        }

        /// <summary>
        /// Integral of ERB from 0 to f
        /// </summary>
        /// <param name="f"></param>
        /// <returns></returns>
        public static double integralofERB(double f)
        {
            double ierb = 2.0766666666666665e-6 * f * (628.3299029068462 + f) * (21857.22386916378 + f);
            return ierb;
        }

        // Inverse function approximation by binary chop (I can' be bothered with the math...)
        public static double invintegralofERB(double i, uint sr)
        {
            return MathUtil.invert(integralofERB, 0, sr / 2, i);
        }


        private const double _bins = 1000;

        public static double f2bin(double f, uint sr, double bins)
        {
            double scale = (bins / ERB.ERBVal(sr / 2));
            return ERB.ERBVal(f) * scale;
        }
        public static double bin2f(double bin, uint sr, double bins)
        {
            double scale = (bins / ERB.ERBVal(sr / 2));
            return ERB.invERBVal(bin / scale);
        }

        public static double f2bin(double f, uint sr)
        {
            double scale = (_bins / ERB.ERBVal(sr / 2));
            return ERB.ERBVal(f) * scale;
        }
        public static double bin2f(double bin, uint sr)
        {
            double scale = (_bins / ERB.ERBVal(sr / 2));
            return ERB.invERBVal(bin / scale);
        }

        public static FilterProfile profile(double[] data, uint sr)
        {
            return profile(data, sr, 1.0);
        }

        /// <summary>
        /// Return a list of freq/gain, at the ERB centers
        /// Assuming you smoothed the data...
        /// </summary>
        /// <param name="data">data from ERB.smooth</param>
        /// <param name="sr">sample rate</param>
        /// <param name="scaleFactor">1.0 for ~38 bands.  0.5 for twice as many...</param>
        /// <returns></returns>
        public static FilterProfile profile(double[] data, uint sr, double scaleFactor)
        {
            FilterProfile pts = new FilterProfile();
            int dl = data.Length;
            for (double j = 1; j < ERB.ERBVal(sr / 2) + 1; j += scaleFactor)
            {
                double f = ERB.invERBVal(j);
                double n = f * 2 * dl / sr;
                if (n < dl)
                {
                    double g = data[(int)n];
                    pts.Add(new FreqGain(f, g));
                }
            }
            return pts;
        }
        
        /// <summary>
        /// Return a list of freq/gain, only including inflection points (where the curve is flat).
        /// </summary>
        /// <param name="data"></param>
        /// <param name="sr"></param>
        /// <returns></returns>
        public static FilterProfile inflections(double[] data, uint sr)
        {
            // Differentiate
            double bins = data.Length + 1;
            double[] diff = new double[data.Length];

            double n = data[0];
            for (int j = 0; j < data.Length; j++)
            {
                double d = data[j];
                diff[j] = d - n;
                n = d;
            }

            // Look for zero-crossings of the first derivative of data[]
            // (always include [0] and [end])
            FilterProfile pts = new FilterProfile();

            int bin;
            double pt = 0.1;
            int last = -1;
            double freq;
            double lastfreq = -1;

            // Always start with points for zero and 10 Hz
            pts.Add(new FreqGain(0, MathUtil.dB(data[0])));

            freq = 10;
            bin = (int)f2bin(freq,sr,bins);
            pts.Add(new FreqGain(freq, MathUtil.dB(data[bin])));
            pt = diff[bin];

            for (int j = bin + 1; j < data.Length; j++)
            {
                if ((pt > 0 && diff[j] <= 0) || (pt < 0 && diff[j] >= 0))
                {
                    freq = bin2f(j,sr,bins);
                    pts.Add(new FreqGain(freq, MathUtil.dB(data[j])));
                    last = j;
                    lastfreq = freq;
                }
                pt = diff[j];
            }
            // Fill in the last few target samples
            if (lastfreq < (sr / 2) - 2050)
            {
                freq = (sr / 2) - 2050;
                bin = (int)f2bin(freq,sr,bins);
                pts.Add(new FreqGain(freq, MathUtil.dB(data[bin])));
            }
            if (lastfreq < (sr / 2) - 1550)
            {
                freq = (sr / 2) - 1550;
                bin = (int)f2bin(freq,sr,bins);
                pts.Add(new FreqGain(freq, MathUtil.dB(data[bin])));
            }
            if (lastfreq < sr / 2)
            {
                freq = sr / 2;
                pts.Add(new FreqGain(freq, MathUtil.dB(data[data.Length - 1])));
            }

            return pts;
        }

        public static double[] smooth(double[] data, int bands)
        {
            int j, k;
            int nn = data.Length;
            double[] smoo = new double[nn];

            // Smooth over the critical band
            // Blackman-Harris window
            double c0 = 0.35875;
            double c1 = 0.48829;
            double c2 = 0.14128;
            double c3 = 0.01168;

            int windowSize;
            windowSize = (int)nn / bands;
            //            windowSize = (int)(bands * MathUtil.ERB(0));

            // Make a little array for the window coefficients
            double[] window = new double[2 * windowSize];
            double windowScale = 0;

            for (k = -windowSize; k < windowSize; k++)
            {
                double frac = (double)(k) / (double)(windowSize);
                double phi = Math.PI * frac;
                double rcos = 1 - (c0 - (c1 * Math.Cos(phi)) + (c2 * Math.Cos(2 * phi)) - (c3 * Math.Cos(3 * phi)));
                window[k + windowSize] = rcos;
                windowScale += rcos;
            }

            for (j = 0; j < nn; j++)
            {
                double v = 0;
                for (k = -windowSize; k < windowSize; k++)
                {
                    double rcos = window[k + windowSize];

                    int i = j + k;
                    if (i < 0) i = -i;
                    if (i >= nn) i = nn - (i % (nn - 1));

                    v += rcos * data[i];
                }
                smoo[j] = v / windowScale;
            }

            // Normalize the smoothed data
            double max = 0.0001;
            for (j = 0; j < nn; j++)
            {
                max = Math.Max(max, smoo[j]);
            }
            for (j = 0; j < nn; j++)
            {
                smoo[j] = smoo[j] / max;
            }

            return smoo;
        }

    }
}
