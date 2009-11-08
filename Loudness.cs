using System;
using System.Collections.Generic;
using System.Text;

namespace DSPUtil
{
    /*
% from iso226.m:
% Generates an Equal Loudness Contour as described in ISO 226
%
% Usage:  [SPL FREQ] = ISO226(PHON);
% 
%         PHON is the phon value in dB SPL that you want the equal
%           loudness curve to represent. (1phon = 1dB @ 1kHz)
%         SPL is the Sound Pressure Level amplitude returned for
%           each of the 29 frequencies evaluated by ISO226.
%         FREQ is the returned vector of frequencies that ISO226
%           evaluates to generate the contour.
%
% Desc:   This function will return the equal loudness contour for
%         your desired phon level.  The frequencies evaulated in this
%         function only span from 20Hz - 12.5kHz, and only 29 selective
%         frequencies are covered.  This is the limitation of the ISO
%         standard.
%
%         In addition the valid phon range should be 0 - 90 dB SPL.
%         Values outside this range do not have experimental values
%         and their contours should be treated as inaccurate.
%
%         If more samples are required you should be able to easily
%         interpolate these values using spline().
%
% Author: Jeff Tackett 03/01/05
     */
    public class Loudness
    {
        static double[] f  = {   20,    25,  31.5,  40.0,  50.0,  63.0,  80.0, 100.0, 125.0, 160.0, 200.0, 250.0, 315.0, 400.0, 500.0, 630.0, 800.0, 1000, 1250, 1600, 2000, 2500, 3150, 4000, 5000, 6300, 8000, 10000, 12500};
        static double[] af = {0.532, 0.506, 0.480, 0.455, 0.432, 0.409, 0.387, 0.367, 0.349, 0.330, 0.315, 0.301, 0.288, 0.276, 0.267, 0.259, 0.253, 0.250, 0.246, 0.244, 0.243, 0.243, 0.243, 0.242, 0.242, 0.245, 0.254, 0.271, 0.301};
        static double[] Lu = {-31.6, -27.2, -23.0, -19.1, -15.9, -13.0, -10.3,  -8.1, -6.2, -4.5, -3.1, -2.0, -1.1, -0.4,  0.0,   0.3,   0.5,   0.0,  -2.7,  -4.1,  -1.0,   1.7,   2.5,   1.2, -2.1, -7.1, -11.2, -10.7, -3.1};
        static double[] Tf = { 78.5, 68.7, 59.5, 51.1, 44.0, 37.5, 31.5, 26.5, 22.1, 17.9, 14.4, 11.4, 8.6, 6.2, 4.4, 3.0, 2.2, 2.4, 3.5, 1.7, -1.3, -4.2, -6.0, -5.4, -1.5, 6.0, 12.6, 13.9, 12.3 };

        /// <summary>
        /// Create a list of dB SPL equal-loudness values for a given 'phon' loudness
        /// (from zero, threshold, to 90)
        /// </summary>
        /// <param name="phon"></param>
        /// <returns>list of {frequency Hz, dB SPL}</returns>
        public static FilterProfile SPL(double phon)
        {
            FilterProfile lfg = new FilterProfile();
            if ((phon < 0) | (phon > 120))
            {
                throw new ArgumentException("Phon value out of bounds!");
            }
            // Setup user-defined values for equation
            double Ln = phon;

            for (int j = 0; j < f.Length; j++)
            {
                // Deriving sound pressure level from loudness level (iso226 sect 4.1)
                double Af = 4.47E-3 * Math.Pow(10, (0.025 * Ln) - 1.15) + Math.Pow(0.4 * Math.Pow(10, (((Tf[j] + Lu[j]) / 10) - 9)), af[j]);
                double Lp = ((10 / af[j]) * Math.Log10(Af)) - Lu[j] + 94;

                // Return user data
                FreqGain fg = new FreqGain(f[j], Lp);
                lfg.Add(fg);
            }
            return lfg;
        }

        public static FilterProfile DifferentialSPL(double phon0, double phon1)
        {
            return DifferentialSPL(phon0, phon1, 1.0);
        }

        public static FilterProfile DifferentialSPL(double phon0, double phon1, double scale)
        {
            FilterProfile spl = new FilterProfile();
            FilterProfile spl0 = Loudness.SPL(phon0);
            FilterProfile spl1 = Loudness.SPL(phon1);
            for (int j = 0; j < spl1.Count; j++)
            {
                FreqGain fg = spl1[j];
                fg.Gain = scale * ( spl0[j].Gain - fg.Gain );
                spl.Add(fg);
            }
            return spl;
        }


        /// <summary>
        /// Calculate the weighted volume of a *single channel* sound source.
        /// NB: this consumes lots of memory for long sources.
        /// </summary>
        /// <param name="src"></param>
        /// <param name="dbSPL"></param>
        /// <returns>Volume (units, not dB)</returns>
        public static double WeightedVolume1(ISoundObj src, double dbSPL, double dbSPLBase)
        {
            if (src.NumChannels != 1)
            {
                throw new ArgumentException("Requires single-channel");
            }

            // Read channel into a buffer
            SoundBuffer buff = new SoundBuffer(src);
            buff.ReadAll();

            // And then double in length to prevent wraparound
            buff.PadTo(buff.Count * 2);
            // Pad to next higher power of two
            buff.PadToPowerOfTwo();
            int n = buff.Count;

            double wvImpulse = WeightedVolume2(buff, dbSPL, dbSPLBase);

            // compare with a Dirac pulse the same length
            CallbackSource dirac = new CallbackSource(1, src.SampleRate, delegate(long j)
            {
                if (j >= n)
                {
                    return null;
                }
                double v = 0;
                if (j == n / 2)
                {
                    v = 1;
                }
                return new Sample(v);
            });
            buff = new SoundBuffer(dirac);
            buff.ReadAll();
            double wvDirac = WeightedVolume2(buff, dbSPL, dbSPLBase);

            buff = null;
            GC.Collect();

            return wvImpulse / wvDirac;
        }

        private static double WeightedVolume2(SoundBuffer src, double dbSPL, double dbSPLBase)
        {
            double v = 0;
            uint sr = src.SampleRate;

            // Read buffer into array of complex
            Complex[][] data = src.ToComplexArray();

            // We only have a single channel
            Complex[] cdata = data[0];

            // FFT in place
            Fourier.FFT(cdata.Length, cdata);

            // Calculate magnitude, weighted by 80-phon loudness, for each loudness band.
            // These are the ISO measured points:
            FilterProfile lfg;
            if (dbSPLBase == 0)
            {
                lfg = SPL(dbSPL);
            }
            else
            {
                lfg = DifferentialSPL(dbSPL, dbSPLBase);
            }
//          lfg.Add(new FreqGain(sr / 2, lfg[lfg.Count - 1].Gain));

            // Cover the ISO measured range (only...)
            int nStart = (int)(lfg[0].Freq * (long)cdata.Length / sr);
            int nEnd = (int)(lfg[lfg.Count - 1].Freq * (long)cdata.Length / sr);

            // Just use linear interpolation (on a dB scale; linear freq scale) of gain between each measured point
            int nfg = 0;

            int startp = nStart;
            int endp = (int)(lfg[nfg + 1].Freq * (long)cdata.Length / sr);     // endpoint of this band
            double dB1 = lfg[nfg].Gain;         // SPL of the ISO223 curve at this freq
            double dB2 = lfg[nfg+1].Gain;       // ...and the next point

            double vThisBand = 0;
            int nThisBand = 0;
            for (int j = nStart; j < nEnd; j++)
            {
                if (j > endp)
                {
                    if (nThisBand > 0) v += Math.Sqrt(vThisBand / nThisBand); // RMS
                    while (j >= endp)
                    {
                        nfg++;
                        startp = j;
                        endp = (int)(lfg[nfg + 1].Freq * (long)cdata.Length / sr);
                        dB1 = lfg[nfg].Gain;
                        dB2 = lfg[nfg + 1].Gain;
                    }
                    vThisBand = 0;
                    nThisBand = 0;
                }
                Complex c = cdata[j];
                double dbHere = dB1 + ((dB2 - dB1) * (double)(j - startp) / (double)(endp - startp));
                vThisBand += (c.Re * c.Re) / MathUtil.gain(dbHere);
                nThisBand++;
            }
            if(nThisBand>0) v += Math.Sqrt(vThisBand / nThisBand);

            return v;
        }

        /// <summary>
        /// Calculate the weighted volume of a sound source.
        /// NB: this consumes lots of memory for long sources.
        /// </summary>
        /// <param name="src"></param>
        /// <param name="dbSPL"></param>
        /// <returns>Volume (units, not dB)</returns>
        public static double WeightedVolume(ISoundObj src, double dbSPL, double dbSPLBase)
        {
            double wv = 0;
            for (ushort c = 0; c < src.NumChannels; c++)
            {
                SingleChannel channel = src.Channel(c);
                wv += Loudness.WeightedVolume1(channel, dbSPL, dbSPLBase);
            }
            src.Reset();
            wv = wv / src.NumChannels;
            return wv;
        }

        /// <summary>
        /// Calculate the weighted volume of a sound source.
        /// NB: this consumes lots of memory for long sources.
        /// </summary>
        /// <param name="src"></param>
        /// <returns>Volume (units, not dB)</returns>
        public static double WeightedVolume(ISoundObj src)
        {
            return WeightedVolume(src, 40, 0);
        }
    }
}
