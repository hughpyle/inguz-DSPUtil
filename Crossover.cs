using System;
using System.Collections.Generic;
using System.Text;

// Copyright (c) 2006, 2007 by Hugh Pyle, inguzaudio.com

namespace DSPUtil
{
    // TEMPORARY
    // until needed for ambisonics

    // Now some synthetic filters, in matched pairs, that approximate the simple IIR responses
    public class FilterPair
    {
        uint _sr;
        double _cornerLo;
        double _cornerHi;
        public FilterPair(double corner1, double corner2, uint sampleRate, int n)
        {
            _sr = sampleRate;
            _cornerLo = corner1;
            _cornerHi = corner2;
            // If cornerLo is at zero frequency: this is a low-pass filter
            // If cornerHi is at Nyquist: this is a high=pass filter
            if (_cornerLo >= _cornerHi && _cornerHi != 0)
            {
                throw new NotSupportedException("Overlapping corner frequencies!");
            }
        }

        /// <summary>
        /// Returns the filter magnitude at a given frequency.
        /// </summary>
        /// <param name="freq">Frequency, Hz</param>
        /// <returns>Magnitude, units (NOT dB)</returns>
        public double FilterMagnitude(double freq)
        {
            // 1st order
            // low pass: f = 1/(1+s)     = 1 / ( 1 + j(w / 2.pi.Fc) )
            // highpass: f = 1/(1+(1/s)) = 1 / ( 1 + j(2.pi.Fc / w) )

            // High pass from _cornerLo
            Complex c1, c2;
            if (freq == 0)
            {
                c1 = new Complex(0, 0);
            }
            else
            {
                c1 = (new Complex(1, 0) / new Complex(1, (_cornerLo) / freq));
            }

            // low pass to cornerHi
            if (_cornerHi == 0)
            {
                c2 = new Complex(0, 0);
            }
            else
            {
                c2 = (new Complex(1, 0) / new Complex(1, freq / (_cornerHi)));
            }

            double mag = (c1 + c2).Magnitude;
            return mag;
        }

        /// <summary>
        /// Returns the inverse filter magnitude at a given frequency.
        /// </summary>
        /// <param name="freq">Frequency, Hz</param>
        /// <returns>Magnitude, units (NOT dB)</returns>
        public double InverseMagnitude(double freq)
        {
            return 1 - FilterMagnitude(freq);
        }

        public ISoundObj Filter()
        {
            return new CallbackSource(1, _sr, delegate(long n)
            {
                return new Sample(0.0);
            });
        }
    }
}
