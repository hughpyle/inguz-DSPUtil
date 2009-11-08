using System;
using System.Collections.Generic;
using System.Text;

namespace DSPUtil
{
    /// <summary>
    /// Channel-shuffle (LR to MS or vice versa)
    /// Requires two-channel input.
    /// </summary>
    public class Shuffler: SoundObj
    {
        public Shuffler()
        {
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
                // Return (a+b)/2 and (a-b)/2
                foreach (ISample sample in _input)
                {
                    yield return _next(sample);
                }
            }
        }

        internal ISample _next(ISample s)
        {
            if (_nc != 2)
            {
                // NOP!
            }
            else
            {
                double L = s[0];
                double R = s[1];

                // Set sample values in-place (quicker than newing another sample)
                s[0] = (L + R) * _sigmaGain;
                s[1] = (L - R) * _deltaGain;
            }
            return s;
        }

        private double _deltaGain = MathUtil.SQRT2;

        public double DeltaGain
        {
            get { return _deltaGain; }
            set { _deltaGain = value; }
        }

        private double _sigmaGain = MathUtil.SQRT2;

        public double SigmaGain
        {
            get { return _sigmaGain; }
            set { _sigmaGain = value; }
        }
	
    }

    // This is symmetrical if you shuffle twice,
    // L, R --> (L+R)/2, (L-R)/2
    // (L+R), (L-R) --> (L/2+R/2+L/2-R/2), (L/2+R/2-L/2+R/2) = L, R
}
