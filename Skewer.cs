using System;
using System.Collections.Generic;
using System.Text;

// Copyright (c) 2006, 2007 by Hugh Pyle, inguzaudio.com

namespace DSPUtil
{
    public class Skewer : SoundObj
    {
        private int _skew;
        private int _s;
        private double[] _data;
        private bool _tail;

        /// <summary>
        /// Create a skew processor for a stereo input
        /// </summary>
        /// <param name="returnTail">True to return all the samples (false to truncate the skew-tail)</param>
        public Skewer(bool returnTail)
        {
            Skew = 0;
            _tail = returnTail;
        }

        public override void Reset()
        {
            Skew = _skew;
            base.Reset();
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

                // Return skewed sample.  Positive skew means delay the right-hand channel.
                // eg if skew=1,
                //  output(t=0) = input(t=0,nul)
                //  output(t=1) = input(t=1,t=0) ...etc...
                // To do this keep a circular buffer _data[] of the channel which is being delayed
                // (right if _skew>0, left if _skew<0)
                // with ptr n
                int n = 0;
                foreach (ISample sample in Input)
                {
                    yield return _nxt(ref n, sample);
                }
                if (_tail)
                {
                    for (int j = 0; j < _s; j++)
                    {
                        yield return _nxt(ref n, new Sample(_nc));
                    }
                }
            }
        }

        ISample _nxt(ref int n, ISample s)
        {
            if (_s == 0 || s.NumChannels != 2)
            {
                return s;
            }
            else
            {
                double L = s[0];
                double R = s[1];
                if (_skew < 0)
                {
                    _data[n % _s] = L;
                    n = (n + 1) % _s;
                    L = _data[n];
                }
                else
                {
                    _data[n % _s] = R;
                    n = (n + 1) % _s;
                    R = _data[n];
                }
                s[0] = L;
                s[1] = R;
                return s;
            }
        }

        /// <summary>
        /// Time difference, in samples, to apply between left and right channels.
        /// Positive values delay the right channel relative to the left channel, and vice versa.
        /// </summary>
        public int Skew
        {
            get { return _skew; }
            set
            {
                _skew = value;
                _s = Math.Abs(value) + 1;
                if (_data==null || _data.Length < _s)
                {
                    // Recreate the skewer cache (even if it means losing data)
                    _data = new double[_s];
                }
            }
        }	
    }
}
