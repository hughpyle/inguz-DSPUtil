using System;
using System.Collections.Generic;
using System.Text;

// Copyright (c) 2006, 2007 by Hugh Pyle, inguzaudio.com

namespace DSPUtil
{
    public class Invertor : SoundObj
    {
        private bool[] _invert;

        public Invertor()
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

                foreach (ISample sample in Input)
                {
                    yield return _nxt(sample);
                }
            }
        }

        internal ISample _nxt(ISample s)
        {
            if (_invert != null)
            {
                for (ushort c = 0; c < s.NumChannels; c++)
                {
                    if (_invert[c])
                    {
                        s[c] = -s[c];
                    }
                }
            }
            return s;
        }

        /// <summary>
        /// Time difference, in samples, to apply between left and right channels.
        /// Positive values delay the right channel relative to the left channel, and vice versa.
        /// </summary>
        public void Invert(ushort nChannel, bool doInvert)
        {
            if (_invert == null)
            {
                // Lazy lazy.  In case input isn't defined yet.  Invertor just handles a max 16 channels.
                _invert = new bool[16];
            }
            if (nChannel > _invert.Length)
            {
                throw new ArgumentOutOfRangeException("nChannel");
            }
            _invert[nChannel] = doInvert;
        }
    }
}
