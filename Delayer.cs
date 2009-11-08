using System;
using System.Collections.Generic;
using System.Text;

// Copyright (c) 2009 by Hugh Pyle, inguzaudio.com

namespace DSPUtil
{
    public class Delayer : SoundObj
    {
        private int _delay;
        private CircularBuffer _buff;
        private bool _tail;

        /// <summary>
        /// Create a delay processor
        /// </summary>
        public Delayer()
        {
            Delay = 0;
            _tail = true;
        }

        /// <summary>
        /// Create a delay processor
        /// </summary>
        /// <param name="returnTail">True to return all the samples (false to truncate the delay-tail)</param>
        public Delayer(bool returnTail)
        {
            Delay = 0;
            _tail = returnTail;
        }

        public override void Reset()
        {
            Delay = _delay;
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

                _buff.Input = _input;

                // Return the delayed
                foreach (ISample sample in _buff)
                {
                    ISample s = _buff[-_delay];
                    yield return s;
                }
                if (_tail)
                {
                    for (int j = 0; j < _delay; j++)
                    {
                        yield return _buff[j-_delay];
                    }
                }
            }
        }


        /// <summary>
        /// Time difference, in samples, to apply between left and right channels.
        /// Positive values delay the right channel relative to the left channel, and vice versa.
        /// </summary>
        public int Delay
        {
            get { return _delay; }
            set
            {
                int d = value;
                if (_buff == null || _delay < d)
                {
                    // Recreate the cache (even if it means losing data)
                    _buff = new CircularBuffer(_input, (uint)(d+1));
                }
                _delay = d;
            }
        }
    }
}
