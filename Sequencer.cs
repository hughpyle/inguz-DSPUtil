using System;
using System.Collections.Generic;
using System.Text;

// Copyright (c) 2006, 2007 by Hugh Pyle, inguzaudio.com

namespace DSPUtil
{
    /// <summary>
    /// Sequencer: join inputs together sequentially.
    /// Create a Sequencer object, then add as many inputs as necessary...
    /// </summary>
    [Serializable]
    public class Sequencer : SoundObj
    {
        List<ISoundObj> _inputs = new List<ISoundObj>();
        List<List<double>> _channelGains = new List<List<double>>();

        public Sequencer()
        {
        }

        // Add another source
        public void Add(ISoundObj input)
        {
            Add(input, new List<double>());
        }

        public void Add(ISoundObj input, List<double> channelGains)
        {
            if (_inputs.Count == 0)
            {
                // Treat this as 'Input'...
                Input = input;
            }
            _inputs.Add(input);
            _channelGains.Add(channelGains);
        }

        public override int Iterations
        {
            get
            {
                int i = 0;
                foreach (ISoundObj input in _inputs)
                {
                    i += input.Iterations;
                }
                return i;
            }
        }

        /// <summary>
        /// Get an iterator for samples
        /// </summary>
        public override IEnumerator<ISample> Samples
        {
            get
            {
                int nIn = 0;
                foreach (ISoundObj input in _inputs)
                {
                    foreach (ISample s in input)
                    {
                        if (_channelGains[nIn].Count > 0)
                        {
                            for(int c=0; c<_channelGains[nIn].Count; c++)
                            {
                                s[c] *= _channelGains[nIn][c];
                            }
                        }
                        yield return s;
                    }
                    nIn++;
                }
            }
        }
    }

}
