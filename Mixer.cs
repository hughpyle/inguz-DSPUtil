using System;
using System.Collections.Generic;
using System.Text;

// Copyright (c) 2006 by Hugh Pyle, inguzaudio.com

namespace DSPUtil
{
    [Serializable]
    public class Mixer : SoundObj
    {
        List<ISoundObj> _inputs = new List<ISoundObj>();
        List<double> _gains = new List<double>();

        public Mixer()
        {
        }

        // Add another source.  Note: gain is units, not db
        public void Add(ISoundObj input, double gainUnits)
        {
            if (_inputs.Count == 0)
            {
                // Treat this as 'Input'...
                Input = input;
            }
            if (_inputs.Count > 0)
            {
                if (input.NumChannels != NumChannels)
                {
                    throw new Exception("Mixer inputs must have the same number of channels.");
                }
            }
            _inputs.Add(input);
            _gains.Add(gainUnits);
        }

        public override int Iterations
        {
            get
            {
                int i = 0;
                foreach (ISoundObj input in _inputs)
                {
                    i = Math.Max(i, input.Iterations);
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
                // Get enumerators for each input
                List<IEnumerator<ISample>> enums = new List<IEnumerator<ISample>>();
                List<bool> mores = new List<bool>();

                foreach (ISoundObj input in _inputs)
                {
                    enums.Add(input.Samples);
                    mores.Add(true);
                }

                bool anymore = true;
                while (anymore)
                {
                    Sample sample = new Sample(NumChannels);
                    int e = 0;
                    anymore = false;
                    foreach (IEnumerator<ISample> src in enums)
                    {
                        ISample s;
                        if (mores[e])
                        {
                            mores[e] = src.MoveNext();
                        }
                        if (mores[e])
                        {
                            s = src.Current;
                        }
                        else
                        {
                            s = new Sample(_inputs[e].NumChannels);
                        }
                        anymore |= mores[e];

                        for (int k = 0; k < s.NumChannels; k++)
                        {
                            sample[k] += (s[k] * _gains[e]);
                        }
                        e++;
                    }
                    yield return sample;
                }
            }
        }
    }
}
