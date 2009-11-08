using System;
using System.Collections.Generic;
using System.Text;

namespace DSPUtil
{
    public delegate ISample SourceCallback(long n);

    /// <summary>
    /// Wraps a callback function to become a SoundObj
    /// Example: new CallbackSource( 1, 44100, delegate(int j) { ... return sample ... } );
    /// </summary>
    public class CallbackSource : SoundObj
    {
        private SourceCallback _fn;

        public CallbackSource()
        {
        }

        public CallbackSource(ushort nChannels, uint sampleRate, SourceCallback fn)
        {
            NumChannels = nChannels;
            SampleRate = sampleRate;
            _fn = fn;
        }

        public override IEnumerator<ISample> Samples
        {
            get
            {
                int n = 0;
                while (true)
                {
                    ISample s = _fn(n);
                    if (s == null)
                    {
                        yield break;
                    }
                    yield return s;
                    n++;
                }
            }
        }

    }

}
