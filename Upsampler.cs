using System;
using System.Collections.Generic;
using System.Text;

/*
namespace DSPUtil
{
    /// <summary>
    /// A very basic up-sampler.
    /// Includes gain.
    /// </summary>
    public class Upsampler:SoundObj
    {
        public Upsampler()
        {
        }
        public override int Iterations
        {
            get
            {
                return (int)(base.Iterations * SampleRate / _input.SampleRate);
            }
        }
        public override IEnumerator<ISample> Samples
        {
            get
            {
                if (_input == null)
                {
                    throw new Exception("Null input");
                }
                uint srIn = _input.SampleRate;
                uint srOut = SampleRate;
                uint lcm = srOut / srIn; //  MathUtil.gcd(srIn, srOut);
                ushort nc = _input.NumChannels;
                foreach (ISample s in _input)
                {
                    for (int j = 0; j < nc; j++)
                    {
                        s[j] *= lcm;
                    }
                    yield return s;
                    for (int j = 1; j < lcm; j++)
                    {
                        yield return (nc == 2) ? new Sample2() : new Sample(nc) as ISample;
                    }
                }
            }
        }
    }


    /// <summary>
    /// A very basic down-sampler.
    /// </summary>
    public class Downsampler : SoundObj
    {
        public Downsampler()
        {
        }
        public override int Iterations
        {
            get
            {
                return (int)(base.Iterations * _input.SampleRate / SampleRate);
            }
        }
        public override IEnumerator<ISample> Samples
        {
            get
            {
                if (_input == null)
                {
                    throw new Exception("Null input");
                }
                uint srIn = _input.SampleRate;
                uint srOut = SampleRate;
                uint lcm = srIn / srOut; // MathUtil.gcd(srIn, srOut);
                uint j = lcm;
                foreach (ISample s in _input)
                {
                    if (j == 0)
                    {
                        yield return s;
                        j = lcm;
                    }
                    j--;
                }
            }
        }
    }
}
*/