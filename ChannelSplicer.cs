using System;
using System.Collections.Generic;
using System.Text;

// Copyright (c) 2006, 2007 by Hugh Pyle, inguzaudio.com

// Various things for dealing with multiple channels

namespace DSPUtil
{
    /// <summary>
    /// SingleChannel: extract a single channel from the input.
    /// </summary>
    [Serializable]
    public class SingleChannel : SoundObj
    {
        ushort _nChannel;
        bool _okIfChannelNotFound = false;

        public SingleChannel(ISoundObj input, ushort nChannel)
        {
            Input = input;
            _nChannel = nChannel;
        }

        public SingleChannel(ISoundObj input, ushort nChannel, bool okIfChannelNotFound)
        {
            Input = input;
            _nChannel = nChannel;
            _okIfChannelNotFound = okIfChannelNotFound;
        }

        public override ushort NumChannels
        {
            get
            {
                return 1;
            }
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
                    throw new Exception("No input");
                }
                ushort nc = _input.NumChannels;
                bool no = false;
                if (nc < _nChannel)
                {
                    no = true;
                    if (!_okIfChannelNotFound)
                    {
                        throw new Exception(String.Format("Input does not have channel {0}", _nChannel));
                    }
                }

                foreach (ISample sample in _input)
                {
                    if (no)
                    {
                        yield return new Sample((ushort)1);
                    }
                    else if (nc == 1)
                    {
                        yield return sample;
                    }
                    else
                    {
                        yield return new Sample(sample[_nChannel]);
                    }
                }
            }
        }
    }


    /// <summary>
    /// TwoChannel: extract a pair of channels from the input.
    /// </summary>
    [Serializable]
    public class TwoChannel : SoundObj
    {
        ushort _nChannel1;
        ushort _nChannel2;

        public TwoChannel(ISoundObj input, ushort nChannel1, ushort nChannel2)
        {
            Input = input;
            _nChannel1 = nChannel1;
            _nChannel2 = nChannel2;
        }

        public override ushort NumChannels
        {
            get
            {
                return 2;
            }
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
                    throw new Exception("No input");
                }
                ushort nc = _input.NumChannels;
                if (nc < _nChannel1)
                {
                    throw new Exception(String.Format("Input does not have channel {0}", _nChannel1));
                }
                if (nc < _nChannel2)
                {
                    throw new Exception(String.Format("Input does not have channel {0}", _nChannel2));
                }

                foreach (ISample sample in _input)
                {
                    yield return new Sample2(sample[_nChannel1], sample[_nChannel2]);
                }
            }
        }
    }


    /// <summary>
    /// ChannelSplicer: join channels together.
    /// Create a ChannelSplicer object, then add as many inputs as necessary...
    /// </summary>
    [Serializable]
    public class ChannelSplicer : SoundObj
    {
        List<ISoundObj> _inputs = new List<ISoundObj>();
        List<double> _gains = new List<double>();
        ushort _nChannels = 0;
        bool _hasInput = false;

        public ChannelSplicer()
        {
        }

        // Add a null source (empty channel)
        public void Add()
        {
            Add(null, 0);
        }

        // Add another source
        public void Add(ISoundObj input)
        {
            Add(input, 1.0);
        }

        // Add another source
        public void Add(ISoundObj input, double gainUnits)
        {
            if (input != null && !_hasInput)
            {
                // Treat this as 'Input'...
                Input = input;
                _hasInput = true;
            }
            _inputs.Add(input);
            _gains.Add(gainUnits);
            _nChannels += (input == null) ? (ushort)1 : input.NumChannels;
        }

        public override ushort NumChannels
        {
            get
            {
                return _nChannels;
            }
        }

        public override int Iterations
        {
            get
            {
                int i = 0;
                foreach (ISoundObj input in _inputs)
                {
                    i = Math.Max(i, (input == null) ? 0 : input.Iterations);
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
                    if (input == null)
                    {
                        enums.Add(null);
                        mores.Add(false);
                    }
                    else
                    {
                        enums.Add(input.Samples);
                        mores.Add(true);
                    }
                }

                bool anymore = true;
                while (anymore)
                {
                    Sample sample = new Sample(NumChannels);
                    int j = 0;
                    int e = 0;
                    anymore = false;
                    foreach (IEnumerator<ISample> src in enums)
                    {
                        Sample s;
                        if (mores[e])
                        {
                            mores[e] = src.MoveNext();
                        }
                        if (mores[e])
                        {
                            s = (Sample)src.Current;
                        }
                        else
                        {
                            s = new Sample(_inputs[e]==null ? (ushort)1 : _inputs[e].NumChannels);
                        }
                        anymore |= mores[e];
                        double gain = _gains[e];
                        e++;

                        ushort nc = s.NumChannels;
                        for(int k=0; k<nc; k++)
                        {
                            sample[j] = s[k] * gain;
                            j++;
                        }
                    }
                    yield return sample;
                }
            }
        }
    }




    /// <summary>
    /// ChannelInvertor: invert any number of the channels.
    /// </summary>
    [Serializable]
    public class ChannelInvertor : SoundObj
    {
        bool[] _invert;

        public ChannelInvertor(ISoundObj input, params bool[] invert)
        {
            Input = input;
            _invert = invert;
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
                    throw new Exception("No input");
                }
                ushort nc = (ushort)Math.Min(_input.NumChannels, _invert.Length);

                foreach (ISample sample in _input)
                {
                    for (ushort c = 0; c < nc; c++)
                    {
                        if (_invert[c])
                        {
                            sample[c] = -sample[c];
                        }
                    }
                    yield return sample;
                }
            }
        }
    }




    /// <summary>
    /// ChannelSwapper: swaps channels in the input.
    /// </summary>
    [Serializable]
    public class ChannelSwapper: SoundObj
    {
        private ushort[] _channels;

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="input">The input</param>
        /// <param name="nChannels">A zero-based list of channels which will be mapped from the input channels.
        /// Use -1 to map null output.
        /// e.g  stereo input: params {1,0}: output is stereo with left and right swapped
        /// e.g. stereo input, params {-1,0,1,-1}: output is 4-channel with left in second channel and right in third.
        /// </param>
        public ChannelSwapper(ISoundObj input, params ushort[] channels)
        {
            Input = input;
            _channels = channels;
        }

        public override ushort NumChannels
        {
            get
            {
                return (ushort)_channels.Length;
            }
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
                    throw new Exception("No input");
                }
                ushort nc = _input.NumChannels;
                ushort nC = (ushort)_channels.Length;

                foreach (ISample sample in _input)
                {
                    ISample s = new Sample(nC);
                    for (ushort c = 0; c < nC; c++)
                    {
                        ushort cc = _channels[c];
                        if (cc >= 0 && cc<=nc)
                        {
                            s[c] = sample[cc];
                        }
                    }
                    yield return s;
                }
            }
        }
    }


}
