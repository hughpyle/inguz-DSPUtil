using System;
using System.Collections.Generic;
using System.Text;

// Copyright (c) 2006, 2007 by Hugh Pyle, inguzaudio.com

namespace DSPUtil
{
    public interface ISampleBuffer
    {
//      void Reset();
        void Skip(int n, out int nn, out bool moreSamples);
        ISample[] Read(int n, out int nn, out bool moreSamples);
        Complex[][] ReadComplex(int n, out int nn, out bool moreSamples);
    }

    // A padder.  Adds silence to the start of the input.
    // If the pad argument is negative, skips that number of samples from input.
    public class Padder : SoundObj
    {
        int _pad;
        public Padder(int pad)
        {
            _pad = pad;
        }
        public Padder(ISoundObj input, int pad)
        {
            Input = input;
            _pad = pad;
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
                Sample nul = new Sample(NumChannels);
                int j = 0;
                while(j < _pad)
                {
                    yield return nul;
                    j++;
                }
                foreach (ISample sample in _input)
                {
                    if (j >= -_pad)
                    {
                        yield return sample;
                    }
                    j++;
                }
            }
        }
    }

    // A trivial buffer which reads its input in batches (and can return in batches)
    // Doesn't support random access.
    [Serializable]
    public class SampleBuffer : SoundObj, ISampleBuffer
    {
        private ISampleBuffer _inputEnum;

        public SampleBuffer()
        {
            Reset();
        }
        public SampleBuffer(ISoundObj input)
        {
            Input = input;
            Reset();
        }

        private ISampleBuffer GetInputEnum()
        {
            if (_inputEnum == null)
            {
                _inputEnum = _input.GetBufferedEnumerator() as ISampleBuffer;
            }
            return _inputEnum;
        }

        public override void Reset()
        {
            _inputEnum = null;
        }

        public ISample[] Read(int n, out int nn, out bool moreSamples)
        {
            return GetInputEnum().Read(n, out nn, out moreSamples);
        }

        public Complex[][] ReadComplex(int n, out int nn, out bool moreSamples)
        {
            throw new NotImplementedException();
        }

        public void Skip(int n, out int nn, out bool moreSamples)
        {
            GetInputEnum().Skip(n, out nn, out moreSamples);
        }

        private int _length = int.MaxValue;
        public int Length
        {
            get { return _length==int.MaxValue ? 0 : _length; }
            set { _length = value; }
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
                int nn;
                int n = Math.Min(DSPUtil.BUFSIZE, _length);
                bool moreSamples = true;
                while (moreSamples)
                {
                    ISample[] tmp = GetInputEnum().Read(n, out nn, out moreSamples);
                    n -= nn;
                    for(int j=0; j<nn; j++)
                    {
                        yield return tmp[j];
                    }
                }
            }
        }

        /// <summary>
        /// Variant of Subset, supporting negative "start" values
        /// </summary>
        /// <param name="start"></param>
        /// <param name="count"></param>
        /// <returns></returns>
        public ISoundObj PaddedSubset(int start, int count)
        {
            bool started = false;
            ISampleBuffer sb = GetInputEnum();
            ushort nc = NumChannels;
            return new CallbackSource(nc, SampleRate, delegate(long n)
            {
                int nn;
                bool more;
                if (start < 0)
                {
                    start++;
                    return new Sample(nc);
                }
                if (!started)
                {
                    sb.Skip(start, out nn, out more);
                    started = true;
                }
                if (n > count)
                {
                    return null;
                }
                ISample[] sa = sb.Read(1, out nn, out more);
                if (sa.Length > 0)
                {
                    return sa[0];
                }
                else
                {
                    return null;
                }
            });
        }

        public ISoundObj Subset(int start)
        {
            return Subset(start, int.MaxValue - start);
        }

        public ISoundObj Subset(int start, int count)
        {
            bool started = false;
            ISampleBuffer sb = GetInputEnum();
            return new CallbackSource(NumChannels, SampleRate, delegate(long n)
            {
                int nn;
                bool more;
                if (!started)
                {
                    sb.Skip(start, out nn, out more);
                    started = true;
                }
                if (n > count)
                {
                    return null;
                }
                ISample[] sa = sb.Read(1, out nn, out more);
                if (sa.Length > 0)
                {
                    return sa[0];
                }
                else
                {
                    return null;
                }
            });
        }
    }


    // A lazy buffer which retains a copy of the entire stream; supports random access.
    // Example uses:
    // - Pull once from input, re-use the stream multiple times.
    //   - More efficient than re-creating the input enum,
    //   - exact same input if e.g. source is noise
    // - Pull only a subset of information from the source, e.g. for re-use or windowing
    //   - we don't need to read the entire thing

    [Serializable]
    public class SoundBuffer : SoundObj, ISampleBuffer
    {
        private List<ISample> _samples;
        private ISampleBuffer _inputEnum = null;
        private bool _moreSamples = true;   // there are samples in the source which we have not yet read
        private int _pos = 0;

        private ISample[] _buff;
        private int _bufflen;

        private Complex[][] _cbuff;
        private int _cbufflen;

        public SoundBuffer()
        {
            _samples = new List<ISample>(DSPUtil.BUFSIZE);
        }
        public SoundBuffer(ISoundObj input)
        {
            _samples = new List<ISample>(DSPUtil.BUFSIZE);
            Input = input;
        }
        public SoundBuffer(List<ISample> samples, ushort numChannels, uint sampleRate)
        {
            _samples = samples;
            Input = null;
            _nc = numChannels;
            _sr = sampleRate;
        }

        /// <summary>
        /// The number of samples currently in the buffer.
        /// This will often be less than the total number of samples available from the input.
        /// </summary>
        public int Count
        {
            get
            {
                return _samples.Count;
            }
        }

        public bool ReadTo(int n)
        {
            if (n < _samples.Count)
            {
                return true;
            }
            if (_input == null)
            {
                _moreSamples = false;
                return _moreSamples;
            }
            if (_inputEnum == null)
            {
                _inputEnum = _input.GetBufferedEnumerator() as ISampleBuffer;
            }
            int nn = n - _samples.Count;
            if (nn > 0 && _moreSamples)
            {
                if (n<int.MaxValue && _samples.Capacity < n) _samples.Capacity = n;
                int tot = 0;
                int nnn;
                while (_moreSamples && tot < nn)
                {
                    ISample[] tmp = _inputEnum.Read(DSPUtil.BUFSIZE, out nnn, out _moreSamples);
                    for (int j = 0; j < nnn; j++)
                    {
                        _samples.Add(tmp[j]);
                    }
                    tot += nnn;
                }
            }
            return nn < 0 || _moreSamples;
        }
        public int ReadAll()
        {
            ReadTo(int.MaxValue);
            return _samples.Count;
        }

        #region ISampleBuffer implementation
        public ISample[] Read(int n, out int nn, out bool moreSamples)
        {
            if (_buff == null || _bufflen < n)
            {
                _buff = new ISample[n];
                _bufflen = n;
            }
            int nnn = _pos;
            int n4 = nnn + n;
            bool more = ReadTo(n4);
            nn = (more ? n : _samples.Count - nnn);
//          Trace.WriteLine("CopyTo {0} {1} {2} {3}", _pos, nnn, nn, _samples.Count);
            _samples.CopyTo(nnn, _buff, 0, (int)nn);
            moreSamples = more;
            _pos += (int)nn;
            return _buff;
        }
        public Complex[][] ReadComplex(int n, out int nn, out bool moreSamples)
        {
            if (_cbuff == null || _cbufflen < n)
            {
                _cbuff = new Complex[_nc][];
                for (ushort c = 0; c < _nc; c++)
                {
                    _cbuff[c] = new Complex[n];
                }
                _cbufflen = n;
            }
            else
            {
                for (ushort c = 0; c < _nc; c++)
                {
                    Array.Clear(_cbuff[c], 0, n);
                }
            }

            int nnn = _pos;
            int n4 = nnn + n;
            bool more = ReadTo(n4);
            nn = (more ? n : _samples.Count - nnn);

            //          Trace.WriteLine("CopyTo {0} {1} {2} {3}", _pos, nnn, nn, _samples.Count);
            _samples.CopyTo(nnn, _buff, 0, (int)nn);
            for (int j = 0; j < nn; j++)
            {
                ISample s = _samples[nnn + j];
                for (ushort c = 0; c < _nc; c++)
                {
                    _cbuff[c][j].Re = s[c];
                }
            }

            moreSamples = more;
            _pos += (int)nn;

            return _cbuff;
        }

//        public void Reset()
//        {
//        }

        public void Skip(int n, out int nn, out bool moreSamples)
        {
            Read(n, out nn, out moreSamples);
        }
        #endregion

        // Operations

        private int _maxPos = int.MaxValue;
        private double _maxVal = 0;
        private void FindMax()
        {
            for (int j = 0; j < _samples.Count; j++)
            {
                ISample s = _samples[j];
                for (int c = 0; c < s.NumChannels; c++)
                {
                    double v = Math.Abs(s[c]);
                    if (v > _maxVal)
                    {
                        _maxVal = v;
                        _maxPos = j;
                    }
                }
            }
        }
        public int MaxPos()
        {
            if (_maxPos == int.MaxValue)
            {
                FindMax();
            }
            return _maxPos;
        }
        public double MaxVal()
        {
            if (_maxPos == int.MaxValue)
            {
                FindMax();
            }
            return _maxVal;
        }

        /// <summary>
        /// Normalize the buffer so maximum value is at (dBfs)
        /// </summary>
        /// <param name="dBfs">dB-of-fullscale for the new peak</param>
        /// <returns></returns>
        public double Normalize(double dBfs)
        {
            return Normalize(dBfs, true);
        }
        public double Normalize(double dBfs, bool doIt)
        {
            // Make sure the whole buffer is scaled
            double gain = 0;
            double max = 0;
            for (int j = 0; j < _samples.Count; j++)
            {
                ISample s = _samples[j];
                for (int c = 0; c < s.NumChannels; c++)
                {
                    max = Math.Max(max, Math.Abs(s[c]));
                }
            }
            if (max == 0)
            {
                return 0;
            }
            gain = MathUtil.gain(dBfs) / max;
            if (doIt)
            {
                ApplyGain(gain);
            }
            return gain;
        }

        /// <summary>
        /// Apply a gain (units) to the whole buffer.
        /// </summary>
        /// <param name="gain">units to scale by</param>
        public void ApplyGain(double gain)
        {
            for (int j = 0; j < _samples.Count; j++)
            {
                ISample s = _samples[j];
                _samples[j] = new Sample(s, gain);
            }
        }

        /// <summary>
        /// Apply a window to the whole buffer.
        /// </summary>
        /// <param name="window">Cosine-based window to apply</param>
        public void ApplyWindow(CosWindow window)
        {
            window.Input = new CallbackSource(NumChannels, SampleRate, delegate(long j)
            {
                if (j >= _samples.Count)
                {
                    return null;
                }
                return _samples[(int)j];
            });
            int n=0;
            foreach (Sample s in window)
            {
                _samples[n++] = s;
            }
        }

        public void PadTo(int n)
        {
            ReadTo(n);
            // Add null samples if we're past eof
            Sample s = new Sample(_input.NumChannels);
            while (_samples.Count < n)
            {
                _samples.Add(s);
            }
        }

        public void PadToPowerOfTwo()
        {
            // Ensure that the buffer's total length is a power of two.
            int n = MathUtil.NextPowerOfTwo(_samples.Count);
            PadTo(n);
        }

        /// <summary>
        /// Get an iterator for samples
        /// </summary>
        public override IEnumerator<ISample> Samples
        {
            get
            {
                _pos = 0;
                while(_pos < _samples.Count || _moreSamples)
                {
                    if (!ReadTo(_pos + 1))
                    {
                        yield break;
                    }
                    yield return _samples[(int)_pos];
                    _pos++;
                }
            }
        }

        /// <summary>
        /// Get an iterator for a subset of samples
        /// </summary>
        /// <param name="startPos">Start position (zero-based)</param>
        /// <param name="count">Number to return</param>
        /// <returns></returns>
        public IEnumerator<ISample> SampleSubset(int start, int count)
        {
            ReadTo(start + count);
            for(int j=start; j<start+count && j<_samples.Count; j++)
            {
                yield return _samples[j];
            }
        }


        public ISoundObj Subset(int start)
        {
            return Subset(start, int.MaxValue-start);
        }

        public ISoundObj Subset(int start, int count)
        {
            ReadTo(start + count);
            return new CallbackSource(NumChannels, SampleRate, delegate(long j)
            {
                if (j >= 0 && j < count && (j + start) < _samples.Count)
                {
                    return _samples[(int)(j + start)];
                }
                return null;
            });
        }

        // Array-index operator
        public ISample this[int arg]
        {
            get
            {
                ReadTo(arg);
                return _samples[arg];
            }
        }


        public ISample[] ToArray()
        {
            return ToArray(0, -1);
        }

        public ISample[] ToArray(int startPos)
        {
            return ToArray(startPos, -1);
        }

        public ISample[] ToArray(int startPos, int count)
        {
            int n = count;
            if (count < 0)
            {
                ReadAll();
                n = _samples.Count - startPos;
            }
            else
            {
                ReadTo(startPos + count);
                n = Math.Min(startPos + count, _samples.Count);
            }
            ISample[] ret = new ISample[n];
            for (int j = startPos; j - startPos < n; j++)
            {
                ret[j - startPos] = _samples[j];
            }
            return ret;
        }

        /// <summary>
        /// Return the buffer contents as an array of doubles
        /// </summary>
        public double[][] ToDoubleArray()
        {
            return ToDoubleArray(0, -1);
        }
        /// <summary>
        /// Return the buffer contents as an array of doubles
        /// </summary>
        /// <param name="startPos">Start position in the buffer</param>
        public double[][] ToDoubleArray(int startPos)
        {
            return ToDoubleArray(startPos, -1);
        }
        /// <summary>
        /// Return the buffer contents as an array of doubles
        /// </summary>
        /// <param name="startPos"></param>
        /// <param name="count">Length of the return array (will be null-padded if count exceeds data size)</param>
        /// <returns></returns>
        public double[][] ToDoubleArray(int startPos, int count)
        {
            int n = count;
            if (count < 0)
            {
                ReadAll();
                n = _samples.Count - startPos;
            }
            else
            {
                ReadTo(startPos + count);
                n = Math.Min(startPos + count, _samples.Count);
            }
            double[][] ret = new double[NumChannels][];
            for (int c = 0; c < NumChannels; c++)
            {
                ret[c] = new double[count<0 ? n : count];
            }
            for (int j = startPos; j - startPos < n; j++)
            {
                ISample s = _samples[j];
                for (int c = 0; c < NumChannels; c++)
                {
                    ret[c][j - startPos] = s[c];
                }
            }
            return ret;
        }


        public Complex[][] ToComplexArray()
        {
            return ToComplexArray(0, -1);
        }
        public Complex[][] ToComplexArray(int startPos)
        {
            return ToComplexArray(startPos, -1);
        }
        public Complex[][] ToComplexArray(int startPos, int count)
        {
            int n = count;
            if (count < 0)
            {
                ReadAll();
                n = _samples.Count - startPos;
            }
            else
            {
                ReadTo(startPos + count);
                n = Math.Min(count, _samples.Count - startPos);
            }
            Complex[][] ret = new Complex[NumChannels][];
            for (int c = 0; c < NumChannels; c++)
            {
                ret[c] = new Complex[n];
            }
            for (int j = startPos; j - startPos < n; j++)
            {
                ISample s = _samples[j];
                for (int c = 0; c < NumChannels; c++)
                {
                    ret[c][j - startPos] = new Complex(s[c],0);
                }
            }
            return ret;
        }
    }


    /// <summary>
    /// A circular buffer.  NB: only have one consumer of the iterator at a time!!
    /// </summary>
    public class CircularBuffer : SoundObj
    {
        private uint _length = 0;
        private uint _pos = 0;
        private ISample[] _data;

        private ISample _peak = null;
        private ISample _mean = null;
        private ISample _meandB = null;
        private ISample _stddev = null;
        private ISample _stddevdB = null;

        /// <summary>
        /// Create a circular buffer with a given size
        /// </summary>
        /// <param name="input">Input stream</param>
        /// <param name="bufsize">Size of the circular buffer</param>
        public CircularBuffer(ISoundObj input, uint bufsize)
        {
            _data = new ISample[bufsize];
            _length = bufsize;
            Input = input;
        }

        /// <summary>
        /// Create a circular buffer with a given size
        /// </summary>
        /// <param name="bufsize">Size of the circular buffer</param>
        public CircularBuffer(uint bufsize)
        {
            _data = new ISample[bufsize];
            _length = bufsize;
        }

        public override IEnumerator<ISample> Samples
        {
            get
            {
                if (_input == null)
                {
                    yield break;
                }

                ushort nc = _input.NumChannels;
                _peak = new Sample(nc);

                foreach (ISample sample in _input)
                {
                    // We'll stash this in our buffer
                    // at position _pos.
                    _data[_pos++] = sample;

                    // Calculate peak
                    for (ushort c = 0; c < nc; c++)
                    {
                        _peak[c] = Math.Max(_peak[c], sample[c]);
                    }

                    // Forget any previous mean & stddev
                    _mean = null;
                    _meandB = null;
                    _stddev = null;
                    _stddevdB = null;

                    // Wrap the buffer
                    if (_pos >= _length)
                    {
                        _pos = 0;
                    }

                    // return it
                    yield return sample;
                }
            }
        }


        /// <summary>
        /// Access the nth sample.
        /// Sample zero is always the "current" sample
        /// (the one most recently returned from the iterator).
        /// Sample[-1] is the previous sample.
        /// Sample [1] is the oldest sample we have in buffer.
        /// </summary>
        /// <param name="arg">n</param>
        /// <returns>sample</returns>
        public ISample this[int arg]
        {
            get
            {
                // uint n = (uint)((arg - _pos) % _length);
                uint n = (uint)((((_pos + arg) % _length) + _length) % _length);
                ISample s = _data[n];
                if (s == null)
                {
                    s = new Sample(_nc);
                }
                return s;
            }
        }

        /// <summary>
        /// Peak value of all samples, ever
        /// </summary>
        /// <returns></returns>
        public ISample Peak()
        {
            return _peak;
        }

        /// <summary>
        /// "Average" - arithmetical mean of all samples in buffer
        /// </summary>
        /// <returns></returns>
        public ISample Mean()
        {
            if (_mean != null)
            {
                return _mean;
            }
            ushort nc = _input.NumChannels;
            _mean = new Sample(nc);
            for (int n = 0; n < _length; n++)
            {
                if (_data[n] != null)
                {
                    for (ushort c = 0; c < nc; c++)
                    {
                        _mean[c] += _data[n][c];
                    }
                }
            }
            for (ushort c = 0; c < nc; c++)
            {
                _mean[c] /= _length;
            }
            return _mean;
        }

        /// <summary>
        /// "Average" - arithmetical mean of dB value of all samples in buffer
        /// </summary>
        /// <returns></returns>
        public ISample MeanDb()
        {
            if (_meandB != null)
            {
                return _meandB;
            }
            ushort nc = _input.NumChannels;
            _meandB = new Sample(nc);
            for (int n = 0; n < _length; n++)
            {
                if (_data[n] != null)
                {
                    for (ushort c = 0; c < nc; c++)
                    {
                        _meandB[c] += MathUtil.dB(_data[n][c]);
                    }
                }
            }
            for (ushort c = 0; c < nc; c++)
            {
                _meandB[c] /= _length;
            }
            return _meandB;
        }

        /// <summary>
        /// Standard deviation of all samples in buffer
        /// </summary>
        /// <returns></returns>
        public ISample StdDev()
        {
            if (_stddev != null)
            {
                return _stddev;
            }
            ISample mean = Mean();
            ushort nc = _input.NumChannels;
            _stddev = new Sample(nc);
            for (int n = 0; n < _length; n++)
            {
                if (_data[n] != null)
                {
                    for (ushort c = 0; c < nc; c++)
                    {
                        double dev = (_data[n][c] - mean[c]);
                        _stddev[c] += (dev * dev);
                    }
                }
            }
            for (ushort c = 0; c < nc; c++)
            {
                _stddev[c] = Math.Sqrt(_stddev[c]/_length);
            }
            return _stddev;
        }

        /// <summary>
        /// Standard deviation in dB of all samples in buffer
        /// </summary>
        /// <returns></returns>
        public ISample StdDevDb()
        {
            if (_stddevdB != null)
            {
                return _stddevdB;
            }
            ISample meandb = MeanDb();
            ushort nc = _input.NumChannels;
            _stddevdB = new Sample(nc);
            for (int n = 0; n < _length; n++)
            {
                if (_data[n] != null)
                {
                    for (ushort c = 0; c < nc; c++)
                    {
                        double devdb = (MathUtil.dB(_data[n][c]) - meandb[c]);
                        _stddevdB[c] += (devdb * devdb);
                    }
                }
            }
            for (ushort c = 0; c < nc; c++)
            {
                _stddevdB[c] = Math.Sqrt(_stddevdB[c] / _length);
            }
            return _stddevdB;
        }
    }


    // A "reader" buffer whose input is a sample array
    public class BufferReader : SoundObj
    {
        private Sample[] _data;

        public BufferReader()
        {
            //_data = null;
        }
        public BufferReader(Sample[] data)
        {
            _data = data;
        }

        /// <summary>
        /// Get an iterator for samples
        /// </summary>
        public override IEnumerator<ISample> Samples
        {
            get
            {
                if (_data == null)
                {
                    yield break;
                }
                foreach (ISample sample in _data)
                {
                    yield return sample;
                }
            }
        }
    }

    // A "reader" buffer whose input is an array of complex arrays (one per input channel)

    public enum ComplexBufferFlags
    {
        Both = 0,
        RealOnly = 1,
        ImaginaryOnly = 2, 
        Magnitude = 3,
        Phase = 4
    }

    public class ComplexBufferReader : SoundObj
    {
        private Complex[][] _data;
        private int _start;
        private int _end;
        private ComplexBufferFlags _flags;

        public ComplexBufferReader()
        {
            // _data = null;
        }

        public ComplexBufferReader(Complex[] data, int start, int end)
        {
            _data = new Complex[1][];
            _data[0] = data;
            _start = Math.Max(start, 0);
            _end = Math.Min(end, data.Length);
            _flags = ComplexBufferFlags.RealOnly;
            NumChannels = 1;
        }
        public ComplexBufferReader(Complex[] data, int start, int end, ComplexBufferFlags flags)
        {
            _data = new Complex[1][];
            _data[0] = data;
            _start = Math.Max(start, 0);
            _end = Math.Min(end, data.Length);
            _flags = flags;
            NumChannels = 1;
        }

        public ComplexBufferReader(Complex[][] data, ushort nChannels, int start, int end, ComplexBufferFlags flags)
        {
            _data = data;
            _start = Math.Max(start, 0);
            _end = Math.Min(end, data[0].Length);
            _flags = flags;
            NumChannels = (ushort)((_flags==ComplexBufferFlags.Both) ? (2 * nChannels) : nChannels);
        }

        public override int Iterations
        {
            get
            {
                return _end - _start;
            }
        }
        /// <summary>
        /// Get an iterator for samples
        /// </summary>
        public override IEnumerator<ISample> Samples
        {
            get
            {
                if (_data == null)
                {
                    yield break;
                }
                for(int j=_start; j<_end; j++)
                {
                    yield return _next(j);
                }
            }
        }

        internal ISample _next(int j)
        {
            ISample s = _nc == 2 ? new Sample2() : new Sample(_nc) as ISample;
            int cc = 0;
            for (int c = 0; c < _nc; c += ((_flags == ComplexBufferFlags.Both) ? 2 : 1))
            {
                if (_flags == ComplexBufferFlags.Both)
                {
                    s[c] = _data[cc][j].Re;
                    s[c + 1] = _data[cc][j].Im;
                }
                else if (_flags == ComplexBufferFlags.RealOnly)
                {
                    s[c] = _data[cc][j].Re;
                }
                else if (_flags == ComplexBufferFlags.ImaginaryOnly)
                {
                    s[c] = _data[cc][j].Im;
                }
                else if (_flags == ComplexBufferFlags.Magnitude)
                {
                    s[c] = _data[cc][j].Magnitude;
                }
                else if (_flags == ComplexBufferFlags.Phase)
                {
                    s[c] = _data[cc][j].Phase;
                }
                cc++;
            }
            return s;
        }
    }
}
