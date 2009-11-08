using System;
using System.Collections.Generic;
using System.Text;
using System.IO;
using System.Threading;

// Copyright (c) 2006, 2007 by Hugh Pyle, inguzaudio.com

using DSPUtil;

namespace DSPUtil
{
    // Convolver for Samples

    // [ServiceContract]
    public interface IConvolver : ISoundObj
    {
        // [OperationContract]
        ISoundObj impulse { set; }

        // [OperationContract]
        // double gain { get; }

        // [OperationContract]
        int partitions { get; set; }

        // [OperationContract]
        bool deconvolve { get; set; }

        string PersistPath { set; }
        string PersistTail { get; set; }
        bool IsPersistTail { get; }
    }

    // Fast (frequency-domain) convolver
    // [Serializable]
    public class FastConvolver : SoundObj, IConvolver
    {
        private Complex[][] _NormalImpulseFFT;
        private Complex[][][] _PartitionedImpulseFFT;

        protected bool _running = false;

        // The impulse
        protected ISoundObj _impulse;
        protected int _impulseLength = DSPUtil.BUFSIZE;
        protected int _impulseLengthOrig = DSPUtil.BUFSIZE;
        protected bool _impulseFFTReady = false;

        #region Constructors
        /// <summary>
        /// Constructor
        /// </summary>
        public FastConvolver()
        {
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="impulseFilter">The impulse</param>
        public FastConvolver(ISoundObj impulseFilter)
        {
            impulse = impulseFilter;
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="input">The input</param>
        /// <param name="impulseFilter">The impulse</param>
        public FastConvolver(ISoundObj input, ISoundObj impulseFilter)
        {
            Input = input;
            impulse = impulseFilter;
        }
        #endregion


        #region IConvolver implementation

        // De-convolution?  Default is normal convolution.
        protected bool _deconvolve;
        public bool deconvolve
        {
            get
            {
                return _deconvolve;
            }
            set
            {
                if (_running) { throw new Exception("Cannot change mode while running."); }
                _deconvolve = value;
            }
        }

        // Persist the convolution tail, and use it to feed the next convolution? Default no (null).
        // Set to a "name of this convolver" (e.g. the squeezebox MAC address it is destined for).
        protected string _persistTail;
        protected string _persistFile;
        protected string _persistPath = Path.GetTempPath();
        public string PersistPath
        {
            set
            {
                _persistPath = value;
            }
        }
        public string PersistTail
        {
            get { return _persistTail; }
            set
            {
                _persistTail = value;
                if (!String.IsNullOrEmpty(_persistTail))
                {
                    string filename = this.GetType() + "." + _persistTail;
                    foreach (char c in System.IO.Path.GetInvalidFileNameChars())
                    {
                        filename = filename.Replace(c, '_');
                    }
                    _persistFile = Path.Combine(_persistPath, filename + ".tail");
                }
//              Trace.WriteLine("Tail {0}", _persistTail);
            }
        }
        public bool IsPersistTail
        {
            get
            {
                // Trace.WriteLine("IsPersistTail? {0}, {1}", _persistTail, (_persistTail != null));
                return (!String.IsNullOrEmpty(_persistTail));
            }
        }

        // Number of partitions -- zero for regular unpartitioned convolution
        protected int _partitions;
        public int partitions
        {
            get
            {
                return _partitions;
            }
            set
            {
                if (_running) { throw new Exception("Cannot change mode while running."); }
                _partitions = value;
                _impulseFFTReady = false;
            }
        }

        #endregion

        #region SoundObj overrides

        public override int Iterations
        {
            get
            {
                if (IsPersistTail)
                {
                    // The tail will be wrapped around.
                    // Expected stream length is same as the input stream.
                    return _input.Iterations;
                }
                return ((_impulse == null) ? 0 : _impulseLengthOrig) + _input.Iterations;
            }
        }

        public override ISoundObj Input
        {
            get
            {
                return base.Input;
            }
            set
            {
                if (value != null && _impulse != null)
                {
                    // If the impulse is one channel, that's always OK; apply the same impulse to each input channel.
                    // If the impulse has the same number of channels as the input, OK: apply each channel of impulse to the corresponding channel of input.
                    // But otherwise, channels mismatch is not allowed (e.g. 2-channel impulse & 1-channel or 4-channel input)
                    if ((_impulse.NumChannels > 1) && (_impulse.NumChannels != value.NumChannels))
                    {
                        throw new ArgumentException(String.Format("Impulse of {0} channels can't be used with an input of {1} channels.", _impulse.NumChannels, value.NumChannels));
                    }
                }
                base.Input = value;
            }
        }
        #endregion

        /// <summary>
        /// Note: setting the impulse reads the entire impulse stream immediately; you can close it after calling.
        /// Note: the convolver creates a buffer of the impulse so you don't need to.
        /// </summary>
        public ISoundObj impulse
        {
            get
            {
                return _impulse;
            }
            set
            {
//                if (_impulse == null && value!=null && _running)
//                {
//                    // Cannot set the impulse while running if we were initialized with no impulse
//                    throw new ArgumentException("Cannot change impulse from null when running.");
//                }
                SoundBuffer imbuff = new SoundBuffer(value);
                int iml = imbuff.ReadAll();
                int npt = MathUtil.NextPowerOfTwo(iml);

                if (value != null)
                {
                    if (_impulse != null && npt != _impulseLength)
                    {
                        throw new ArgumentException("Cannot change impulse size when running.");
                    }
                }
                _impulse = imbuff;
                if (value == null)
                {
                    // Still run, just pass the input thru without any processing.
                    _enabled = false;
                    return;
                }
                _impulseLengthOrig = iml;
                _impulseLength = npt;
                if (_input!=null)
                {
                    // If the impulse is one channel, that's always OK; apply the same impulse to each input channel.
                    // If the impulse has the same number of channels as the input, OK: apply each channel of impulse to the corresponding channel of input.
                    // But otherwise, channels mismatch is not allowed (e.g. 2-channel impulse & 1-channel or 4-channel input)
                    if ((_impulse.NumChannels > 1) && (_impulse.NumChannels != _input.NumChannels))
                    {
                        throw new ArgumentException(String.Format("Impulse of {0} channels can't be used with an input of {1} channels.", _impulse.NumChannels, _input.NumChannels));
                    }
                }
                // Compute FFT of the new impulse
                _enabled = true;
                _impulseFFTReady = false;
                ComputeImpulseFFT();
            }
        }


        /// <summary>
        /// Get an iterator for convolved samples
        /// </summary>
        public override IEnumerator<ISample> Samples
        {
            get
            {
                // Unpartitioned overlap-and-save algorithm:
                // - Impulse is length N
                //   - Null-padded to length Nh = N*2
                //
                // OR:
                //
                // Uniformly partitioned overlap-and-save, from Stockham1966, Kulp1998, and
                // Armelloni/Giottoli/Farina2003 "implementation of real-time partitioned convolution on a DSP board"
                // - Impulse is broken into P blocks length K = N/_partitions
                //   - each block is null-padded to length L = nextpowerof2(K*2)
                //   - each block FFTd, so we have a collection of impulse FFTs S
                // - Input is read in overlapping blocks length L (each block begins L-K points after the previous)
                //   - FFT each, multiply with the S
                //   - Sum in a collection of freq-domain accumulators
                // - At the end, IFFT the first accumulator, and yield the latest L-K points from output data.

                if (_input == null)
                {
                    Trace.WriteLine("Null input");
                    yield break;
                }

                if (_impulse.NumChannels == 0 || !_enabled)
                {
                    // Pass through; just return the samples one at a time
                    foreach (ISample sample in _input)
                    {
                        yield return sample;
                    }
                    yield break;
                }

                init();
                p = 0;
                inputSamples = _input.GetBufferedEnumerator() as ISampleBuffer;
                bool moreSamples = true;
                bool tail = false;
                int tailLen = 0;

                // Storage for the tail of this convolution
                List<ISample> thisTailSamples = new List<ISample>(_impulseLengthOrig>0 ? _impulseLengthOrig : DSPUtil.BUFSIZE);

                while (moreSamples || tail)
                {
                    int nn;
                    if (_partitions <= 0)
                    {
                        compute0(out nn, out moreSamples);
                        // Yield the first half of the output
                        for (int j = 0; j < nn; j++)
                        {
                            yield return _next(j);
                        }
                    }
                    else
                    {
                        computeP(out nn, out moreSamples);
                        // Yield the first half of the output
                        for (int j = 0; j < nn; j++)
                        {
                            yield return _next(j);
                        }
                        // Move to the next partition
                        p = (p + 1) % P;
                    }


                    if (tail)
                    {
                        // There's more tail
                        for (int j = 0; j < Math.Min(K,tailLen); j++)
                        {
                            thisTailSamples.Add(_next(j));
                        }
                        if (K > tailLen)
                        {
                            // Done collecting the tail
                            // Trace.WriteLine("Write tail? ");
                            if (IsPersistTail)
                            {
                                // Write the tail to disk.
                                // Tail becomes input to the next convolution operation.
                                Trace.WriteLine("Write {0} ({1})", _persistFile.Replace(_persistPath, ""), thisTailSamples.Count);
                                try
                                {
                                    WaveWriter tailWriter = new WaveWriter(_persistFile);
                                    //tailWriter.Input = new ComplexBufferReader(output, nChannels, nn, nn + K, ComplexBufferFlags.Both);
                                    tailWriter.Input = new CallbackSource(nChannels, SampleRate, delegate(long n)
                                    {
                                        if (n >= thisTailSamples.Count)
                                        {
                                            return null;
                                        }
                                        return thisTailSamples[(int)n];
                                    });
                                    tailWriter.SampleRate = SampleRate;
                                    tailWriter.Format = WaveFormat.IEEE_FLOAT;
                                    tailWriter.BitsPerSample = 64;
                                    tailWriter.Run();
                                    tailWriter.Close();
                                }
                                catch (Exception e)
                                {
                                    Trace.WriteLine("Failed, {0}", e.Message);
                                    // Carry on regardless
                                }
                            }
                            else
                            {
                                // yield it
                                foreach (ISample s in thisTailSamples)
                                {
                                    yield return s;
                                }
                            }
                            yield break;
                        }
                        else
                        {
                            // There's yet more tail...
                            tailLen -= K;
                        }
                    }

                    if (!moreSamples && !tail)
                    {
                        int nK = Math.Min(_impulseLengthOrig, K - nn);
                        for (int j = 0; j < nK; j++)
                        {
                            thisTailSamples.Add(_next(j + nn));
                        }
                        tail = true;
                        tailLen = _impulseLengthOrig - nK;
                    }
                }

                done();
            }
        }


        ushort nChannels;
        int N;
        int Nh;
        int K, L, P;
        int p;
        Complex[][] src;
        Complex[][] output;
        Complex[][][] accum;
        ISampleBuffer inputSamples;

        List<ISample> prevTailSamples;
        IEnumerator<ISample> prevTailEnum = null;

        ushort nImpulseChannels;

        private ISample _next(int j)
        {
            // Yield the output sample (real-only)
            if (prevTailEnum != null)
            {
                if (prevTailEnum.MoveNext())
                {
                    for (ushort c = 0; c < nChannels; c++)
                    {
                        output[c][j].Re += prevTailEnum.Current[c];
                    }
                }
                else
                {
                    prevTailEnum = null;
                }
            }
            if (nChannels == 2)
            {
                return new Sample2(output[0][j].Re, output[1][j].Re);
            }
            ISample ret = new Sample(nChannels);
            for (ushort c = 0; c < nChannels; c++)
            {
                ret[c] = output[c][j].Re;
            }
            return ret;
        }

        private void init()
        {
            _running = true;

            ComputeImpulseFFT();

            nChannels = NumChannels;
            // Size for work area
            N = _impulseLength;
            nImpulseChannels = _impulse.NumChannels;
            Nh = N << 1;

            K = (_partitions <= 0) ? N : (N / _partitions);
            L = MathUtil.NextPowerOfTwo(K << 1);
            P = (int)Math.Ceiling((double)N / (double)K);

            //              Trace.WriteLine("{0} partitions of size {1}, chunk size {2}", P, L, K);

            // Allocate buffers for the source and output
            src = new Complex[nChannels][];
            output = new Complex[nChannels][];
            for (ushort c = 0; c < nChannels; c++)
            {
                src[c] = new Complex[Nh];
                output[c] = new Complex[Nh];
            }

            // For partitioned convolution, allocate arrays for the fft accumulators
            if (_partitions > 0)
            {
                accum = new Complex[nChannels][][];
                for (ushort c = 0; c < nChannels; c++)
                {
                    accum[c] = new Complex[P][];
                    for (p = 0; p < P; p++)
                    {
                        accum[c][p] = new Complex[L];
                    }
                }
            }

            if (IsPersistTail)
            {
                // If there's any *unprocessed* leftover from a previous convolution,
                // load it into prevTailSamples before we begin
                if (System.IO.File.Exists(_persistFile))
                {
                    WaveReader tailReader = null;
                    try
                    {
                        tailReader = new WaveReader(_persistFile);
                        prevTailSamples = new List<ISample>(tailReader.Iterations);
                        //                            Trace.WriteLine("Read {0} tail samples", tailReader.Iterations);
                        if (tailReader.NumChannels == NumChannels)
                        {
                            foreach (ISample s in tailReader)
                            {
                                prevTailSamples.Add(s);
                            }
                        }
                        prevTailEnum = prevTailSamples.GetEnumerator();
                    }
                    catch (Exception e)
                    {
                        Trace.WriteLine("Could not read tail {0}: {1}", _persistFile.Replace(_persistPath, ""), e.Message);
                    }
                    // finally...
//                    ThreadPool.QueueUserWorkItem(delegate(object o)
//                    {
                        try
                        {
                            if (tailReader != null)
                            {
                                tailReader.Close();
                                tailReader = null;
                            }
                            System.IO.File.Delete(_persistFile);
                        }
                        catch (Exception)
                        {
                            // ignore, we're done
                        }
//                    });
                }
            }
//              Trace.WriteLine("{0} partitions of size {1}, chunk size {2}", P, L, K);
        }

        struct mre
        {
            public ManualResetEvent e;
            public int p;
            public ushort c;
        }

        // Unpartitioned compute
        unsafe int compute0(out int nn, out bool moreSamples)
        {
            ushort nc = nChannels;

            // Input is always buffered
            // Read the next Nh/2 data points into the second half of the src buffer
            // NB ReadComplex guarantees the buffer to be padded with zeros past nn
            Complex[][] buf = inputSamples.ReadComplex(K, out nn, out moreSamples);

            mre[] ev = new mre[nc];
            ManualResetEvent[] ee = new ManualResetEvent[nc];
            for (ushort c = 0; c < nc; c++)
            {
                ManualResetEvent mr = new ManualResetEvent(false);
                ev[c].e = mr;
                ev[c].c = c;
                ee[c] = mr;
                ThreadPool.QueueUserWorkItem(delegate(object e)
                {
                    mre m = (mre)e;
                    ushort _c = m.c;

                    Complex[] impulseFFT = _NormalImpulseFFT[(nImpulseChannels == 1) ? (ushort)0 : _c];

                    Array.Copy(buf[_c], 0, src[_c], K, K);
                    // if (kj != 0) Array.Clear(src[_c], jj + K, K - jj);

                    Fourier.Convolve(Nh, impulseFFT, src[_c], _deconvolve);

                    Array.Copy(src[_c], output[_c], K);
                    Array.Copy(buf[_c], src[_c], K);
                    // if (kj != 0) Array.Clear(src[_c], jj, K - jj);

                    m.e.Set();
                }, ev[c]);
            }
            WaitHandle.WaitAll(ee);
            //Thread.Sleep(1); // (1);
            return nn;
        }

        // Partitioned compute
        unsafe int computeP(out int nn, out bool moreSamples)
        {
            ushort nc = nChannels;

            // Input is always buffered
            // Read the next Nh/2 data points into the second half of the src buffer
            // NB ReadComplex guarantees the buffer to be padded with zeros past nn
            Complex[][] buf = inputSamples.ReadComplex(K, out nn, out moreSamples);

            mre[] ev = new mre[nc];
            ManualResetEvent[] ee = new ManualResetEvent[nc];
            for (ushort c = 0; c < nc; c++)
            {
                ManualResetEvent mr = new ManualResetEvent(false);
                ev[c].e = mr;
                ev[c].c = c;
                ev[c].p = p;
                ee[c] = mr;
                ThreadPool.QueueUserWorkItem(delegate(object e)
                {
                    mre m = (mre)e;
                    ushort _c = m.c;

                    Complex[][] impulseFFT = _PartitionedImpulseFFT[(nImpulseChannels == 1) ? (ushort)0 : _c];

                    Array.Copy(buf[_c], 0, src[_c], K, K);

                    Fourier.ConvolvePart(L, m.p, P, impulseFFT, src[_c], accum[_c], output[_c], _deconvolve);

                    Array.Copy(buf[_c], src[_c], K);

                    m.e.Set();
                }, ev[c]);
            }
            WaitHandle.WaitAll(ee);
            //Thread.Sleep(1); // (1);
            return nn;
        }

        void done()
        {
            // Done
            _running = false;
            src = null;
            output = null;
            accum = null;
            inputSamples = null;
        }

        #region Private methods

        private void ComputeImpulseFFT()
        {
            if (!_impulseFFTReady)
            {
                // Compute FFT of the new impulse
                if (_partitions == 0)
                {
                    _impulseFFTReady = ComputeNormalImpulseFFT();
                }
                else
                {
                    _impulseFFTReady = ComputePartitionedImpulseFFT();
                }
            }
        }

        private bool ComputeNormalImpulseFFT()
        {
            if (_impulse == null)
            {
                return false;
            }
            ushort nChannels;
            if (_input == null)
            {
                nChannels = _impulse.NumChannels;
            }
            else
            {
                nChannels = _input.NumChannels;
            }
            int N = _impulseLength;
            int Nh = N << 1;

            // Initialize the array of impulse-FFTs
            _NormalImpulseFFT = new Complex[nChannels][];
            for (ushort c = 0; c < nChannels; c++)
            {
                _NormalImpulseFFT[c] = new Complex[Nh];
            }

            // Read all samples from the impulse
            // into the second half of the src buffer (the first half is all zeros)
            int i = 0;
            foreach (ISample sample in _impulse)
            {
                for (ushort c = 0; c < _impulse.NumChannels; c++)
                {
                    _NormalImpulseFFT[c][i + N] = new Complex(sample[c],0);
                }
                i++;
                if (i >= _impulseLength) { break; } // source gave us more data than we'd expected
            }
            _impulseLengthOrig = i - 1;

            // Do the fft of impulse, in-place
            for (ushort c = 0; c < nChannels; c++)
            {
                if (c >= _impulse.NumChannels)
                {
                    // Handle the case where impulse has only one channel, but input (hence nChannels) is wider;
                    // we already computed the fft in _NormalImpulseFFT[0], so just duplicate it.
                    // In fact we can just ref the whole buffer -- no need even to make a deep copy
                    _NormalImpulseFFT[c] = _NormalImpulseFFT[0];
                }
                else
                {
                    Fourier.FFT(Nh, _NormalImpulseFFT[c]);
                }
            }
            return true;
        }


        private bool ComputePartitionedImpulseFFT()
        {
            if (_impulse == null)
            {
                return false;
            }
            ushort nChannels;
            if (_input == null)
            {
                nChannels = _impulse.NumChannels;
            }
            else
            {
                nChannels = _input.NumChannels;
            }
            ushort p;
            int N = _impulseLength;
            int Nh = N << 1;
            int K = (int)(N / _partitions);
            int L = MathUtil.NextPowerOfTwo(K << 1);
            int P = (int)Math.Ceiling((double)N / (double)K);

            // Initialize the arrays of impulse-FFTs
            _PartitionedImpulseFFT = new Complex[nChannels][][];
            Complex[][] tmp = new Complex[nChannels][];
            for (ushort c = 0; c < nChannels; c++)
            {
                tmp[c] = new Complex[Nh];
                _PartitionedImpulseFFT[c] = new Complex[P][];
                for (p = 0; p < P; p++)
                {
                    _PartitionedImpulseFFT[c][p] = new Complex[L];
                }
            }

            // Read all samples from the impulse, & fft in blocks of size K (padded to L, L~=2K)
            int i = 0;
            int n = 0;
            p = 0;
            foreach (ISample sample in _impulse)
            {
                // Reading a segment of the impulse (into src[])
                n++;
                if (n > _impulseLength) { break; } // source gave us more data than we'd expected
                for (ushort c = 0; c < _impulse.NumChannels; c++)
                {
                    tmp[c][i + K] = new Complex(sample[c],0); // TBD
                }
                i++;
                if (i >= K)
                {
                    // Do the fft of this segment of the impulse, into fftImpulse[]
                    for (ushort c = 0; c < nChannels; c++)
                    {
                        if (c >= _impulse.NumChannels)
                        {
                            // Handle the case where impulse has only one channel, but input (hence nChannels) is wider;
                            // we already computed the fft in _PartitionedImpulseFFT[0][p], so just duplicate it.
                            // In fact we can just ref the whole buffer -- no need even to make a deep copy
                            _PartitionedImpulseFFT[c][p] = _PartitionedImpulseFFT[0][p];
                        }
                        else
                        {
                            Array.Copy(tmp[c], _PartitionedImpulseFFT[c][p], L);
                            Fourier.FFT(L, _PartitionedImpulseFFT[c][p]);
                        }
                    }
                    i = 0; p++;
                }
            }
            _impulseLengthOrig = n - 1;
            if (p < P)
            {
                // Fill any extra space with nulls
                for (int ii = i; ii < K; ii++)
                {
                    for (ushort c = 0; c < _impulse.NumChannels; c++)
                    {
                        tmp[c][ii + K] = new Complex(0, 0); // TBD
                    }
                }
                // FFT the last segment
                for (ushort c = 0; c < nChannels; c++)
                {
                    if (c >= _impulse.NumChannels)
                    {
                        // Handle the case where impulse has only one channel, but input (hence nChannels) is wider;
                        // we already computed the fft in _PartitionedImpulseFFT[0][p], so just duplicate it.
                        // In fact we can just ref the whole buffer -- no need even to make a deep copy
                        _PartitionedImpulseFFT[c][p] = _PartitionedImpulseFFT[0][p];
                    }
                    else
                    {
                        Array.Copy(tmp[c], _PartitionedImpulseFFT[c][p], L);
                        Fourier.FFT(L, _PartitionedImpulseFFT[c][p]);
                    }
                }
            }
            return true;
        }
        #endregion
    }


}
