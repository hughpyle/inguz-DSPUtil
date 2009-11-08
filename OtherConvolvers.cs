using System;
using System.Collections.Generic;
using System.Text;
using System.IO;
using System.Threading;

// Copyright (c) 2006 by Hugh Pyle, inguzaudio.com

using DSPUtil;

namespace DSPUtil
{
    // Slower-than-fast convolvers for Samples


    // Slow (time-domain) convolver
    // Overlap-and-add method
    [Serializable]
    public class SlowConvolver : SoundObj, IConvolver
    {
        #region Instance members
        // A convolver for each channel
        TimeDomainConvolver[] _convolver;

        // Impulse length rounded up to power of two
        int _length;

        // Buffers for each channel of the two data inputs
        double[][] _impbuff;
        double[][] _databuff;

        // The convolver's output buffer (one for each channel)
        double[][] _outbuff;

        // A carry-forward buffer
        double[][] _copybuff;
        #endregion


        #region IConvolver impl

        // The impulse
        protected ISoundObj _impulse;
        protected int _impulseLength = DSPUtil.BUFSIZE;
        public virtual ISoundObj impulse
        {
            get
            {
                return _impulse;
            }
            set
            {
                // Note: it's OK to reset the impulse while running
                _impulse = value;
                _impulseLength = (_impulse == null) ? 0 : MathUtil.NextPowerOfTwo(_impulse.Iterations);

                // Read the impulse into a buffer, separating the channels for performance
                // Pad the impulse buffer to 2^n
                SoundBuffer buff = new SoundBuffer(_impulse);
                _impbuff = buff.ToDoubleArray(0, ImpulseLength);

            }
        }

        // De-convolution?  Default is normal convolution.
        protected bool _deconvolve;
        public bool deconvolve
        {
            get
            {
                return false;
            }
            set
            {
                if (value)
                {
                    throw new Exception("Time domain deconvolution not implemented.");
                }
            }
        }

        // Persist the convolution tail, and use it to feed the next convolution? Default no (null).
        // Set to a "name of this convolver" (e.g. the squeezebox MAC address it is destined for).
        // NOT IMPLEMENTED FOR TIME_DOMAIN CONVOLVER
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
                    _persistFile = Path.Combine(Path.GetTempPath(), filename + ".tail");
                }
            }
        }
        /// <summary>
        /// This returns a bogus value. Tail persistence is not implemented in SlowConvolver.
        /// </summary>
        public bool IsPersistTail
        {
            get
            {
                return !String.IsNullOrEmpty(_persistTail);
            }
        }

        // Number of partitions -- zero for regular unpartitioned convolution
        public int partitions
        {
            get
            {
                return 0;
            }
            set
            {
                // ignore
            }
        }


        int ImpulseLength
        {
            get
            {
                // For now only iterate the impulse once
                if (_length == 0)
                {
                    _length = MathUtil.NextPowerOfTwo(_impulse.Iterations);
                }
                return _length;
            }
        }

        #endregion

        // Gain, zero if not specified
        protected double _gain;
        public double gain
        {
            get
            {
                return _gain;
            }
        }

        public override int Iterations
        {
            get
            {
                return ((_impulse == null) ? 0 : _impulseLength) + _input.Iterations;
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
                    yield break;
                }

                // Verify the input etc
                ushort nChannels = _input.NumChannels;

                // If the impulse is one channel, that's always OK; apply the same impulse to each input channel.
                // If the impulse has the same number of channels as the input, OK: apply each channel of impulse to the corresponding channel of input.
                // But otherwise, channels mismatch is not allowed (e.g. 2-channel impulse & 1-channel or 4-channel input)
                if ((_impulse.NumChannels > 1) && (_impulse.NumChannels != _input.NumChannels))
                {
                    throw new ArgumentException(String.Format("Impulse of {0} channels can't be used with an input of {1} channels.", _impulse.NumChannels, _input.NumChannels));
                }


                // Allocate an array to hold the output buffers
                // (the buffers themselves are allocated by the convolver)
                _outbuff = new double[nChannels][];

                // Allocate buffers for data & copyforward
                _databuff = new double[nChannels][];
                _copybuff = new double[nChannels][];
                for (int c = 0; c < nChannels; c++)
                {
                    _databuff[c] = new double[ImpulseLength];
                    _copybuff[c] = new double[ImpulseLength];
                }

                // Read the impulse into a buffer, separating the channels for performance
                // Pad the impulse buffer to 2^n
//                SoundBuffer buff = new SoundBuffer(_impulse);
//                _impbuff = buff.ToDoubleArray(0, ImpulseLength);

                // Create the convolvers, one per channel
                _convolver = new TimeDomainConvolver[nChannels];
                for (int c = 0; c < nChannels; c++)
                {
                    double[] imp = _impbuff[_impulse.NumChannels == 1 ? 0 : c];
                    _convolver[c] = new TimeDomainConvolver(imp, _deconvolve);
                }

                IEnumerator<ISample> inputSamples = _input.Samples;
                bool moreSamples = true;
                bool tail = false;
                while (moreSamples) // || tail)
                {
                    // Read input, up to size of the impulse
                    int x = 0;
                    for (int j = 0; j < ImpulseLength; j++)
                    {
                        ISample s;
                        if (moreSamples)
                        {
                            moreSamples = inputSamples.MoveNext();
                        }
                        if (moreSamples)
                        {
                            s = inputSamples.Current;
                            x = j;
                        }
                        else
                        {
                            s = new Sample(nChannels);  // null data
                        }
                        for (ushort c = 0; c < nChannels; c++)
                        {
                            _databuff[c][j] = s[c];
                        }
                    }
                    tail = (!moreSamples && !tail);

                    // Convolve input with the impulse
                    for (int k = 0; k < nChannels; k++)
                    {
                        if (_outbuff[k] != null)
                        {
                            // Save the latter half of the output buffer for next time overlap-add
                            for (int j = 0; j < ImpulseLength - 1; j++)
                            {
                                _copybuff[k][j] = _outbuff[k][j + ImpulseLength];
                            }
                        }

                        _outbuff[k] = _convolver[k].Convolve(_databuff[k], 1.0f);
                        _gain = _convolver[k].Gain;
                    }

                    if (!tail)
                    {
                        x++;
                    }

                    // Yield the samples from convolver's output, added to any copy-buffer data
                    for (int n = 0; n < x; n++)
                    {
                        ISample ret = nChannels == 2 ? new Sample2() : new Sample(nChannels) as ISample;
                        for (int c = 0; c < nChannels; c++)
                        {
                            ret[c] = _copybuff[c][n] + _outbuff[c][n];
                        }
                        yield return ret;
                    }
                }
            }
        }

    }


    // Direct (and reasonably efficient) time domain convolver, expects equal-size 2^n inputs only,
    // See http://www.musicdsp.org/showone.php?id=65,
    //     http://www.musicdsp.org/showone.php?id=66,
    //     http://mathworld.wolfram.com/KaratsubaMultiplication.html

    [Serializable]
    class TimeDomainConvolver
    {
        double[] _impulse;
        double[] _buffer;
        double[] _output;
        double _gain;
        int _size;
        bool _deconvolve;

        // Util for convolution

        public TimeDomainConvolver(double[] impulse, bool deConvolve)
        {
            _impulse = impulse;
            _deconvolve = deConvolve;
            _size = _impulse.Length;
            if (!MathUtil.IsPowerOfTwo(_size))
            {
                throw new ArgumentException("Convolve: input array size must be power of two");
            }

            _buffer = new double[_size * 2];
            _output = new double[_size * 2];
        }

        public unsafe double[] Convolve(double[] data, double gain)
        {
            if (_size != data.Length)
            {
                throw new ArgumentException("Convolve: input arrays must be the same size");
            }

            // arr_mul_knuth(_output, 0, _impulse, 0, data, 0, _buffer, 0, _size);

            fixed (double* o = _output, a = data, b = _impulse, tmp = _buffer)
            {
                mul_knuth(o, a, b, tmp, (uint)_size);
            }

            _gain = gain;
            if (_gain == 0)
            {
                // Gain not specified, calculate it from peak of in2 & output
                // and reduce further by 3dB
                double peak1 = 0;
                for (int j = 0; j < _impulse.Length; j++)
                {
                    peak1 = Math.Max(peak1, Math.Abs(data[j]));
                }
                if (peak1 != 0)
                {
                    double peak2 = 0;
                    for (int j = 0; j < _output.Length; j++)
                    {
                        peak2 = Math.Max(peak2, Math.Abs(_output[j]));
                    }
                    if (!double.IsNaN(peak2))
                    {
                        _gain = peak1 / peak2;
                    }
                    else
                    {
                        _gain = 1.0;
                    }
                }
            }
            // Apply the specified gain
            //            for (int j = 0; j < _output.Length; j++)
            //            {
            //                _output[j] *= _gain;
            //            }

            return _output;
        }

        public double Gain
        {
            get
            {
                return _gain;
            }
        }

        // Array methods, safe

        private void arr_mul_brute(double[] r, uint ro, double[] a, uint ao, double[] b, uint bo, uint w)
        {
            for (uint i = 0; i < w + w; i++)
                r[ro + i] = 0;
            if (_deconvolve)
            {
                for (uint i = 0; i < w; i++)
                {
                    for (uint j = 0; j < w; j++)
                        r[ro + i + j] += a[ao + i] / b[bo + j];
                }
            }
            else
            {
                for (uint i = 0; i < w; i++)
                {
                    for (uint j = 0; j < w; j++)
                        r[ro + i + j] += a[ao + i] * b[bo + j];
                }
            }
        }

        private void arr_mul_knuth(double[] r, uint ro, double[] a, uint ao, double[] b, uint bo, double[] tmp, uint tmpo, uint w)
        {
            if (w < 25)
            {
                arr_mul_brute(r, ro, a, ao, b, bo, w);
            }
            else
            {
                uint m = w >> 1;

                for (uint i = 0; i < m; i++)
                {
                    r[ro + i] = a[ao + m + i] - a[ao + i];
                    r[ro + i + m] = b[bo + i] - b[bo + m + i];
                }

                arr_mul_knuth(tmp, tmpo, r, ro, r, ro + m, tmp, tmpo + w, m);
                arr_mul_knuth(r, ro, a, ao, b, bo, tmp, tmpo + w, m);
                arr_mul_knuth(r, ro + w, a, ao + m, b, bo + m, tmp, tmpo + w, m);

                for (uint i = 0; i < m; i++)
                {
                    double bla = r[ro + m + i] + r[ro + w + i];
                    r[ro + m + i] = bla + r[ro + i] + tmp[tmpo + i];
                    r[ro + w + i] = bla + r[ro + w + m + i] + tmp[tmpo + m + i];
                }
            }
        }

        // Pointer methods, unsafe, approx twice as fast as the safe methods

        private unsafe void mul_brute(double* r, double* a, double* b, uint w)
        {
            for (uint i = 0; i < w + w; i++)
                r[i] = 0;
            if (_deconvolve)
            {
                for (uint i = 0; i < w; i++)
                {
                    for (uint j = 0; j < w; j++)
                        r[i + j] += a[i] / b[j];
                }
            }
            else
            {
                for (uint i = 0; i < w; i++)
                {
                    for (uint j = 0; j < w; j++)
                        r[i + j] += a[i] * b[j];
                }
            }
        }


        // tmp must be of length 2*w
        private unsafe void mul_knuth(double* r, double* a, double* b, double* tmp, uint w)
        {
            if (w < 25)
            {
                mul_brute(r, a, b, w);
            }
            else
            {
                uint m = w >> 1;

                for (uint i = 0; i < m; i++)
                {
                    r[i] = a[m + i] - a[i];
                    r[i + m] = b[i] - b[m + i];
                }

                mul_knuth(tmp, r, r + m, tmp + w, m);
                mul_knuth(r, a, b, tmp + w, m);
                mul_knuth(r + w, a + m, b + m, tmp + w, m);

                for (uint i = 0; i < m; i++)
                {
                    double bla = r[m + i] + r[w + i];
                    r[m + i] = bla + r[i] + tmp[i];
                    r[w + i] = bla + r[w + m + i] + tmp[m + i];
                }
            }
        }
    }

}
