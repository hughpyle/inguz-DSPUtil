using System;
using System.Collections.Generic;
using System.Text;

// Copyright (c) 2006, 2007 by Hugh Pyle, inguzaudio.com

namespace DSPUtil
{
    public enum DitherType
    {
        NONE = 0,
        TRIANGULAR = 1,
        SHAPED = 2,
        SBM = 3
    }

    public class Dither
    {
        // F-weighted

        // These from "Psychoacoustically Optimal Noise Shaping," by Robert A. Wannamaker
        // (From AES Journal Volume 40, No. 7/8, 1992 July/August)
        // !!!NOTE: This dither FIR is only good for 44.1kHz sampling rate!!!
        private static double[] _Wannamaker = {
                                      2.412f,
                                     -3.370f,
                                      3.937f,
                                     -4.174f,
                                      3.353f,
                                     -2.205f,
                                      1.281f,
                                     -0.569f,
                                      0.0847f };

        // These from the Sony Super Bit Mapping (SBM)
        private static double[] _SBM = {
                                      1.47933,
                                     -1.59032,
                                      1.64436,
                                     -1.36613,
                                      9.26704e-1,
                                     -5.57931e-1,
                                      2.67859e-1,
                                     -1.06726e-1,
                                      2.85161e-2,
                                      1.23066e-3,
                                     -6.16555e-3,
                                      3.06700e-3};

        private const double _scale8 = 128f;
        private const double _scale16 = 32768f;
        private const double _scale24 = 8388608f;
        private const double _scale32 = 2147483648f;

        private DitherType _type = DitherType.NONE;
        private uint _sampleRate;
        private ushort _bitsPerSample;
        private double[] _filter;

        private double _minv;
        private double _maxv;
        private double _peak;

        // Error history array for noise shaping
        private double[] _EH;
        private int _HistPos;

        // Cached random-number generator
        private Random _random;

        // Previously tried to cache a block of random numbers.  But fail:
        // with any reasonably small cache the LF behavior is audible pumping
        // at (samplerate/cachesize).  So instead just use the proper rand gen
        // even though it's a bit slower.
        /*
        private const int _cachelen = 4096;
        private static double[] _randArr;
        private static double[] _tpdfArr;
        private static int _randPos;
        private static int _tpdfPos;
        */

        private bool _clipping = false;

        public Dither(DitherType type, uint sampleRate, ushort bitsPerSample)
        {
            _type = type;
            _sampleRate = sampleRate;
            _bitsPerSample = bitsPerSample;

            _minv = MinValue(_bitsPerSample);
            _maxv = MaxValue(_bitsPerSample);

            _filter = FilterArray(type);
            _EH = new double[2 * _filter.Length];
            _HistPos = _filter.Length - 1;

            // Tweak the random number source
            _random = new Random();
        }

        static private double[] FilterArray(DitherType type)
        {
            switch (type)
            {
                case DitherType.SHAPED:
                    return _Wannamaker;
                case DitherType.SBM:
                    return _SBM;
                case DitherType.NONE:
                case DitherType.TRIANGULAR:
                default:
                    return new double[0];
            }
        }

        private static double MaxValue(ushort bitsPerSample)
        {
            switch (bitsPerSample)
            {
                case 8:
                    return _scale8; // -1;
                case 16:
                    return _scale16; // -1;
                case 24:
                    return _scale24; // -1;
                case 32:
                    return _scale32; // -1;
                default:
                    return Math.Pow(2, bitsPerSample - 1) - 1;
            }
        }

        private static double MinValue(ushort bitsPerSample)
        {
            switch (bitsPerSample)
            {
                case 8:
                    return -_scale8;
                case 16:
                    return -_scale16;
                case 24:
                    return -_scale24;
                case 32:
                    return -_scale32;
                default:
                    return -(Math.Pow(2, bitsPerSample - 1));
            }
        }

        // Random value from -1 to +1 with rectangular PDF
        private double NextRandom()
        {
            return (_random.NextDouble() - 0.5);
            /*
            lock (_random)
            {
                if (_randArr == null)
                {
                    _randArr = new double[_cachelen];
                    for (int j = 0; j < _cachelen; j++)
                    {
                        _randArr[j] = (_random.NextDouble() - 0.5) * 2;
                    }
                    _randPos = 0;
                }
                _randPos++;
                if (_randPos >= _cachelen)
                {
                    _randPos = 0;
                }
                return _randArr[_randPos];
            }
            */
        }
        // Random value from -1 to +1 with triangular PDF
        private double NextRandom2()
        {
            return (_random.NextDouble() + _random.NextDouble() - 1);
            /*
            lock (_random)
            {
                if (_tpdfArr == null)
                {
                    _tpdfArr = new double[_cachelen];
                    for (int j = 0; j < _cachelen; j++)
                    {
                        _tpdfArr[j] = (_random.NextDouble() + _random.NextDouble() - 1) * 2;
                    }
                    _tpdfPos = 0;
                }
                _tpdfPos++;
                if (_tpdfPos >= _cachelen)
                {
                    _tpdfPos = 0;
                }
                return _tpdfArr[_tpdfPos];
            }
            */
        }

        public bool clipping
        {
            get { return _clipping; }
            set { _clipping = value; }
        }

        public double dbfsPeak
        {
            get
            {
                // _peak is from 0 to 1 (unless we clipped, in which case it might be more than 1)
                return MathUtil.dB(_peak);
            }
        }

        public double processDouble(double samp)
        {
            int quantized = process(samp);
            return (double)quantized / _maxv;
        }

        /// <summary>
        /// Process a ISample (in-place; the original object is returned, with modified data)
        /// </summary>
        /// <param name="samp"></param>
        public ISample process(ISample samp)
        {
            for (int c = 0; c < samp.NumChannels; c++)
            {
                samp[c] = processDouble(samp[c]);
            }
            return samp;
        }

        // Dither double-precision float to arbitrary bits-per-sample (up to 32)
        public int process(double samp)
        {
            // peak is absolute, normal range should be 0 thru 1
            if (double.IsNaN(samp) || double.IsInfinity(samp))
            {
                samp = 0;
            }
            _peak = Math.Max(_peak, Math.Abs(samp));

            int output = 0;
            double noise;
            double error;
            try
            {
                // Scale the output first,
                // - we're still working with floats
                // - max scale is dependent on the target bit-depth (eg. 32768 for 16-bit)
                // - min scale is +/- 1, any fractional component will be rounded
                samp *= _maxv;
                switch (_type)
                {
                    case DitherType.NONE:
                        output = (int)Math.Round(samp); //, MidpointRounding.AwayFromZero);
                        break;

                    case DitherType.TRIANGULAR:
                        // Triangular dither
                        // Make +/- 1LSB white randomness
                        noise = NextRandom2();
                        // Add this noise
                        samp += noise;
                        // Round the sample
                        output = (int)Math.Round(samp); //, MidpointRounding.AwayFromZero);
                        break;

                    case DitherType.SHAPED:
                    case DitherType.SBM:
                        // Make +/- 1LSB white randomness
                        noise = NextRandom();
                        // Add this noise and subtract the previous noise (filtered)
                        samp += noise;
                        for (int x = 0; x < _filter.Length; x++)
                        {
                            samp -= _filter[x] * _EH[_HistPos + x];
                        }
                        // Round the sample
                        output = (int)Math.Round(samp); //, MidpointRounding.AwayFromZero);
                        break;
                }

                if (_filter.Length > 0)
                {
                    // Find the error of this output sample (before we clip...)
                    error = (double)output - samp;

                    if (_HistPos < 1) _HistPos += _filter.Length;
                    _HistPos--;

                    // Update buffer (both copies)
                    _EH[_HistPos + _filter.Length] = _EH[_HistPos] = error;
                }
            }
            catch (OverflowException)
            {
                if (samp > (_maxv-1))
                    output = (int)(_maxv-1);
                else
                    output = (int)_minv;
            }
            if (output < _minv)
            {
                if (!_clipping)
                {
                    _clipping = true;
                    Trace.WriteLine("CLIPPING-");
                }
                output = (int)_minv;
            }
            else if (output > (_maxv-1))
            {
                if (!clipping)
                {
                    _clipping = true;
                    Trace.WriteLine("CLIPPING+");
                }
                output = (int)(_maxv-1);
            }

            return output;
        }

        public int process(float samp)
        {
            return process((double)samp);
        }
        
        public void reset()
        {
            for (int x = 0; x < _filter.Length; x++)
            {
                _EH[x] = 0.0;
            }
            _HistPos = _filter.Length - 1;
        }
    }
}