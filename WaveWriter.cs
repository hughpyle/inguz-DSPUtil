using System;
using System.Collections.Generic;
using System.Text;
using System.IO;

// Copyright (c) 2006 by Hugh Pyle, inguzaudio.com
// WAV file writer based originally on Garbe.Sound

namespace DSPUtil
{
    public enum NormalizationType
    {
        PEAK_DBFS = 1,
        PEAK_GAIN,
        RMS_DBFS,
        RMS_GAIN
    }

    /// <summary> Write the samples processed in a wave file </summary>
    [Serializable]
    public sealed class WaveWriter : SoundObj
    {
        // From WinBase.h
        internal const int FILE_TYPE_DISK = 0x0001;
        internal const int FILE_TYPE_CHAR = 0x0002;
        internal const int FILE_TYPE_PIPE = 0x0003;

        // Note, these are #defines used to extract handles, and are NOT handles.
        internal const int STD_INPUT_HANDLE = -10;
        internal const int STD_OUTPUT_HANDLE = -11;
        internal const int STD_ERROR_HANDLE = -12;

        [System.Runtime.InteropServices.DllImport("Kernel32.dll")]
        internal static extern int GetFileType(IntPtr i_Handle);

        [System.Runtime.InteropServices.DllImport("Kernel32.dll", SetLastError = true)]
        internal static extern IntPtr GetStdHandle(int i_Handle);  // param is NOT a handle, but it returns one!

        // Members

        private FileStream _fs;
        private BinaryWriter _w;
        private BufferedStream _bs;

        private bool _raw = false;
        private WaveFormat _audioFormat;
        private WaveFormatEx _formatEx;
        private SpeakerChannelMask _channelMask;
        private DitherType _dither;
        private ushort _bitsPerSample;
        private int _sampleCount;
        private int _iterations = 0;
        private double _gain = 0;
        private bool _ignoreclipping = false;
        private double _normalization = double.NaN;
        private NormalizationType _normType = NormalizationType.PEAK_DBFS;

        private double[] _gains = null;

        private bool _doneHeader;
        private bool _isConsole;

        //        private double _scale8 = 128f;
        //        private double _scale16 = 32768f;
        //        private double _scale24 = 8388608f;
        //        private double _scale32 = 2147483648f;

        private Dither[] _ditherProcessors;

        #region Constructors

        // //////////////////////////////////////////////////////////////////////////////////////////////////
        // Constructors
        // //////////////////////////////////////////////////////////////////////////////////////////////////

        public WaveWriter(string fileName, ushort numChannels, uint sampleRate, ushort bitsPerSample, DitherType dither, WaveFormat format, bool rewrite)
        {
            Initialize(fileName, numChannels, sampleRate, bitsPerSample, dither, format, rewrite);
        }
        public WaveWriter(string fileName, ushort numChannels, uint sampleRate, ushort bitsPerSample, DitherType dither, WaveFormat format)
        {
            Initialize(fileName, numChannels, sampleRate, bitsPerSample, dither, format, true);
        }
        public WaveWriter(string fileName, ushort numChannels, uint sampleRate, ushort bitsPerSample, DitherType dither)
        {
            Initialize(fileName, numChannels, sampleRate, bitsPerSample, dither, WaveFormat.ANY, true);
        }
        public WaveWriter(string fileName, ushort numChannels, uint sampleRate, ushort bitsPerSample)
        {
            Initialize(fileName, numChannels, sampleRate, bitsPerSample, DitherType.NONE, WaveFormat.ANY, true);
        }
        public WaveWriter(string fileName)
        {
            Initialize(fileName, 0, 0, 0, DitherType.NONE, WaveFormat.ANY, true);
        }

        public WaveWriter(Stream output, ushort numChannels, uint sampleRate, ushort bitsPerSample, DitherType dither, WaveFormat format)
        {
            Initialize(output, numChannels, sampleRate, bitsPerSample, dither, format);
        }
        public WaveWriter(Stream output, ushort numChannels, uint sampleRate, ushort bitsPerSample, DitherType dither)
        {
            Initialize(output, numChannels, sampleRate, bitsPerSample, dither, WaveFormat.ANY);
        }
        public WaveWriter(Stream output, ushort numChannels, uint sampleRate, ushort bitsPerSample)
        {
            Initialize(output, numChannels, sampleRate, bitsPerSample, DitherType.NONE, WaveFormat.ANY);
        }
        public WaveWriter(Stream output)
        {
            Initialize(output, 0, 0, 0, DitherType.NONE, WaveFormat.ANY);
        }

        public WaveWriter()
        {
            Stream output = System.Console.OpenStandardOutput();
            _isConsole = true;
            Initialize(output, 0, 0, 0, DitherType.NONE, WaveFormat.ANY);
        }

        // //////////////////////////////////////////////////////////////////////////////////////////////////
        // Helpers
        // //////////////////////////////////////////////////////////////////////////////////////////////////

        private void Initialize(string fileName, ushort numChannels, uint sampleRate, ushort bitsPerSample, DitherType dither, WaveFormat format, bool rewrite)
        {
            NumChannels = numChannels;
            SampleRate = sampleRate;
            BitsPerSample = bitsPerSample;
            _audioFormat = format;
            _dither = dither;
            _sampleCount = 0;
            _doneHeader = false;

            if (File.Exists(fileName))
            {
                if (rewrite == false)
                    throw (new Exception("File already exists: " + fileName));
            }

            _fs = new FileStream(fileName, FileMode.Create);
            _bs = new BufferedStream(_fs);
            _w = new BinaryWriter(_bs);
        }
        private void Initialize(Stream output, ushort numChannels, uint sampleRate, ushort bitsPerSample, DitherType dither, WaveFormat format)
        {
            NumChannels = numChannels;
            SampleRate = sampleRate;
            BitsPerSample = bitsPerSample;
            _audioFormat = format;
            _dither = dither;
            _sampleCount = 0;
            _doneHeader = false;

            _fs = null;
            _bs = new BufferedStream(output);
            _w = new BinaryWriter(_bs);
        }

        /// <summary>
        /// Gain, units
        /// </summary>
        public double Gain
        {
            get
            {
                return _gain;
            }
            set
            {
                _gain = value;
                // If you set gain, we stop applying normalization
                _normalization = double.NaN;
            }
        }

        /// <summary>
        /// Set gain for individual channels.
        /// This is multiplied by the global Gain if applicable.
        /// </summary>
        /// <param name="channel">Channel number</param>
        /// <param name="gain">Gain (units), or double.NaN to reset</param>
        public void SetChannelGain(ushort channel, double gain)
        {
            if (_gains==null)
            {
                _gains = new double[_nc];
                for (ushort c = 0; c < _nc; c++)
                {
                    _gains[c] = double.NaN;
                }
            }
            _gains[channel] = gain;
            _normalization = double.NaN;
        }

        public bool IgnoreClipping
        {
            get
            {
                return _ignoreclipping;
            }
            set
            {
                _ignoreclipping = value;
            }
        }

        private void WriteWaveHeader()
        {
            if (_raw)
            {
                // Don't write any header for RAW files...
                _doneHeader = true;
                return;
            }
            ushort nChannels = NumChannels;
            if (nChannels == 0)
            {
                throw new NotSupportedException("Number of channels cannot be zero");
            }
            ushort bPerSample = BitsPerSample;
            if (bPerSample == 0)
            {
                throw new NotSupportedException("Bits per sample cannot be zero");
            }
            uint sRate = SampleRate;
            if (sRate == 0)
            {
                throw new NotSupportedException("Sample rate cannot be zero");
            }

            ushort blockSize = (ushort)((nChannels * bPerSample) / 8);
            uint dataSize = (uint)(Iterations * blockSize);
            uint fmtSize = (uint)(_audioFormat == WaveFormat.EXTENSIBLE ? 40 : 16);

            // Write Riff ///////////////////////////////////////////////

            _w.Write('R');
            _w.Write('I');
            _w.Write('F');
            _w.Write('F');

            // RIFF size: size of the rest of the riff chunk following this
            // = size of "WAVE" + size of 'fmt' + size of "DATA" + datasize
            // = size of the entire file (bytes) - 8
            uint riffSize = 4 + (8 + fmtSize) + (8 + dataSize);
            _w.Write((uint)riffSize);

            // Write Wave //////////////////////////////////////////////

            _w.Write('W');
            _w.Write('A');
            _w.Write('V');
            _w.Write('E');

            // Write Format ////////////////////////////////////////////

            _w.Write('f');
            _w.Write('m');
            _w.Write('t');
            _w.Write(' ');

            _w.Write((uint)fmtSize);         // size of the fmt block
            _w.Write((ushort)_audioFormat);  // wave format

            _w.Write((ushort)nChannels);              // Number of channels
            _w.Write((uint)sRate);                 // SampleRate (Hz)
            _w.Write((uint)(blockSize * sRate));    //ByteRate
            _w.Write((ushort)blockSize);                //BlockAlign
            _w.Write((ushort)bPerSample);            //BitsPerSample

            if (_audioFormat == WaveFormat.EXTENSIBLE)
            {
                _w.Write((UInt16)22);                       // size of this block
                _w.Write((UInt16)bPerSample);               // union{} = valid bits per sample, in this case
                _w.Write((UInt32)_channelMask);             // channel mask
                _w.Write(_formatEx.guid.ToByteArray());     // GUID, 16 bytes
            }

            // Write Data ///////////////////////////////////////////////

            _w.Write('d');
            _w.Write('a');
            _w.Write('t');
            _w.Write('a');

            // Write Data Size //////////////////////////////////////////

            _w.Write((uint)dataSize);

            _doneHeader = true;
        }

        #endregion


        
        public override IEnumerator<ISample> Samples
        {
            get
            {
                bool err;
                if (_audioFormat == WaveFormat.ANY)
                {
                    throw new Exception("WaveWriter: format not specified");
                }

                MakeDither();
                WriteWaveHeader();

                foreach (ISample sample in _buff())
                {
                    ISample s = _next(sample, out err);
                    if (err)
                    {
                        yield break;
                    }
                    yield return s;
                }
            }
        }

        private ISoundObj _buff()
        {
            if (double.IsNaN(_normalization))
            {
                return _input;
            }
            // We've been asked to normalize
            SoundBuffer b = new SoundBuffer(_input);
            b.ReadAll();
            _gain = b.Normalize(_normalization, false);
            return b;
        }

        private ISample _next(ISample sample, out bool err)
        {
            try
            {
                //            if (_gain == 0) Gain = 1;
                _sampleCount++;

                if (!_ignoreclipping && (_sampleCount % _sr == 0))
                {
                    // Every second, check for clipping and wind back the gain by 0.5dB if so
                    bool clipped = false;
                    for (int n = 0; n < _nc; n++)
                    {
                        if (_ditherProcessors[n].clipping)
                        {
                            clipped = true;
                        }
                    }
                    if (clipped)
                    {
                        // Reduce gain by 0.5dB to back off from clipping
                        Gain = _gain * 0.94406087628592338036438049660227;
                        //                            Trace.WriteLine("Gain {0} dB", MathUtil.dB(_gain));
                        for (int n = 0; n < _nc; n++)
                        {
                            _ditherProcessors[n].clipping = false;
                        }
                    }
                }

                for (int n = 0; n < _nc; n++)
                {
                    // dither processor does the float-to-PCM conversion
                    // (dither is not applied to 32-f output, only to PCM)
                    double val = sample[n];

                    if (_gain != 0 && !double.IsNaN(_gain))
                    {
                        val *= _gain;
                    }
                    if (_gains != null && !double.IsNaN(_gains[n]))
                    {
                        val *= _gains[n];
                    }

                    switch (_bitsPerSample)
                    {
                        case 8:
                            _w.Write((byte)(_ditherProcessors[n].process(val) + 127));
                            break;
                        case 16:
                            _w.Write((short)_ditherProcessors[n].process(val));
                            break;
                        case 24:
                            // Little-endian, signed 24-bit
                            int c = _ditherProcessors[n].process(val);
                            _w.Write((ushort)(c & 0xFFFF));
                            _w.Write((sbyte)(c >> 16 & 0xFF));
                            break;
                        case 32:
                            if (_audioFormat == WaveFormat.PCM || _audioFormat == WaveFormat.EXTENSIBLE)
                            {
                                _w.Write((int)_ditherProcessors[n].process(val));
                            }
                            else if (_audioFormat == WaveFormat.IEEE_FLOAT)
                            {
                                // Internals are double; just cast to float and discard any extra resolution
                                // (really we should dither here too, to approx 24 bits?)
                                _w.Write((float)val);
                            }
                            break;
                        case 64:
                            // we only do float, not PCM, 64-bits
                            _w.Write((double)val);
                            break;
                        default:
                            throw new Exception(String.Format("Bits per sample cannot be {0}", BitsPerSample));
                    }
                }

                err = false;
                if (_isConsole)
                {
                    // Check the stdout stream is still alive
                    int Err = System.Runtime.InteropServices.Marshal.GetLastWin32Error();
                    if (Err != 0)
                    {
                        // Err 997: "Overlapped I/O is in progress" (don't know cause)
                        // Err 183: cannot create a file... caused in Trace
                        // Err 2:   cannot find a file... caused in Trace
                        if (Err != 997 && Err != 183 && Err != 2)
                        {
                            System.ComponentModel.Win32Exception e = new System.ComponentModel.Win32Exception(Err);
                            Trace.WriteLine("Write fault " + Err + ": " + e.Message);
                            err = true;// yield break;
                        }
                    }
                }
            }
            catch (Exception e)
            {
                if (DSPUtil.IsMono() && e.Message.Contains("Write fault on path"))
                {
                    // This is the usual end-of-stream error on Mono
                    Trace.WriteLine("Write finished; " + e.Message);
                    err = true; // yield break
                }
                else if (e.GetHashCode() == 33574638)
                {
                    // "The specified network name is no longer available", from softsqueeze
                    Trace.WriteLine("Write finished; " + e.Message);
                    err = true; // yield break
                }
                else
                {
                    // Trace.WriteLine("Interrupted: " + e.Message);
                    throw e;
                }
            }
            return sample;
        }

        /// <summary> Number of iterations expected to do the signal processing </summary>
        public override int Iterations
        {
            get
            {
                if (_iterations == 0 && base._input!=null)
                {
                    return (base._input.Iterations);
                }
                return _iterations;
            }
        }

        /// <summary> Gets the number of bits per sample of the signal </summary>
        public ushort BitsPerSample
        {
            get { return _bitsPerSample; }
            set { _bitsPerSample = value; }
        }

        private void MakeDither()
        {
            _ditherProcessors = new Dither[NumChannels];
            for (int j = 0; j < NumChannels; j++)
            {
                _ditherProcessors[j] = new Dither(_dither, SampleRate, BitsPerSample);
            }
        }

        /// <summary>
        /// Set or get the dither type
        /// </summary>
        public DitherType Dither
        {
            get
            {
                return _dither;
            }
            set
            {
                // Dither can be changed at runtime...
                bool recreate = (value != _dither);
                _dither = value;
                if (recreate)
                {
                    // Create new dither-processors.
                    MakeDither();
                }
            }
        }

        /// <summary>
        /// Set or get the wave format
        /// </summary>
        public WaveFormat Format
        {
            get
            {
                return _audioFormat;
            }
            set
            {
                _audioFormat = value;
//                if (_audioFormat == WaveFormat.IEEE_FLOAT)
//                {
//                    // Float is always 32-bit
//                    BitsPerSample = 32;
//                }
                /* else */
                if (_audioFormat == WaveFormat.INTERNAL_DOUBLE)
                {
                    // Double is always 64-bit
                    BitsPerSample = 64;
                }
            }
        }

        /// <summary> Set or get the extensible Wave subtype </summary>
        public WaveFormatEx FormatEx
        {
            get { return _formatEx; }
            set
            {
                _audioFormat = WaveFormat.EXTENSIBLE;
                _formatEx = value;
            }
        }

        public SpeakerChannelMask ChannelMask
        {
            get { return _channelMask; }
            set { _channelMask = value; }
        }

        public bool Raw
        {
            get
            {
                return _raw;
            }
            set
            {
                _raw = value;
            }
        }
        // Normalization (default: peak dBFS)
        /// <summary>
        /// Normalization (default: peak dBFS)
        /// (EXPENSIVE)
        /// </summary>
        public double Normalization
        {
            get
            {
                return _normalization;
            }
            set
            {
                _normalization = value;
            }
        }
        /// <summary>
        /// Type of normalization
        /// </summary>
        public NormalizationType NormalizationType
        {
            get
            {
                return _normType;
            }
            set
            {
                _normType = value;
            }
        }

        /// <summary>
        /// After writing some stuff - what's the peak level we wrote?
        /// </summary>
        public double dbfsPeak
        {
            get
            {
                double pk = -200;
                for (int j = 0; j < NumChannels; j++)
                {
                    pk = Math.Max(pk, _ditherProcessors[j].dbfsPeak);
                }
                return pk;
            }
        }

        /// <summary> Close the wave file </summary>
        public void Close()
        {
            if (_sampleCount != Iterations || !_doneHeader)
            {
                // Try rewind to write correct length into the header
                try
                {
                    if (_w.BaseStream.CanSeek)
                    {
                        _iterations = _sampleCount;
                        _w.Seek(0, SeekOrigin.Begin);
                        WriteWaveHeader();
                    }
                }
                catch (Exception)
                {
                }
            }

            // Flush and close silently
            try
            {
                _bs.Flush();
                _w.Close(); _w = null;
                _bs.Close(); _bs = null;
                if (_fs != null) { _fs.Close(); _fs = null; }
            }
            catch (Exception)
            {
            }
        }

    }
}
