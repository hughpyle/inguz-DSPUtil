using System;
using System.Collections.Generic;
using System.Text;
using System.IO;

// Copyright (c) 2006, 2007 by Hugh Pyle, inguzaudio.com
// WAV file reader based originally on Garbe.Sound

namespace DSPUtil
{
    public enum WaveFormat
    {
        // The only types of data we understand are...
        ANY = 0,        // Pseudo-type, means we don't know

        PCM = 1,
        ADPCM = 2,
        IEEE_FLOAT = 3,
        EXTENSIBLE = 65534,
        INTERNAL_DOUBLE = 65533     // Our invention, 64-bit doubles
    }


    public enum SpeakerChannelMask
    {
        none = 0,
        stereo = 0x3,
        itu51 = 0x3F,

        FRONT_LEFT = 0x1,
        FRONT_RIGHT = 0x2,
        FRONT_CENTER = 0x4,
        LOW_FREQUENCY = 0x8,
        BACK_LEFT = 0x10,
        BACK_RIGHT = 0x20,
        FRONT_LEFT_OF_CENTER = 0x40,
        FRONT_RIGHT_OF_CENTER = 0x80,
        BACK_CENTER = 0x100,
        SIDE_LEFT = 0x200,
        SIDE_RIGHT = 0x400,
        TOP_CENTER = 0x800,
        TOP_FRONT_LEFT = 0x1000,
        TOP_FRONT_CENTER = 0x2000,
        TOP_FRONT_RIGHT = 0x4000,
        TOP_BACK_LEFT = 0x8000,
        TOP_BACK_CENTER = 0x10000,
        TOP_BACK_RIGHT = 0x20000,
        /*RESERVED = 0x80000000*/
    }

    public class WaveFormatEx
    {
        public static WaveFormatEx PCM = new WaveFormatEx("00000001-0000-0010-8000-00aa00389b71");
        public static WaveFormatEx IEEE_FLOAT = new WaveFormatEx("00000003-0000-0010-8000-00aa00389b71");
        public static WaveFormatEx AMBISONIC_B_FORMAT_PCM = new WaveFormatEx("00000001-0721-11d3-8644-C8C1CA000000");
        public static WaveFormatEx AMBISONIC_B_FORMAT_IEEE_FLOAT = new WaveFormatEx("00000003-0721-11d3-8644-C8C1CA000000");
        public Guid guid;
        public WaveFormatEx(Guid g)
        {
            guid = g;
        }
        public WaveFormatEx(byte[] g)
        {
            guid = new Guid(g);
        }
        public WaveFormatEx(String g)
        {
            guid = new Guid(g);
        }
        public override bool Equals(object obj)
        {
            if (obj == null || GetType() != obj.GetType())
                return false;
            return (((WaveFormatEx)obj).guid.Equals(guid));
        }
        public override int GetHashCode()
        {
            return guid.GetHashCode();
        }
        public override string ToString()
        {
            return guid.ToString();
        }
    }


    /// <summary> Read a wave file </summary>
    [Serializable]
    public sealed class WaveReader : SoundObj /*, ISampleBuffer */
    {
        private bool _ok;
        private FileStream fs;
        private BinaryReader _rdr;
        private Stream bs;

        private string _filename;
        private string _riff;
        private uint _length;
        private string _wave;
        private string _format;

        private uint _size;
        private WaveFormat _audioFormat;
        private WaveFormatEx _formatEx;
        private uint _byteRate;
        private ushort _blockAlign;
        private ushort _bitsPerSample;
        private bool _bigEndian = false;
        private bool _isSPDIF = false;

        private string _data;
        private uint _dataSize;
        private uint _channelMask = 0x3; // default to stereo, unless WAVEFORMATEXTENSIBLE says otherwise

        private long _pos; // position relative to start of data
        private long _max; // maximum position (number of samples in stream), or uint.MaxValue for an endless stream
        private long _seekpos;  // byte position of the first data
        private ISample _first;
        private bool _moreThanFirst;
        private ISample _current;

        // Buffers for ISampleBuffer
        /*
        private ISample[] _buff;
        private int _bufflen;
        private Complex[][] _cbuff;
        private int _cbufflen;
        */

        private const double _scale8 = 1 / 128f;
        private const double _scale16 = 1 / 32768f;
        private const double _scale24 = 1 / 8388608f;
        private const double _scale32 = 1 / 2147483648f;

        #region Constructors

        /// <summary> Read a wave file </summary>
        /// <param name="fileName">Name of the wave file</param>
        public WaveReader(string fileName)
        {
            OpenFile(fileName);
            ReadWaveHeader(WaveFormat.ANY, true);
            ReadSPDIF();
        }
        public WaveReader(string fileName, TimeSpan startTime)
        {
            OpenFile(fileName);
            ReadWaveHeader(WaveFormat.ANY, true);
            ReadSPDIF();
            SkipToStart(startTime);
        }
        public WaveReader(string fileName, WaveFormat format)
        {
            OpenFile(fileName);
            _audioFormat = format;
            ReadWaveHeader(format, true);
            ReadSPDIF();
        }
        public WaveReader(string fileName, WaveFormat format, ushort bitsPerSample, ushort numChannels)
        {
            // To read raw
            OpenFile(fileName);
            _audioFormat = format;
            ReadWaveHeader(format, false);
            NumChannels = numChannels;
            _bitsPerSample = bitsPerSample;
            ReadSPDIF();
        }
        public WaveReader(string fileName, WaveFormat format, ushort bitsPerSample, ushort numChannels, TimeSpan startTime)
        {
            // To read raw
            OpenFile(fileName);
            _audioFormat = format;
            ReadWaveHeader(format, false);
            NumChannels = numChannels;
            _bitsPerSample = bitsPerSample;
            ReadSPDIF();
            SkipToStart(startTime);
        }

        public WaveReader(Stream input)
        {
            //fs = null;
            //bs = null;
            _rdr = new BinaryReader(input);
            ReadWaveHeader(WaveFormat.ANY, true);
            ReadSPDIF();
        }
        public WaveReader(Stream input, WaveFormat format)
        {
            _audioFormat = format;
            //fs = null;
            //bs = null;
            _rdr = new BinaryReader(input);
            ReadWaveHeader(format, true);
            ReadSPDIF();
        }

        private void OpenFile(string fileName)
        {
            _filename = fileName;
            if (fileName == null || fileName=="-")
            {
                // use stdin
                //fs = null;
                //bs = null;
                //              Trace.WriteLine("Read stdin");
                Stream stdin = System.Console.OpenStandardInput();
                bs = new BufferedStream(stdin);
                _rdr = new BinaryReader(bs);
            }
            else if (File.Exists(fileName))
            {
                //              Trace.WriteLine("Read {0}", fileName);
                fs = new FileStream(fileName, FileMode.Open, FileAccess.Read, FileShare.Read, 65536, true);
                bs = new BufferedStream(fs);
                _rdr = new BinaryReader(bs);
            }
            else
            {
                throw (new FileNotFoundException("File " + fileName + " not found"));
            }
        }


        private void ReadWaveHeader(WaveFormat format, bool expectHeader)
        {
            _ok = false;
            _pos = 0;

            //          Trace.WriteLine("ReadWaveHeader {0}", format);
            if (!expectHeader)
            {
                // Input file is raw data
                _riff = "(no header)";
                _length = 0;
                _wave = "(no header)";
                _format = "(no header)";
                _size = 16;
                _data = "data";
                _dataSize = 0;
                _max = uint.MaxValue;
                _audioFormat = format;

                // Assume defaults, they can be overridden later
                NumChannels = 2;
                SampleRate = 44100;

                if (format == WaveFormat.PCM || format == WaveFormat.EXTENSIBLE)
                {
                    // Raw PCM, assume 16 bit 44k1 stereo PCM
                    _bitsPerSample = 16;
                    _ok = true;
                }
                else if (format == WaveFormat.IEEE_FLOAT)
                {
                    // IEEE Float; assume 32 bit 44k1 stereo IEEE_FLOAT
                    _bitsPerSample = 32;
                    _ok = true;
                }
                else if (format == WaveFormat.INTERNAL_DOUBLE)
                {
                    // our 64 bit 44k1 stereo 'double-precision'
                    _bitsPerSample = 32;
                    _ok = true;
                }

                _blockAlign = (ushort)((NumChannels * _bitsPerSample) >> 3);
                _byteRate = (uint)(_blockAlign * SampleRate);

                return;
            }

            // Read 'RIFF' (WAV) or 'FORM' (AIFF) tag ///////////////////////////////////////////////

            char[] hdr = _rdr.ReadChars(4);
            _riff = new string(hdr);
            if (_riff != "RIFF" && _riff != "FORM")
            {
                if (hdr.Length == 0)
                {
                    throw (new Exception("File could not be read: no data."));
                }
                string x = "";
                for (int j = 0; j < hdr.Length; j++)
                {
                    x += String.Format("{0:X} ", (int)hdr[j]);
                }
                throw (new Exception(String.Format("File is not WAV: no 'RIFF' tag found, instead '{0}'.", x)));
            }

            // File length
            int fileLen = _rdr.ReadInt32();
            _length = (uint)fileLen;

            // Read Wave //////////////////////////////////////////////

            _wave = new string(_rdr.ReadChars(4));
            if (_wave != "WAVE" && _wave != "AIFF")
                throw (new Exception(String.Format("File is not WAV: no 'WAVE' tag found, instead {0}", _wave)));
            if (_wave == "AIFF")
            {
                // The whole file is big-endian, including lengths in the header
                BigEndian = true;
                _length = (uint)System.Net.IPAddress.NetworkToHostOrder(fileLen);
            }
            // Read Format ////////////////////////////////////////////

            _format = new string(_rdr.ReadChars(4));

            // WMP11-ripped WAV files begin with 'LIST' metadata - even before the format tag.
            // AIFF files could have this too (haven't seen it yet).
            // Skip any chunks up to the format header.
            while (_format.Length > 0 && _format != "fmt " && _format != "COMM")
            {
                int chunkSize = _rdr.ReadInt32();
                if (BigEndian)
                    chunkSize = System.Net.IPAddress.NetworkToHostOrder(chunkSize);
                Trace.WriteLine("Skipping {0} ({1} bytes)", _format, chunkSize);
                _rdr.ReadBytes(chunkSize);
                _format = new string(_rdr.ReadChars(4));
            }

            if (_format == "fmt ")
            {
                // WAV file-format chunk
                _size = _rdr.ReadUInt32();
                if (_size < 16)
                    throw (new Exception("File could not be read: don't know how to read 'fmt' size " + _size));

                _audioFormat = (WaveFormat)_rdr.ReadUInt16();

                if (_audioFormat == WaveFormat.PCM ||
                    _audioFormat == WaveFormat.ADPCM ||
                    _audioFormat == WaveFormat.IEEE_FLOAT ||
                    _audioFormat == WaveFormat.INTERNAL_DOUBLE ||
                    _audioFormat == WaveFormat.EXTENSIBLE)
                {
                    //                                      // WAVEFORMATEX wFormatTag          2 bytes
                    NumChannels = _rdr.ReadUInt16();        // WAVEFORMATEX nChannels           2
                    SampleRate = _rdr.ReadUInt32();         // WAVEFORMATEX nSamplesPerSec      4
                    _byteRate = _rdr.ReadUInt32();          // WAVEFORMATEX nAvgBytesPerSec     4
                    _blockAlign = _rdr.ReadUInt16();        // WAVEFORMATEX nBlockAlign         2 (channels * bitspersample / 8)
                    _bitsPerSample = _rdr.ReadUInt16();     // WAVEFORMATEX wBitsPerSample      2 (the *container* size)
                    if (_size > 16)
                    {
                        uint skip = 16;
                        if (_audioFormat == WaveFormat.EXTENSIBLE)
                        {
                            UInt16 kip = _rdr.ReadUInt16();
                            UInt16 union = _rdr.ReadUInt16();       // the Samples union, wdc
                            _channelMask = _rdr.ReadUInt32();       // channel mask
                            // then the GUID, 16 bytes
                            _formatEx = new WaveFormatEx(_rdr.ReadBytes(16));
                            skip = 40;
                            if (_formatEx.Equals(WaveFormatEx.PCM) || _formatEx.Equals(WaveFormatEx.AMBISONIC_B_FORMAT_PCM))
                            {
                                _audioFormat = WaveFormat.PCM;
                            }
                            else if (_formatEx.Equals(WaveFormatEx.IEEE_FLOAT) || _formatEx.Equals(WaveFormatEx.AMBISONIC_B_FORMAT_IEEE_FLOAT))
                            {
                                _audioFormat = WaveFormat.IEEE_FLOAT;
                            }
                        }
                        // Read and discard the rest of the 'fmt' structure
                        _rdr.ReadBytes((int)(_size - skip));
                    }
                }
                else
                {
                    throw (new Exception("File could not be read: don't know how to read audio format " + _audioFormat));
                }
            }
            else if (_format == "COMM")
            {
                // AIFF file-format chunk
                _size = (uint)System.Net.IPAddress.NetworkToHostOrder(_rdr.ReadInt32());
                if (_size < 18)
                    throw (new Exception("File could not be read: don't know how to read 'COMM' size " + _size));

                _audioFormat = WaveFormat.PCM;
                NumChannels = (ushort)System.Net.IPAddress.NetworkToHostOrder(_rdr.ReadInt16());
                uint numFrames = (uint)System.Net.IPAddress.NetworkToHostOrder(_rdr.ReadInt32()); // number of sample frames
                _bitsPerSample = (ushort)System.Net.IPAddress.NetworkToHostOrder(_rdr.ReadInt16());

                // SampleRate is 10-byte IEEE_extended format.  Don't bother converting that
                // properly, just check for good known values (yuk!)
                byte[] ext = _rdr.ReadBytes(10);
                if (ext[0] == 64 && ext[1] == 14 && ext[2] == 172 && ext[3] == 68)
                {
                    SampleRate = 44100;
                }
                else if (ext[0] == 64 && ext[1] == 14 && ext[2] == 187 && ext[3] == 128)
                {
                    SampleRate = 48000;
                }
                else if (ext[0] == 64 && ext[1] == 15 && ext[2] == 187 && ext[3] == 128)
                {
                    SampleRate = 96000;
                }
                else
                {
                    throw (new Exception("File could not be read: don't know how to interpret sample rate."));
                }

                _blockAlign = (ushort)((NumChannels * _bitsPerSample) / 8);
                _byteRate = (uint)(_blockAlign * SampleRate);

                if (_size > 18)
                {
                    // Read and discard the rest of the 'fmt' structure
                    _rdr.ReadBytes((int)(_size - 18));
                }
            }
            else
            {
                throw (new Exception(String.Format("File could not be read: no 'fmt' tag found, instead {0}", _format)));
            }


            // Read Data ///////////////////////////////////////////////

            _data = new string(_rdr.ReadChars(4));
            while (_data.Length > 0 && _data != "data" && _data != "SSND")
            {
                // Not a data chunk, ignore
                int miscSize = _rdr.ReadInt32();
                if (BigEndian)
                    miscSize = System.Net.IPAddress.NetworkToHostOrder(miscSize);
                _rdr.ReadBytes(miscSize);
                _data = new string(_rdr.ReadChars(4));
            }

            // Read the data size
            if (BigEndian)
            {
                _dataSize = (uint)System.Net.IPAddress.NetworkToHostOrder(_rdr.ReadInt32());
            }
            else
            {
                _dataSize = _rdr.ReadUInt32();
            }

            // See if we can read this
            if (NumChannels > 0)
            {
                switch (_audioFormat)
                {
                    case WaveFormat.PCM:
                    case WaveFormat.EXTENSIBLE:
                        switch (_bitsPerSample)
                        {
                            case 8:
                            case 16:
                            case 24:
                            case 32:
                                _ok = true;
                                break;
                            default:
                                break;
                        }
                        break;
                    case WaveFormat.IEEE_FLOAT:
                        switch (_bitsPerSample)
                        {
                            case 32:
                            case 64:
                                _ok = true;
                                break;
                        }
                        break;
                    case WaveFormat.INTERNAL_DOUBLE:
                        switch (_bitsPerSample)
                        {
                            case 64:
                                _ok = true;
                                break;
                        }
                        break;
                    default:
                        break;
                }
            }
            if (_ok)
            {
                _max = (uint)((_dataSize / (_bitsPerSample / 8)) / NumChannels);
                if (_dataSize==4294967292)
                {
                    Trace.WriteLine("Wave file: unknown size from header");
                    _max = uint.MaxValue;
                }
                if (_max == 0)
                {
                    Trace.WriteLine("Wave file: zero size from header");
                    _max = uint.MaxValue;
                }
            }
        }

        private void ReadSPDIF()
        {
            if (_ok)
            {
                byte[] spdif = { 0x72, 0xf8, 0x1f, 0x4e };
                if (_rdr.BaseStream.CanSeek)
                {
                    _seekpos = _rdr.BaseStream.Position;
                }

                // Read the first sample "manually"
                // so we can check whether there's a SPDIF or other magic number
                // at the beginning of the stream.
                int nFirst = (_nc * _bitsPerSample / 8);
                byte[] firstBytes = _rdr.ReadBytes(nFirst);
                MemoryStream ms = new MemoryStream(firstBytes);
                BinaryReader mr = new BinaryReader(ms);

                // Save the sample for later use
                _first = Next(mr, out _moreThanFirst);

                if (firstBytes.Length >= nFirst)
                {
                    // Check whether it's SPDIF-wrapped
                    _isSPDIF = true;
                    for (int b = 0; _isSPDIF && b < spdif.Length && b < nFirst; b++)
                    {
                        _isSPDIF &= (spdif[b] == firstBytes[b]);
                    }
                }
            }
        }

        private void SkipToStart(TimeSpan ts)
        {
            if (_ok)
            {
                // Number of samples = seconds * samplerate
                long pos = (long)(ts.TotalSeconds * _sr);
                if (pos > 0)
                {
                    if (_rdr.BaseStream.CanSeek)
                    {
                        Trace.WriteLine("Skip to time {0} (sample {1})", ts, pos);
                        _rdr.BaseStream.Seek(_seekpos + (pos * _nc * (_bitsPerSample / 8)), SeekOrigin.Begin);

                        // remember our new base-position, so Seek() works relative to this
                        _seekpos = _rdr.BaseStream.Position;

                        // Read the first sample again (see ReadSPDIF)
                        int nFirst = (_nc * _bitsPerSample / 8);
                        byte[] firstBytes = _rdr.ReadBytes(nFirst);
                        MemoryStream ms = new MemoryStream(firstBytes);
                        BinaryReader mr = new BinaryReader(ms);
                        _first = Next(mr, out _moreThanFirst);
                    }
                    else
                    {
                        // Laboriously skip samples until we reach the right place?
                        // Bah.  TODO.  For now we only care about file-input sources, which are seekable.
                    }
                }
            }
        }

        #endregion

//        public override void Reset()
//        {
//           Seek(1);
//        }

        /// <summary>
        /// Seek to position in the data stream
        /// </summary>
        /// <param name="pos">Sample position, 0-based</param>
        private void Seek(long pos)
        {
            if (_pos != pos)
            {
                if (!_rdr.BaseStream.CanSeek)
                {
                    string msg = String.Format("Cannot rewind this input stream (from {0} to {1}).", _pos, pos);
                    throw new NotSupportedException(msg);
                }
                _rdr.BaseStream.Seek(_seekpos + (pos * _nc * (_bitsPerSample / 8)), SeekOrigin.Begin);
                _pos = pos;
            }
        }

        private double NextDouble(BinaryReader rdr)
        {
            double val = 0;
            if (_bigEndian)
            {
                // For now only handle 16-bit big-endian data   
                if (_bitsPerSample == 16)
                {
                    short beword = rdr.ReadInt16();
                    beword = System.Net.IPAddress.NetworkToHostOrder(beword);
                    val = ((double)beword * _scale16);
                }
            }
            else
            {
                switch (_bitsPerSample)
                {
                    case 16:
                        val = ((double)rdr.ReadInt16() * _scale16);
                        break;
                    case 8:
                        val = ((double)rdr.ReadByte() - 128) * _scale8;        // 8-bit PCM uses unsigned bytes
                        break;
                    case 24:
                        // Little-endian, signed 24-bit
                        int a = (int)rdr.ReadUInt16();
                        int b = (int)rdr.ReadSByte();
                        int c = (b << 16) + a;
                        val = ((double)c) * _scale24;
                        break;
                    case 32:
                        if (_audioFormat == WaveFormat.IEEE_FLOAT)
                        {
                            val = (double)rdr.ReadSingle();
                        }
                        else
                        {
                            val = ((double)rdr.ReadInt32()) * _scale32;
                        }
                        break;
                    case 64:
                        if ((_audioFormat == WaveFormat.IEEE_FLOAT) || (_audioFormat == WaveFormat.INTERNAL_DOUBLE))
                        {
                            val = rdr.ReadDouble();
                        }
                        else
                        {
                            // throw new Exception("64-bit PCM not handled");
                            val = 0;
                        }
                        break;
                }
            }
            return val;
        }

        private ISample First(out bool more)
        {
            // Looking for spdif, we cached the first sample (and whether there were more samples after that)
            more = _moreThanFirst;
            return _first;
        }

        private ISample This(out bool more)
        {
            // The current sample
            more = (_current != null);
            return _current;
        }

        private ISample Next(BinaryReader rdr, out bool more)
        {
            try
            {
                if (_nc == 2)
                {
                    if (_pos < _max)
                    {
                        _pos++;
                        double a = NextDouble(rdr);
                        double b = NextDouble(rdr);
                        more = true;
                        _current = new Sample2(a, b);
                    }
                    else
                    {
                        more = false;
                        _current = null;
                    }
                }
                else
                {
                    if (_pos < _max)
                    {
                        _pos++;
                        ISample sample = new Sample(_nc);
                        for (int n = 0; n < _nc; n++)
                        {
                            sample[n] = NextDouble(rdr);
                        }
                        more = true;
                        _current = sample;
                    }
                    else
                    {
                        more = false;
                        _current = null;
                    }
                }
            }
            catch (EndOfStreamException)
            {
                Trace.WriteLine("End of input ({0}) at {1}.", streamName, _pos);
                more = false;
                _current = null;
            }
            catch (IOException e)
            {
                Trace.WriteLine("IO error ({0}) at {1}: {2}", streamName, _pos, e.Message);
                more = false;
                _current = null;
            }
            return _current;
        }

        private string streamName
        {
            get
            {
                if (String.IsNullOrEmpty(_filename))
                {
                    return "stream";
                }
                return Path.GetFileName(_filename);
            }
        }

        /// <summary>
        /// Get an iterator for samples
        /// </summary>
        public override IEnumerator<ISample> Samples
        {
            get
            {
                if (!_ok)
                {
                    Trace.WriteLine("WaveReader: Cannot process");
                    yield break;
                }
                // Reset();
                bool more = true;
                long pos = 0;
                ISample s = First(out more);
                if (s!=null)
                {
                    pos++;
                    yield return s;
                }
                while (more)
                {
                    if (_pos != pos)
                    {
                        // Someone else moved the stream forward!
                        // Seek to the right place, please
                        if (_pos == pos + 1)
                        {
                            s = This(out more);
                        }
                        else
                        {
                            Seek(pos);
                            s = Next(_rdr, out more);
                        }
                    }
                    else
                    {
                        s = Next(_rdr, out more);
                    }
                    if (s != null)
                    {
                        pos++;
                        yield return s;
                    }
                }
                yield break;
            }
        }

        #region ISampleBuffer implementation
        /*
        public void Skip(int n, out int nn, out bool moreSamples)
        {
            Read(n, out nn, out moreSamples);
        }
        public ISample[] Read(int n, out int nn, out bool moreSamples)
        {
            moreSamples = true;
            int j;
            if (_buff == null || _bufflen < n)
            {
                _buff = new ISample[n];
                _bufflen = n;
            }
            for (j = 0; j < n && moreSamples; j++)
            {
                _buff[j] = Next(_rdr, out moreSamples);
            }
            nn = moreSamples ? j : (j - 1);
            return _buff;
        }

        /// <summary>
        /// Read, into an array of Complex.
        /// Unused elements in the array are guaranteed null.
        /// </summary>
        /// <param name="n"></param>
        /// <param name="nn"></param>
        /// <param name="moreSamples"></param>
        /// <returns></returns>
        public Complex[][] ReadComplex(int n, out int nn, out bool moreSamples)
        {
            moreSamples = true;
            int j = 0;
            ushort nc = NumChannels;
            if (_cbuff == null || _cbufflen < n)
            {
                _cbuff = new Complex[nc][];
                for (ushort c = 0; c < nc; c++)
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
            for (j = 0; j < n && moreSamples; j++)
            {
                ISample s = Next(_rdr, out moreSamples);
                if (s != null)
                {
                    for (ushort c = 0; c < nc; c++)
                    {
                        _cbuff[c][j].Re = s[c];
                    }
                }
            }
            nn = moreSamples ? j : (j - 1);
            return _cbuff;
        }
        */
        #endregion

        /// <summary> Number of iterations expected to do the signal processing == number of samples </summary>
        public override int Iterations
        {
            get { return ((int)_max); }
        }

        #region Methods

        /// <summary> Gets the number of bits per sample of the signal </summary>
        public ushort BitsPerSample
        {
            get { return _bitsPerSample; }
        }

        /// <summary> Close the wave file </summary>
        public void Close()
        {
            if (_rdr != null) { _rdr.Close(); }
            if (bs != null) { bs.Close(); bs = null; }
            if (fs != null) { fs.Close(); fs = null; }
        }

        #endregion


        #region Miscellaneous Properties

        public bool IsSPDIF
        {
            get { return _isSPDIF; }
        }

        // Set big-endian-ness (only useful for raw, since we always read the header as little-endian)
        public bool BigEndian
        {
            set { _bigEndian = value; }
            get { return _bigEndian; }
        }

        /// <summary> Get the RIFF tag </summary>
        public string RiffTag
        {
            get { return _riff; }
        }

        /// <summary> Get the length tag </summary>
        public uint RiffLength
        {
            get { return _length; }
        }

        /// <summary> Get the Wave tag </summary>
        public string WaveTag
        {
            get
            { return _wave; }
        }

        /// <summary> Get the Format tag </summary>
        public string FormatTag
        {
            get
            { return _format; }
        }

        /// <summary> Get the Size tag </summary>
        public uint FormatSize
        {
            get { return _size; }
        }

        /// <summary> Get the Audio Format tag </summary>
        public WaveFormat Format
        {
            get { return _audioFormat; }
        }

        /// <summary> Get the extensible Wave subtype </summary>
        public WaveFormatEx FormatEx
        {
            get { return _formatEx; }
        }

        /// <summary> Get the Byte Rate tag </summary>
        public uint ByteRate
        {
            get { return _byteRate; }
        }

        /// <summary> Get the Block Align tag </summary>
        public ushort BlockAlign
        {
            get { return _blockAlign; }
        }

        /// <summary> Get the Data tag </summary>
        public string Data
        {
            get { return _data; }
        }

        /// <summary> Get the Data Size tag </summary>
        public uint DataSize
        {
            get { return _dataSize; }
        }

        /// <summary> Get the channel mask (default to stereo) </summary>
        public uint ChannelMask
        {
            get { return _channelMask; }
        }
        #endregion
    }
}