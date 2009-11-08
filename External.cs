using System;
using System.Collections.Generic;
using System.Text;
using System.IO;

// Copyright (c) 2006, 2007 by Hugh Pyle, inguzaudio.com

namespace DSPUtil
{
    public class SPDIFWrappedExternal : SoundObj
    {
        string _exeName;
        string _exeArgs;
        WaveWriter _writer = null;
        System.Diagnostics.Process _process = null;
        StreamWriter _in = null;
        StreamReader _out = null;

        double _gain;
        double[] _gains = null;
        ushort _bitsPerSample;
        DitherType _dither;
        WaveFormat _audioFormat;
        WaveFormatEx _formatEx;

        public SPDIFWrappedExternal(string exeName, string exeArgs)
        {
            _exeName = exeName;
            _exeArgs = exeArgs;
        }

        /// <summary>
        /// Gain, units
        /// </summary>
        public double Gain
        {
            set
            {
                _gain = value;
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
            if (_gains == null)
            {
                _gains = new double[_nc];
                for (ushort c = 0; c < _nc; c++)
                {
                    _gains[c] = double.NaN;
                }
            }
            _gains[channel] = gain;
        }

        /// <summary> Set the number of bits per sample of the signal </summary>
        public ushort BitsPerSample
        {
            set { _bitsPerSample = value; }
        }


        /// <summary>
        /// Set the dither type
        /// </summary>
        public DitherType Dither
        {
            set
            {
                _dither = value;
            }
        }

        /// <summary>
        /// Set the wave format
        /// </summary>
        public WaveFormat Format
        {
            set
            {
                _audioFormat = value;
            }
        }

        /// <summary> Set the extensible Wave subtype </summary>
        public WaveFormatEx FormatEx
        {
            set
            {
                _audioFormat = WaveFormat.EXTENSIBLE;
                _formatEx = value;
            }
        }

        /// <summary>
        /// Number of channels output is always 2, but beware: these are not
        /// audio samples anymore, they're SPDIF chunk pairs.
        /// </summary>
        public override ushort NumChannels
        {
            get
            {
                return 2;
            }
        }

        private bool valid()
        {
            bool ok = false;
            if (_input != null)
            {
                ok = true;
                if (_sr != 48000 && _sr != 44100 && _sr != 32000)
                {
                    throw new Exception("Invalid sample rate for SPDIF");
                }
            }
            return ok;
        }

        /// <summary>
        /// Get an iterator for samples
        /// </summary>
        public override IEnumerator<ISample> Samples
        {
            get
            {
                if (!valid())
                {
                    yield break;
                }
                // Spin up the external process and start feeding it data
                Start();
                // Read the external process' output and wrap in SPDIF
                foreach (ISample sample in _input)
                {
                    if (_process.HasExited)
                    {
                        break;
                    }
                    yield return sample;
                }
                Stop();
            }
        }

        private void Start()
        {
            _process = new System.Diagnostics.Process();
            System.Diagnostics.ProcessStartInfo processInfo = new System.Diagnostics.ProcessStartInfo();
            processInfo.Arguments = _exeArgs;
            processInfo.FileName = _exeName;
            processInfo.UseShellExecute = false;
            processInfo.RedirectStandardInput = true;
            processInfo.RedirectStandardOutput = true;
            processInfo.RedirectStandardError = false;
            _in = _process.StandardInput;
            _out = _process.StandardOutput;

            _writer = new WaveWriter(_in.BaseStream);
            _writer.Input = _input;
            _writer.BitsPerSample = 16;

            _process.StartInfo = processInfo;
            try
            {
                _process.Start();
            }
            catch (Exception e)
            {
                Trace.WriteLine("That didn't seem to work: {0}", e.Message);
                _process = null;
            }
        }

        private void Stop()
        {
        }
    }
}
