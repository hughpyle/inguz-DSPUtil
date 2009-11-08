using System;
using System.Collections.Generic;
using System.Text;

// Copyright (c) 2006 by Hugh Pyle, inguzaudio.com

namespace DSPUtil
{
    public interface ISweepGenerator : ISoundObj
    {
        ISoundObj Inverse { get; }
    }

    [Serializable]
    public class SweepGenerator : SoundObj, ISweepGenerator
    {
        double _lengthSecs;       // length of just the sweep, not the whole file
        double _startFreq;
        double _endFreq;
        double _gap;
        bool _pulse;
        double _gain;
        bool _repeat;

        double _phi;

        public SweepGenerator(ushort numChannels, double lengthSecs, uint startFreq, uint endFreq, uint sampleRate, double gapSecs, bool pulse, double gain, bool repeat)
        {
            _lengthSecs = lengthSecs;
            _lengthSamples = (int)(sampleRate * ( lengthSecs + 2*gapSecs ) + (pulse ? (1.025f*sampleRate) : 0 ));
            _startFreq = startFreq;
            _endFreq = endFreq;
            NumChannels = numChannels;
            SampleRate = sampleRate;
            _gap = gapSecs;
            _pulse = pulse;
            _gain = gain;
            _repeat = repeat;
            _phi = 0;
        }

        public SweepGenerator(ushort numChannels, int lengthSamples, uint startFreq, uint endFreq, uint sampleRate, double gain, bool repeat)
        {
            _lengthSecs = lengthSamples / sampleRate;
            _lengthSamples = lengthSamples;
            _startFreq = startFreq;
            _endFreq = endFreq;
            NumChannels = numChannels;
            SampleRate = sampleRate;
            _gap = 0;
            _pulse = false;
            _gain = gain;
            _repeat = repeat;
            _phi = 0;
        }

        /// <summary>
        /// Get an iterator for samples
        /// </summary>
        public override IEnumerator<ISample> Samples
        {
            get
            {
                if (_pulse)
                {
                    // First a half-second gap before the pulse
                    int j = 0;
                    for (; j < SampleRate / 2; j++)
                    {
                        yield return new Sample(NumChannels);
                    }
                    // A reversed sweep, sub-second, 10kHz down to 100Hz
                    SoundObj sweep = new SweepGenerator(NumChannels, 0.025f, 100, 10000, SampleRate, 0.0f, false, _gain, false);
                    SoundObj rever = new Reverser();
                    rever.Input = sweep;
                    foreach (ISample sample in rever)
                    {
                        yield return sample;
                        j++;
                    }
                    // Then half-second gap from the *start* of the pulse
                    for (; j < SampleRate; j++)
                    {
                        yield return new Sample(NumChannels);
                    }
                }
                if (_gap > 0)
                {
                    // The gap before and after the sweep
                    for (int j = 0; j < (_gap * SampleRate); j++)
                    {
                        yield return new Sample(NumChannels);
                    }
                }

                // The sweep itself
                double logStart = Math.Log(_startFreq);
                double logEnd = Math.Log(_endFreq);
                bool more = true;
                while (more)
                {
                    for (int j = 0; j < _lengthSecs * SampleRate; j++)
                    {
                        Sample sample = new Sample(NumChannels);
                        double value = (Math.Sin(_phi));

                        // instantaneous frequency, radians/sec
                        double f = Math.Exp(logStart + j * (logEnd - logStart) / (_lengthSecs * SampleRate));

                        double delta = 2 * Math.PI * f / SampleRate;
                        _phi += delta;

                        value *= _gain;
                        for (int c = 0; c < NumChannels; c++)
                        {
                            sample[c] = value;
                        }
                        yield return sample;
                    }
                    more = _repeat;
                }
            }
        }

        /// <summary> Number of iterations expected to do the signal processing </summary>
        int _lengthSamples;
        public override int Iterations
        {
            get { return (_lengthSamples); }
        }

                /// <summary>
        /// Get an iterator for samples
        /// </summary>
        public ISoundObj Inverse
        {
            get
            {
                // Build an inverse filter, for deconvolution of the recorded sweep.
                // We only build inverse of the main sweep, and don't include pulse;
                // assume those will be removed by post-processing the recorded sweep.
                //
                // Per Farina swept sine distortion paper,
                // inverse of a log sweep is the same sweep, time-reversed, then
                // amplitude reduced by 6dB/oct from 0 to -6.log2(w2/w1)

                double dbEnd = -6 * Math.Log(_endFreq / _startFreq, 2);

                SoundObj revsw = new Reverser();
                revsw.Input = this;

                SoundObj adjust = new LinearDbEnvelope(0, dbEnd, this.Iterations);
                adjust.Input = revsw;

                return adjust;
            }
        }
    }

    public unsafe class FFTSweepGenerator : SoundObj, ISweepGenerator
    {
        double _lengthSecs;       // length of just the sweep, not the whole file
        double _startFreq;
        double _endFreq;
        double _gain;
        double _A;
        double _B;
        bool _repeat;
        Complex[] _data;
        bool _gotdata;

        public FFTSweepGenerator(ushort numChannels, int lengthSamples, uint startFreq, uint endFreq, uint sampleRate, double gain, bool repeat)
        {
            _lengthSecs = lengthSamples / sampleRate;
            _lengthSamples = lengthSamples;
            _startFreq = startFreq;
            _endFreq = endFreq;
            NumChannels = numChannels;
            SampleRate = sampleRate;
            _gain = gain;
            _repeat = repeat;
        }

        /// <summary> Number of iterations expected to do the signal processing </summary>
        int _lengthSamples;
        public override int Iterations
        {
            get { return (_lengthSamples); }
        }

        public override IEnumerator<ISample> Samples
        {
            get
            {
                ushort nc = NumChannels;
                ISoundObj sweep = CalculateSweep();
                bool more = true;
                while (more)
                {
                    foreach (Sample s in sweep)
                    {
                        if (nc == 1)
                        {
                            yield return s;
                        }
                        else if (nc == 2)
                        {
                            yield return new Sample2(s[0], s[0]);
                        }
                    }
                    more = _repeat;
                }
            }
        }

        public ISoundObj Inverse
        {
            get
            {
                CalculateSweep();

                int fftSize = _data.Length;
                int N = fftSize / 2;

                // Copy the data buffer (in case someone's reading from it)
                Complex[] data1 = new Complex[fftSize];

                for (int j = 0; j < N; j++)
                {
                    data1[j] = new Complex(_data[j].Re, 0);
                }

                // FFT the sweep data
                Fourier.FFT(fftSize, data1);
                

                Complex unity = new Complex(1, 0);
                Complex[] data2 = new Complex[fftSize];
                data2[N-1] = unity;
                Fourier.FFT(fftSize, data2);

                for (int j = 0; j < fftSize; j++)
                {
                    data1[j].idiv(data2[j]);
                }

                // IFFT
                Fourier.IFFT(fftSize, data1); //, 1/Math.Sqrt(n));
                ComplexBufferReader cbr = new ComplexBufferReader(data1, 0, fftSize );
                cbr.SampleRate = _sr;
                return cbr;
            }
        }

        double T(double f)
        {
            return _A + (_B * Math.Log(f));
        }
        double phi(double f)
        {
            return (2 * Math.PI * f * (T(f) - _B)) % (2 * Math.PI);
        }
        double mag(double f)
        {
            double ff = f - _startFreq;
            if (ff <= 0)
            {
                return (Math.Cos(Math.PI*ff/_startFreq)+1)/2;
            }
            double dbNow = 1 - (10 * Math.Log10(f / _startFreq));
            double gainNow = MathUtil.gain(dbNow);
            return gainNow;
        }
        ISoundObj CalculateSweep()
        {
            // Per http://www.anselmgoertz.de/Page10383/Monkey_Forest_dt/Manual_dt/aes-swp-english.pdf
            int fftSize = MathUtil.NextPowerOfTwo(_lengthSamples * 2);
            int N = fftSize / 2;

            // Center the sweep in time to reduce impact of its extremities
            double fNyq = _sr/2;
            double FStart = _startFreq;
            double FEnd = _endFreq;
            double SStart = (N - _lengthSamples)/2;
            double TStart = SStart / _sr;
            double TEnd = TStart + _lengthSecs;
            _B = (TEnd - TStart) / Math.Log(FEnd / FStart);
            _A = TStart - _B * Math.Log(FStart);

            // Make the complex spectrum
            //double ph = 0;
            double df = (double)_sr / N;
            double phiNyq = phi(fNyq);
            double phiAdj = phiNyq % (2 * Math.PI);

            if (!_gotdata)
            {
                _data = new Complex[fftSize];
                fixed (Complex* cdata = _data)
                {
                    for (int j = 0; j < N; j++)
                    {
                        int m = j + 1;
                        double f = (double)m * _sr / fftSize;
                        double ph = phi(f) - (f / fNyq) * phiAdj;
                        double v = mag(f);
                        double Re = Math.Cos(ph) * v;
                        double Im = Math.Sin(ph) * v;
                        _data[j] = new Complex(Re, Im);
                    }
                    Fourier.IFFT(fftSize, _data, Math.Sqrt(fftSize) * _gain * MathUtil.gain(20));
                }

                // Look for values beyond the end
                // whose magnitude greater than our allowed threshold;
                // if present, window then ifft then start again.

                // Below doesn't seem to converge well
                // so just window and be done
                /*
                double threshold = MathUtil.gain(-90);
                bool iterate = true;
                while (iterate)
                {
                    iterate = false;
                    for (n = (int)(TEnd * sr + SStart * 2); n < fftSize; n++)
                    {
                        if (_data[n].Magnitude > threshold)
                        {
                            iterate = true;
                            break;
                        }
                    }
                    if (iterate)
                    {
                        Blackman bh = new Blackman((int)(_lengthSamples / 2 + SStart), _lengthSamples / 200, _lengthSamples / 2);
                        bh.Input = cbr;
                        n=0;
                        foreach (ISample s in bh)
                        {
                            _data[n++] = new Complex(s[0],0);
                        }
                        Fourier.FFT(fftSize, _data);
                        for (n = 0; n < N; n++)
                        {
                            int m = n + 1;
                            double f = (double)m * sr / fftSize;
                            double ph = _data[n].Phase;
                            double v = mag(f);
                            double Re = Math.Cos(ph) * v;
                            double Im = Math.Sin(ph) * v;
                            _data[n] = new Complex(Re, Im);
                        }
                        Fourier.IFFT(fftSize, _data);
                    }
                }
                */
                
                CosWindow bh = new Hamming((int)(_lengthSamples / 2 + SStart), (int)(SStart), _lengthSamples / 2);
                IEnumerator<double> gains = bh.Gains;
                for (int j = 0; j < N; j++)
                {
                    gains.MoveNext();
                    double g = gains.Current;
                    _data[j].mul(g);
                    _data[fftSize - j - 1].mul(g);
                }
                
                _gotdata = true;
            }

            ComplexBufferReader cbr = new ComplexBufferReader(_data, 0, N);

            //bh.Input = cbr;

            return cbr;
        }
        
    }


    /*
[Serializable]
public class EQSweepGenerator : SoundObj
{
    // Generates a log-sweep, with amplitude controlled by parameters
    // When this sweep is deconvolved with a log-sweep, you get the impulse response
    // of the parametric filter defined by those parameters (and linear phase).

    private double _lengthSecs;
    private uint _startFreq;
    private uint _endFreq;

    public EQSweepGenerator(double lengthSecs)
    {
        _lengthSecs = lengthSecs;
        _startFreq = 20;
        _endFreq = 20480;
    }

    private FilterProfile _coeffs;
    public FilterProfile Coefficients
    {
        get { return _coeffs; }
        set { _coeffs = value; }
    }

    #region Overrides
    /// <summary> Number of iterations expected to do the signal processing </summary>
    public override int Iterations
    {
        get { return (int)(SampleRate * _lengthSecs); }
    }

    #endregion

    public SoundObj RawSweep
    {
        get
        {
            SoundObj source;
            source = new NoiseGenerator(NoiseType.WHITE_FLAT, NumChannels, _lengthSecs, SampleRate, 1.0, false);
//                source = new SweepGenerator("S", _lengthSecs, _startFreq, _endFreq, _sampleRate, 0, false, 1.0);
            SoundBuffer buff = new SoundBuffer(source);
            return buff;
        }
    }

    public SoundObj InverseSweep
    {
        get
        {
            SoundObj sweep = RawSweep;

            SoundObj revsw = new Reverser();
            revsw.Input = sweep;

            return revsw;
        }
    }

    /// <summary>
    /// Get an iterator for samples
    /// </summary>
    public override IEnumerator<ISample> Samples
    {
        get
        {
            double lengthSecs = 2;

            // Construct a sweep from 20Hz to 20kHz
            // (10 octaves; 20, 40, 80, 160, 320, 640, 1280, 2560, 5120, 10240, 20480)
            // over lengthSecs
            SoundObj sweep = RawSweep;
            double logStart = Math.Log(_startFreq);
            double logEnd = Math.Log(_endFreq);

            // Frequency (Hz) at sample j is
            // double omega = Math.Exp(logStart + j * (logEnd - logStart) / (lengthSecs * SampleRate));
            // f = 2.pi.omega
            //
            // ln(f/2pi) = logStart + j*(...)
            // j*(...) = ln(f/2pi) - logStart
            //
            // Sample j for frequency f is
            // j = (ln(f/2pi) - logStart) / ( (logEnd - logStart) / (lengthSecs * SampleRate) )

            // Construct shelf envelopes
            ShelfEnvelope rcsePrev = null;
            for(int n=1; n<_coeffs.Count; n++)
            {
                double freq1 = _coeffs[n-1].Freq;
                double freq2 = _coeffs[n].Freq;

                // gains in dB
                double gain1 = _coeffs[n - 1].Gain;
                double gain2 = _coeffs[n].Gain;

                // Since we're chaining these filters together,
                // we actually want to adjust for pre-sweep == 0dB, then end-sweep=(difference)
                double gainEnd = gain2 - gain1;

                double startSample = (Math.Log(freq1 * (2 * Math.PI)) - logStart) / ((logEnd - logStart) / (lengthSecs * SampleRate));
                double endSample = (Math.Log(freq2 * (2 * Math.PI)) - logStart) / ((logEnd - logStart) / (lengthSecs * SampleRate));

                ShelfEnvelope rcse = new ShelfEnvelope(0, gainEnd, (int)startSample, (int)endSample);

                // Input of this shelf is output of the previous shelf
                // (or, for the first shelf, the sweep)
                if (rcsePrev == null)
                {
                    rcse.Input = sweep;
                }
                else
                {
                    rcse.Input = rcsePrev;
                }
                rcsePrev = rcse;
            }

            // Ready
            foreach (ISample sample in rcsePrev)
            {
                yield return sample;
            }
        }
    }
}
*/
}
