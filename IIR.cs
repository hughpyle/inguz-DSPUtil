using System;
using System.Collections.Generic;
using System.Text;

// Copyright (c) 2006, 2007 by Hugh Pyle, inguzaudio.com

namespace DSPUtil
{
    /// <summary>
    /// Base class for first-order IIR filters
    /// If input is null, produces the Dirac response.
    /// </summary>
    public abstract class IIR1 : SoundObj
    {
        protected int _n;
        protected double _a1, _b0, _b1;

        public override ushort NumChannels
        {
            get
            {
                return _input == null ? (ushort)1 : _input.NumChannels;
            }
        }

        public override int Iterations
        {
            get { return _n; }
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
                    double a = 1.0;
                    double z0 = 0;
                    double z1 = 0;
                    for (int j = 0; j < _n; j++)
                    {
                        double v = a * _b0 + z0 * _b1 - z1 * _a1;
                        z1 = v;
                        z0 = a;
                        a = 0;
                        yield return new Sample(v);
                    }
                }
                else
                {
                    int nc = _input.NumChannels;
                    double[] z0 = new double[nc];
                    double[] z1 = new double[nc];
                    foreach (ISample sample in _input)
                    {
                        ISample s = sample;
                        for (int c = 0; c < nc; c++)
                        {
                            double a = s[c];
                            double v = a * _b0 + z0[c] * _b1 - z1[c] * _a1;
                            z1[c] = v;
                            z0[c] = a;
                            s[c] = v;
                        }
                        yield return sample;
                    }
                }
            }
        }
    }

    /// <summary>
    /// First-order high pass IIR
    /// </summary>
    public class IIR1HP : IIR1
    {
        /// <summary>
        /// Construct a generator for first-order high pass IIR
        /// </summary>
        /// <param name="fS">Sample rate, Hz</param>
        /// <param name="fC">Cutoff frequency, Hz</param>
        /// <param name="n">Length of impulse (samples)</param>
        public IIR1HP(uint fS, double fC, int n)
        {
            base.SampleRate = fS;
            _n = n;
            double w = Math.Tan(Math.PI * fC / fS);
            double a = 1 / (1 + w);
            _b0 = a;
            _b1 = -_b0;
            _a1 = a * (w - 1);
        }

        /// <summary>
        /// Apply a first-order high pass IIR
        /// </summary>
        /// <param name="fC">Cutoff frequency, Hz</param>
        /// <param name="input">Input</param>
        public IIR1HP(double fC, ISoundObj input)
        {
            uint fS = input.SampleRate;
            base.SampleRate = fS;
            base.Input = input;
            double w = Math.Tan(Math.PI * fC / fS);
            double a = 1 / (1 + w);
            _b0 = a;
            _b1 = -_b0;
            _a1 = a * (w - 1);
        }
    }

    /// <summary>
    /// First-order low pass IIR
    /// </summary>
    public class IIR1LP : IIR1
    {
        /// <summary>
        /// Construct a generator for first-order low pass IIR
        /// </summary>
        /// <param name="fS">Sample rate, Hz</param>
        /// <param name="fC">Cutoff frequency, Hz</param>
        /// <param name="n">Length of impulse (samples)</param>
        public IIR1LP(uint fS, double fC, int n)
        {
            base.SampleRate = fS;
            _n = n;
            double w = Math.Tan(Math.PI * fC / fS);
            double a = 1 / (1 + w);
            _b0 = w * a;
            _b1 = _b0;
            _a1 = a * (w - 1);
        }

        /// <summary>
        /// Apply a first-order low pass IIR
        /// </summary>
        /// <param name="fC">Cutoff frequency, Hz</param>
        /// <param name="input">Input</param>
        public IIR1LP(double fC, ISoundObj input)
        {
            uint fS = input.SampleRate; 
            base.SampleRate = fS;
            base.Input = input;
            double w = Math.Tan(Math.PI * fC / fS);
            double a = 1 / (1 + w);
            _b0 = w * a;
            _b1 = _b0;
            _a1 = a * (w - 1);
        }
    }
}
