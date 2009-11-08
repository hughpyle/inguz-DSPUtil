using System;
using System.Collections;
using System.Collections.Generic;
using System.Runtime.Remoting;
using System.Runtime.Remoting.Channels;
using System.Threading;

using System.Runtime.Serialization;

// Copyright (c) 2006 by Hugh Pyle, inguzaudio.com
// SoundObj interface structure based originally on Garbe.Sound

namespace DSPUtil
{
    // [ServiceContract]
    public interface ISoundObj : IEnumerable<ISample>
    {
        ISoundObj Input { get; set;}
        IEnumerator<ISample> Samples { get; }
        IEnumerator<ISample> GetBufferedEnumerator();
        int Run();
        void Reset();

        int Iterations { get; }
        ushort NumChannels { get; }
        uint SampleRate { get; set; }
        bool Enabled { get; set; }

        SingleChannel Channel(ushort n);
    }

    /// <summary> Basic implementation of ISoundObj interface </summary>
    // This is not abstract, but it's not very useful unless subclassed... :-)
    //
    // Note: to enable remoting, derive from MarshalByRefObject
    // but performance is almost halved (guess due mostly to inlining?)
    // - Future may be worthwhile optimizing by hand for that
    // - Current remoting implementation not likely to be worthwhile

    public class SoundObj : /* MarshalByRefObject, */ ISoundObj
    {
        // input to this object
        protected ISoundObj _input;

        protected bool _enabled = true;
        protected bool _useBuff = true;
        protected uint _sr;
        protected ushort _nc;

        public SoundObj()
        {
            // _input = null;
        }

        public SoundObj(ISoundObj input)
        {
            Input = input;
        }

//        public void Dispose()
        //{
        //}

        /// <summary> Get or Set the input object </summary>
        public virtual ISoundObj Input
        {
            get { return (_input); }
            set
            {
                _input = value;
                _nc = (_input == null) ? (ushort)0 : _input.NumChannels;
                _sr = (_input == null) ? (uint)0 : _input.SampleRate;
            }
        }

        /// <summary> Get or Set whether this processor is enabled. False means this processor simply passes its input unchanged. </summary>
        public bool Enabled
        {
            get { return (_enabled); }
            set { _enabled = value; }
        }

        // Run() with no arguments: we don't care about the "output" from this,
        public int Run()
        {
            int j = 0;
            SampleEnum e = new SampleEnum(this);
            e.SkipAll(out j);
            return j;
        }

        public int Run(int nSamples)
        {
            int j = 0;
            bool bMore = false;
            SampleEnum e = new SampleEnum(this);
            e.Skip(nSamples, out j, out bMore);
            return j;
        }

        /// <summary>
        /// Default iterator for samples.  Most implementations will override this.
        /// </summary>
        public virtual IEnumerator<ISample> Samples
        {
            get
            {
                if (_input == null)
                {
                    yield break;
                }
                foreach (ISample sample in Input)
                {
                    yield return sample;
                }
            }
        }
        
        IEnumerator IEnumerable.GetEnumerator()
        {
            return this.Samples as IEnumerator;
        }
        IEnumerator<ISample> IEnumerable<ISample>.GetEnumerator()
        {
            return this.Samples as IEnumerator<ISample>;
        }
        
        public IEnumerator<ISample> GetBufferedEnumerator()
        {
            return new SampleEnum(this);
        }

        public virtual void Reset()
        {
            if (_input == null)
            {
                return;
            }
            _input.Reset();
        }

        // Each channel is accessible as its own SoundObj
        public virtual SingleChannel Channel(ushort n)
        {
            return new SingleChannel(this, n);
        }


        /// <summary> Number of iterations expected to do the signal processing </summary>
        public virtual int Iterations
        {
            get { return (_input==null ? 0 : _input.Iterations); }
        }

        /// <summary> Gets the sample rate of the signal </summary>
        public virtual uint SampleRate
        {
            get { return ((_sr==0) ? (_input == null ? 0 : _input.SampleRate) : _sr); }
            set { _sr = value; }
        }

        /// <summary> Gets the number of channels of the signal </summary>
        public virtual ushort NumChannels
        {
            get { return ((_nc==0) ? (_input == null ? (ushort)0 : _input.NumChannels) : _nc); }
            set { _nc = value; }
        }
    }


    // SampleEnum: a buffering enumerator.
    // Used as the default enumerator for ISoundObj for a couple reasons:
    // - Buffering generally more efficient
    // - Dramatically so if you're remote
    // - This can be made serializable (unlike the default enumerator provided by 'yield')
    public class SampleEnum : /* MarshalByRefObject, */ ISampleBuffer, IEnumerator<ISample>
    {
        ISoundObj _obj;
        IEnumerator<ISample> _enum;
        ISample[] _cache;
        int _cacheCount;
        int _cachePos;
        ushort _nc;

        private ISample[] _buff;
        private int _bufflen;

        private Complex[][] _cbuff;
        private int _cbufflen;

        private bool _moreSamples = true;   // there are samples in the source which we have not yet read

        public SampleEnum(ISoundObj obj)
        {
            _obj = obj;
            _nc = _obj.NumChannels;
            Reset();
        }

        object IEnumerator.Current { get { return _cache[_cachePos]; } }
        public ISample Current
        {
            get
            {
                return _cache[_cachePos];
            }
        }
  
        public bool MoveNext()
        {
            _cachePos++;
            if (_cache!=null && _cachePos >= _cacheCount)
            {
                _cache = null;
            }
            if (_cache==null && _moreSamples)
            {
                _cache = Read(DSPUtil.BUFSIZE, out _cacheCount, out _moreSamples);
                _cachePos = 0;
            }
            return (_cachePos < _cacheCount);
        }

        public void Reset()
        {
            _obj.Reset();
            _enum = _obj.Samples;
        }

        public void Dispose()
        {
            _enum.Dispose();
            _enum = null;
            _obj = null;
        }

        public void Skip(int n, out int nn, out bool moreRemaining)
        {
            bool more = true;
            int j = 0;
            for (j = 0; j < n && more; j++)
            {
                more = _enum.MoveNext();
            }
            nn = more ? j : (j - 1);
            moreRemaining = more;
        }

        public void SkipAll(out int nn)
        {
            bool more = true;
            int j = 0;
            for (j = 0; more; j++)
            {
                more = _enum.MoveNext();
            }
            nn = j - 1;
        }

        public ISample[] Read(int n, out int nn, out bool moreRemaining)
        {
            ISample[] buff;
            if (_obj is ISampleBuffer)
            {
                buff = (_obj as ISampleBuffer).Read(n, out nn, out moreRemaining);
            }
            else
            {
                if (_buff == null || _bufflen < n)
                {
                    _buff = new ISample[n];
                    _bufflen = n;
                }
                buff = _buff;
                moreRemaining = true;
                int j;
                for (j = 0; j < n && moreRemaining; j++)
                {
                    if (_enum.MoveNext())
                    {
                        buff[j] = _enum.Current;
                    }
                    else
                    {
                        moreRemaining = false;
                    }
                }
                nn = moreRemaining ? j : (j - 1);
            }
            return buff;
        }

        /// <summary>
        /// Read, into an array of Complex.
        /// Unused elements in the array are guaranteed null.
        /// </summary>
        /// <param name="n"></param>
        /// <param name="nn"></param>
        /// <param name="moreRemaining"></param>
        /// <returns></returns>
        public Complex[][] ReadComplex(int n, out int nn, out bool moreRemaining)
        {
            if (_obj is ISampleBuffer)
            {
                return (_obj as ISampleBuffer).ReadComplex(n, out nn, out moreRemaining);
            }
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
            moreRemaining = true;
            int j;
            for (j = 0; j < n && moreRemaining; j++)
            {
                if (_enum.MoveNext())
                {
                    ISample s = _enum.Current;
                    for (ushort c = 0; c < _nc; c++)
                    {
                        _cbuff[c][j].Re = s[c];
                    }
                }
                else
                {
                    moreRemaining = false;
                }
            }
            nn = moreRemaining ? j : (j - 1);
            return _cbuff;
        }
    }

}
