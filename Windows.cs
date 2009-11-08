using System;
using System.Collections.Generic;
using System.Text;

// Copyright (c) 2006, 2007 by Hugh Pyle, inguzaudio.com

namespace DSPUtil
{
    public class WindowedBuffer : SoundObj
    {
        SoundBuffer _buff;

        public WindowedBuffer(ISoundObj input, CosWindow window, int start, int count)
        {
            _buff = new SoundBuffer(new SampleBuffer(input).Subset(start, count));
            _buff.ApplyWindow(window);
        }

        /// <summary>
        /// Get an iterator for samples
        /// </summary>
        public override IEnumerator<ISample> Samples
        {
            get
            {
                return _buff.Samples;
            }
        }
    }

    public class OverlappingWindows
    {
        int _length;
        ISoundObj _input;
        double _overlap;
        CosWindow _window;

        public OverlappingWindows(ISoundObj input, int length)
        {
            _input = input;
            _length = length;
            _overlap = 0.5;
            _window = new Hamming(length / 2, length / 2);
        }

        public OverlappingWindows(ISoundObj input, int length, double overlap, CosWindow window)
        {
            _input = input;
            _length = length;
            _overlap = overlap;
            _window = window;
        }

        /// <summary>
        /// Return an endless enumeration of ISoundObj
        /// representing overlapped windowed slices of the input.
        /// You can tell when the input is empty when your slices are empty.
        /// </summary>
        public IEnumerator<ISoundObj> Items
        {
            get
            {
                int start = 0;
                while (true)
                {
                    yield return new WindowedBuffer(_input, _window, start, _length);
                    start += _length;
                }
            }
        }
    }
}
