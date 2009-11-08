using System;
using System.Collections.Generic;
using System.Text;

// Copyright (c) 2006, 2007 by Hugh Pyle, inguzaudio.com

namespace DSPUtil
{
    public class MorseCode : SoundObj
    {
        // Apply a morse-code envelope

        static string[] _code = { " ", ".-", "-...", "-.-.", "-..", ".", "..-.", "--.", "....", "..", ".---", "-.-", ".-..", "--", "-.", "---", ".--.", "--.-", ".-.", "...", "-", "..-", "...-", ".--", "-..-", "-.--", "--.." };

        string _text;
        double _speed;
        bool _repeat;
        List<bool> _envelope;

        public MorseCode(string text)
        {
            _speed = 5;
            _repeat = true;
            Text = text;
        }

        public MorseCode(string text, double speed, bool repeat)
        {
            _speed = speed;
            _repeat = repeat;
            Text = text;
        }

        public string Text
        {
            set
            {
                _text = value;
                _envelope = Code(_text);
            }
        }

        public static List<bool> Code(string s)
        {
            List<bool> list = new List<bool>();
            foreach (char c in s.ToUpperInvariant())
            {
                int i = c - 'A' + 1;
                if (i < 0) i = 0;
                if (i > _code.Length) i = 0;
                foreach (char d in _code[i])
                {
                    switch (d)
                    {
                        case ' ':
                            // word space 7 dits, but we'll have the three-dit letter-space too
                            list.Add(false);
                            list.Add(false);
                            list.Add(false);
                            list.Add(false);
                            break;
                        case '.':
                            list.Add(true);
                            break;
                        case '-':
                            list.Add(true);
                            list.Add(true);
                            list.Add(true);
                            break;
                        default:
                            break;
                    }
                    // symbol-space
                    list.Add(false);
                }
                // letter-space
                list.Add(false);
                list.Add(false);
            }
            return list;
        }

        /// <summary>
        /// Length of a dit, in seconds
        /// </summary>
        public double DitSeconds
        {
            get
            {
                // Speed is PARIS words (50 dits) per minute.
                // One dit = 50/speed minutes = (50/speed)/60 seconds
                return 50 / (_speed * 60);
            }
        }

        public double DitSamples
        {
            get
            {
                return SampleRate * DitSeconds;
            }
        }

        // Length of the text, in seconds
        public double LengthSeconds
        {
            get
            {
                return DitSeconds * _envelope.Count;
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

                int n = 0;
                foreach (ISample sample in Input)
                {
                    int tick = (int)(n / DitSamples);
                    if (_repeat)
                    {
                        tick = tick % _envelope.Count;
                    }

                    double gain = 0;
                    if (tick < _envelope.Count)
                    {
                        gain = _envelope[tick] ? 1 : 0;
                    }

                    Sample s = new Sample(sample, gain);
                    n++;
                    yield return s;
                }
            }
        }
    }
}
