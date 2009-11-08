using System;
using System.Collections.Generic;
using System.Text;

// Copyright (c) 2006 by Hugh Pyle, inguzaudio.com

namespace DSPUtil
{
    [Serializable]
    public class Reverser : SoundObj
    {
        private List<ISample> _data;

        public Reverser()
        {
            _data = new List<ISample>();
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
                // Read the whole input into our buffer
                foreach (ISample sample in Input)
                {
                    _data.Add(sample);
                }
                // Write them back out in reverse order
                _data.Reverse();
                foreach (ISample sample in _data)
                {
                    yield return sample;
                }
            }
        }
    }
}
