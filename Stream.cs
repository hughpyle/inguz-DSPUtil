using System;
using System.Collections.Generic;
using System.Text;
using System.IO;

// Copyright (c) 2006, 2007 by Hugh Pyle, inguzaudio.com

namespace DSPUtil
{
    /// <summary>
    /// Slightly seekable stream (circular buffer)
    /// </summary>
    class ButteredStream : Stream
    {
        Stream _s;
        public ButteredStream(Stream s)
        {
            _s = s;
        }

        public override bool CanRead
        {
            get { return _s.CanRead; }
        }

        public override bool CanSeek
        {
            get { return _s.CanSeek; }
        }

        public override bool CanWrite
        {
            get { return _s.CanWrite; }
        }

        public override void Flush()
        {
            _s.Flush();
        }

        public override long Length
        {
            get { return _s.Length; }
        }

        public override long Position
        {
            get
            {
                return _s.Position;
            }
            set
            {
                _s.Position = value;
            }
        }

        public override int Read(byte[] buffer, int offset, int count)
        {
            return _s.Read(buffer, offset, count);
        }

        public override long Seek(long offset, SeekOrigin origin)
        {
            return _s.Seek(offset, origin);
        }

        public override void SetLength(long value)
        {
            _s.SetLength(value);
        }

        public override void Write(byte[] buffer, int offset, int count)
        {
            _s.Write(buffer, offset, count);
        }
    }
}
