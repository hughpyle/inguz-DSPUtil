using System;
using System.Collections.Generic;
using System.Text;
using System.IO;

namespace DSPUtil
{
    public class WaveChunkFactory
    {
        public static WaveChunk ReadChunk(BinaryReader rdr, bool bigEndian)
        {
            WaveChunk chunk = null;
            string id = new string(rdr.ReadChars(4));
            if (id.Length > 0)
            {
                int chunkSize = rdr.ReadInt32();
                if (bigEndian)
                    chunkSize = System.Net.IPAddress.NetworkToHostOrder(chunkSize);

                switch (id)
                {
                    case "fmt ":
                        chunk = new WaveFmtChunk(chunkSize);
                        break;
                    default:
                        chunk = new WaveChunk(id, chunkSize);
                        break;
                }
            }
            return chunk;
        }
    }

    /// <summary>
    /// An undifferentiated Wave chunk
    /// </summary>
    public class WaveChunk
    {
        protected string _id = "";
        protected int _chunkSize = 0;
        protected int _structSize = 0;
        protected int _dataSize;

        public WaveChunk(string id, int chunkSize)
        {
            _id = id;
            _chunkSize = chunkSize;
        }
        public string ID { get { return _id; } }
        public virtual int DataSize { get { return _dataSize; } set { _dataSize = value; } }
        public int TotalSize { get { return _dataSize + _structSize; } }
        public virtual void Skip(BinaryReader rdr)
        {
            rdr.ReadBytes(_chunkSize);
        }
        public virtual void Write(BinaryWriter w)
        {
        }
    }

    class WaveFmtChunk : WaveChunk
    {
        public WaveFmtChunk(int chunkSize) : base("fmt ", chunkSize)
        {
        }

        public override int DataSize
        {
            set
            {
                throw new InvalidOperationException();
            }
        }
    }


}
