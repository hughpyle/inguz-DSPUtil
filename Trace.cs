using System;
using System.Collections.Generic;
using System.Text;
using System.IO;

// Copyright (c) 2006 by Hugh Pyle, inguzaudio.com

namespace DSPUtil
{
    public class Trace
    {
        private static bool _init; // = false;
        private static readonly object _lock = new Object();

        private static bool _useConsole; // = false;
        public static bool UseConsole
        {
            get { return _useConsole; }
            set { _useConsole = value; }
        }

        private static TextWriter _useTextWriter = null;
        public static TextWriter UseTextWriter
        {
            get { return _useTextWriter; }
            set { _useTextWriter = value; }
        }

        private static string _tracefile; // = null;
        public static string FilePath
        {
            set
            {
                _tracefile = value;
                _init = true;
            }
        }

        private static string _prefix = "";
        public static string Prefix
        {
            set
            {
                _prefix = value;
            }
        }

        private static bool TraceFile()
        {
            lock (_lock)
            {
                if (!_useConsole && !_init)
                {
                    try
                    {
                        System.Configuration.AppSettingsReader rdr = new System.Configuration.AppSettingsReader();
                        _tracefile = (string)rdr.GetValue("trace", typeof(string));
                    }
                    catch (Exception)
                    {
                    }
                }
                _init = true;
            }
            return (_init && _tracefile!=null);
        }

        public static void Write(string format, params object[] args)
        {
            //          System.Diagnostics.Trace.WriteLine(String.Format(format, args));
            if (_useTextWriter != null)
            {
                _useTextWriter.WriteLine(format, args);
            }
            else if (_useConsole)
            {
                Console.Write(format, args);
            }
            else if (TraceFile())
            {
                lock (_lock)
                {
                    try
                    {
                        string msg = String.Format(format, args);
                        File.AppendAllText(_tracefile, DateTime.Now.ToString("yyyyMMddHHmmss: ") + _prefix + msg);
                    }
                    catch (Exception)
                    {
                    }
                }
            }
        }

        public static void WriteLine()
        {
            Trace.Write(System.Environment.NewLine);
        }

        public static void WriteLine(string format, params object[] args)
        {
            Trace.Write(format + System.Environment.NewLine, args);
        }

        public static void Clear()
        {
            if (TraceFile())
            {
                try
                {
                    System.IO.File.Delete(_tracefile);
                }
                catch (Exception)
                {
                }
            }
        }
    }

}
