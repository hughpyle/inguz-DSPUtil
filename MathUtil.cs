using System;
using System.Collections.Generic;
using System.Text;

// Copyright (c) 2006, 2009 by Hugh Pyle, inguzaudio.com

namespace DSPUtil
{
    // Package constants
    public class DSPUtil
    {
        public const int BUFSIZE = 8192; // 512;
        public static Version VERSION { get { return new Version("0.9.32"); } }
        public static DateTime EXPIRY { get { return DateTime.MaxValue; } }  // expire just about never
        public static String GetVersionInfo()
        {
            string s = String.Format("version {0}{1}", VERSION, IsMono() ? " on Mono" : "");
            return s;
        }
        public static bool IsMono()
        {
            Type t = Type.GetType("Mono.Runtime");
            return (t != null) ? true : false;
        }
    }


    // flags (bit)
    public enum ChannelFlag
    {
        NONE  = 0x00,
        LEFT  = 0x01,
        RIGHT = 0x02,
        BOTH  = 0x03
    }

    // Math utils
    public sealed class MathUtil
    {
        public static double SQRT2 = Math.Sqrt(2);
        public static double INVSQRT2 = 1/Math.Sqrt(2);

        public static bool IsPowerOfTwo(int n)
        {
            // http://graphics.stanford.edu/~seander/bithacks.html#DetermineIfPowerOf2
            return ((n & (n - 1)) == 0) && n>0;
        }
        public static int NextPowerOfTwo(int n)
        {
            // http://graphics.stanford.edu/~seander/bithacks.html#RoundUpPowerOf2
            n--;
            n |= n >> 1;
            n |= n >> 2;
            n |= n >> 4;
            n |= n >> 8;
            n |= n >> 16;
            n++;
            return n;
        }

        public static double dB(double gain)
        {
            return (20 * Math.Log(gain,10));
        }
        public static double gain(double dB)
        {
            return (Math.Pow(10, dB / 20));
        }


        public static double Feet(int samples, uint sampleRate)
        {
            return 1116.43701 * samples / sampleRate;
        }
        public static int FSamples(double feet, uint sampleRate)
        {
            return (int)((feet / 1116.43701) * sampleRate);
        }

        public static double Metres(int samples, uint sampleRate)
        {
            return 340.29 * samples / sampleRate;
        }
        public static int MSamples(double metres, uint sampleRate)
        {
            return (int)((metres / 340.29) * sampleRate);
        }

        public static double FcFromFeet(double r)
        {
            return 1116.43701 / (2 * Math.PI * r);
        }

        public static double FcFromMetres(double r)
        {
            return 340.29 / (2 * Math.PI * r);
        }

        /// <summary>
        /// Compute the Bark critical band rate z at frequency f (Hz)
        /// </summary>
        /// <param name="f">Frequency (Hz)</param>
        /// <returns>z</returns>
        public static double Bark(double f)
        {
            // from http://www.ling.su.se/STAFF/hartmut/bark.htm
            double z = (26.81 / (1 + 1960 / f)) - 0.53;
            return z;
        }

        /// <summary>
        /// Compute the Bark critical bandwidth Cb at frequency f (Hz)
        /// </summary>
        /// <param name="f">Frequency (Hz)</param>
        /// <returns>Cb</returns>
        public static double BarkCb(double f)
        {
            // from http://www.ling.su.se/STAFF/hartmut/bark.htm
            double z = Bark(f);
            double Cb = 52548 / ((z * z) - (52.56 * z) + 690.39);
            return Cb;
        }

        public delegate double InvertDelegate(double f);
        public static double invert(InvertDelegate fn, double low, double high, double i)
        {
            double probe;
            while (high - low > 0.001)
            {
                probe = (low + high) / 2;
                if (fn(probe) > i)
                    high = probe;
                else
                    low = probe;
            }
            return low;
        }


        // Greatest common divisor
        public static uint gcd(uint a, uint b)
        {
            while (b != 0)
            {
                uint t = b;
                b = a % b;
                a = t;
            }
            return a;
        }

        // Least common multiple
        public static uint lcm(uint a, uint b)
        {
            long g = gcd(a, b);
            uint l = (uint)(a / g) * b;
            return l;
        }

        public static double Radians(double degrees)
        {
            return degrees * Math.PI / 180;
        }

        public static double Degrees(double radians)
        {
            return radians * 180 / Math.PI;
        }
    }
}
