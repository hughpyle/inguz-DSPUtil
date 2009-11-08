using System;
using System.Collections.Generic;
using System.Text;

// Copyright (c) 2009 by Hugh Pyle, inguzaudio.com

namespace DSPUtil
{
    abstract public class Matrix
    {
        // how many channels
        protected ushort _dimension;

        // square matrix [rows][columns]
        protected double[,] _matrix;

        // parameter (usually degrees, sometimes something else)
        protected double _param = double.NaN;
        virtual public double Param
        {
            get { return _param; }
            set { _param = value; }
        }

        public ISample Process(ISample sample)
        {
            if (sample.NumChannels != _dimension)
            {
                throw new InvalidOperationException();
            }
            Sample s = new Sample(_dimension);
            // row of matrix
            // channel of output
            for (ushort j = 0; j < _dimension; j++)
            {
                // channel of original sample
                // column of matrix
                for (ushort c = 0; c < _dimension; c++)
                {
                    s[j] += sample[c] * _matrix[j,c];
                }
            }
            return s;
        }
    }

    abstract public class Amb1Matrix : Matrix
    {
        public Amb1Matrix()
        {
            _dimension = 4;
            _matrix = new double[4,4];
        }
    }


    // Classic rotation matrices for first-order
    // http://www.muse.demon.co.uk/fmhrotat.html

    public class Amb1RotAboutX : Amb1Matrix
    {
        public Amb1RotAboutX(double degrees)
            : base()
        {
            Param = degrees;
        }
        override public double Param
        {
            set
            {
                if (_param != value)
                {
                    _param = value;
                    double r = MathUtil.Radians(value);
                    _matrix[0, 0] = 1;
                    _matrix[1, 1] = 1;
                    _matrix[2, 2] = Math.Cos(r);
                    _matrix[3, 2] = Math.Sin(r);
                    _matrix[2, 3] = -Math.Sin(r);
                    _matrix[3, 3] = Math.Cos(r);
                }
            }
        }
    }

    public class Amb1RotAboutY : Amb1Matrix
    {
        public Amb1RotAboutY(double degrees)
            : base()
        {
            Param = degrees;
        }
        override public double Param
        {
            set
            {
                if (_param != value)
                {
                    _param = value;
                    double r = MathUtil.Radians(value);
                    _matrix[0, 0] = 1;
                    _matrix[1, 1] = Math.Cos(r);
                    _matrix[3, 1] = Math.Sin(r);
                    _matrix[2, 2] = 1;
                    _matrix[1, 3] = -Math.Sin(r);
                    _matrix[3, 3] = Math.Cos(r);
                }
            }
        }
    }

    public class Amb1RotAboutZ : Amb1Matrix
    {
        public Amb1RotAboutZ(double degrees)
            : base()
        {
            Param = degrees;
        }
        override public double Param
        {
            set
            {
                if (_param != value)
                {
                    _param = value;
                    double r = MathUtil.Radians(value);
                    _matrix[0, 0] = 1;
                    _matrix[1, 1] = Math.Cos(r);
                    _matrix[2, 1] = Math.Sin(r);
                    _matrix[1, 2] = -Math.Sin(r);
                    _matrix[2, 2] = Math.Cos(r);
                    _matrix[3, 3] = 1;
                }
            }
        }
    }

    // Classic dominance and
    // Joseph Anderson's push etc
    // http://ambisonics.iem.at/symposium2009/proceedings/ambisym09-josephanderson-ambitk-poster.pdf/at_download/file

    public class Amb1Dominance : Amb1Matrix
    {
        private double _az, _el;
        public Amb1Dominance(double lambda, double az, double el)
            : base()
        {
            Param = lambda;
            _az = az;
            _el = el;
        }
        override public double Param
        {
            set
            {
                if (_param != value)
                {
                    _param = value;
                    double lambda = value;
                    _matrix[0, 0] = 0.5 * (lambda + (1 / lambda));
                    _matrix[1, 0] = (1 / Math.Sqrt(2)) * (lambda - (1 / lambda));
                    _matrix[0, 1] = (1 / Math.Sqrt(8)) * (lambda - (1 / lambda));
                    _matrix[1, 1] = 0.5 * (lambda + (1 / lambda));
                    _matrix[2, 2] = 1;
                    _matrix[3, 3] = 1;
                }
            }
        }
    }

    public class Amb1Focus : Amb1Matrix
    {
        private double _az, _el;
        public Amb1Focus(double degrees, double az, double el)
            : base()
        {
            Param = degrees;
            _az = az;
            _el = el;
        }
        override public double Param
        {
            set
            {
                if (_param != value)
                {
                    _param = value;
                    double r = MathUtil.Radians(value);
                    double sinr = Math.Sin(r);
                    double sqr2 = Math.Sqrt(2);
                    _matrix[0, 0] = (1 / (1 + Math.Abs(sinr)));
                    _matrix[1, 0] = sqr2 * (sinr / (1 + Math.Abs(sinr)));
                    _matrix[0, 1] = (1 / sqr2) * (sinr / (1 + Math.Abs(sinr)));
                    _matrix[1, 1] = (1 / (1 + Math.Abs(sinr)));
                    _matrix[2, 2] = Math.Sqrt((1 - Math.Abs(sinr)) / (1 + Math.Abs(sinr)));
                    _matrix[3, 3] = Math.Sqrt((1 - Math.Abs(sinr)) / (1 + Math.Abs(sinr)));
                }
            }
        }
    }

    public class Amb1Push : Amb1Matrix
    {
        private double _az, _el;
        public Amb1Push(double degrees, double az, double el)
            : base()
        {
            Param = degrees;
            _az = az;
            _el = el;
        }
        override public double Param
        {
            set
            {
                if (_param != value)
                {
                    _param = value;
                    double r = MathUtil.Radians(value);
                    _matrix[0, 0] = 1;
                    _matrix[1, 0] = Math.Sqrt(2) * Math.Abs(Math.Sin(r)) * Math.Sin(r);
                    _matrix[1, 1] = Math.Cos(r) * Math.Cos(r);
                    _matrix[2, 2] = Math.Cos(r) * Math.Cos(r);
                    _matrix[3, 3] = Math.Cos(r) * Math.Cos(r);
                }
            }
        }
    }

    public class Amb1Press : Amb1Matrix
    {
        private double _az, _el;
        public Amb1Press(double degrees, double az, double el)
            : base()
        {
            Param = degrees;
            _az = az;
            _el = el;
        }
        override public double Param
        {
            set
            {
                if (_param != value)
                {
                    _param = value;
                    double r = MathUtil.Radians(value);
                    _matrix[0, 0] = 1;
                    _matrix[1, 0] = Math.Sqrt(2) * Math.Abs(Math.Sin(r)) * Math.Sin(r);
                    _matrix[1, 1] = Math.Cos(r) * Math.Cos(r);
                    _matrix[2, 2] = Math.Cos(r);
                    _matrix[3, 3] = Math.Cos(r);
                }
            }
        }
    }

}
