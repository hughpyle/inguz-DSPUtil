using System;
using System.Collections.Generic;
using System.Text;

// Copyright (c) 2006, 2007 by Hugh Pyle, inguzaudio.com

namespace DSPUtil
{
    /// <summary>
    /// Implementation of a double-precision (64-bit float) complex number
    /// </summary>
    //[Serializable]
    public struct Complex
    {
        // Real and Imaginary parts of a Complex number
        public double Re;
        public double Im;

        public Complex(Complex c)
        {
            this.Re = c.Re;
            this.Im = c.Im;
        }

        public Complex(double real, double imaginary)
        {
            this.Re = real;
            this.Im = imaginary;
        }

        public Complex(double phase)
        {
            this.Re = Math.Cos(phase);
            this.Im = Math.Sin(phase);
        }

        // Accessor methods for accessing/setting private variables
        /// <summary>
        /// Real component
        /// </summary>
        public double Real
        {
            get { return Re; }
            set { Re = value; }
        }

        /// <summary>
        /// Imaginary component
        /// </summary>
        public double Imaginary
        {
            get { return Im; }
            set { Im = value; }
        }

        /// <summary>
        /// Magnitude (always positive)
        /// </summary>
        public double Magnitude
        {
            get
            {
                double abs = Math.Sqrt(Im * Im + Re * Re);
//                if (Re > 0) abs = -abs;
                return abs;
            }
        }

        /// <summary>
        /// Phase, radians
        /// </summary>
        public double Phase
        {
            get
            {
                return Math.Atan2(Im, Re);
            }
        }

        /// <summary>
        /// Set value
        /// </summary>
        /// <param name="re">Real</param>
        /// <param name="im">Imaginary</param>
        public void Set(double re, double im)
        {
            Re = re;
            Im = im;
        }

        /// <summary>
        /// In-place multiply
        /// </summary>
        /// <param name="c">Complex</param>
        public void mul(Complex c)
        {
            double gem = c.Re * Re - c.Im * Im;
            Im = c.Im * Re + c.Re * Im;
            Re = gem;
        }

        /// <summary>
        /// In-place multiply
        /// </summary>
        /// <param name="d">double</param>
        public void mul(double d)
        {
            Im *= d;
            Re *= d;
        }

        /// <summary>
        /// In-place divide: this=this/arg
        /// </summary>
        /// <param name="c">Complex</param>
        public void div(Complex c)
        {
            double d = c.Re * c.Re + c.Im * c.Im;
            //if (d == 0) throw new DivideByZeroException();
            double r = (Re * c.Re + Im * c.Im) / d;
            double i = (Im * c.Re - Re * c.Im) / d;
            Re = r;
            Im = i;
        }

        /// <summary>
        /// In-place inverted divide: this = arg/this
        /// </summary>
        /// <param name="c">Complex</param>
        public void idiv(Complex c)
        {
            double d = Re * Re + Im * Im;
            //if (d == 0) throw new DivideByZeroException();
            double r = (c.Re * Re + c.Im * Im) / d;
            double i = (c.Im * Re - c.Re * Im) / d;
            Re = r;
            Im = i;
        }

        /// <summary>
        /// Raises to a power
        /// </summary>
        /// <param name="c">Complex, power</param>
        public void Pow(Complex y)
        {
            // http://home.att.net/~srschmitt/complexnumbers.html
            double r = Magnitude;
            double t = Phase;
            double c = y.Re;
            double d = y.Im;
            Re = Math.Pow(r, c) * Math.Exp(-d * t) * Math.Cos(c * t + d * Math.Log(r));
            Im = Math.Pow(r, c) * Math.Exp(-d * t) * Math.Sin(c * t + d * Math.Log(r));
        }

        /// <summary>
        /// Implicit conversion of Complex to double
        /// </summary>
        /// <param name="c"></param>
        /// <returns></returns>
        public static implicit operator double(Complex c)
        {
            return c.Re;
        }

        /// <summary>
        /// Explicit cast of double to Complex
        /// </summary>
        /// <param name="f"></param>
        /// <returns></returns>
        public static explicit operator Complex(double f)
        {
            return new Complex(f, 0);
        }

        /// <summary>
        /// operator add
        /// </summary>
        /// <param name="c"></param>
        /// <returns></returns>
        public static Complex operator +(Complex c)
        {
            return c;
        }

        /// <summary>
        /// operator subtract
        /// </summary>
        /// <param name="c"></param>
        /// <returns></returns>
        public static Complex operator -(Complex c)
        {
            return new Complex(-c.Re, -c.Im);
        }

        /// <summary>
        /// Complex Conjugate operator
        /// </summary>
        /// <param name="c"></param>
        /// <returns></returns>
        public static Complex operator ~(Complex c)
        {
            return new Complex(c.Re, -c.Im);
        }

        /// <summary>
        /// operator add
        /// </summary>
        /// <param name="c1"></param>
        /// <param name="c2"></param>
        /// <returns></returns>
        public static Complex operator +(Complex c1, Complex c2)
        {
            return new Complex(c1.Re + c2.Re, c1.Im + c2.Im);
        }

        public static Complex operator +(Complex c1, double num)
        {
            return new Complex(c1.Re + num, c1.Im);
        }

        public static Complex operator +(double num, Complex c1)
        {
            return new Complex(c1.Re + num, c1.Im);
        }

        /// <summary>
        /// Operator subtract
        /// </summary>
        /// <param name="c1"></param>
        /// <param name="num"></param>
        /// <returns></returns>
        public static Complex operator -(Complex c1, double num)
        {
            return new Complex(c1.Re - num, c1.Im);
        }

        public static Complex operator -(double num, Complex c1)
        {
            return new Complex(num - c1.Re, -c1.Im);
        }

        public static Complex operator -(Complex c1, Complex c2)
        {
            return new Complex(c1.Re - c2.Re, c1.Im - c2.Im);
        }

        /// <summary>
        /// Operator multiply
        /// </summary>
        /// <param name="c1"></param>
        /// <param name="c2"></param>
        /// <returns></returns>
        public static Complex operator *(Complex c1, Complex c2)
        {
            return new Complex((c1.Re * c2.Re) - (c1.Im * c2.Im),  (c1.Re * c2.Im) + (c1.Im * c2.Re));
        }

        public static Complex operator *(Complex c1, double num)
        {
            return new Complex(c1.Re * num, c1.Im * num);
        }

        public static Complex operator *(double num, Complex c1)
        {
            return new Complex(c1.Re * num, c1.Im * num);
        }

        /// <summary>
        /// Operator divide
        /// </summary>
        /// <param name="c1"></param>
        /// <param name="c2"></param>
        /// <returns></returns>
        public static Complex operator /(Complex c1, Complex c2)
        {
            double div = c2.Re * c2.Re + c2.Im * c2.Im;
            //if (div == 0) throw new DivideByZeroException();

            return new Complex((c1.Re * c2.Re + c1.Im * c2.Im) / div,
                                (c1.Im * c2.Re - c1.Re * c2.Im) / div);
        }

        public static Complex operator /(double num, Complex c1)
        {
            return new Complex(num / c1.Re, num / c1.Im);
        }

        public static Complex operator /(Complex c1, double num)
        {
            return new Complex(c1.Re / num, c1.Im / num);
        }

        /// <summary>
        /// Operator equalitytest
        /// </summary>
        /// <param name="c1"></param>
        /// <param name="c2"></param>
        /// <returns></returns>
        public static bool operator ==(Complex c1, Complex c2)
        {
            return (c1.Re == c2.Re) && (c1.Im == c2.Im);
        }

        /// <summary>
        /// Operator negative equalitytest
        /// </summary>
        /// <param name="c1"></param>
        /// <param name="c2"></param>
        /// <returns></returns>
        public static bool operator !=(Complex c1, Complex c2)
        {
            return (c1.Re != c2.Re) || (c1.Im != c2.Im);
        }

        public override int GetHashCode()
        {
            return (Re.GetHashCode() ^ Im.GetHashCode());
        }

        public override bool Equals(object obj)
        {
            return (obj is Complex) ? (this == (Complex)obj) : false;
        }

        public override string ToString()
        {
            return (String.Format("{0} + {1}i", Real, Imaginary));
        }
    }
}
