using System;
using System.Collections.Generic;
using System.Text;
using System.Xml;
using System.Globalization;

// Copyright (c) 2007 by Hugh Pyle, inguzaudio.com


namespace DSPUtil
{
    public class UHJCoeffs
    {
        public Complex LX;
        public Complex LW;
        public Complex LY;
        public Complex RX;
        public Complex RW;
        public Complex RY;
        public UHJCoeffs()
        {
            // Proper UHJ

            // Transformation filters
            // Gerzon 1985 "Ambisonics in Multichannel Broadcasting and Video"
            //
            // Coefficients: http://www.york.ac.uk/inst/mustech/3d_audio/ambis2.htm
            // Left = (0.0928 + 0.255j)X + (0.4699 - 0.171j)W + (0.3277)Y
            // Right= (0.0928 - 0.255j)X + (0.4699 + 0.171j)W - (0.3277)Y
            // Higher precision versions on http://en.wikipedia.org/wiki/Ambisonic_UHJ_format
            // S = 0.9396926*W + 0.1855740*X
            // D = j(-0.3420201*W + 0.5098604*X) + 0.6554516*Y
            // Left = (S + D)/2.0  = ( (0.94 - 0.342j)W + (0.1856 + 0.51j)X + 0.655Y )/2
            // Right = (S - D)/2.0 = ( (0.94 + 0.342j)W + (0.1856 - 0.51j)X - 0.655Y )/2

            // The Mid-Side versions are simpler
            // L+R = (0.0928 + 0.255j)X + (0.4699 - 0.171j)W + (0.3277)Y + ((0.0928 - 0.255j)X + (0.4699 + 0.171j)W - (0.3277)Y)
            //     = (0.1856)X          + (0.9398)W
            // L-R = (0.0928 + 0.255j)X + (0.4699 - 0.171j)W + (0.3277)Y - ((0.0928 - 0.255j)X + (0.4699 + 0.171j)W - (0.3277)Y)
            //     =          (0.510j)X +         (-0.342j)W + (0.6554)Y
            // but since we're delaying signal via convolution anyway, not *too* much extra processing to do in LR mode...

            LX = new Complex( 0.0927870,  0.25493020);
            LW = new Complex( 0.4698463, -0.17101005);
            LY = new Complex( 0.3277258,  0.00000000);
            RX = new Complex( 0.0927870, -0.25493020);
            RW = new Complex( 0.4698463,  0.17101005);
            RY = new Complex(-0.3277258,  0.00000000);
        }
        public UHJCoeffs(double j)
        {
            // Blumlein with a j*W mix (mucks up the low-mid timbre a little but adds really nice spaciousness)

            // For default mic angle of 90 degrees, mulX=mulY=sqrt(2)/2
            double mulX = Math.Sqrt(2); // Math.Cos(MathUtil.Radians(_ambiMicAngle / 2));
            double mulY = Math.Sqrt(2); // Math.Sin(MathUtil.Radians(_ambiMicAngle / 2));

            // Left = (mulX)X + ( - 0.171j)W + (mulY)Y
            // Right= (mulX)X + ( + 0.171j)W - (mulY)Y
            LW = new Complex(0, -j);
            LX = new Complex(mulX, 0);
            LY = new Complex(mulY, 0);
            RW = new Complex(0, j);
            RX = new Complex(mulX, 0);
            RY = new Complex(-mulY, 0);
        }
        public UHJCoeffs(double angleDegrees, double w, double j)
        {
            // Blumlein with a j*W mix

            // For default mic angle of 90 degrees, mulX=mulY=sqrt(2)/2
            double mulX = Math.Cos(MathUtil.Radians(angleDegrees / 2));
            double mulY = Math.Sin(MathUtil.Radians(angleDegrees / 2));

            // Left = (mulX)X + ( - 0.171j)W + (mulY)Y
            // Right= (mulX)X + ( + 0.171j)W - (mulY)Y
            LW = new Complex(w * MathUtil.INVSQRT2, -j * MathUtil.INVSQRT2);
            LX = new Complex(mulX, 0);
            LY = new Complex(mulY, 0);
            RW = new Complex(w * MathUtil.INVSQRT2, j * MathUtil.INVSQRT2);
            RX = new Complex(mulX, 0);
            RY = new Complex(-mulY, 0);
        }
    }

    
    /*
     * this intended to become a general multi-band amb decoder
     * design is incomplete
     * 
      class Transform
      {
          public Transform()
          {
          }
      }

      struct Band
      {
          double from;
          double to;
      }

      public class AmbisonicDecoder : SoundObj
      {
          string _name;
          bool _useShelf;
          bool _useDistance;
          int _inputs;
          int _outputs;
          bool[] _usesInput;
          List<Transform[]> _frequencyBands;

          public AmbisonicDecoder(uint sampleRate)
          {
              _sr = sampleRate;
              _usesInput = new bool[4];   // W, X, Y, Z (and so on)
              _frequencyBands = new List<Transform[]>();
          }
            
          public void Load(XmlDocument config)
          {
              XmlNode node;

              node = config.SelectSingleNode("//Matrix/@Name");
              _name = (node == null) ? null : node.Value;

              node = config.SelectSingleNode("//Matrix/@Outputs");
              _outputs = (node == null) ? 0 : int.Parse(node.Value, CultureInfo.InvariantCulture);

              node = config.SelectSingleNode("//Matrix/@Shelf");
              _useShelf = (node == null) ? false : bool.Parse(node.Value, CultureInfo.InvariantCulture);

              node = config.SelectSingleNode("//Matrix/@Distance");
              _useDistance = (node == null) ? false : bool.Parse(node.Value, CultureInfo.InvariantCulture);

              if (_useShelf)
              {
              }
              else
              {
                  // Start with a "null-band" transform
                  // (no band, meaning full-frequency-range or shelved)
              }

              // Count the inputs
              _inputs = 0;
              XmlNodeList inputs = config.SelectNodes("//In");
              foreach (XmlNode input in inputs)
              {
                  // Check each iput channel is actually sensible
                  int c = int.Parse(input.Attributes["Channel"].Value, CultureInfo.InvariantCulture);
                  if (c >= _inputs)
                  {
                      _inputs = c+1;
                  }
              }
              _usesInput = new bool[_inputs];

              // Count the outputs
              XmlNodeList outputs = config.SelectNodes("//Matrix/Out");
              foreach (XmlNode output in outputs)
              {
                  // Check each output channel is actually sensible
                  int c = int.Parse(output.Attributes["Channel"].Value, CultureInfo.InvariantCulture);
                  if (c >= _outputs)
                  {
                      throw new NotSupportedException("Output channel " + c + " must be in the range 0 to " + (_outputs-1));
                  }
              }

  //            _ambiChannels = new Complex[_ambiOutputs][];
  //            for (int cOut = 0; cOut < _ambiOutputs; cOut++)
  //            {
  //                _ambiChannels[cOut] = new Complex[4];
  //                for (int cIn = 0; cIn < 4; cIn++)
  //                {
  //                    double re = nodeValueDouble(ambiMatrixDoc, String.Format("//Matrix/Out[@Channel={0}]/In[@Channel={1}]/Re", cOut, cIn));
  //                    double im = nodeValueDouble(ambiMatrixDoc, String.Format("//Matrix/Out[@Channel={0}]/In[@Channel={1}]/Im", cOut, cIn));
  //                    _ambiChannels[cOut][cIn] = new Complex(re, im);
  //                }
  //            }
          }
      }
     */
}
