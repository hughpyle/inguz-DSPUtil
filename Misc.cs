using System;
using System.Collections.Generic;
using System.Text;
using System.Globalization;

// Support classes for SlimCLI

namespace DSPUtil
{
    class Value
    {
        private string[] _value;
        public Value(string s)
        {
            char[] c = { ':' };
            _value = s.Split(c, 2);
        }
        private string Val
        {
            get
            {
                if (_value == null)
                {
                    return "";
                }
                if (_value.Length > 1)
                {
                    return _value[1];
                }
                else if (_value.Length == 1)
                {
                    return _value[0];
                }
                return "";

            }
        }
        public string Tag
        {
            get
            {
                if (_value.Length > 1)
                {
                    return _value[0];
                }
                return "";
            }
        }
        public string StringValue
        {
            get
            {
                return Val;
            }
        }
        public bool BoolValue
        {
            get
            {
                switch (Val)
                {
                    case "0":
                        return false;
                    case "1":
                        return true;
                    default:
                        throw new Exception("Invalid response, expected 0 or 1");
                }
            }
        }
        public int IntValue
        {
            get
            {
                return int.Parse(Val, CultureInfo.InvariantCulture);
            }
        }
        public float FloatValue
        {
            get
            {
                return float.Parse(Val, CultureInfo.InvariantCulture);
            }
        }
    }
}
