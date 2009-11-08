using System;
using System.Collections.Generic;
using System.Text;

namespace DSPUtil
{

    /// <summary>
    /// SlimPlayer class provides methods to access a single player's functionality.
    /// </summary>
    public class SlimPlayer
    {
        private SlimCLI _server;
        private Dictionary<string, string> _attributes;

        public SlimPlayer(SlimCLI server, string playerID)
        {
            _server = server;
            _attributes = new Dictionary<string, string>();
            _attributes["playerid"] = playerID;
        }

        public SlimPlayer(SlimCLI server, Dictionary<string, string> attributes)
        {
            _server = server;
            _attributes = attributes;
        }

        public override string ToString()
        {
            if (_attributes.ContainsKey("name"))
            {
                return _attributes["name"];
            }
            return base.ToString();
        }

        public IEnumerable<string> Attributes
        {
            get
            {
                foreach (string key in _attributes.Keys)
                {
                    yield return key;
                }
            }
        }
        public string Attribute(string name)
        {
            return _attributes[name];
        }

        public string PlayerID
        {
            get
            {
                return _attributes["playerid"];
            }
        }

        public bool Power
        {
            get
            {
                // command is: playerid power ?
                // response is: playerid power <poweron>
                // <poweron> can be 0 or 1
                string[] response = _server.SendCommand(PlayerID, "power", "?");
                return new Value(response[2]).BoolValue;
            }
            set
            {
                _server.SendCommand(PlayerID, "power", value ? "1" : "0");
            }
        }

        public int SignalStrength
        {
            get
            {
                string[] response = _server.SendCommand(PlayerID, "signalstrength", "?");
                return new Value(response[2]).IntValue;
            }
        }

        public bool Connected
        {
            get
            {
                string[] response = _server.SendCommand(PlayerID, "connected", "?");
                return new Value(response[2]).BoolValue;
            }
        }

        public int LinesPerScreen
        {
            get
            {
                string[] response = _server.SendCommand(PlayerID, "linesperscreen", "?");
                return new Value(response[2]).IntValue;
            }
        }

        /// <summary>
        /// Volume from 0 to 100 (fractional values are possible)
        /// </summary>
        public float Volume
        {
            get
            {
                string[] response = _server.SendCommand(PlayerID, "mixer volume", "?");
                return new Value(response[2]).FloatValue;
            }
            set
            {
                _server.SendCommand(PlayerID, "mixer volume", value.ToString());
            }
        }

        public int Rate
        {
            get
            {
                string[] response = _server.SendCommand(PlayerID, "rate", "?");
                return new Value(response[2]).IntValue;
            }
            set
            {
                _server.SendCommand(PlayerID, "rate", value.ToString());
            }
        }

        public int Sleep
        {
            get
            {
                string[] response = _server.SendCommand(PlayerID, "sleep", "?");
                return new Value(response[2]).IntValue;
            }
            set
            {
                _server.SendCommand(PlayerID, "sleep", value.ToString());
            }
        }

        public string[] Lines
        {
            get
            {
                string[] response = _server.SendCommand(PlayerID, "display", "?", "?");
                string[] lines = new string[2];
                lines[0] = response[2];
                lines[1] = response[3];
                return lines;
            }
            set
            {
                // Lines may have up to three parameters.  Third is number of seconds (defaults to 1, I think)
                string[] lines = new string[3];
                lines[0] = value.Length > 0 ? value[0] : "";
                lines[1] = value.Length > 1 ? value[1] : "";
                lines[1] = value.Length > 2 ? value[2] : "";
                _server.SendCommand(PlayerID, "display", lines[0], lines[1], lines[2]);
            }
        }

        /// <summary>
        /// Show a message briefly on the player screen.
        /// </summary>
        /// <param name="message">Message to display</param>
        public void ShowBriefly(string message)
        {
            _server.SendCommand(PlayerID, "display", "", message, "");
        }

        /// <summary>
        /// Show a message on the player screen.
        /// </summary>
        /// <param name="message">Message to display</param>
        public void Show(string header, string message, int secs)
        {
            _server.SendCommand(PlayerID, "display", header, message, secs.ToString());
        }

        //

        public void Button(string button)
        {
            _server.SendCommand(PlayerID, "button", button);
        }

        // Playlist etcetera

        public string Mode
        {
            get
            {
                string[] response = _server.SendCommand(PlayerID, "mode", "?");
                return new Value(response[2]).StringValue;
            }
            set
            {
                _server.SendCommand(PlayerID, "mode", value);
            }
        }

        public bool Paused
        {
            set
            {
                _server.SendCommand(PlayerID, "pause", value ? "1" : "0");
            }
        }

        public float Time
        {
            get
            {
                string[] response = _server.SendCommand(PlayerID, "time", "?");
                return new Value(response[2]).FloatValue;
            }
            set
            {
                _server.SendCommand(PlayerID, "time", value.ToString());
            }
        }
    }

}
