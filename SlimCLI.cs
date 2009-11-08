using System;
using System.Collections.Generic;
using System.Text;
using System.Net;
using System.Net.Sockets;
using System.Web;

// SlimCLI package wraps the SlimServer command-line interface RPC
// Basic connectivity lifted from http://www.phred.org/~alex/SlimserverConnector.cs

namespace DSPUtil
{
    public class SlimCLI
    {
        private const ushort DEFAULT_PORT = 9090;
        private Socket _socket;
        private byte[] _receiveBuffer;

        public SlimCLI()
        {
        }

        ~SlimCLI()
        {
            Close();
        }

        #region Internal methods

        internal bool IsOpen()
        {
            return (_socket != null);
        }

        /// <summary>
        /// Send a single command to the slimserver and wait for a response.  
        /// The command shouldn't contain any newlines.  They are added 
        /// as defined by the protocol.  The response is put into a 
        /// tokenized array.
        /// </summary>
        /// <param name="player">the ID of the player, or "player"</param>
        /// <param name="parameters">up to 4 command parameters</param>
        /// <returns>The response from the server</returns>
        internal string[] SendCommand(string player, params string[] parameters)
        {
            lock (this)
            {
                if (!IsOpen())
                {
                    throw new InvalidOperationException("SlimCLI must be opened before use.");
                }
                // build our command string.  This consists of the player id
                // and each parameter which has been URL encoded.  Spaces
                // seperate all parameters.  The command terminates with CR.
                StringBuilder sb = new StringBuilder();
                sb.Append(player);
                sb.Append(" ");
                for (int i = 0; i < parameters.Length; i++)
                {
                    string s = parameters[i];
                    s = HttpUtility.UrlEncode(s);
                    s = s.Replace("+", "%20");
                    sb.Append(s);
                    // if this is the last parameter then add a space, 
                    // otherwise add a CR
                    sb.Append((i == (parameters.Length - 1)) ? "\r\n" : " ");
                }
                // send our command
                byte[] commandBytes = Encoding.ASCII.GetBytes(sb.ToString());
                _socket.Send(commandBytes);

                // start waiting for a response
                string response = null;		// set to a response when CR is found
                int iReceiveBuffer = 0;		// offset into m_receiveBuffer
                while (response == null)
                {
                    // read some more bytes from the socket
                    int cNewBytes = _socket.Receive(
                        _receiveBuffer,
                        iReceiveBuffer,
                        _receiveBuffer.Length - iReceiveBuffer,
                        SocketFlags.None);

                    // slimserver never sends more than one line to us at
                    // a time.  We should find the CR, CRLF, or LF at the
                    // end of the last data received.  If we do not then
                    // we need to read some more data.
                    int chompBytes = 0;
                    // look at the last two bytes for a newline
                    byte[] lastBytes = new byte[2] 
				{
					_receiveBuffer[iReceiveBuffer + cNewBytes - 2],
					_receiveBuffer[iReceiveBuffer + cNewBytes - 1]
				};

                    if (lastBytes[0] == '\r' && lastBytes[1] == '\n')
                    {
                        chompBytes = 2;
                    }
                    else if (lastBytes[1] == '\r' || lastBytes[1] == '\n')
                    {
                        chompBytes = 1;
                    }

                    if (chompBytes == 0)
                    {
                        // we need more data
                        iReceiveBuffer += cNewBytes;
                        if (iReceiveBuffer >= _receiveBuffer.Length - 1)
                        {
                            throw new Exception("Invalid protocol");
                        }
                    }
                    else
                    {
                        // the receive buffer contains our reponse.  
                        // chompBytes tells us how many bytes are the tail for
                        // the newline.  We encode everything but the newline
                        response = Encoding.ASCII.GetString(
                            _receiveBuffer,
                            0,
                            iReceiveBuffer + cNewBytes - chompBytes);
                    }
                }

                // now we have the response as a string.  break it into
                // individual words
                string[] responseParams = response.Split(new char[] { ' ' });

                // URL decode each of the strings
                for (int i = 0; i < responseParams.Length; i++)
                {
                    responseParams[i] = HttpUtility.UrlDecode(responseParams[i]);
                }

                // return the response array
                return responseParams;
            }
        }

        internal IEnumerable<string[]> ChunkedCommand(string command)
        {
            yield break;
        }
        #endregion

        /// <summary>
        /// Open a connection to the local SlimServer
        /// </summary>
        public void Open()
        {
            Open("127.0.0.1", DEFAULT_PORT, null, null);
        }

        /// <summary>
        /// Open a connection to the specified SlimServer
        /// </summary>
        /// <param name="hostname"></param>
        public void Open(string hostname)
        {
            Open(hostname, DEFAULT_PORT, null, null);
        }

        /// <summary>
        /// Open a connection to the specified SlimServer
        /// </summary>
        /// <param name="hostname"></param>
        /// <param name="port"></param>
        public void Open(string hostname, ushort port)
        {
            Open(hostname, port, null, null);
        }
        public void Open(string hostname, ushort port, string username, string password)
        {
            IPHostEntry hostInfo = Dns.GetHostEntry(hostname);
            IPAddress[] ipAddresses = hostInfo.AddressList;

            // use the default port if they didn't specify one
            if (port == 0) port = DEFAULT_PORT;

            // try to connect to each host in turn
            for (int i = 0; i < ipAddresses.Length && _socket == null; i++)
            {
                IPEndPoint endPoint =
                    new IPEndPoint(ipAddresses[i], port);
                _socket = new Socket(
                    AddressFamily.InterNetwork,
                    SocketType.Stream,
                    ProtocolType.Tcp);
                _socket.Connect(endPoint);
                if (!_socket.Connected)
                {
                    _socket = null;
                }
            }
            if (_socket == null)
            {
                throw new Exception("Unable to open connection to server.");
            }
            _receiveBuffer = new byte[1024 * 4];

            if (username != null)
            {
                SendCommand("login", username, password);
                // Server will disconnect if login is invalid
            }
        }

        /// <summary>
        /// Close this connection
        /// </summary>
        public void Close()
        {
            if (_socket != null)
            {
                _socket.Close();
            }
            _socket = null;
            _receiveBuffer = null;
        }


        // Player access

        /// <summary>
        /// Get the number of players connected to the server
        /// </summary>
        /// <returns>Number of players</returns>
        public int PlayerCount()
        {
            string[] response = SendCommand("player", "count", "?");
            return new Value(response[2]).IntValue;
        }

        /// <summary>
        /// Get an enumerable list of all players
        /// </summary>
        public IEnumerable<SlimPlayer> Players
        {
            get
            {
                // Query for up to 10 players at a time
                List<Dictionary<string, string>> responses = new List<Dictionary<string, string>>();
                int start = 0;
                int chunk = 10;
                int total = 1;
                while (total > start)
                {
                    string[] resp = SendCommand("players", start.ToString(), chunk.ToString());
                    Value count = new Value(resp[3]);
                    if (count.Tag != "count")
                    {
                        throw new Exception("no count");
                    }
                    total = count.IntValue;

                    Dictionary<string, string> attributes = new Dictionary<string, string>();
                    for (int j = 4; j < resp.Length; j++)
                    {
                        Value val = new Value(resp[j]);
                        if (val.Tag == "playerindex" && (attributes.Count > 0))
                        {
                            yield return new SlimPlayer(this, attributes);
                            attributes = new Dictionary<string, string>();
                        }
                        attributes.Add(val.Tag, val.StringValue);
                    }
                    if (attributes.Count > 0)
                    {
                        yield return new SlimPlayer(this, attributes);
                    }
                    start += chunk;
                }
                yield break;
            }
        }


        /// <summary>
        /// Get a particular player.
        /// </summary>
        /// <param name="playerID">The player ID</param>
        /// <returns>The SlimPlayer, or null if not known</returns>
        public SlimPlayer OpenPlayer(string playerID)
        {
            foreach (SlimPlayer player in Players)
            {
                if (player.PlayerID == playerID)
                {
                    return player;
                }
            }
            return null;
        }


    }

}
