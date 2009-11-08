using System;
using System.Collections.Generic;
using System.Text;
using System.Xml;
using System.Security.Cryptography;
using System.Security.Cryptography.X509Certificates;
using System.Net;
using System.IO;

// Copyright (c) 2006, 2007

namespace DSPUtil
{
    public class License
    {
        // The X.509 certificate of a software publisher
        private const string CERT =
"-----BEGIN CERTIFICATE-----" +
"MIICjjCCAfegAwIBAgIBADANBgkqhkiG9w0BAQQFADCBjDELMAkGA1UEBhMCVVMx" +
"CzAJBgNVBAgTAk1BMRIwEAYDVQQHEwlDYW1icmlkZ2UxEzARBgNVBAoTClJ1bWkg" +
"U291bmQxFTATBgNVBAsTDEFwcGxpY2F0aW9uczENMAsGA1UEAxMEUnVtaTEhMB8G" +
"CSqGSIb3DQEJARYScnVtaUBydW1pc291bmQuY29tMB4XDTA3MDUwNDAxMzgwOFoX" +
"DTcxMDMwMzE5MDk1MlowgYwxCzAJBgNVBAYTAlVTMQswCQYDVQQIEwJNQTESMBAG" +
"A1UEBxMJQ2FtYnJpZGdlMRMwEQYDVQQKEwpSdW1pIFNvdW5kMRUwEwYDVQQLEwxB" +
"cHBsaWNhdGlvbnMxDTALBgNVBAMTBFJ1bWkxITAfBgkqhkiG9w0BCQEWEnJ1bWlA" +
"cnVtaXNvdW5kLmNvbTCBnzANBgkqhkiG9w0BAQEFAAOBjQAwgYkCgYEA1Xky6QkC" +
"kSCDf7ZJpBfyr5iAGSGNHyZZ7H2PuMSEFhTNAW/SyAhBE3/Nv2j2YqQht7BF9gyn" +
"I9uHjey1KcjIhaGKxtSRi5Qebr8kh3/H7dZcWRuH4C3g4a9t+x9OggJO1y6zJH3t" +
"HfU/8IeTZtD6DlNAQMbkDBIGXf9jPNGUfasCAwEAATANBgkqhkiG9w0BAQQFAAOB" +
"gQCoaJKz6Q5W+dJ2Ny8IZx0mHPebRIzK0ob+MzZqnVKY6hOr8RPpkbhrNW+cBrsQ" +
"UwlczrPYjQsvKL5MachEXRm/6d2mk/f0/3yP+IG9Ji5FwU87obOVsHjQ8sBwPl+T" +
"yDd7IxAmXwKD1zZu5E20qfvZ2e+AuB8x/0Htb5NEyeqDQA==" +
"-----END CERTIFICATE-----";

        private const string URL = "http://rumi-sound.com/lic/upd.php?guid={0}&mac={1}";
        private const string SIG = "ok,guid={0},mac={1}";

        /// <summary>
        /// Verify a signature.
        /// Throws only on unexpected failures.
        /// </summary>
        /// <param name="data">Data, utf8</param>
        /// <param name="signature">Signature, base64</param>
        /// <returns>True or false (or throws)</returns>
        public static bool Verify(string data, string signature)
        {
            // Construct the cert
            byte[] certdata = new System.Text.ASCIIEncoding().GetBytes(CERT);
            byte[] rawdata = Encoding.UTF8.GetBytes(data);

            // Compute hash of the raw data.
            SHA1 sha = new SHA1CryptoServiceProvider();
            byte[] contenthash = sha.ComputeHash(rawdata);

            // Base64-decode the signature
            int sl = (signature.Length % 4);
            if(sl>0)
            for (int j = 0; j < 4-sl; j++) { signature += "="; }
            byte[] sigdata = Convert.FromBase64String(signature);

            bool ok = false;
            X509Certificate2 cert = new X509Certificate2(certdata);
            RSAPKCS1SignatureDeformatter RSADeformatter = new RSAPKCS1SignatureDeformatter();
            RSADeformatter.SetHashAlgorithm("SHA1");
            RSADeformatter.SetKey(cert.PublicKey.Key);
            ok=RSADeformatter.VerifySignature(contenthash, sigdata);
            return ok;
        }


        /// <summary>
        /// Provide key to the server, check response.
        /// Throws on failure.
        /// </summary>
        /// <param name="mac">MAC</param>
        /// <param name="guid">GUID</param>
        public static void Update(string mac, string guid)
        {
            string url = String.Format(URL, guid, mac);
            string sig = String.Format(SIG, guid, mac);
            string val;

            HttpWebRequest req = (HttpWebRequest)WebRequest.Create(url);
            req.Timeout = 5000; // milliseconds
            HttpWebResponse rsp = (HttpWebResponse)req.GetResponse();
            XmlDocument doc = new XmlDocument();
            doc.Load(rsp.GetResponseStream());
            rsp.Close();
            val = doc.SelectSingleNode("//sig").InnerText;

            bool ok = Verify(sig, val);
            if (!ok)
            {
                throw new Exception("Signature could not be verified.");
            }
        }

    }
}
