using System;
using System.Collections.Generic;
using System.IO;
using System.Text;
using System.Threading.Tasks;
using System.Linq;
using System.Diagnostics;
using SurvivalCS;
using SurvivalCS.Distributions;
using SurvivalCS.Functions;
using SurvivalCSTest;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;


namespace SurvivalCSConsole
{
    class Program
    {
        static void Main(string[] args)
        {
            var ut = new UnitTests();
            ut.TestLogLogisticWFeatures();
        }        
    }
}
