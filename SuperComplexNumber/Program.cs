using System.Numerics;

var real = SuperComplex.Real;
var i = SuperComplex.Imaginary;
var th = SuperComplex.Theta;
var mu = SuperComplex.Mu;

SuperComplex v1 = 3 * real + 4 * i + 2 * th + 4 * mu;

System.Console.WriteLine(v1);
System.Console.WriteLine((1/v1+2*v1*v1).ToString("0.000"));
System.Console.WriteLine(v1.Apply(x=>1/x+2*x*x).ToString("0.000"));
