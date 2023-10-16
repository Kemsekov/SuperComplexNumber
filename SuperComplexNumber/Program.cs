using System.Numerics;

var real = SuperComplex.Real;
var i = SuperComplex.Imaginary;
var th = SuperComplex.Theta;
var mu = SuperComplex.Mu;

SuperComplex v1 = 6 * real - 19 * i - 6 * th - 24 * mu;

System.Console.WriteLine(v1);
System.Console.WriteLine((1/v1+2*v1*v1).ToString("0.000"));
System.Console.WriteLine(v1.Apply(x=>1/x+2*x*x).ToString("0.000"));
System.Console.WriteLine(v1.Apply(Complex.Sin).ToString("0.000"));
System.Console.WriteLine(v1.Apply(x=>Math.Sin(x.Real+x.Imaginary)).ToString("0.000"));
