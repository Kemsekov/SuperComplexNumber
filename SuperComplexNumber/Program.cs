using System.Numerics;

var real = SuperComplex.Real;
var i = SuperComplex.Imaginary;
var th = SuperComplex.Theta;
var mu = SuperComplex.Mu;

SuperComplex v = 1;
v = real+i;
System.Console.WriteLine(v);
System.Console.WriteLine((v*v).ToString("0.00"));
System.Console.WriteLine(v.ApplyComplex(x=>x*x).ToString("0.00"));