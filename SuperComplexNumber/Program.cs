using System.Numerics;

var real = SuperComplex.Real;
var i = SuperComplex.Imaginary;
var th = SuperComplex.Theta;
var mu = SuperComplex.Mu;

SuperComplex v = 1;
v = 1+3*i-th+1.5*mu;
System.Console.WriteLine(v);
System.Console.WriteLine((v*v).ToString("0.00"));
System.Console.WriteLine(v.ApplyReal(x=>x*x).ToString("0.00"));
System.Console.WriteLine(v.ApplyComplex(x=>x*x).ToString("0.00"));