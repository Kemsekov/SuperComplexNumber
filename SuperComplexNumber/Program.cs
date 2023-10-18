using System.Numerics;
using Accord;

var real = SuperComplex.Real;
var i = SuperComplex.Imaginary;
var th = SuperComplex.Theta;
var mu = SuperComplex.Mu;

SuperComplex v = 1;
v = 8*real+i-th*i+mu;
v/=v.Magnitude;
var transform = v.ApplyComplex(Complex.Sqrt);
var actual = transform.ApplyComplex(v=>Complex.Pow(v,2));

// var transform = v.ApplyComplex(Complex.Sin);
// var actual = transform.ApplyComplex(Complex.Asin);

System.Console.WriteLine("Exptected \n"+v.ToMatrix());
System.Console.WriteLine("Actual \n"+actual.ToMatrix());
System.Console.WriteLine("Exptected \n"+v.ToString("0.00"));
System.Console.WriteLine("Actual \n"+actual.ToString("0.00"));



