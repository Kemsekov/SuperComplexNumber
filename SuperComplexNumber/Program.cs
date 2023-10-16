using System.Numerics;

var real = SuperComplex.Real;
var i = SuperComplex.Imaginary;
var th = SuperComplex.Theta;
var mu = SuperComplex.Mu;

SuperComplex v = real+th/2+3*i-mu;
System.Console.WriteLine(v);
// System.Console.WriteLine((v-1)*(v-1));
System.Console.WriteLine(v.Apply(x=>x*x));
System.Console.WriteLine(v.ApplyOnSuperComplex(x=>x*x));
