using System.Numerics;

var real = SuperComplex.Real;
var i = SuperComplex.Imaginary;
var th = SuperComplex.Theta;
var mu = SuperComplex.Mu;

SuperComplex v1 = 3 * real + 4 * i + 2 * th + 4 * mu;
SuperComplex v2 = 6 * real + 1.5 * i - 6 * th - mu;

System.Console.WriteLine(v1);
System.Console.WriteLine(v2);
System.Console.WriteLine(v1*v2);
System.Console.WriteLine(v2.Apply(x=>x+1));
System.Console.WriteLine(v2.Apply(Complex.Sin));
