double func(double x)=>Math.Sqrt(x);
double approximationPoint = 1;
int order = 10;

var s2 = TaylorSeries.BuildSeriesEtalon(func,approximationPoint,order);
var s3 = TaylorSeries.BuildSeriesCustom(func,approximationPoint,order,0.01);


foreach(var s in s2)
    System.Console.WriteLine(s);
System.Console.WriteLine("------------");
foreach(var s in s3)
    System.Console.WriteLine(s);

