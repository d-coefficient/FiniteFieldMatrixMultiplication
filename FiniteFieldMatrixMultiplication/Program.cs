using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using static System.Math;

namespace FiniteFieldMatrixMultiplication
{
    internal class Program
    {
        static void Main(string[] args)
        {
            const byte p = 2;
            const byte n = 2;
            var ffield = FiniteField(p, n);
            SquareMatrix.ConfigureAlgebra(2, (int)Pow(p, n), ffield.sum, ffield.product);
            var list = new int[] { 20, 66, 129, 69, 73, 77, 81, 97, 113 };
            var matrices = list.Select((x) => new SquareMatrix(x)).ToArray();
            for (int i = 0; i < list.Length; i++)
            {
                for (int j = 0; j < list.Length; j++)
                {
                    Display(SquareMatrix.Product(matrices[i], matrices[j]));
                }
            }
            Console.ReadKey(true);
        }

        static void Display(SquareMatrix matrix)
        {
            char openbracket = '[';
            char closebracket = ']';
            for (int i = 0; i < SquareMatrix.Dimension; i++)
            {
                //if (SquareMatrix.Dimension != 1)
                //{
                //    if (i == 0)
                //    {
                //        openbracket = '\u23a1';
                //        closebracket = '\u23a4';
                //    }
                //    else if (i == SquareMatrix.Dimension - 1)
                //    {
                //        openbracket = '\u23a3';
                //        closebracket = '\u23a6';
                //    }
                //    else if (i == 1)
                //    {
                //        openbracket = '\u23a2';
                //        closebracket = '\u23a5';
                //    }
                //}
                Console.Write(openbracket);
                for (int j = 0; j < SquareMatrix.Dimension; j++)
                {
                    var term = matrix.Value[j, i];
                    if (j != 0)
                    {
                        Console.Write(" ");
                    }
                    if (term < 10)
                    {
                        Console.Write(" ");
                    }
                    Console.Write(term);
                }
                Console.WriteLine(closebracket);
            }
            Console.WriteLine();
        }

        static (Func<int, int, int> sum, Func<int, int, int> product) FiniteField(byte prime, byte exponent)
        {
            Func<int, byte[]> vectorize = (int x) =>
            {
                x = (int)(x % Pow(prime, exponent));
                var v = new byte[exponent];
                for (int i = 0; i < exponent; i++)
                {
                    v[i] = (byte)(x % prime);
                    x /= prime;
                }
                return v;
            };

            Func<byte[], int> devectorize = (byte[] v) =>
            {
                var x = 0;
                for (int i = 0; i < exponent; i++)
                {
                    x += (int)(v[i] * Pow(prime, i));
                }
                return x;
            };

            Func<int, byte[]> loosevectorize = (int x) =>
            {
                var v = new List<byte>();
                while (x > 0)
                {
                    v.Add((byte)(x % prime));
                    x /= prime;
                }
                return v.ToArray();
            };

            Func<byte[], byte[], byte[]> loosesum = (byte[] a, byte[] b) =>
            {
                var aplusb = new byte[2 * exponent];
                for (int i = 0; i < 2 * exponent; i++)
                {
                    byte x;
                    try { x = a[i]; }
                    catch (IndexOutOfRangeException) { x = 0; }
                    byte y;
                    try { y = b[i]; }
                    catch (IndexOutOfRangeException) { y = 0; }
                    aplusb[i] = (byte)((x + y) % prime);
                }
                return aplusb;
            };

            Func<byte[], byte[], byte[]> sum = (byte[] a, byte[] b) =>
            {
                return loosesum(a, b).Take(exponent).ToArray();
            };

            Func<int, byte, byte[], byte[]> subproduct = (int i, byte alpha, byte[] b) =>
            {
                var alphatimesb = new byte[exponent + i];
                for (int j = 0; j < exponent; j++)
                {
                    alphatimesb[i + j] = (byte)(alpha * b[j] % prime);
                }
                return alphatimesb;
            };

            Func<byte[], byte, byte> apply = (byte[] polynomial, byte argument) =>
            {
                byte value = 0;
                for (int i = 0; i < polynomial.Length; i++)
                {
                    value = (byte)((value + polynomial[i] * Pow(argument, i % (prime - 1))) % prime);
                }
                return value;
            };

            var modulus = new byte[2];

            if (exponent > 1)
            {
                var i = (int)Pow(prime, exponent);
                var test = 0;
                while (test == 0)
                {
                    i++;
                    if (i % prime == 0)
                    {
                        continue;
                    }
                    test = 1;
                    for (byte j = 0; j < prime; j++)
                    {
                        test *= apply(loosevectorize(i), j);
                    }
                }
                modulus = loosevectorize(i);
            }

            Func<byte[], byte[], byte[]> product = (byte[] a, byte[] b) =>
            {
                var atimesb = new byte[2 * exponent];
                for (int i = 0; i < exponent; i++)
                {
                    atimesb = loosesum(atimesb, subproduct(i, a[i], b));
                }
                for (int i = exponent - 1; i >= 0; i--)
                {
                    var leadingdigit = atimesb[exponent + i];
                    if (leadingdigit != 0)
                    {
                        atimesb = loosesum(atimesb, subproduct(i, (byte)(prime - leadingdigit), modulus));
                    }
                }
                return atimesb.Take(exponent).ToArray();
            };

            return ((int a, int b) => devectorize(sum(vectorize(a), vectorize(b))), (int a, int b) => devectorize(product(vectorize(a), vectorize(b))));
        }
    }

    class SquareMatrix
    {
        static int dimension = 1;
        public static int Dimension { get { return dimension; } }
        static int order = 1;
        static Func<int, int, int> sum = (int a, int b) => a + b;
        static Func<int, int, int> product = (int a, int b) => a * b;
        int[,] value = new int[dimension, dimension];
        public int[,] Value { get { return value; } }
        public SquareMatrix(long code)
        {
            for (int i = 0; i < dimension; i++)
            {
                for (int j = 0; j < dimension; j++)
                {
                    value[i, j] = (int)(code % order);
                    code /= order;
                }
            }
        }
        public SquareMatrix(int[,] value)
        {
            this.value = value;
        }
        public SquareMatrix(bool identity)
        {
            if (identity)
            {
                for (int i = 0; i < dimension; i++)
                {
                    value[i, i] = 1;
                }
            }
        }
        public static void ConfigureAlgebra(int dimension, int order, Func<int, int, int> sum, Func<int, int, int> product)
        {
            SquareMatrix.dimension = dimension;
            SquareMatrix.order = order;
            SquareMatrix.sum = sum;
            SquareMatrix.product = product;
        }
        public static SquareMatrix Sum(SquareMatrix a, SquareMatrix b)
        {
            var result = new int[dimension, dimension];
            for (int i = 0; i < dimension; i++)
            {
                for (int j = 0; j < dimension; j++)
                {
                    result[i, j] = sum(a.value[i, j], b.value[i, j]);
                }
            }
            return new SquareMatrix(result);
        }
        public static SquareMatrix Product(SquareMatrix a, SquareMatrix b)
        {
            var result = new int[dimension, dimension];
            for (int i = 0; i < dimension; i++)
            {
                for (int j = 0; j < dimension; j++)
                {
                    var term = 0;
                    for (int k = 0; k < dimension; k++)
                    {
                        term = sum(term, product(a.value[i, k], b.value[k, j]));
                    }
                    result[i, j] = term;
                }
            }
            return new SquareMatrix(result);
        }
        public long Encode()
        {
            long result = 0;
            for (int i = 0; i < dimension; i++)
            {
                for (int j = 0; j < dimension; j++)
                {
                    result += value[i, j] * (long)Pow(order, i * dimension + j);
                }
            }
            return result;
        }
    }
}
