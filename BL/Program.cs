using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace BL
{
    class Program
    {
        public class BLClass
        {
            public double Swi = 0.363;
            public double Sor = 0.205;
            public double a1 = 1;
            public double a2 = 0.78;
            public double visco = 2; // cp
            public double viscw = 1; // cp
            public double m = 2.56;
            public double n = 3.72;

            public double GetSwd(double Sw)
            {
                return (Sw - Swi) / (1 - Sor - Swi);
            }

            public double GetKro(double Swd)
            {
                return a1 * Math.Pow((1 - Swd), m);
            }

            public double GetKrw(double Swd)
            {
                return a2 * Math.Pow(Swd, n);
            }

            public double Getfw(double Swd) // Функция fw в точке Swd
            {
                double A = (a1 * viscw) / (a2 * visco);
                return Math.Pow(Swd, n) / (Math.Pow(Swd, n) + A * Math.Pow((1 - Swd), m));
            }

            public double GetdfwdSw(double Swd) // Производная fw` в точке Swd
            {
                double A = (a1 * viscw) / (a2 * visco);
                return (1 / (1 - Sor - Swi)) * A * (n * Math.Pow(Swd, (n - 1)) * Math.Pow((1 - Swd), m) + m * Math.Pow(Swd, n) * Math.Pow((1 - Swd), (m - 1))) / (Math.Pow((Math.Pow(Swd, n) + A * Math.Pow((1 - Swd), m)), 2));
            }

            public double GetSwf(double Sw)
            {
                double L = 0;
                double R = 1;
                double C = 0;

                double Swd; // Нормированная насыщенность
                double dfwdSw;
                double Y = 1;
                double eps = 1e-5;

                while (Math.Abs(Y) > eps)
                {
                    C = (L + R) * 0.5;

                    Swd = GetSwd(C);
                    dfwdSw = GetdfwdSw(Swd);

                    // Проводим линию касательной к функции fw, если в точке Swi она не пересекает ноль, значит решение пока не найдено

                    Y = Getfw(Swd) + dfwdSw * (Sw - C);

                    if (Y < -eps)
                    {
                        L = C;
                    }

                    if (Y > +eps)
                    {
                        R = C;
                    }
                }

                return C;
            } // Получить насыщенность на фронте вытеснения

            public double GetLiquidRate(double La)
            {
                return (1.127 * 0.200 * 20 * 300 * 500) / (La * 1000);
            }

            public double GetLaAverage(double Sw) // Определение средней вязкости в интревале насыщенности от Sw до (1 - Sor)
            {
                // Решение интеграла методом трапеции

                int N = 49;
                double dSw = (1 - Sor - Sw) / N;
                double Swd;
                double Kro, Krw;
                double La1, La2;
                double dfwdSw1, dfwdSw2;
                double dfwdSwi;
                double LaSum;

                double h;

                Swd = GetSwd(Sw);
                Kro = GetKro(Swd);
                Krw = GetKrw(Swd);

                La1 = 1 / (Kro / visco + Krw / viscw);

                if (dSw == 0) return La1;

                dfwdSw1 = GetdfwdSw(Swd);

                dfwdSwi = dfwdSw1;

                LaSum = 0;

                for (int i = 0; i < (N); ++i)
                {
                    Sw = Sw + dSw;

                    if (Sw > (1 - Sor)) Sw = 1 - Sor;

                    Swd = GetSwd(Sw);
                    Kro = GetKro(Swd);
                    Krw = GetKrw(Swd);
                    La2 = 1 / (Kro / visco + Krw / viscw);
                    dfwdSw2 = GetdfwdSw(Swd);

                    h = dfwdSw1 - dfwdSw2;

                    LaSum = LaSum + 0.5 * (La1 + La2) * h;

                    La1 = La2;
                    dfwdSw1 = dfwdSw2;
                }

                return LaSum / dfwdSwi;
            }
        }

        static void Main(string[] args)
        {
            BLClass BL = new BLClass();

            using (TextWriter text = new StreamWriter("BL.out"))
            {
                text.WriteLine("  Buckley-Leverett Equation");
                text.WriteLine();
                text.WriteLine("------------------------------------");
                text.WriteLine();

                text.WriteLine("\tSwi = " + BL.Swi);
                text.WriteLine("\tSor = " + BL.Sor);
                text.WriteLine("\ta1 = " + BL.a1);
                text.WriteLine("\ta2 = " + BL.a2);
                text.WriteLine("\tm = " + BL.m);
                text.WriteLine("\tn = " + BL.n);
                text.WriteLine("\tvisco = " + BL.visco + " cp");
                text.WriteLine("\tviscw = " + BL.viscw + " cp");
                text.WriteLine();

                text.WriteLine(" Relative Perms");
                text.WriteLine("------------------------------------");
                text.WriteLine();

                int N = 10;
                double dSw = (1 - BL.Sor - BL.Swi) / N; // Шаг по насыщенности

                text.WriteLine("Sw\tSwd\tKro\tKrw\tfw\tdfw/dSw");
                text.WriteLine();

                double Sw;
                double Swd;

                for (int J = 0; J <= N; ++J)
                {
                    Sw = BL.Swi + dSw * J;
                    Swd = BL.GetSwd(Sw);

                    text.WriteLine($"{(Sw):N4}\t{Swd:N4}\t{BL.GetKro(Swd):N4}\t{BL.GetKrw(Swd):N4}\t{BL.Getfw(Swd):N4}\t{BL.GetdfwdSw(Swd):N4}");
                }

                double Swf = BL.GetSwf(BL.Swi); // Насыщение на фронте вытеснения

                text.WriteLine();
                text.WriteLine($" Front Water Saturation {Swf:N4}");
                text.WriteLine();

                N = 51;
                dSw = (1 - BL.Sor - Swf) / (N - 1);
                Sw = 1 - BL.Sor;
                text.WriteLine($"Sw\tvisc.av\tQi\n");

                double Qi;
                for (int i = 0; i < N; ++i)
                {
                    Qi = BL.GetdfwdSw(BL.GetSwd(Sw));

                    if (Qi == 0)
                        Qi = Double.PositiveInfinity ;
                    else
                        Qi = 1 / Qi;

                    text.WriteLine($"{Sw:N4}\t{BL.GetLaAverage(Sw):N4}\t{Qi:N4}");

                    Sw = Sw - dSw;
                }

                text.WriteLine();
                text.WriteLine(" Computed perfomance at linear waterflooding");
                text.WriteLine(" at constant pressure drop of 500 psi");
                text.WriteLine("--------------------------------------------------");
                text.WriteLine();
                text.WriteLine("Time\tqo\tqw\tqt\tWCT\tQi\tNp\tNp\tLa");
                text.WriteLine("days\tB/D\tB/D\tB/D\t\tPV\tPV\tbbl\tcp");

                // Нулевой момент времени
                double Qibt = 1 / BL.GetdfwdSw(BL.GetSwd(Swf));
                double dQibt = Qibt / 10;
                double Labt = BL.GetLaAverage(Swf);

                double time = 0;
                double Qi_t_prev = 0;
                double qt_t_prev = 0;
                double qw = 0;
                double qo = 0;
                double WCT = 0;

                for (int J = 0; J <= 10; ++J) // Десять расчетных шагов до начала обводнения
                {
                    double Qi_t = dQibt * J;
                    double La_t = BL.visco + (Labt - BL.visco) * J / 10;
                    double qt_t = BL.GetLiquidRate(La_t);

                    if (J > 0)
                    {
                        time = time + 2 * (Qi_t - Qi_t_prev) * 160285 / (qt_t + qt_t_prev);
                    }

                    qo = qt_t;

                    if (J == 10) // Прорыв в добывающую скважину
                    {
                        var fw = BL.Getfw(BL.GetSwd(Swf));
                        qw = qt_t * fw;
                        qo = qt_t - qw;
                        WCT = fw;
                    }

                    text.WriteLine($"{time:F1}\t{qo:N2}\t{qw:N2}\t{qt_t:N2}\t{WCT:N3}\t{Qi_t:N3}\t{Qi_t:N3}\t{Qi_t*160.285:N3}\t{La_t:N4}");

                    Qi_t_prev = Qi_t;
                    qt_t_prev = qt_t;
                }


                N = 51;
                dSw = (1 - BL.Sor - Swf) / (N - 1);

                Sw = Swf;

                for (int i = 0; i < N - 2; ++i)
                {
                    Qi = 1 / BL.GetdfwdSw(BL.GetSwd(Sw));
                    var La = BL.GetLaAverage(Sw);
                    var qt_t = BL.GetLiquidRate(La);
                    
                    var fw = BL.Getfw(BL.GetSwd(Sw));

                    var Swav = Sw + (1 - fw) / BL.GetdfwdSw(BL.GetSwd(Sw));

                    qw = qt_t * fw;
                    qo = qt_t - qw;

                    WCT = fw;

                    if (WCT > 0.999) break;

                    if (i > 0)
                    {
                        time = time + 2 * (Qi - Qi_t_prev) * 160285 / (qt_t + qt_t_prev);
                        text.WriteLine($"{time:F1}\t{qo:N2}\t{qw:N2}\t{qt_t:N2}\t{WCT:N3}\t{Qi:N3}\t{Qi:N3}\t{160.285 * (Swav - BL.Swi):N3}\t{La:N4}");
                    }

                    Qi_t_prev = Qi;
                    qt_t_prev = qt_t;

                    Sw = Sw + dSw;
                }

                text.Close();
            }
        }
    }
}