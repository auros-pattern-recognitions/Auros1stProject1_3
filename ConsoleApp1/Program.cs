using System;
using System.Collections;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Numerics;

using static System.Console;
using static System.Math;

namespace ConsoleApp1
{
    class Program
    {
        static double Degree2Radian(double angle) { return PI * angle / 180.0; }
        static double Radian2Degree(double angle) { return angle * (180.0 / PI); }
        static void Main(string[] args)
        {
            // Si_new.txt 파일 로딩(물성값)
            string[] Si_data = File.ReadAllLines(@"..\SiN.txt");
            int dataNum = Si_data.Length;

            // 파장, 굴절률, 소광계수 각각의 데이터를 담을 배열 선언
            List<double> waveLength = new List<double>();
            List<double> refractiveIndex = new List<double>();
            List<double> extinctionCoefficient = new List<double>();

            for (int i = 2; i < dataNum; i++)
            {
                string[] temp = Si_data[i].Split((char)0x09);
                waveLength.Add(double.Parse(temp[0]));
                refractiveIndex.Add(double.Parse(temp[1]));
                extinctionCoefficient.Add(double.Parse(temp[2]));
            }

            // 입사각에 대한 반사계수 비 (40 ~ 85, 2도 간격)

            double Nj = 1.0; // 공기의 굴절률 = 1

            // 파장이 700nm인 데이터의 인덱스 검색
            int IndexOf700nm = 0;
            for (int i = 0; i < waveLength.Count; i++)
            {
                if (waveLength[i] >= 700.0)
                {
                    IndexOf700nm = i;
                    break;
                }
            }

            // 파장이 700nm일때 굴절률과 소광계수를 통해 복소 굴절률을 구한다.
            Complex Nk = new Complex(refractiveIndex[IndexOf700nm], -extinctionCoefficient[IndexOf700nm]);

            double Nk_length = Sqrt(Pow(Nk.Real, 2) + Pow(Nk.Imaginary, 2));
            double Nk_angle = Atan(Nk.Imaginary / Nk.Real);
            double Nk_exp = Nk_length * Math.Exp(Nk_angle);

            List<double> refractiveCoefficient_p = new List<double>();
            List<double> refractiveCoefficient_s = new List<double>();
            List<double> incidenceAngle = new List<double>();

            for (int theta_j = 40; theta_j <= 85; theta_j += 2)
            {
                double radianTheta_j = Degree2Radian((double)theta_j);

                // 스넬 법칙
                double theta_k = Asin((Nj * Sin(radianTheta_j)) / Nk_exp);

                // p파 반사계수
                double r01p = ((Nk_exp * Cos(radianTheta_j)) - Nj * Cos(theta_k)) /
                              ((Nk_exp * Cos(radianTheta_j)) + Nj * Cos(theta_k));

                // s파 반사계수
                double r01s = ((Nj * Cos(radianTheta_j)) - Nk_exp * Cos(theta_k)) /
                               ((Nj * Cos(radianTheta_j)) + Nk_exp * Cos(theta_k));

                incidenceAngle.Add(theta_j);
                refractiveCoefficient_p.Add(Pow(r01p, 2));
                refractiveCoefficient_s.Add(Pow(r01s, 2));
            }
            #region 반사계수 출력

            //WriteLine("=================입사각===========================");
            //foreach (var item in incidenceAngle)
            //{
            //    WriteLine(item);
            //}

            //WriteLine("===================p파 반사계수=========================");
            //foreach (var item in refractiveCoefficient_p)
            //{
            //    WriteLine(item);
            //}

            //WriteLine("===================s파 반사계수=========================");
            //foreach (var item in refractiveCoefficient_s)
            //{
            //    WriteLine(item);
            //}
            #endregion

            #region 입사각에 따른 반사율 파일 저장

            //int LoopNum = incidenceAngle.Count();
            //// 파일 쓰기.
            //using (StreamWriter NewSpectrumOutputFile = new StreamWriter("r01p, r01s_test.dat"))
            //{
            //    // 컬럼 명 쓰기.
            //    NewSpectrumOutputFile.WriteLine(
            //        "AOI" + "\t"
            //        + "r01p" + "\t"
            //        + "r01s");    // 컬럼명 쓰기.
            //    // WriteLine(Columns);

            //    // 스펙트럼 데이터 쓰기.
            //    for (int i = 0; i < LoopNum; i++)
            //    {
            //        // tsv 데이터 형식으로 데이터를 쓴다.
            //        NewSpectrumOutputFile.WriteLine(
            //            incidenceAngle[i] + "\t"
            //            + refractiveCoefficient_p[i] + "\t"
            //            + refractiveCoefficient_s[i]);

            //    }
            //}
            #endregion
            //WriteLine(Complex.Log(E));

            List<double> alpha = new List<double>();
            List<double> beta = new List<double>();

            for (int theta_j = 40; theta_j <= 85; theta_j += 5)
            {
                string filename = theta_j.ToString();
                using (StreamWriter NewSpectrumOutputFile = new StreamWriter("AOI_" + filename + ".dat"))
                {
                    // 컬럼 명 쓰기.
                    NewSpectrumOutputFile.WriteLine(
                        "waveLength" + "\t"
                        + "Psi" + "\t"
                        + "Delta");    // 컬럼명 쓰기.

                    double radianTheta_j = Degree2Radian((double)theta_j);
                    //WriteLine("======================================================");
                    for (int i = 0; i < waveLength.Count; i++)
                    {
                        Complex Nk2 = new Complex(refractiveIndex[i], -extinctionCoefficient[i]);
                        //WriteLine(waveLength[i] + " " + Nk2);

                        double Nk2_length = Sqrt(Pow(Nk2.Real, 2) + Pow(Nk2.Imaginary, 2));
                        double Nk2_angle = Atan(Nk2.Imaginary / Nk2.Real);
                        double Nk2_exp = Nk2_length * Math.Exp(Nk2_angle);


                        // 스넬 법칙
                        double theta_k = Asin((Nj * Sin(radianTheta_j)) / Nk2_exp);

                        // p파 반사계수
                        double r01p = ((Nk2_exp * Cos(radianTheta_j)) - Nj * Cos(theta_k)) /
                                      ((Nk2_exp * Cos(radianTheta_j)) + Nj * Cos(theta_k));

                        // s파 반사계수
                        double r01s = ((Nj * Cos(radianTheta_j)) - Nk2_exp * Cos(theta_k)) /
                                       ((Nj * Cos(radianTheta_j)) + Nk2_exp * Cos(theta_k));

                        double ratio = r01p / r01s;

                        //Complex ratioComp = Complex.Exp(Complex.Log(ratio));

                        double Psi = Atan(ratio);
                        double Delta = Log(Abs(ratio));

                        NewSpectrumOutputFile.WriteLine(
                            waveLength[i] + "\t"
                            + Radian2Degree(Psi) + "\t"
                            + Abs(Radian2Degree(Delta)));
                    }
                    //WriteLine("======================================================");
                }
            }
        }
    }
}



