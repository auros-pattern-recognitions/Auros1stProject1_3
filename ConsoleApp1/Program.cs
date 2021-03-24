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
            string[] Si_data = File.ReadAllLines(@"..\Si_new.txt");
            int dataNum = Si_data.Length;

            // 파장, 굴절률, 소광계수 각각의 데이터를 담을 배열 선언
            List<double> waveLength = new List<double>();
            List<double> refractiveIndex = new List<double>();
            List<double> extinctionCoefficient = new List<double>();

            for (int i = 1; i < dataNum; i++)
            {
                string[] temp = Si_data[i].Split((char)0x09);
                waveLength.Add(double.Parse(temp[0]));
                refractiveIndex.Add(double.Parse(temp[1]));
                extinctionCoefficient.Add(double.Parse(temp[2]));
            }

            // 입사각에 대한 반사계수 비 (40 ~ 85, 2도 간격)



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
            Complex Nj = new Complex(1, 0); // 공기의 굴절률 = 1

            List<double> refractiveCoefficient_p = new List<double>();
            List<double> refractiveCoefficient_s = new List<double>();
            List<double> incidenceAngle = new List<double>();

            for (int theta_j = 40; theta_j <= 85; theta_j += 2)
            {

                Complex Sintheta_j = new Complex(Sin(Degree2Radian((double)theta_j)), 0);
                Complex Costheta_j = new Complex(Cos(Degree2Radian((double)theta_j)), 0);

                // 스넬 법칙

                Complex Sintheta_k = (Nj / Nk) * Sintheta_j;
                double theta_k = Asin(Sintheta_k.Real);
                Complex Costheta_k = new Complex(Cos(theta_k), 0);

                // p파 반사계수
                Complex r01p = ((Nk * Costheta_j) - (Nj * Costheta_k)) /
                                ((Nk * Costheta_j) + (Nj * Costheta_k));

                // s파 반사계수
                Complex r01s = ((Nj * Costheta_j) - (Nk * Costheta_k)) /
                                ((Nj * Costheta_j) + (Nk * Costheta_k));


                incidenceAngle.Add(theta_j);
                refractiveCoefficient_p.Add(Pow(r01p.Magnitude, 2));
                refractiveCoefficient_s.Add(Pow(r01s.Magnitude, 2));
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

            int LoopNum = incidenceAngle.Count();
            // 파일 쓰기.
            using (StreamWriter NewSpectrumOutputFile = new StreamWriter("r01p, r01s_test.dat"))
            {
                // 컬럼 명 쓰기.
                NewSpectrumOutputFile.WriteLine(
                    "AOI" + "\t"
                    + "r01p" + "\t"
                    + "r01s");    // 컬럼명 쓰기.
                // WriteLine(Columns);

                // 스펙트럼 데이터 쓰기.
                for (int i = 0; i < LoopNum; i++)
                {
                    // tsv 데이터 형식으로 데이터를 쓴다.
                    NewSpectrumOutputFile.WriteLine(
                        incidenceAngle[i] + "\t"
                        + refractiveCoefficient_p[i] + "\t"
                        + refractiveCoefficient_s[i]);

                }
            }
            #endregion

            #region 반사계수를 이용한 alpha, beta 계산

            double polarizerAngle = Degree2Radian(45.0);

            for (int theta_j = 40; theta_j <= 85; theta_j += 5)
            {
                string filename = theta_j.ToString();
                using (StreamWriter NewSpectrumOutputFile = new StreamWriter("AOI_" + filename + "_ab.dat"))
                {
                    // 컬럼 명 쓰기.
                    NewSpectrumOutputFile.WriteLine(
                        "waveLength" + "\t"
                        + "alpha" + "\t"
                        + "beta");    // 컬럼명 쓰기.

                    double radianTheta_j = Degree2Radian((double)theta_j);
                    Complex Sintheta_j2 = new Complex(Sin(radianTheta_j), 0);
                    Complex Costheta_j2 = new Complex(Cos(radianTheta_j), 0);

                    for (int i = 0; i < waveLength.Count; i++)
                    {


                        Complex Nk2 = new Complex(refractiveIndex[i], -extinctionCoefficient[i]);
                        Complex Nj2 = new Complex(1, 0); // 공기의 굴절률 = 1

                        // 스넬 법칙
                        Complex Sintheta_k2 = (Nj2 / Nk2) * Sintheta_j2;
                        double theta_k2 = Asin(Sintheta_k2.Real);
                        Complex Costheta_k2 = new Complex(Cos(theta_k2), 0);

                        // p파 반사계수
                        Complex r01p = ((Nk2 * Costheta_j2) - (Nj2 * Costheta_k2)) /
                                        ((Nk2 * Costheta_j2) + (Nj2 * Costheta_k2));

                        // s파 반사계수
                        Complex r01s = ((Nj2 * Costheta_j2) - (Nk2 * Costheta_k2)) /
                                        ((Nj2 * Costheta_j2) + (Nk2 * Costheta_k2));

                        Complex ratio = r01p / r01s;

                        //Complex ratioComp = Complex.Exp(Complex.Log(ratio));

                        double Psi = Atan(ratio.Magnitude);
                        double Delta = ratio.Phase;

                        double alpha = (Pow(Tan(Psi), 2) - Pow(Tan(polarizerAngle), 2)) /
                                        (Pow(Tan(Psi), 2) + Pow(Tan(polarizerAngle), 2));

                        double beta = (2 * Tan(Psi) * Cos(Delta) * Tan(polarizerAngle)) /
                                        (Pow(Tan(Psi), 2) + Pow(Tan(polarizerAngle), 2));

                        NewSpectrumOutputFile.WriteLine(
                            waveLength[i] + "\t"
                            + alpha + "\t"
                            + beta);


                    }
                    //WriteLine("======================================================");
                }
            }
            #endregion

            #region 측정값과 계산값의 MSE 계산

            List<double> alpha_cal = new List<double>();
            List<double> beta_cal = new List<double>();


            waveLength.Clear();


            string[] Si2nm_data = File.ReadAllLines(@"SiO2 2nm_on_Si_new.dat");
            int DataNum = Si2nm_data.Length;
            for (int i = 1; i < DataNum; i++)
            {
                string[] temp = Si2nm_data[i].Split((char)0x09);
                waveLength.Add(double.Parse(temp[0]));
                alpha_cal.Add(double.Parse(temp[2]));
                beta_cal.Add(double.Parse(temp[3]));
            }


            string[] Si_exp_data = File.ReadAllLines(@"AOI_65_ab.dat");
            DataNum = Si_exp_data.Length;

            List<double> alpha_exp = new List<double>();
            List<double> beta_exp = new List<double>();

            for (int i = 1; i < DataNum; i++)
            {
                string[] temp = Si_exp_data[i].Split((char)0x09);
                waveLength.Add(double.Parse(temp[0]));
                alpha_exp.Add(double.Parse(temp[1]));
                beta_exp.Add(double.Parse(temp[2]));
            }

            double sum = 0;

            using (StreamWriter NewSpectrumOutputFile = new StreamWriter("MSE_ab.dat"))
            {

                DataNum = alpha_exp.Count();
                for (int i = 0; i < DataNum; i++)
                {
                    double difference_MSE = Pow((alpha_exp[i] - alpha_cal[i]), 2) + Pow((beta_exp[i] - beta_cal[i]), 2);
                    sum += difference_MSE;

                    double diff_alpha = (alpha_exp[i] - alpha_cal[i]);
                    double diff_beta = (beta_exp[i] - beta_cal[i]);
                    WriteLine(waveLength[i] + " " + "alpha 의 차이 : " + diff_alpha + "    beta 의 차이 : " + diff_beta);

                    NewSpectrumOutputFile.WriteLine(waveLength[i] + "\t" + diff_alpha + "\t" + diff_beta);
                }
            }

            double MSE = sum / DataNum;
            WriteLine("최종 MSE : " + MSE);

            #endregion
        }
    }
}



