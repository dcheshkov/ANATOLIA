/*****************************************************/
/***           ANATOLIA-X based on V1.2            ***/
/***      Free open-source software for total      ***/
/***       lineshape analysis of NMR spectra       ***/
/*****************************************************/
/***       Copyright 2025, Dmitry Cheshkov,        ***/
/***      Dmitry Sinitsyn, Kirill Sheberstov       ***/
/*****************************************************/
/***             dcheshkov@gmail.com               ***/
/***          http://anatolia.nmrclub.ru           ***/
/***     https://github.com/dcheshkov/ANATOLIA     ***/
/*****************************************************/
/***  This program is free software: you can       ***/
/***  redistribute it and/or modify it under the   ***/
/***  terms of the GNU General Public License as   ***/
/***  published by the Free Software Foundation,   ***/
/***  either version 3 of the License, or (at      ***/
/***  your option) any later version.              ***/
/*****************************************************/
/***  This program is distributed in the hope that ***/
/***  it will be useful,but WITHOUT ANY WARRANTY;  ***/
/***  without even the implied warranty of         ***/
/***  MERCHANTABILITY or FITNESS FOR A PARTICULAR  ***/
/***  PURPOSE. See the GNU General Public License  ***/
/***  for more details.                            ***/
/*****************************************************/
/***  You should have received a copy of the GNU   ***/
/***  General Public License along with this       ***/
/***  program.                                     ***/
/***  If not, see <http://www.gnu.org/licenses/>.  ***/
/*****************************************************/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstring>
#include <cstdint>
#include <cmath>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#ifdef _WIN32
#include <direct.h>
#else
#include <unistd.h>
#endif

using namespace std;

char Title[256];
int nSpins = 0;
int* Offs = NULL;
int** JCoups = NULL;
int nSSParams = 0;
double* SSParams = NULL;
int nFreqsFiltered = 0;
double* FreqsFiltered = NULL;
double* IntensFiltered = NULL;

#define defaultprecision 6
#define IntensityThreshold 0.005

int nBroadenings = 0;
double* LBs = NULL;

char textline[256];

#define uint_t uint64_t
#define BF_t uint64_t

#ifdef _WIN32
#define chdir _chdir
#define exit_ exit_with_pause()
bool ExitWithPause = true;
void exit_with_pause(void)
{

	if (!ExitWithPause) exit(0);
	char textline[8];
	cout << "Press any key..." << endl;
	cin.getline(textline, 8);
	exit(0);

}
#else
#define exit_ exit(0)
#endif

void print_logo(ostream& ostr)
{

	ostr
		<< "*****************************************************" << endl
		<< "***           ANATOLIA-X based on V1.2            ***" << endl
		<< "***      Free open-source software for total      ***" << endl
		<< "***       lineshape analysis of NMR spectra       ***" << endl
		<< "*****************************************************" << endl
		<< "***       Copyright 2025, Dmitry Cheshkov,        ***" << endl
		<< "***      Dmitry Sinitsyn, Kirill Sheberstov       ***" << endl
		<< "*****************************************************" << endl
		<< endl;

}

void print_citation(ostream& ostr)
{

	ostr
		<< "When publishing results obtained with this program, please cite:" << endl
		<< "D.A. Cheshkov, K.F. Sheberstov, D.O. Sinitsyn, V.A. Chertkov, ANATOLIA: NMR" << endl
		<< "software for spectral analysis of total lineshape. Magn. Reson. Chem., 2018," << endl
		<< "56, 449, DOI: 10.1002/mrc.4689." << endl;

}

bool isunsignint(char* text)
{

	bool res = int(text[0]) != 0;
	for (int i = 0; int(text[i]) != 0; i++)
		res = res && isdigit(text[i]);
	return(res);

}

bool isunsignreal(char* text)
{

	bool res = int(text[0]) != 0;
	int point = 0;
	for (int i = 0; int(text[i]) != 0; i++)
	{
		if (text[i] == '.') point++;
		res = res && (point < 2) && (isdigit(text[i]) || text[i] == '.');
	}
	return(res);

}

bool isreal(char* text)
{

	int j = (text[0] == '-') || (text[0] == '+') ? 1 : 0;
	return isunsignreal(text + j);

}

bool isemptyline(char* text)
{

	bool res = true;
	for (int i = 0; int(text[i]) != 0; i++)
		res = res && isspace(text[i]);
	return(res);

}

inline int getbitsum(BF_t number)
{
    //int result = 0;
    //BF_t mask = (BF_t)1 << (nSpins - 1);
    //for (int i = 1; i <= nSpins; i++)
    //{
    //    if (number & mask) result++;
    //    mask /= 2;
    //}
    //return result;
#ifdef _MSC_VER
	return (int)__popcnt64(number);
#else
	return __builtin_popcountll(number);
#endif

}

#define BOBYQA_SUCCESS           0 // algorithm converged
#define BOBYQA_BAD_NPT          -1 // NPT is not in the required interval
#define BOBYQA_TOO_CLOSE        -2 // insufficient space between the bounds
#define BOBYQA_ROUNDING_ERRORS  -3 // too much cancellation in a denominator
#define BOBYQA_STEP_FAILED      -5 // a trust region step has failed to reduce Q

int bobyqa(long n, long npt, double (*objfun)(double*), double* x,
	double* xl, double* xu, double rhobeg, double rhoend, double* w);

class Data
{

public:
	int nPoints;
	double Magnitude;
	double* Points;
	char* Filename;

	Data(void)
	{

		Filename = new char[256];
		strcpy(Filename, "");
		nPoints = 0;
		Magnitude = 0;
		Points = NULL;

	}

	void Create(void)
	{

		if (Points != NULL) { delete[] Points; Points = NULL; }

		if (nPoints > 0)
		{
			Points = new double[(uint_t)nPoints + 1];
			for (int i = 0; i <= nPoints; i++)
				Points[i] = 0;
		}

	}

	void Clean(void)
	{

		if (Points != NULL) { delete[] Points; Points = NULL; }
		nPoints = 0;
		Magnitude = 0;

	}

	void Rescale(int PivotPointIndex)
	{

		double rescaling = Magnitude / (PivotPointIndex == 0 ? CalcMagnitude() : Points[PivotPointIndex]);
		for (int i = 1; i <= nPoints; i++)
			Points[i] = rescaling * Points[i];

	}

	double CalcMagnitude(void)
	{

		double MaxIntens = 0;
		for (int i = 1; i <= nPoints; i++)
			if (MaxIntens < abs(Points[i])) MaxIntens = abs(Points[i]);
		return(MaxIntens);

	}

	void SaveSpecToFile(void)
	{

		if (Filename[0] == '-') return;
		int point = 0;
		ofstream ostr(Filename, ios::out | ios::binary);
		for (int i = 1; i <= nPoints; i++)
		{
			point = int(Points[i]);
			ostr.write((char*)&point, 4);
		}
		ostr.close();

	}

	void LoadSpecFromFile(void)
	{

		int point = 0;

		if (int(Filename[0]) == 0) { cout << "Can't load spectrum, filename is empty!" << endl; exit_; }
		ifstream istr(Filename, ios::in | ios::binary | ios::ate);
		if (!istr) { cout << "File " << Filename << " does not exists!" << endl; exit_; }

		nPoints = (int)istr.tellg() / 4;
		istr.seekg(0);
		if (Points != NULL) delete[] Points;
		Points = new double[(uint_t)nPoints + 1];
		Points[0] = 0;

		Magnitude = 0;
		for (int i = 1; i <= nPoints; i++)
		{
			istr.read((char*)&point, 4);
			Points[i] = double(point);
			if (Magnitude < abs(Points[i])) Magnitude = abs(Points[i]);
		}
		istr.close();

	}

	void ZeroData(void)
	{

		for (int i = 0; i <= nPoints; i++)
			Points[i] = 0;

	}

};

class Spectrum
{

public:
	char* DatasetPath;
	int ExpProcNo;
	int BrExpProcNo;
	int TheorProcNo;
	Data ExperimentalSpec;
	Data ExperimentalSpecWithBroadening;
	Data TheoreticalSpec;
	//double LineWidth; => Global LineWidth
	double Offset;
	double FreqStep;
	double SF;
	double LB;
	int nIntervals;
	int* StartPoint;
	int* EndPoint;
	int nPointsRated;
	double* ExpSpecPointsOnIntervals;
	double* TheorSpecPointsOnIntervals;
	double* FreqsOnIntervals;
	double* Lorentz;
	double SumOfExpSquaresOnIntervals;
	char* SpectraTextOutputFilename;
	bool ScaleOpt;
	bool FineCalc;

	Spectrum(ifstream& istr)
	{

		DatasetPath = new char[256];
		SpectraTextOutputFilename = new char[256];
		SF = 0;
		LB = 0;
		ScaleOpt = false;
		FineCalc = true;
		SumOfExpSquaresOnIntervals = 0;
		Lorentz = NULL;

		istr.getline(textline, 256); // Spectra parameters
		istr >> textline >> DatasetPath;
		for (int i = 0; int(DatasetPath[i]) != 0; i++) if (DatasetPath[i] == '\\') DatasetPath[i] = char('/');

		// Reading data from ACQUS
		double BF = 0, O1_h = 0, SW_h = 0;
		sprintf(textline, "%s/acqus", DatasetPath);
		ifstream dataset(textline);
		if (!dataset) { cout << "File ACQUS does not exists!" << endl; exit_; }
		while (dataset)
		{
			dataset >> textline;
			if (strstr(textline, "$BF1=") != NULL)  dataset >> BF;
			if (strstr(textline, "$O1=") != NULL)   dataset >> O1_h;
			if (strstr(textline, "$SW_h=") != NULL) dataset >> SW_h;
			if (strstr(textline, "$SW_H=") != NULL) dataset >> SW_h;
		}
		dataset.close();
		if(BF == 0 || SW_h == 0) { cout << "File ACQUS corrupted!" << endl; exit_; }

		istr >> textline >> textline;
		if (!isunsignint(textline)) { cout << "Wrong experimental spectrum processing number." << endl; exit_; }
		ExpProcNo = atoi(textline);

		// Reading data from PROCS of experimental spectrum
		sprintf(ExperimentalSpec.Filename, "%s/pdata/%i/procs", DatasetPath, ExpProcNo);
		dataset.open(ExperimentalSpec.Filename);
		if (!dataset) { cout << "File PROCS in processing number " << ExpProcNo << " does not exists!" << endl; exit_; }
		sprintf(ExperimentalSpec.Filename, "%s/pdata/%i/1r", DatasetPath, ExpProcNo);
		while (dataset)
		{
			dataset >> textline;
			if (strstr(textline, "$SF=") != NULL) dataset >> SF;
			if (strstr(textline, "$SI=") != NULL) dataset >> TheoreticalSpec.nPoints;
		}
		dataset.close();

		if(SF == 0 || TheoreticalSpec.nPoints == 0) { cout << "File PROCS in processing number " << ExpProcNo << " corrupted!" << endl; exit_; }
		ExperimentalSpecWithBroadening.nPoints = TheoreticalSpec.nPoints;

		Offset = (BF - SF) * 1e6 + O1_h + 0.5 * SW_h; // SR = SF - BF; SR *= 1e6; Offset = O1_h - SR + 0.5 * SW_h;
		FreqStep = SW_h / TheoreticalSpec.nPoints;

		istr >> textline >> textline;
		if (textline[0] == '-')
			strcpy(ExperimentalSpecWithBroadening.Filename, "-");
		else
		{
			if (!isunsignint(textline)) { cout << "Check the processing number for experimental spectrum with broadening." << endl; exit_; }
			BrExpProcNo = atoi(textline);
			sprintf(ExperimentalSpecWithBroadening.Filename, "%s/pdata/%i/1r", DatasetPath, BrExpProcNo);
		}

		istr >> textline >> textline;
		if (textline[0] == '-')
			strcpy(TheoreticalSpec.Filename, "-");
		else
		{
			if (!isunsignint(textline)) { cout << "Check the processing number for theoretical spectrum." << endl; exit_; }
			TheorProcNo = atoi(textline);
			sprintf(TheoreticalSpec.Filename, "%s/pdata/%i/1r", DatasetPath, TheorProcNo);
		}
		istr.getline(textline, 256); // Rest of proc. no. line
		istr.getline(textline, 256); // Empty line
		if (!isemptyline(textline)) { cout << "Empty line should follow the section with spectra parameters!" << endl; exit_; }

		LoadIntervals();

		ExperimentalSpec.LoadSpecFromFile();
		if (ExperimentalSpec.nPoints != TheoreticalSpec.nPoints) { cout << "File with experimental spectrum is corrupted." << endl; exit_; }

		//CalcExpSpecMagnOnIntervals
		ExperimentalSpec.Magnitude = 0;
		for (int i = 1; i <= nIntervals; i++)
			for (int j = StartPoint[i]; j <= EndPoint[i]; j++)
				if (ExperimentalSpec.Magnitude < abs(ExperimentalSpec.Points[j])) ExperimentalSpec.Magnitude = abs(ExperimentalSpec.Points[j]);

		ExperimentalSpecWithBroadening.Magnitude = ExperimentalSpec.Magnitude;

		TheoreticalSpec.Create();
		ExperimentalSpecWithBroadening.Create();

	}

	void LoadIntervals(void)
	{

		stringstream parse;
		int position = 0;
		bool ok = false;
		sprintf(textline, "%s/pdata/%i/integrals.txt", DatasetPath, ExpProcNo);
		ifstream istr(textline);
		if (!istr) { cout << "File integrals.txt does not exist in the processing folder with experimental spectrum!" << endl; exit_; }

		nIntervals = 0;

		while (istr)
		{
			istr.getline(textline, 256); position++;
			if (strstr(textline, "Integral") != NULL) { ok = true; break; }
		}
		if (!ok) { cout << "File integrals.txt in the processing folder with experimental spectrum is corrupted!" << endl; exit_; }

		int i = 0;
		istr.getline(textline, 256);
		while (istr)
		{
			i++;
			parse << textline;
			parse >> nIntervals;
			if (nIntervals != i) { cout << "File integrals.txt in the processing folder with experimental spectrum is corrupted!" << endl; exit_; }
			parse.str(""); parse.clear();
			istr.getline(textline, 256);
		}

		if (nIntervals < 1) { cout << "There is no defined intervals in the integrals.txt file for the experimental spectrum!" << endl; exit_; }

		istr.clear(); istr.seekg(0);
		for (i = 1; i <= position; i++) istr.getline(textline, 256);

		int j = 0;
		double StartPPM = 0, EndPPM = 0, IntVal = 0;
		double PPMOffset = Offset / SF;
		double PPMStep = FreqStep / SF;

		StartPoint = new int[(uint_t)nIntervals + 1];
		EndPoint = new int[(uint_t)nIntervals + 1];
		StartPoint[0] = 0;
		EndPoint[0] = 0;
		nPointsRated = 0;

		for (i = 1; i <= nIntervals; i++)
		{
			istr.getline(textline, 256);
			parse << textline;
			parse >> j >> StartPPM >> EndPPM >> IntVal;
			parse.str(""); parse.clear();
			StartPoint[i] = int((PPMOffset - StartPPM) / PPMStep) - 1;
			if (StartPoint[i] < 1) StartPoint[i] = 1;
			if (StartPoint[i] > TheoreticalSpec.nPoints) StartPoint[i] = TheoreticalSpec.nPoints;
			EndPoint[i] = int((PPMOffset - EndPPM) / PPMStep) + 1;
			if (EndPoint[i] < 1) EndPoint[i] = 1;
			if (EndPoint[i] > TheoreticalSpec.nPoints) EndPoint[i] = TheoreticalSpec.nPoints;
			if (EndPoint[i] - StartPoint[i] <= 2) { cout << "Spectrum interval " << i << " is too small. Please increase the number of points." << endl; exit_; }
			nPointsRated += (EndPoint[i] - StartPoint[i] + 1);
		}

		istr.close();

		ExpSpecPointsOnIntervals = new double[(uint_t)nPointsRated + 1];
		TheorSpecPointsOnIntervals = new double[(uint_t)nPointsRated + 1];
		FreqsOnIntervals = new double[(uint_t)nPointsRated + 1];
		for (i = 0; i <= nPointsRated; i++)
		{
			ExpSpecPointsOnIntervals[i] = 0;
			TheorSpecPointsOnIntervals[i] = 0;
			FreqsOnIntervals[i] = 0;
		}

		int n = 1;
		for (int i = 1; i <= nIntervals; i++)
			for (int j = StartPoint[i]; j <= EndPoint[i]; j++)
				FreqsOnIntervals[n++] = Offset - FreqStep * ((double)j - 1);

	}

	void CheckSpinOffsets(ostream& ostr)
	{

		int MaxOffs = Offs[nSpins];

		bool check = false;

		for (int i = 1; i <= MaxOffs; i++)
		{
			double CurrOffset = SSParams[i];
			check = false;
			for (int j = 1; j <= nIntervals; j++)
				if (CurrOffset <= (Offset - FreqStep * (StartPoint[j] - 1)) && CurrOffset >= (Offset - FreqStep * (EndPoint[j] - 1))) check = true;
			if (!check) ostr << "Warning! Chemical shift no. " << i << " (" << CurrOffset << ") does not fall into any of defined spectral intervals." << endl;
		}

		for (int i = 1; i <= nIntervals; i++)
		{
			double StrartFreq = Offset - FreqStep * (StartPoint[i] - 1);
			double EndFreq = Offset - FreqStep * (EndPoint[i] - 1);
			check = false;
			for (int j = 1; j <= MaxOffs; j++)
				if (SSParams[j] <= StrartFreq && SSParams[j] >= EndFreq) check = true;
			if (!check) ostr << "Warning! Spectral interval no. " << i << " does not contain any chemical shift." << endl;
		}

		if (!check) ostr << endl;

	}

	void InitBroadening(double MaxLb)
	{

		int nPoints = int(ceill(sqrt(199) * MaxLb / (2 * FreqStep)));
		if (nPoints > 3) Lorentz = new double[(uint_t)nPoints + 1];

	}

	void BroadOnIntervals(double LWBroad)
	{

		LB = LWBroad;
		ExperimentalSpecWithBroadening.ZeroData();
		int step = 1;

		int nPoints = int(ceill(sqrt(199) * LB / (2 * FreqStep)));
		if (nPoints > 3)
		{
			double Betta = 2 * FreqStep / LB; Betta *= Betta;
			for (int i = 0; i <= nPoints; i++)
				Lorentz[i] = 1 / (1 + Betta * i * i);

			for (int i = 1; i <= nIntervals; i++)
				for (int j = StartPoint[i]; j <= EndPoint[i]; j++)
				{
					double* Point = ExperimentalSpecWithBroadening.Points + j;
					int startk = j - nPoints; if (startk <= 0) startk = 1;
					int stopk = j + nPoints;  if (stopk > ExperimentalSpec.nPoints) stopk = ExperimentalSpec.nPoints;
					for (int k = startk; k <= stopk; k++)
						*Point += ExperimentalSpec.Points[k] * Lorentz[abs(j - k)];
				}
			ExperimentalSpecWithBroadening.Rescale(0);
			step = int(LB / (2 * FreqStep)); if ((step == 0) || FineCalc) step = 1;
		}
		else
		{
			LB = 0;
			for (int i = 1; i <= nIntervals; i++)
				for (int j = StartPoint[i]; j <= EndPoint[i]; j++)
					ExperimentalSpecWithBroadening.Points[j] = ExperimentalSpec.Points[j];
		}

		SumOfExpSquaresOnIntervals = 0;
		int n = 1;
		double tmp = 0;
		for (int i = 1; i <= nIntervals; i++)
			for (int j = StartPoint[i]; j <= EndPoint[i]; j += step)
			{
				FreqsOnIntervals[n] = Offset - FreqStep * ((double)j - 1);
				tmp = ExperimentalSpecWithBroadening.Points[j];
				ExpSpecPointsOnIntervals[n++] = tmp; tmp *= tmp;
				SumOfExpSquaresOnIntervals += tmp;
			}
		nPointsRated = n - 1;

	}

	void CalcSpecOnIntervals(void)
	{
		if ((!FineCalc) && ScaleOpt) { CalcSpecOnIntervalsRoughWithScOpt(); return; }
		if (!FineCalc) { CalcSpecOnIntervalsRoughFixedSc(); return; }
		if (ScaleOpt) { CalcSpecOnIntervalsFineWithScOpt(); return; }
		else { CalcSpecOnIntervalsFineFixedSc(); return; }
	}

	void CalcSpecOnIntervalsRoughWithScOpt(void)
	{

		double LW = (abs(SSParams[nSSParams]) + LB) / 2, tmp = 0, tmp1 = 0, tmp2 = 0, tmp3 = 0, tmp4 = 0, CurrFreq = 0;
		double sqLW = LW * LW;

		tmp2 = 25 * LW;
		for (int i = 1; i <= nPointsRated; i++)
		{
			CurrFreq = FreqsOnIntervals[i];
			tmp1 = 0;
			for (int j = 1; j <= nFreqsFiltered; j++)
			{
				tmp = abs(CurrFreq - FreqsFiltered[j]);
				if (tmp < tmp2)
				{
					tmp *= tmp;
					tmp1 += IntensFiltered[j] / (tmp + sqLW);
				}
			}
			TheorSpecPointsOnIntervals[i] = tmp1;
			tmp3 += tmp1 * tmp1;
			tmp4 += ExpSpecPointsOnIntervals[i] * tmp1;
		}

		// Rescale TheorSpec on intervals
		tmp = tmp4 / tmp3;
		for (int i = 1; i <= nPointsRated; i++)
			TheorSpecPointsOnIntervals[i] *= tmp;

	}

	void CalcSpecOnIntervalsRoughFixedSc(void)
	{

		double LW = (abs(SSParams[nSSParams]) + LB) / 2, tmp = 0, tmp1 = 0, tmp2 = 0, CurrFreq = 0;
		double sqLW = LW * LW;
		double MaxIntense = 0;

		tmp2 = 25 * LW;
		for (int i = 1; i <= nPointsRated; i++)
		{
			CurrFreq = FreqsOnIntervals[i];
			tmp1 = 0;
			for (int j = 1; j <= nFreqsFiltered; j++)
			{
				tmp = abs(CurrFreq - FreqsFiltered[j]);
				if (tmp < tmp2)
				{
					tmp *= tmp;
					tmp1 += IntensFiltered[j] / (tmp + sqLW);
				}
			}
			TheorSpecPointsOnIntervals[i] = tmp1;
			if (MaxIntense < tmp1) MaxIntense = tmp1;
		}

		// Rescale TheorSpec on intervals
		tmp = abs(TheoreticalSpec.Magnitude / MaxIntense);
		for (int i = 1; i <= nPointsRated; i++)
			TheorSpecPointsOnIntervals[i] *= tmp;

	}

	void CalcSpecOnIntervalsFineWithScOpt(void)
	{

		double LW = (abs(SSParams[nSSParams]) + LB) / 2, tmp = 0, tmp1 = 0, tmp2 = 0, tmp3 = 0, CurrFreq = 0;
		double sqLW = LW * LW;

		for (int i = 1; i <= nPointsRated; i++)
		{
			CurrFreq = FreqsOnIntervals[i];
			tmp1 = 0;
			for (int j = 1; j <= nFreqsFiltered; j++)
			{
				tmp = CurrFreq - FreqsFiltered[j]; tmp *= tmp;
				tmp1 += IntensFiltered[j] / (tmp + sqLW);
			}
			TheorSpecPointsOnIntervals[i] = tmp1;
			tmp2 += tmp1 * tmp1;
			tmp3 += ExpSpecPointsOnIntervals[i] * tmp1;
		}

		// Rescale TheorSpec on intervals
		tmp = tmp3 / tmp2;
		for (int i = 1; i <= nPointsRated; i++)
			TheorSpecPointsOnIntervals[i] *= tmp;

	}

	void CalcSpecOnIntervalsFineFixedSc(void)
	{

		double LW = (abs(SSParams[nSSParams]) + LB) / 2, tmp = 0, tmp1 = 0, CurrFreq = 0;
		double sqLW = LW * LW;
		double MaxIntense = 0;

			for (int i = 1; i <= nPointsRated; i++)
			{
				CurrFreq = FreqsOnIntervals[i];
				tmp1 = 0;
				for (int j = 1; j <= nFreqsFiltered; j++)
				{
					tmp = CurrFreq - FreqsFiltered[j]; tmp *= tmp;
					tmp1 += IntensFiltered[j] / (tmp + sqLW);
				}
				TheorSpecPointsOnIntervals[i] = tmp1;
				if (MaxIntense < tmp1) MaxIntense = tmp1;
			}

		// Rescale TheorSpec on intervals
		tmp = abs(TheoreticalSpec.Magnitude / MaxIntense);
		for (int i = 1; i <= nPointsRated; i++)
			TheorSpecPointsOnIntervals[i] *= tmp;

	}

	void UpdateTheorSpecMagnitude(void)
	{

		double MaxIntense = 0;
		for (int i = 1; i <= nPointsRated; i++)
			if (MaxIntense < TheorSpecPointsOnIntervals[i]) MaxIntense = TheorSpecPointsOnIntervals[i];
		TheoreticalSpec.Magnitude = MaxIntense;

	}

	double Badness(void)
	{
		CalcSpecOnIntervals();
		return SumOfSquareDeviationsOnIntervals();
	}

	double SumOfSquareDeviationsOnIntervals(void)
	{
		double tmp = 0, res = 0;

		for (int i = 1; i <= nPointsRated; i++)
		{
			tmp = ExpSpecPointsOnIntervals[i] - TheorSpecPointsOnIntervals[i]; tmp *= tmp;
			res += tmp;
		}
		return res;
	}

	void CalcFullSpectrum(void)
	{

		double LW = (abs(SSParams[nSSParams]) + LB) / 2;
		double tmp = 0, tmp1 = 0, sqLW = LW * LW, CurrFreq = 0;

		for (int i = 1; i <= TheoreticalSpec.nPoints; i++)
		{
			CurrFreq = Offset - FreqStep * ((double)i - 1);
			tmp1 = 0;
			for (int j = 1; j <= nFreqsFiltered; j++)
			{
				tmp = CurrFreq - FreqsFiltered[j]; tmp *= tmp;
				tmp1 += IntensFiltered[j] / (tmp + sqLW);
			}
			TheoreticalSpec.Points[i] = tmp1;
		}

		// Rescaling of the full spectrum taking into account spectral intervals.
		int index = 0;
		tmp = 0;
		for (int i = 1; i <= nIntervals; i++)
			for (int j = StartPoint[i]; j <= EndPoint[i]; j++)
				if (tmp < TheoreticalSpec.Points[j]) { tmp = TheoreticalSpec.Points[j]; index = j; }
		TheoreticalSpec.Rescale(index);

	}

	double CalcRFactor(void)
	{

		double rfactor = 1;
		if (SumOfExpSquaresOnIntervals > 0)
			rfactor = SumOfSquareDeviationsOnIntervals() / SumOfExpSquaresOnIntervals;
		return(100 * sqrt(rfactor));

	}

	double CalcSpectraColleration(void)
	{

		if (SumOfExpSquaresOnIntervals > 0)
		{
			double SumOfTheorSquaresOnIntervals = 0;
			double tmp = 0;
			for (int i = 1; i <= nPointsRated; i++)
			{
				SumOfTheorSquaresOnIntervals += TheorSpecPointsOnIntervals[i] * TheorSpecPointsOnIntervals[i];
				tmp += TheorSpecPointsOnIntervals[i] * ExpSpecPointsOnIntervals[i];
			}
			return tmp / sqrt(SumOfTheorSquaresOnIntervals * SumOfExpSquaresOnIntervals);
		}

		return 0;

	}

	void SaveSpecsOnIntervalsTXT(void)
	{

		if (SpectraTextOutputFilename[0] == '-') return;

		ofstream ostr(SpectraTextOutputFilename);
		ostr << setw(16) << left << "Freq(Hz)";
		ostr << setw(14) << right << "Exp.Intens.";
		ostr << setw(16) << right << "Theor.Intens.";
		ostr << setw(13) << right << "Diff." << endl;

		ostr.precision(defaultprecision);
		ostr << fixed;
		for (int i = 1; i <= nPointsRated; i++)
		{
			ostr << setw(11) << left << FreqsOnIntervals[i];
			ostr << setw(16) << right << int(ExpSpecPointsOnIntervals[i]);
			ostr << setw(16) << right << int(TheorSpecPointsOnIntervals[i]);
			ostr << setw(16) << right << int(ExpSpecPointsOnIntervals[i]) - int(TheorSpecPointsOnIntervals[i]) << endl;
		}

		ostr.close();

	}

};

class Hamiltonian
{
public:
	int nSpinGroups;
	int* SpinsGroup;
	int* nSpinsInGroup = NULL;
	BF_t* GroupMasks = NULL;
	int nBl = 0;
	int* BlockSize = NULL;
	BF_t** BFInBlock = NULL;
	gsl_matrix** Ham;
	gsl_matrix** EVec;
	gsl_vector** EVal;
	gsl_eigen_symmv_workspace** W;
	int**** Perturb = NULL;
	int* nBlockPairs = NULL;
	int** BlockFrom = NULL;
	int** BlockTo = NULL;
	int* nFreqs = NULL;
	double** Freqs = NULL;
	double** Intens = NULL;
	double* Perturbik = NULL;

	Hamiltonian()
	{

		// --- 1. Determination of strongly coupled spin groups ---
		nSpinGroups = 0;
		SpinsGroup = new int[nSpins + 1]();

		for (int i = 1; i <= nSpins; i++)
		{
			if (SpinsGroup[i] == 0) SpinsGroup[i] = ++nSpinGroups;
			for (int j = i + 1; j <= nSpins; j++)
				if ((abs(SSParams[JCoups[i][j]]) > 1e-6 ? abs(((SSParams[Offs[i]] - SSParams[Offs[j]]) / SSParams[JCoups[i][j]])) : abs(1e6 * (SSParams[Offs[i]] - SSParams[Offs[j]]))) < 10) SpinsGroup[j] = SpinsGroup[i];
		}

		int Counter = 1;
		for (int i = 2; i <= nSpins; i++)
		{
			if (SpinsGroup[i] < Counter || SpinsGroup[i] > Counter + 1) { cout << "Error: Strongly coupled spins must be ordered sequentially." << endl; exit(1); }
			Counter = SpinsGroup[i];
		}

		nSpinsInGroup = new int[nSpinGroups + 1]();
		GroupMasks = new BF_t[nSpinGroups + 1]();

		for (int i = 1; i <= nSpins; i++) {
			int g = SpinsGroup[i];
			nSpinsInGroup[g]++;
			GroupMasks[g] |= ((BF_t)1 << (i - 1));
		}

		BF_t nBF = (BF_t)1 << nSpins;
		nBl = 1;
		for (int i = 1; i <= nSpinGroups; i++) nBl *= (nSpinsInGroup[i] + 1);

		// --- 2. Generate Block Patterns ---

		int** BlockFz = new int* [nBl + 1]();
		int** BlockFzSorted = new int* [nBl + 1]();
		for (int i = 1; i <= nBl; i++)
		{
				BlockFz[i] = new int[nSpinGroups + 1]();
				BlockFzSorted[i] = new int[nSpinGroups + 1]();
		}

		int Multiplier = 1;
		for (int i = 1; i <= nSpinGroups; i++)
		{
			Counter = 1;
			int maxFz = nSpinsInGroup[i];

			while (Counter <= nBl)
				for (int k = maxFz; k >= 0; k--)
					for (int m = 1; m <= Multiplier; m++)
					{
						BlockFz[Counter][i] = k;
						Counter++;
					}
			Multiplier *= (maxFz + 1);
		}

		for (int i = 1; i <= nBl; i++) {
			int sum = 0;
			for (int j = 1; j <= nSpinGroups; j++) sum += BlockFz[i][j];
			BlockFz[i][0] = sum;
		}

		Counter = 1;
		for (int i = nSpins; i >= 0; i--)
			for (int j = 1; j <= nBl; j++)
				if (BlockFz[j][0] == i)
				{
					memcpy(BlockFzSorted[Counter], BlockFz[j], (nSpinGroups + 1) * sizeof(int));
					Counter++;
				}

		for (int i = 1; i <= nBl; i++) delete[] BlockFz[i];
		delete[] BlockFz;

		// --- 3. Distribute Basis Functions into Blocks ---

		// Pass 1: Calculate sizes
		BlockSize = new int[nBl + 1]();
		for (int i = 1; i <= nBl; i++)
		{
			int Counter = 0;
			for (BF_t j = 0; j < nBF; j++)
			{
				bool match = true;
				for (int k = 1; (k <= nSpinGroups) && match; k++)
					match = (getbitsum(j & GroupMasks[k]) == BlockFzSorted[i][k]);
				if (match) Counter++;
			}
			BlockSize[i] = Counter;
		}

		// Pass 2: Fill arrays
		BFInBlock = new BF_t * [nBl + 1]();
		for (int i = 1; i <= nBl; i++)
		{
			BFInBlock[i] = new BF_t[BlockSize[i] + 1]();
			BF_t* BFInBlockI = BFInBlock[i];

			int Counter = 1;
			for (BF_t j = 0; j < nBF; j++)
			{
				bool match = true;
				for (int k = 1; (k <= nSpinGroups) && match; k++)
					match = (getbitsum(j & GroupMasks[k]) == BlockFzSorted[i][k]);
				if (match)
				{
					BFInBlockI[Counter] = j;
					Counter++;
				}
			}
		}

		// --- 4. GSL Allocations ---
		Ham = new gsl_matrix * [(size_t)nBl + 1]();
		EVec = new gsl_matrix * [(size_t)nBl + 1]();
		EVal = new gsl_vector * [(size_t)nBl + 1]();
		W = new gsl_eigen_symmv_workspace * [(size_t)nBl + 1]();

		int maxBlockSize = 0;
		for (int i = 1; i <= nBl; i++)
		{
			int bs = BlockSize[i];
			if (bs > maxBlockSize) maxBlockSize = bs;
			Ham[i] = gsl_matrix_alloc(bs, bs);
			EVec[i] = gsl_matrix_alloc(bs, bs);
			EVal[i] = gsl_vector_alloc(bs);
			if (bs > 1) W[i] = gsl_eigen_symmv_alloc(bs);
		}

		Perturbik = new double[maxBlockSize + 1]();

		// --- 5. Selection Rules ---
		nBlockPairs = new int[nSpinGroups + 1]();
		BlockFrom = new int* [nSpinGroups + 1]();
		BlockTo = new int* [nSpinGroups + 1]();

		for (int i = 1; i <= nSpinGroups; i++) {
			nBlockPairs[i] = nBl * nSpinsInGroup[i] / (nSpinsInGroup[i] + 1);
			BlockFrom[i] = new int[nBlockPairs[i] + 1]();
			BlockTo[i] = new int[nBlockPairs[i] + 1]();
		}

		for (int i = 1; i <= nSpinGroups; i++)
		{
			Counter = 1;
			for (int j = 1; j <= nBl - 1; j++)
			{
				int* sortedJ = BlockFzSorted[j];
				for (int k = j + 1; k <= nBl; k++)
				{
					int* sortedK = BlockFzSorted[k];
					bool match = true;

					if (abs(sortedJ[i] - sortedK[i]) != 1) match = false;
					else 
						for (int l = 1; (l <= nSpinGroups) && match; l++)
							if (l != i) match = (sortedJ[l] == sortedK[l]);

					if (match)
					{
						BlockFrom[i][Counter] = j;
						BlockTo[i][Counter] = k;
						Counter++;
					}
				}
			}
		}

		for (int i = 1; i <= nBl; i++) delete[] BlockFzSorted[i];
		delete[] BlockFzSorted;

		// --- 6. Sparse Perturbation Matrices ---
		Perturb = new int*** [nSpinGroups + 1]();
		nFreqs = new int[nSpinGroups + 1]();
		int nFreqsTotal = 0;

		for (int i = 1; i <= nSpinGroups; i++)
		{
			Perturb[i] = new int** [nBlockPairs[i] + 1]();
			//nFreqs[i] = 0;

			for (int j = 1; j <= nBlockPairs[i]; j++)
			{
				int bFrom = BlockFrom[i][j];
				int bTo = BlockTo[i][j];
				int sizeFrom = BlockSize[bFrom];
				int sizeTo = BlockSize[bTo];

				Perturb[i][j] = new int* [sizeFrom + 1]();
				int** PerturbIJ = Perturb[i][j];
				nFreqs[i] += sizeFrom * sizeTo;

				BF_t* BFFrom = BFInBlock[bFrom];
				BF_t* BFTo = BFInBlock[bTo];

				for (int k = 1; k <= sizeFrom; k++)
				{
					PerturbIJ[k] = new int[sizeTo + 1]();
					int* PerturbIJK = PerturbIJ[k];
					BF_t stateK = BFFrom[k];
					for (int l = 1; l <= sizeTo; l++)
						if (getbitsum(stateK ^ BFTo[l]) == 1) PerturbIJK[l] = 1;
				}
			}
			nFreqsTotal += nFreqs[i];
		}

		Freqs = new double* [nSpinGroups + 1]();
		Intens = new double* [nSpinGroups + 1]();

		for (int i = 1; i <= nSpinGroups; i++)
		{
			int nf = nFreqs[i];
			Freqs[i] = new double[nf + 1]();
			Intens[i] = new double[nf + 1]();
		}

		FreqsFiltered = new double[nFreqsTotal + 1]();
		IntensFiltered = new double[nFreqsTotal + 1]();

	}

	void Build(void)
	{

		for (int i = 1; i <= nBl; i++)
		{
			int bs = BlockSize[i];
			double* Hami = Ham[i]->data;
			BF_t* BF_block = BFInBlock[i];

			for (int j = 1; j <= bs; j++)
			{
				int row_offset = (j - 1) * bs - 1;
				BF_t state_j = BF_block[j];

				for (int k = 1; k <= bs; k++)
				{
					double matrix_element = 0;
					BF_t state_k = BF_block[k];

					if (j == k)
						for (int l = 1; l <= nSpins; l++)
						{
							int val_l = (state_j & ((BF_t)1 << (l - 1))) ? 1 : -1;
							matrix_element += val_l * SSParams[Offs[l]];
							for (int m = l + 1; m <= nSpins; m++)
							{
								int val_m = (state_j & ((BF_t)1 << (m - 1))) ? 1 : -1;
								double term = 0.5 * SSParams[JCoups[l][m]];
								if (val_l != val_m) term = -term;
								matrix_element += term;
							}
						}
					else
					{
						BF_t diff = state_j ^ state_k;
						if (getbitsum(diff) == 2)
						{
							int idx[] = { 0, 0 };
							int found = 0;
							for (int bit = 1; bit <= nSpins; bit++)
								if (diff & ((BF_t)1 << (bit - 1)))
								{
									idx[found++] = bit;
									if (found == 2) break;
								}
							matrix_element = SSParams[JCoups[idx[0]][idx[1]]];
						}
					}
					Hami[row_offset + k] = matrix_element;
				}
			}
		}
	};

	void FindEigensystem(void)
	{
		for (int i = 1; i <= nBl; i++)
		{
			if (BlockSize[i] > 1)
			{
				gsl_eigen_symmv(Ham[i], EVal[i], EVec[i], W[i]);
				//gsl_eigen_symmv_sort(EVal[i], EVec[i], GSL_EIGEN_SORT_VAL_ASC);

				double* matrix = EVec[i]->data;
				int bs = BlockSize[i];

				for (int j = 0; j < bs; j++)
				{
					for (int k = j + 1; k < bs; k++)
					{
						double* a = matrix + j * bs + k;
						double* b = matrix + k * bs + j;
						double tmp = *a;
						*a = *b;
						*b = tmp;
					}
				}
			}
			else
			{
				EVal[i]->data[0] = Ham[i]->data[0];
				EVec[i]->data[0] = 1.0;
			}
		}
	}

	void CalcFreqIntens(void)
	{

		// 1. Calculate Frequencies and Intensities for all groups
		double threshold = 0;
		for (int i = 1; i <= nSpinGroups; i++)
		{

			double* FreqsI = Freqs[i];
			double* IntensI = Intens[i];
			int nPairs = nBlockPairs[i];

			int Counter = 1;
			for (int j = 1; j <= nPairs; j++)
			{
				int blFrom = BlockFrom[i][j];
				int blTo = BlockTo[i][j];
				int sizeFrom = BlockSize[blFrom];
				int sizeTo = BlockSize[blTo];

				double* eval_from = EVal[blFrom]->data;
				double* eval_to = EVal[blTo]->data;
				double* vec_from = EVec[blFrom]->data;
				double* vec_to = EVec[blTo]->data;

				int** P_matrix = Perturb[i][j];

				for (int k = 0; k < sizeFrom; k++)
				{
					double* vec_from_row = vec_from + k * sizeFrom;
					memset(Perturbik, 0, (sizeTo + 1) * sizeof(double));

					for (int m = 1; m <= sizeFrom; m++)
					{
						int* P_row_m = P_matrix[m];
						double val_vec = vec_from_row[m - 1];
						for (int l = 1; l <= sizeTo; l++)
							if (P_row_m[l]) Perturbik[l] += val_vec;
					}

					double val_k = eval_from[k];
					for (int l = 0; l < sizeTo; l++)
					{
						double* vec_to_row = vec_to + l * sizeTo;
						double currentIntens = 0;

						for (int m = 1; m <= sizeTo; m++)
							currentIntens += Perturbik[m] * vec_to_row[m - 1];

						currentIntens *= currentIntens;
						IntensI[Counter] = currentIntens;
						FreqsI[Counter] = val_k - eval_to[l];
						if (threshold < currentIntens) threshold = currentIntens;

						Counter++;
					}
				}
			}
		}

		// 2. Filtering with sequential storage
		threshold *= IntensityThreshold;
		int Counter = 0;
		for (int i = 1; i <= nSpinGroups; i++)
		{
			int nf = nFreqs[i];
			double* curIntens = Intens[i];
			double* curFreqs = Freqs[i];

			for (int j = 1; j <= nf; j++)
				if (curIntens[j] > threshold)
				{
					Counter++;
					IntensFiltered[Counter] = curIntens[j];
					FreqsFiltered[Counter] = curFreqs[j] * 0.5;
				}
		}
		nFreqsFiltered = Counter;
	}

	void ComputeFreqIntens(void)
	{
		Build();
		FindEigensystem();
		CalcFreqIntens();
		return;
	}

	void PrintXApproximationStructure(ostream& ostr)
	{
		ostr << "Spectrum computed under weak coupling X-approximation!" << endl;
		int strongGroups = 0;
		int LastManySpinGroup = 0;
		for (int i = 1; i <= nSpinGroups; i++)
			if (nSpinsInGroup[i] > 1) { strongGroups++; LastManySpinGroup = i; }
		if (strongGroups > 1)
		{
			ostr << "Strong coupling accounted for ";
			bool firstGroup = true;
			for (int i = 1; i <= nSpins; ) {
				int group = SpinsGroup[i];
				if (nSpinsInGroup[group] > 1) {
					if (!firstGroup) ostr << ", ";
					firstGroup = false;
					ostr << '(';
					for (int j = i; j < i + nSpinsInGroup[group]; j++) { if (j > i) ostr << ", "; ostr << j; }
					ostr << ')';
					i += nSpinsInGroup[group];
				}
				else i++;
			}
			ostr << " spin groups." << endl;
		}
		else if (strongGroups == 1)
		{
			ostr << "Strong coupling accounted for spins ";
			int tmp = 1;
			while (tmp <= nSpins && SpinsGroup[tmp] != LastManySpinGroup) tmp++;
			for (int i = tmp; i < tmp + nSpinsInGroup[LastManySpinGroup]; i++) { if (i > tmp) ostr << ", "; ostr << i; }
			ostr << '.' << endl;
		}
	}

};

Spectrum* Spec = NULL;
Hamiltonian* Hami = NULL;

double GlobalBadnessLW(double* x)
{
	SSParams[nSSParams] = x[1];
	return Spec->Badness();
}

class OptHamiltonian
{
public:
	//int nParams; => Global nSSParams
	int nVarParams;
	int* VarParamsIndx;
	double* VarParams;
	double* VarParamsUB;
	double* VarParamsLB;
	//int nBroadenings = 0; => Global nBroadenings
	//double* LBs; => Global LBs
	char* InputParameters;
	char* OutputParameters;
	//double* Parameters; => Global SSParams
	double* ParameterErrors;
	double** ParameterCorrelations;
	char** ParNames;
	bool ErrorsComputed;
	bool LWPreOpt;
	double* WS;

	OptHamiltonian(ifstream& istr)
	{

		nVarParams = 0;
		LWPreOpt = false;
		ErrorsComputed = false;

		// All Spin Systems parameters including Magnitude - nSSParams + 1
		ParNames = new char* [(uint_t)nSSParams + 2];
		ParNames[0] = NULL;
		ParameterErrors = new double[(uint_t)nSSParams + 2]; // + Magnitude
		ParameterCorrelations = new double* [(uint_t)nSSParams + 2]; // + Magnitude
		ParameterErrors[0] = 0; ParameterCorrelations[0] = NULL;
		VarParamsIndx = new int[(uint_t)nSSParams + 2]; VarParamsIndx[0] = 0;

		for (int i = 1; i <= nSSParams + 1; i++)
		{
			ParNames[i] = new char[256];
			ParameterErrors[i] = 0;
			ParameterCorrelations[i] = new double[(uint_t)nSSParams + 2];
			for (int j = 0; j <= nSSParams + 1; j++)
				ParameterCorrelations[i][j] = 0;
			strcpy(ParNames[i], "");
			VarParamsIndx[i] = 0;
		}

		VarParams = new double[(uint_t)nSSParams + 1];
		VarParamsUB = new double[(uint_t)nSSParams + 1];
		VarParamsLB = new double[(uint_t)nSSParams + 1];
		for (int i = 0; i <= nSSParams; i++)
		{
			VarParams[i] = 0;
			VarParamsUB[i] = 0;
			VarParamsLB[i] = 0;
		}

		InputParameters = new char[256];
		OutputParameters = new char[256];
		istr.getline(textline, 256);                         // Optimization parameters
		istr >> textline >> InputParameters;                 // Input parameters filename
		istr >> textline >> OutputParameters;                // Output parameters filename
		istr >> textline >> Spec->SpectraTextOutputFilename; // Filename for spectra in ASCII text format
		istr.getline(textline, 256);                         // Rest of SpectraTextOutputFilename line

		// Parsing and loading of Broadenings sequence
		char textline1[256];
		stringstream parse;
		istr.getline(textline, 256);   //  LBs
		parse << textline;
		while (parse) { parse >> textline1; nBroadenings++; }
		nBroadenings -= 2;
		if (nBroadenings < 1) { cout << "Wrong number of broadenings (" << nBroadenings << "), should be at least 1." << endl; exit_; }
		LBs = new double[(uint_t)nBroadenings + 1];
		LBs[0] = 0;
		parse.clear(); parse << textline;
		parse >> textline;

		double MaxLb = 0;
		for (int i = 1; i <= nBroadenings; i++)
		{
			parse >> textline;
			if (!isunsignreal(textline)) { cout << "Wrong value of broadening number " << i << "." << endl; exit_; }
			LBs[i] = atof(textline);
			if(MaxLb < LBs[i]) MaxLb = LBs[i];
		}

		Spec->InitBroadening(MaxLb);

		CheckBroadSequence(cout);

		LoadParameters();

		Spec->CheckSpinOffsets(cout);

		// Magnitude from exp. spectrum
		bool tmp;
		istr >> textline >> tmp; // Magnitude from exp. spectrum
		if (istr.fail()) { cout << "Wrong MagnitudeFromExpSpec flag value, sould be 0 or 1." << endl; exit_; }
		if (tmp) Spec->TheoreticalSpec.Magnitude = Spec->ExperimentalSpec.Magnitude;

	}

	void CheckBroadSequence(ostream& ostr)
	{

		bool printwarn = false;
		for (int i = 2; i <= nBroadenings; i++)
			if (LBs[i] >= LBs[i - 1]) printwarn = true;
		if (printwarn) ostr << "Warning! The sequence of additional broadenings should monotonically decrease!" << endl;

	}

	void LoadParameters(void)
	{

		int ParNum = 0;
		ifstream istr(InputParameters);
		if (!istr) { cout << "File " << InputParameters << " does not exists!" << endl; exit_; }

		for (int i = 1; i <= nSSParams - 1; i++)
		{
			istr >> textline;
			if (!isunsignint(textline)) { cout << "Incorrect number of parameter " << i << "." << endl; exit_; }
			ParNum = atoi(textline);
			if (ParNum != i) { cout << "Incorrect number of parameter " << i << "." << endl; exit_; }
			istr >> ParNames[i];
			istr >> textline;
			if (!isreal(textline)) { cout << "Incorrect value of parameter " << i << "." << endl; exit_; }
			SSParams[i] = atof(textline);
			istr.getline(textline, 256);
		}

		istr >> textline;
		if (!isunsignint(textline)) { cout << "Incorrect number of parameter " << nSSParams << " (linewidth)." << endl; exit_; }
		ParNum = atoi(textline);
		if (ParNum != nSSParams) { cout << "Incorrect number of parameter " << nSSParams << " (linewidth)." << endl; exit_; }
		istr >> ParNames[nSSParams];
		istr >> textline;
		if (!isunsignreal(textline)) { cout << "Incorrect value of parameter " << nSSParams << " (linewidth)." << endl; exit_; }
		SSParams[nSSParams] = atof(textline);
		if (SSParams[nSSParams] == 0) { cout << "Incorrect value of parameter " << nSSParams << " (linewidth)." << endl; exit_; }
		istr.getline(textline, 256);

		istr >> textline;
		if (!isunsignint(textline)) { cout << "Incorrect number of parameter " << nSSParams + 1 << " (spectrum magnitude)." << endl; exit_; }
		ParNum = atoi(textline);
		if (ParNum != nSSParams + 1) { cout << "Incorrect number of parameter " << nSSParams + 1 << " (spectrum magnitude)." << endl; exit_; }
		istr >> ParNames[nSSParams + 1];

		istr >> Spec->TheoreticalSpec.Magnitude;
		if (istr.fail() || Spec->TheoreticalSpec.Magnitude <= 0) { cout << "Incorrect value of parameter " << nSSParams + 1 << " (spectrum magnitude)." << endl; exit_; }

		istr.close();

	}

	void LoadVarParameters(ifstream& istr)
	{

		istr.getline(textline, 256); // Rest of MagnitudeFromExpSpec line
		istr.getline(textline, 256); // Empty line
		if (!isemptyline(textline)) { cout << "Empty line should follow the section with optimization parameters!" << endl; exit_; }
		istr.getline(textline, 256); // List of optimized parameters

		istr >> textline;
		while (istr && nVarParams < nSSParams + 1)
		{
			if (!isunsignint(textline)) { cout << "Wrong index of varied parameter " << nVarParams << "." << endl; exit_; }
			VarParamsIndx[++nVarParams] = atoi(textline);
			if (VarParamsIndx[nVarParams] > nSSParams + 1) { cout << "Wrong index of varied parameter " << nVarParams << "." << endl; exit_; }
			istr >> textline;
		}

		if (nVarParams == 0) { cout << "Nothing to optimize!" << endl; exit_; }

		for (int i = 2; i <= nVarParams; i++)
			if (VarParamsIndx[i] <= VarParamsIndx[i - 1]) { cout << "Wrong index of varied parameter " << i << "." << endl; exit_; }

		if (VarParamsIndx[nVarParams] == nSSParams + 1) { Spec->ScaleOpt = true; nVarParams--; } // By default Spec->ScaleOpt is false.

		if (nVarParams == 0) { cout << "Nothing to optimize!" << endl; exit_; }

		LWPreOpt = VarParamsIndx[nVarParams] == nSSParams && nVarParams > 1;

		for (int i = 1; i <= nVarParams; i++) VarParams[i] = SSParams[VarParamsIndx[i]];

		WS = new double[(uint_t)(5 * nVarParams * (3 * nVarParams + 11) / 2 + 6)];

	}

	double Badness(double* VarPars)
	{

		for (int i = 1; i <= nVarParams; i++)
			SSParams[VarParamsIndx[i]] = VarPars[i];
		Hami->ComputeFreqIntens();
		return Spec->Badness();

	}

	double Badness(void)
	{

		for (int i = 1; i <= nVarParams; i++)
			SSParams[VarParamsIndx[i]] = VarParams[i];
		Hami->ComputeFreqIntens();
		return Spec->Badness();

	}

	void PowellOpt(void);

	void UpdateParamsAfterOpt(void)
	{
		for (int i = 1; i <= nVarParams; i++) SSParams[VarParamsIndx[i]] = VarParams[i];
	}

	void ComputeErrors(void)
	{

		ErrorsComputed = false;

		double** TheorSpecDerivativesOnIntervals = new double* [(uint_t)nSSParams + 2];
		TheorSpecDerivativesOnIntervals[0] = NULL;
		int nPoints = Spec->nPointsRated;
		for (int i = 1; i <= nSSParams + 1; i++)
		{
			TheorSpecDerivativesOnIntervals[i] = new double[(uint_t)nPoints + 1];
			for (int j = 0; j <= nPoints; j++)
				TheorSpecDerivativesOnIntervals[i][j] = 0;
		}

		double* Data = new double[(uint_t)nPoints + 1];
		Data[0] = 0;

		double Badn = 0;
		for (int i = 1; i <= nPoints; i++)
		{
			double tmp = Spec->TheorSpecPointsOnIntervals[i];
			Badn += (tmp - Spec->ExpSpecPointsOnIntervals[i]) * (tmp - Spec->ExpSpecPointsOnIntervals[i]);
			Data[i] = tmp;
			TheorSpecDerivativesOnIntervals[nSSParams + 1][i] = tmp;
		}

		for (int i = 1; i <= nSSParams; i++)
		{
			double step = 1.0e-4;
			double Par = SSParams[i];
			SSParams[i] += step;
			Hami->ComputeFreqIntens();
			Spec->CalcSpecOnIntervals();
			SSParams[i] = Par;
			for (int j = 1; j <= nPoints; j++)
				TheorSpecDerivativesOnIntervals[i][j] = (Spec->TheorSpecPointsOnIntervals[j] - Data[j]) / step;
		}

		for (int i = 1; i <= nPoints; i++)
			Spec->TheorSpecPointsOnIntervals[i] = Data[i];
		delete[] Data;

		gsl_matrix* DTD = gsl_matrix_alloc((uint_t)nSSParams + 1, (uint_t)nSSParams + 1);

		for (int i = 1; i <= nSSParams + 1; i++)
			for (int j = i; j <= nSSParams + 1; j++)
			{
				double tmp = 0;
				for (int k = 1; k <= nPoints; k++)
					tmp += TheorSpecDerivativesOnIntervals[i][k] * TheorSpecDerivativesOnIntervals[j][k];
				DTD->data[(nSSParams + 1) * (i - 1) + j - 1] = tmp;
			}

		for (int i = 0; i < nSSParams + 1; i++)
			for (int j = i + 1; j < nSSParams + 1; j++)
				DTD->data[(nSSParams + 1) * j + i] = DTD->data[(nSSParams + 1) * i + j];

		for (int i = 1; i <= nSSParams + 1; i++)
			delete[] TheorSpecDerivativesOnIntervals[i];
		delete[] TheorSpecDerivativesOnIntervals;

		gsl_matrix* evec = gsl_matrix_alloc((uint_t)nSSParams + 1, (uint_t)nSSParams + 1);
		gsl_vector* eval = gsl_vector_alloc((uint_t)nSSParams + 1);
		gsl_eigen_symmv_workspace* w = gsl_eigen_symmv_alloc((uint_t)nSSParams + 1);

		gsl_eigen_symmv(DTD, eval, evec, w);
		gsl_eigen_symmv_free(w);

		gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_ABS_ASC);

		if (eval->data[0] > 1e-4)
		{
			gsl_matrix_set_all(DTD, 0);

			for (int i = 0; i < nSSParams + 1; i++)
				for (int j = 0; j < nSSParams + 1; j++)
					for (int k = 0; k < nSSParams + 1; k++)
						DTD->data[(nSSParams + 1) * i + j] += evec->data[(nSSParams + 1) * i + k] * evec->data[(nSSParams + 1) * j + k] / eval->data[k];

			for (int i = 0; i < nSSParams + 1; i++)
			{
				ParameterErrors[i + 1] = sqrt(DTD->data[i * (nSSParams + 1) + i] * Badn / (nPoints - nSSParams - 1));
				for (int j = 0; j <= i; j++)
					ParameterCorrelations[i + 1][j + 1] = DTD->data[i * (nSSParams + 1) + j] / sqrt(DTD->data[i * (nSSParams + 1) + i] * DTD->data[j * (nSSParams + 1) + j]);
			}

			ParameterErrors[nSSParams + 1] *= Spec->TheoreticalSpec.Magnitude;

			ErrorsComputed = true;
		}

		gsl_matrix_free(evec);
		gsl_vector_free(eval);
		gsl_matrix_free(DTD);

	}

	void SaveParameters(void)
	{

		int parnamelen = 0;
		ofstream ostr(OutputParameters);

		for (int i = 1; i <= nSSParams + 1; i++)
		{
			int tmp = (int)strlen(ParNames[i]);
			if (parnamelen < tmp) parnamelen = tmp;
		}
		parnamelen += 4;

		if (SSParams[nSSParams] < 0) SSParams[nSSParams] *= -1;

		ostr << fixed;
		for (int i = 1; i <= nSSParams; i++)
		{
			ostr << setw(4) << left << i << setw(parnamelen) << left << ParNames[i];
			ostr << setw(13) << right << SSParams[i];
			if (ErrorsComputed) ostr << setw(4) << ' ' << "+-" << ParameterErrors[i];
			ostr << endl;
		}

		ostr << scientific;
		ostr << setw(4) << left << nSSParams + 1 << setw(parnamelen) << left << ParNames[nSSParams + 1];
		ostr << setw(13) << right << Spec->TheoreticalSpec.Magnitude;
		if (ErrorsComputed) ostr << setw(4) << ' ' << "+-" << ParameterErrors[nSSParams + 1];
		ostr << fixed << endl;

		if (!ErrorsComputed) ostr << "Errors are not calculated due to singularity of normal equations matrix!" << endl;
		ostr << endl;

		ostr << Title << endl << endl;
		Hami->PrintXApproximationStructure(ostr); //ANATOLIA-X addon
		Spec->CheckSpinOffsets(ostr);
		CheckBroadSequence(ostr);

		ostr << "Line Broadening: " << fixed << Spec->LB << endl;
		ostr << "Theoretical spectrum linewidth: " << Spec->LB + SSParams[nSSParams] << endl;
		ostr << "RSS Value: " << scientific << Badness() << endl;
		ostr.precision(2);
		ostr << "R-Factor: " << fixed << Spec->CalcRFactor() << " %" << endl;
		ostr.precision(defaultprecision);
		ostr << "Spectra correlation coefficient: " << Spec->CalcSpectraColleration() << endl;
		ostr << endl;

		// Print Spin System
		ostr.precision(3);
		ostr << "Chemical shifts (ppm):" << endl;
		for (int i = 1; i <= nSpins; i++)
			ostr << setw(10) << right << SSParams[Offs[i]] / Spec->SF;
		ostr << endl;
		ostr.precision(4);
		ostr << "J-coupling constants (Hz):" << endl;
		for (int i = 1; i < nSpins; i++)
		{
			for (int j = 1; j <= nSpins; j++)
				if (j <= i) ostr << setw(10) << ' ';
				else ostr << setw(10) << right << SSParams[JCoups[i][j]];
			ostr << endl;
		}
		ostr.precision(defaultprecision);
		ostr << endl;

		if (ErrorsComputed)
		{
			ostr << "Parameters correlation coefficients:" << endl;
			for (int i = 1; i <= nSSParams + 1; i++)
			{
				ostr << setw(3) << right << i << ' ';
				for (int j = 1; j <= i; j++)
					ostr << setw(10) << right << ParameterCorrelations[i][j];
				ostr << endl;

			}
			ostr << endl;
		}

		print_logo(ostr);
		print_citation(ostr);
		ostr.close();

	}

};

OptHamiltonian* HamOpt = NULL;

double GlobalBadness(double* x)
{
	return HamOpt->Badness((double*)x);
}

void OptHamiltonian :: PowellOpt(void)
{

	for (int i = 1; i <= nVarParams; i++)
	{
		VarParamsUB[i] = VarParams[i] + 100;
		VarParamsLB[i] = VarParams[i] - 100;
	}

	if (LWPreOpt)
	{
		cout << "Spectrum linewidth preoptimization" << endl;

		Hami->ComputeFreqIntens();

		int offset = nVarParams - 1;
		bobyqa(1, 3, GlobalBadnessLW, VarParams + offset, VarParamsLB + offset, VarParamsUB + offset, 10.0, 1e-10, WS);

		VarParamsLB[nVarParams] = VarParams[nVarParams] - 100;
		VarParamsUB[nVarParams] = VarParams[nVarParams] + 100;

		cout << "Main optimization" << endl;
	}

	bobyqa(nVarParams, 2 * nVarParams + 1, GlobalBadness, VarParams, VarParamsLB, VarParamsUB, 10.0, 1e-10, WS);

}

void LoadSpinSystem(ifstream& istr)
{

	// Reading SpinSystem => nSpins, Offs, JCoups
	istr.getline(textline, 256); // Spin System
	istr >> textline >> textline;
	if (!isunsignint(textline)) { cout << "Check the number of spins (Nspins)!" << endl; exit_; }
	nSpins = atoi(textline);
	if (nSpins > 8 * (int)sizeof(unsigned int)) { cout << "Number of spins (Nspins) exceeds the maximum value (" << 8 * sizeof(unsigned int) << ")!" << endl; exit_; }
	istr.getline(textline, 256);
	istr.getline(textline, 256); // Shifts indices

	Offs = new int[(uint_t)nSpins + 1]; Offs[0] = 0;
	JCoups = new int* [(uint_t)nSpins + 1]; JCoups[0] = NULL;

	for (int i = 1; i <= nSpins; i++)
	{
		Offs[i] = 0;
		JCoups[i] = new int[(uint_t)nSpins + 1];
		for (int j = 0; j <= nSpins; j++)
			JCoups[i][j] = 0;
	}

	// Spin offset indices reading
	for (int i = 1; i <= nSpins; i++)
	{
		istr >> textline;
		if (!isunsignint(textline)) { cout << "Wrong input of spin " << i << " chemical shift index." << endl; exit_; }
		Offs[i] = atoi(textline);
	}

	istr.getline(textline, 256); // Rest of the line
	istr.getline(textline, 256); // Coupling Indices

	// Spin offset indices checking
	int tmp = 1;
	if (Offs[1] != 1) { cout << "Spin 1 should have chemical shift index 1." << endl; exit_; }
	for (int i = 2; i <= nSpins; i++)
	{
		if (Offs[i] < tmp) { cout << "Wrong input of spin " << i << " chemical shift index." << endl; exit_; }
		if (Offs[i] > tmp + 1) { cout << "Wrong input of spin " << i << " chemical shift index." << endl; exit_; }
		tmp = Offs[i];
	}
	tmp++;

	// Coupling constant indices reading
	for (int i = 1; i <= nSpins; i++)
		for (int j = i + 1; j <= nSpins; j++)
		{
			istr >> textline;
			if (!isunsignint(textline)) { cout << "Wrong index for J-coupling constant " << i << "," << j << "." << endl; exit_; }
			JCoups[i][j] = atoi(textline);
		}
	istr.getline(textline, 256); // Rest of line
	istr.getline(textline, 256); // Empty line
	if (!isemptyline(textline)) { cout << "Empty line should follow the section with spin system description!" << endl; exit_; }

	// Coupling constant indices checking
	if (JCoups[1][2] != tmp) { cout << "J-coupling constant 1,2 should have index " << tmp << " instead of " << JCoups[1][2] << "." << endl; exit_; }
	for (int i = 1; i <= nSpins; i++)
		for (int j = i + 1; j <= nSpins; j++)
		{
			if (JCoups[i][j] <= Offs[nSpins]) { cout << "Wrong index for J-coupling constant " << i << "," << j << "." << endl; exit_; }
			if (JCoups[i][j] > tmp + 1) { cout << "Wrong index for J-coupling constant " << i << "," << j << "." << endl; exit_; }
			if (JCoups[i][j] > tmp) tmp = JCoups[i][j];
		}

	nSSParams = tmp + 1;

	SSParams = new double[(uint_t)nSSParams + 1];
	for (int i = 0; i <= nSSParams; i++)
		SSParams[i] = 0;

}

int main(int argc, char* argv[])
{

	print_logo(cout);

// If running under Windows Platform
#ifdef _WIN32
	char* env = NULL;
	if (getenv("XWINNMRHOME") != NULL) ExitWithPause = false; // Running under TopSpin
	if ((env = getenv("ExitWithoutPause")) != NULL)           // Checking ExitWithoutPause enviroment variable
		if (strcmp(env, "1") == 0) ExitWithPause = false;
#endif

	if (argc == 2)
	{
		strcpy(textline, argv[1]);
		for (int i = 0; int(textline[i]) != 0; i++) if (textline[i] == '\\') textline[i] = '/';
		int n = (int)strlen(textline);
		if ((textline[n - 1] != '/') && (n < 255)) { textline[n] = '/'; textline[n + 1] = 0; }
		if (chdir(textline) != 0) { cout << "Failed to change working directory!" << endl; exit_; }
	}
	if (argc > 2) { cout << "Wrong command line argument!" << endl; exit_; }

	bool SimMode = true;

	// Parsing input control file and initialize relevant data structures.
	ifstream input("Input_Data.txt");
	if (!input) { cout << "File Input_Data.txt does not exists!" << endl; exit_; }
	input.getline(Title, 256);    // Title
	input.getline(textline, 256); // Empty line
	if (!isemptyline(textline)) { cout << "Empty line should follow the title!" << endl; exit_; }
	input >> textline >> SimMode; // Sim mode
	if (input.fail()) { cout << "Wrong simulation mode value, sould be 0 or 1." << endl; exit_; }
	input.getline(textline, 256); // Rest of sim mode line
	input.getline(textline, 256); // Empty line
	if (!isemptyline(textline)) { cout << "Empty line should follow the SimMode line!" << endl; exit_; }

	LoadSpinSystem(input);

	Spec = new Spectrum(input);

	HamOpt = new OptHamiltonian(input);

	Hami = new Hamiltonian();

	Hami->PrintXApproximationStructure(cout);
	cout << endl;

	if (SimMode)
	{
		Spec->BroadOnIntervals(LBs[1]);
		Spec->ExperimentalSpecWithBroadening.SaveSpecToFile();
		Hami->ComputeFreqIntens();
		Spec->CalcSpecOnIntervals();
		Spec->SaveSpecsOnIntervalsTXT();
		Spec->LB = 0;
		Spec->CalcFullSpectrum();
		Spec->TheoreticalSpec.SaveSpecToFile();
		print_citation(cout);
		input.close();
		exit_;
	}

	HamOpt->LoadVarParameters(input);

	input.close(); // End of INPUT file parsing

	Spec->FineCalc = false;
	for (int i = 1; i <= nBroadenings; i++)
	{

		if (i == nBroadenings) Spec->FineCalc = true;

		cout << "Broadening:  " << LBs[i] << endl;

		Spec->BroadOnIntervals(LBs[i]);

		HamOpt->PowellOpt();

		cout << endl;

	}

	HamOpt->UpdateParamsAfterOpt();
	Hami->ComputeFreqIntens();
	Spec->CalcSpecOnIntervals();
	if (Spec->ScaleOpt)
	{
		Spec->UpdateTheorSpecMagnitude();
		Spec->ScaleOpt = false;
	}

	cout << setprecision(2) << fixed
		<< "R-Factor: " << Spec->CalcRFactor() << " %" << endl
		<< setprecision(defaultprecision);

	Spec->SaveSpecsOnIntervalsTXT();
	Spec->CalcFullSpectrum();
	Spec->TheoreticalSpec.SaveSpecToFile();
	Spec->ExperimentalSpecWithBroadening.SaveSpecToFile();
	HamOpt->ComputeErrors();
	HamOpt->SaveParameters();

	cout << endl;
	print_logo(cout);
	print_citation(cout);

	exit_;

}

/*
* BOBYQA
* The MIT License (MIT) Copyright (c) 2015
*
* Copyright (c) 2009, Mike Powell (FORTRAN version).
* Copyright (c) 2015, Eric Thiebaut (C version).
*
* Permission is hereby granted, free of charge, to any person obtaining a copy
* of this software and associated documentation files (the "Software"), to deal
* in the Software without restriction, including without limitation the rights
* to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the Software is
* furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included in all
* copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
* OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
* SOFTWARE.
*
* Implementation of Mike Powell's BOBYQA algorithm for minimizing a function
* of many variables.  The method is "derivatives free" (only the function
* values are needed) and accounts for bound constraints on the variables.  The
* algorithm is described in:
*
*   M.J.D. Powell, "The BOBYQA Algorithm for Bound Constrained Optimization
*   Without Derivatives."  Technical report, Department of Applied Mathematics
*   and Theoretical Physics, University of Cambridge (2009).
*
* The present code is based on the original FORTRAN version written by Mike
* Powell who kindly provides his code on demand (at mjdp@cam.ac.uk) and has
* been converted to C by E. Thiebaut.
*
*/

#define LOOP(var,num)    for (var = 1; var <= num; ++var)
#define MAX(a,b)         ((a) >= (b) ? (a) : (b))
#define MIN(a,b)         ((a) <= (b) ? (a) : (b))
#define HOW_MANY(a,b)    ((((b) - 1) + (a))/(b))
#define ROUND_UP(a,b)    (HOW_MANY(a,b)*(b))

#define XPT(a1,a2)       xpt[(a2)*npt + a1]
#define BMAT(a1,a2)      bmat[(a2)*ndim + a1]
#define ZMAT(a1,a2)      zmat[(a2)*npt + a1]
#define PTSAUX(a1,a2)    ptsaux[(a2)*2 + a1]

void prelim(long n, long npt,
	double (*objfun)(double*),
	double* x, double* xl, double* xu,
	double rhobeg, double* xbase, double* xpt,
	double* fval, double* gopt, double* hq, double* pq,
	double* bmat, double* zmat, long ndim,
	double* sl, double* su, long* nf, long* kopt)
{
	double half = 0.5;
	double one = 1.0;
	double two = 2.0;
	double zero = 0.0;
	double diff, f, fbeg, recip, rhosq, stepa, stepb, temp;
	long i, i1, ih, ipt, itemp, j, jpt, k, nfm, nfx, np;

	x -= 1;
	xl -= 1;
	xu -= 1;
	xbase -= 1;
	xpt -= 1 + npt;
	fval -= 1;
	gopt -= 1;
	hq -= 1;
	pq -= 1;
	bmat -= 1 + ndim;
	zmat -= 1 + npt;
	sl -= 1;
	su -= 1;
	stepa = zero;
	stepb = zero;
	fbeg = zero;
	ipt = 0;
	jpt = 0;
	rhosq = rhobeg * rhobeg;
	recip = one / rhosq;
	np = n + 1;

	LOOP(j, n) {
		xbase[j] = x[j];
		LOOP(k, npt) {
			XPT(k, j) = zero;
		}
		LOOP(i, ndim) {
			BMAT(i, j) = zero;
		}
	}
	i1 = n * np / 2;
	LOOP(ih, i1) {
		hq[ih] = zero;
	}
	LOOP(k, npt) {
		pq[k] = zero;
		i1 = npt - np;
		LOOP(j, i1) {
			ZMAT(k, j) = zero;
		}
	}

	*nf = 0;
	do {
		nfm = *nf;
		nfx = *nf - n;
		++(*nf);
		if (nfm <= 2 * n) {
			if (nfm >= 1 && nfm <= n) {
				stepa = rhobeg;
				if (su[nfm] == zero) {
					stepa = -stepa;
				}
				XPT(*nf, nfm) = stepa;
			}
			else if (nfm > n) {
				stepa = XPT(*nf - n, nfx);
				stepb = -rhobeg;
				if (sl[nfx] == zero) {
					stepb = two * rhobeg;
					stepb = MIN(stepb, su[nfx]);
				}
				if (su[nfx] == zero) {
					stepb = -two * rhobeg;
					stepb = MAX(stepb, sl[nfx]);
				}
				XPT(*nf, nfx) = stepb;
			}
		}
		else {
			itemp = (nfm - np) / n;
			jpt = nfm - itemp * n - n;
			ipt = jpt + itemp;
			if (ipt > n) {
				itemp = jpt;
				jpt = ipt - n;
				ipt = itemp;
			}
			XPT(*nf, ipt) = XPT(ipt + 1, ipt);
			XPT(*nf, jpt) = XPT(jpt + 1, jpt);
		}

		LOOP(j, n) {
			temp = xbase[j] + XPT(*nf, j);
			temp = MAX(temp, xl[j]);
			x[j] = MIN(temp, xu[j]);
			if (XPT(*nf, j) == sl[j]) {
				x[j] = xl[j];
			}
			if (XPT(*nf, j) == su[j]) {
				x[j] = xu[j];
			}
		}
		f = objfun(x); // objfun(n, &x[1], data);
		fval[*nf] = f;
		if (*nf == 1) {
			fbeg = f;
			*kopt = 1;
		}
		else if (f < fval[*kopt]) {
			*kopt = *nf;
		}

		if (*nf <= 2 * n + 1) {
			if (*nf >= 2 && *nf <= n + 1) {
				gopt[nfm] = (f - fbeg) / stepa;
				if (npt < *nf + n) {
					BMAT(1, nfm) = -one / stepa;
					BMAT(*nf, nfm) = one / stepa;
					BMAT(npt + nfm, nfm) = -half * rhosq;
				}
			}
			else if (*nf >= n + 2) {
				ih = nfx * (nfx + 1) / 2;
				temp = (f - fbeg) / stepb;
				diff = stepb - stepa;
				hq[ih] = two * (temp - gopt[nfx]) / diff;
				gopt[nfx] = (gopt[nfx] * stepb - temp * stepa) / diff;
				if (stepa * stepb < zero) {
					if (f < fval[*nf - n]) {
						fval[*nf] = fval[*nf - n];
						fval[*nf - n] = f;
						if (*kopt == *nf) {
							*kopt = *nf - n;
						}
						XPT(*nf - n, nfx) = stepb;
						XPT(*nf, nfx) = stepa;
					}
				}
				BMAT(1, nfx) = -(stepa + stepb) / (stepa * stepb);
				BMAT(*nf, nfx) = -half / XPT(*nf - n, nfx);
				BMAT(*nf - n, nfx) = -BMAT(1, nfx) - BMAT(*nf, nfx);
				ZMAT(1, nfx) = sqrt(two) / (stepa * stepb);
				ZMAT(*nf, nfx) = sqrt(half) / rhosq;
				ZMAT(*nf - n, nfx) = -ZMAT(1, nfx) - ZMAT(*nf, nfx);
			}
		}
		else {
			ih = ipt * (ipt - 1) / 2 + jpt;
			ZMAT(1, nfx) = recip;
			ZMAT(*nf, nfx) = recip;
			ZMAT(ipt + 1, nfx) = -recip;
			ZMAT(jpt + 1, nfx) = -recip;
			temp = XPT(*nf, ipt) * XPT(*nf, jpt);
			hq[ih] = (fbeg - fval[ipt + 1] - fval[jpt + 1] + f) / temp;
		}
	} while (*nf < npt);
}

void altmov(long n, long npt, double xpt[],
	double* xopt, double* bmat, double* zmat, long ndim,
	double* sl, double* su, long kopt, long knew,
	double adelt, double* xnew, double* xalt, double* alpha,
	double* cauchy, double* glag, double* hcol, double* w)
{

	double half = 0.5;
	double one = 1.0;
	double zero = 0.0;
	double one_plus_sqrt2 = one + M_SQRT2;
	double bigstp, csave, curv, dderiv, diff, distsq, ggfree, gw, ha,
		predsq, presav, scale, slbd, step, stpsav, subd, sumin,
		temp, tempa, tempb, tempd, vlag, wfixsq, wsqsav;
	long i, i1, ibdsav, iflag, ilbd, isbd, iubd, j, k, ksav;

	xpt -= 1 + npt;
	xopt -= 1;
	bmat -= 1 + ndim;
	zmat -= 1 + npt;
	sl -= 1;
	su -= 1;
	xnew -= 1;
	xalt -= 1;
	glag -= 1;
	hcol -= 1;
	w -= 1;
	csave = zero;
	stpsav = zero;
	step = zero;
	ksav = 0;
	ibdsav = 0;
	wsqsav = 0;

	LOOP(k, npt) {
		hcol[k] = zero;
	}
	i1 = npt - n - 1;
	LOOP(j, i1) {
		temp = ZMAT(knew, j);
		LOOP(k, npt) {
			hcol[k] += temp * ZMAT(k, j);
		}
	}
	*alpha = hcol[knew];
	ha = half * (*alpha);

	LOOP(i, n) {
		glag[i] = BMAT(knew, i);
	}
	LOOP(k, npt) {
		temp = zero;
		LOOP(j, n) {
			temp += XPT(k, j) * xopt[j];
		}
		temp = hcol[k] * temp;
		LOOP(i, n) {
			glag[i] += temp * XPT(k, i);
		}
	}

	presav = zero;
	LOOP(k, npt) {
		if (k == kopt) {
			continue;
		}
		dderiv = zero;
		distsq = zero;
		LOOP(i, n) {
			temp = XPT(k, i) - xopt[i];
			dderiv += glag[i] * temp;
			distsq += temp * temp;
		}
		subd = adelt / sqrt(distsq);
		slbd = -subd;
		ilbd = 0;
		iubd = 0;
		sumin = MIN(one, subd);

		LOOP(i, n) {
			temp = XPT(k, i) - xopt[i];
			if (temp > zero) {
				if (slbd * temp < sl[i] - xopt[i]) {
					slbd = (sl[i] - xopt[i]) / temp;
					ilbd = -i;
				}
				if (subd * temp > su[i] - xopt[i]) {
					subd = (su[i] - xopt[i]) / temp;
					subd = MAX(subd, sumin);
					iubd = i;
				}
			}
			else if (temp < zero) {
				if (slbd * temp > su[i] - xopt[i]) {
					slbd = (su[i] - xopt[i]) / temp;
					ilbd = i;
				}
				if (subd * temp < sl[i] - xopt[i]) {
					subd = (sl[i] - xopt[i]) / temp;
					subd = MAX(subd, sumin);
					iubd = -i;
				}
			}
		}

		if (k == knew) {
			diff = dderiv - one;
			step = slbd;
			vlag = slbd * (dderiv - slbd * diff);
			isbd = ilbd;
			temp = subd * (dderiv - subd * diff);
			if (abs(temp) > abs(vlag)) {
				step = subd;
				vlag = temp;
				isbd = iubd;
			}
			tempd = half * dderiv;
			tempa = tempd - diff * slbd;
			tempb = tempd - diff * subd;
			if (tempa * tempb < zero) {
				temp = tempd * tempd / diff;
				if (abs(temp) > abs(vlag)) {
					step = tempd / diff;
					vlag = temp;
					isbd = 0;
				}
			}
		}
		else {
			step = slbd;
			vlag = slbd * (one - slbd);
			isbd = ilbd;
			temp = subd * (one - subd);
			if (abs(temp) > abs(vlag)) {
				step = subd;
				vlag = temp;
				isbd = iubd;
			}
			if (subd > half) {
				if (abs(vlag) < 0.25) {
					step = half;
					vlag = 0.25;
					isbd = 0;
				}
			}
			vlag *= dderiv;
		}

		temp = step * (one - step) * distsq;
		predsq = vlag * vlag * (vlag * vlag + ha * temp * temp);
		if (predsq > presav) {
			presav = predsq;
			ksav = k;
			stpsav = step;
			ibdsav = isbd;
		}
	}

	LOOP(i, n) {
		temp = xopt[i] + stpsav * (XPT(ksav, i) - xopt[i]);
		temp = MIN(temp, su[i]);
		xnew[i] = MAX(temp, sl[i]);
	}
	if (ibdsav < 0) {
		xnew[-ibdsav] = sl[-ibdsav];
	}
	if (ibdsav > 0) {
		xnew[ibdsav] = su[ibdsav];
	}

	bigstp = adelt + adelt;
	iflag = 0;
	for (;;) {
		wfixsq = zero;
		ggfree = zero;
		LOOP(i, n) {
			w[i] = zero;
			tempa = xopt[i] - sl[i];
			tempa = MIN(tempa, glag[i]);
			tempb = xopt[i] - su[i];
			tempb = MAX(tempb, glag[i]);
			if (tempa > zero || tempb < zero) {
				w[i] = bigstp;
				ggfree += glag[i] * glag[i];
			}
		}
		if (ggfree == zero) {
			*cauchy = zero;
			return;
		}

		do {
			temp = adelt * adelt - wfixsq;
			if (temp <= zero) break;
			wsqsav = wfixsq;
			step = sqrt(temp / ggfree);
			ggfree = zero;
			LOOP(i, n) {
				if (w[i] == bigstp) {
					temp = xopt[i] - step * glag[i];
					if (temp <= sl[i]) {
						w[i] = sl[i] - xopt[i];
						wfixsq += w[i] * w[i];
					}
					else if (temp >= su[i]) {
						w[i] = su[i] - xopt[i];
						wfixsq += w[i] * w[i];
					}
					else {
						ggfree += glag[i] * glag[i];
					}
				}
			}
		} while (wfixsq > wsqsav && ggfree > zero);

		gw = zero;
		LOOP(i, n) {
			if (w[i] == bigstp) {
				w[i] = -step * glag[i];
				temp = xopt[i] + w[i];
				temp = MIN(temp, su[i]);
				xalt[i] = MAX(temp, sl[i]);
			}
			else if (w[i] == zero) {
				xalt[i] = xopt[i];
			}
			else if (glag[i] > zero) {
				xalt[i] = sl[i];
			}
			else {
				xalt[i] = su[i];
			}
			gw += glag[i] * w[i];
		}

		curv = zero;
		LOOP(k, npt) {
			temp = zero;
			LOOP(j, n) {
				temp += XPT(k, j) * w[j];
			}
			curv += hcol[k] * temp * temp;
		}
		if (iflag == 1) {
			curv = -curv;
		}
		if (curv > -gw && curv < -one_plus_sqrt2 * gw) {
			scale = -gw / curv;
			LOOP(i, n) {
				temp = xopt[i] + scale * w[i];
				temp = MIN(temp, su[i]);
				xalt[i] = MAX(temp, sl[i]);
			}
			temp = half * gw * scale;
			*cauchy = temp * temp;
		}
		else {
			temp = gw + half * curv;
			*cauchy = temp * temp;
		}

		if (iflag != 0) {
			break;
		}
		LOOP(i, n) {
			glag[i] = -glag[i];
			w[n + i] = xalt[i];
		}
		csave = *cauchy;
		iflag = 1;
	}
	if (csave > *cauchy) {
		LOOP(i, n) {
			xalt[i] = w[n + i];
		}
		*cauchy = csave;
	}
}

void trsbox(long n, long npt, double* xpt,
	double* xopt, double* gopt, double* hq, double* pq,
	double* sl, double* su, double delta, double* xnew,
	double* d, double* gnew, double* xbdi, double* s,
	double* hs, double* hred, double* dsq, double* crvmin)
{

	double half = 0.5;
	double one = 1.0;
	double onemin = -1.0;
	double zero = 0.0;
	double angbd, angt, beta, blen, cth, delsq, dhd, dhs, dredg, dredsq, ds,
		ggsav, gredsq, qred, rdnext, rdprev, redmax, rednew, redsav, resid,
		sdec, shs, sredg, ssq, stepsq, sth, stplen, temp, tempa, tempb, xsav, xsum;
	long i, iact, ih, isav, itcsav, iterc, itermax, iu, j, k, nact;

	xpt -= 1 + npt;
	xopt -= 1;
	gopt -= 1;
	hq -= 1;
	pq -= 1;
	sl -= 1;
	su -= 1;
	xnew -= 1;
	d -= 1;
	gnew -= 1;
	xbdi -= 1;
	s -= 1;
	hs -= 1;
	hred -= 1;

	angbd = zero;
	dredg = zero;
	dredsq = zero;
	ggsav = zero;
	gredsq = zero;
	rdnext = zero;
	sredg = zero;
	xsav = zero;
	iact = 0;
	itcsav = 0;
	itermax = 0;

	iterc = 0;
	nact = 0;
	LOOP(i, n) {
		xbdi[i] = zero;
		if (xopt[i] <= sl[i]) {
			if (gopt[i] >= zero) {
				xbdi[i] = onemin;
			}
		}
		else if (xopt[i] >= su[i]) {
			if (gopt[i] <= zero) {
				xbdi[i] = one;
			}
		}
		if (xbdi[i] != zero) {
			++nact;
		}
		d[i] = zero;
		gnew[i] = gopt[i];
	}
	delsq = delta * delta;
	qred = zero;
	*crvmin = onemin;

L20:
	beta = zero;
L30:
	stepsq = zero;
	LOOP(i, n) {
		if (xbdi[i] != zero) {
			s[i] = zero;
		}
		else if (beta == zero) {
			s[i] = -gnew[i];
		}
		else {
			s[i] = beta * s[i] - gnew[i];
		}
		stepsq += s[i] * s[i];
	}
	if (stepsq == zero) {
		goto L190;
	}
	if (beta == zero) {
		gredsq = stepsq;
		itermax = iterc + n - nact;
	}
	if (gredsq * delsq <= qred * 1e-4 * qred) {
		goto L190;
	}
	goto L210;
L50:
	resid = delsq;
	ds = zero;
	shs = zero;
	LOOP(i, n) {
		if (xbdi[i] == zero) {
			resid -= d[i] * d[i];
			ds += s[i] * d[i];
			shs += s[i] * hs[i];
		}
	}
	if (resid <= zero) {
		goto L90;
	}
	temp = sqrt(stepsq * resid + ds * ds);
	if (ds < zero) {
		blen = (temp - ds) / stepsq;
	}
	else {
		blen = resid / (temp + ds);
	}
	if (shs > zero) {
		stplen = gredsq / shs;
		stplen = MIN(blen, stplen);
	}
	else {
		stplen = blen;
	}
	iact = 0;
	LOOP(i, n) {
		if (s[i] != zero) {
			xsum = xopt[i] + d[i];
			if (s[i] > zero) {
				temp = (su[i] - xsum) / s[i];
			}
			else {
				temp = (sl[i] - xsum) / s[i];
			}
			if (temp < stplen) {
				stplen = temp;
				iact = i;
			}
		}
	}
	sdec = zero;
	if (stplen > zero) {
		++iterc;
		temp = shs / stepsq;
		if (iact == 0 && temp > zero) {
			*crvmin = MIN(*crvmin, temp);
			if (*crvmin == onemin) {
				*crvmin = temp;
			}
		}
		ggsav = gredsq;
		gredsq = zero;
		LOOP(i, n) {
			gnew[i] += stplen * hs[i];
			if (xbdi[i] == zero) {
				gredsq += gnew[i] * gnew[i];
			}
			d[i] += stplen * s[i];
		}
		sdec = stplen * (ggsav - half * stplen * shs);
		sdec = MAX(sdec, zero);
		qred += sdec;
	}
	if (iact > 0) {
		++nact;
		xbdi[iact] = one;
		if (s[iact] < zero) {
			xbdi[iact] = onemin;
		}
		delsq -= d[iact] * d[iact];
		if (delsq <= zero) {
			goto L90;
		}
		goto L20;
	}
	if (stplen < blen) {
		if (iterc == itermax) {
			goto L190;
		}
		if (sdec <= qred * 0.01) {
			goto L190;
		}
		beta = gredsq / ggsav;
		goto L30;
	}
L90:
	*crvmin = zero;
L100:
	if (nact >= n - 1) {
		goto L190;
	}
	dredsq = zero;
	dredg = zero;
	gredsq = zero;
	LOOP(i, n) {
		if (xbdi[i] == zero) {
			dredsq += d[i] * d[i];
			dredg += d[i] * gnew[i];
			gredsq += gnew[i] * gnew[i];
			s[i] = d[i];
		}
		else {
			s[i] = zero;
		}
	}
	itcsav = iterc;
	goto L210;
L120:
	++iterc;
	temp = gredsq * dredsq - dredg * dredg;
	if (temp <= qred * 1e-4 * qred) {
		goto L190;
	}
	temp = sqrt(temp);
	LOOP(i, n) {
		if (xbdi[i] == zero) {
			s[i] = (dredg * d[i] - dredsq * gnew[i]) / temp;
		}
		else {
			s[i] = zero;
		}
	}
	sredg = -temp;
	angbd = one;
	iact = 0;
	LOOP(i, n) {
		if (xbdi[i] == zero) {
			tempa = xopt[i] + d[i] - sl[i];
			tempb = su[i] - xopt[i] - d[i];
			if (tempa <= zero) {
				++nact;
				xbdi[i] = onemin;
				goto L100;
			}
			else if (tempb <= zero) {
				++nact;
				xbdi[i] = one;
				goto L100;
			}
			ssq = d[i] * d[i] + s[i] * s[i];
			temp = xopt[i] - sl[i];
			temp = ssq - temp * temp;
			if (temp > zero) {
				temp = sqrt(temp) - s[i];
				if (angbd * temp > tempa) {
					angbd = tempa / temp;
					iact = i;
					xsav = onemin;
				}
			}
			temp = su[i] - xopt[i];
			temp = ssq - temp * temp;
			if (temp > zero) {
				temp = sqrt(temp) + s[i];
				if (angbd * temp > tempb) {
					angbd = tempb / temp;
					iact = i;
					xsav = one;
				}
			}
		}
	}
	goto L210;
L150:
	shs = zero;
	dhs = zero;
	dhd = zero;
	LOOP(i, n) {
		if (xbdi[i] == zero) {
			shs += s[i] * hs[i];
			dhs += d[i] * hs[i];
			dhd += d[i] * hred[i];
		}
	}
	redmax = zero;
	isav = 0;
	redsav = zero;
	iu = (long)(angbd * 17.0 + 3.1);
	LOOP(i, iu) {
		angt = angbd * (double)i / (double)iu;
		sth = (angt + angt) / (one + angt * angt);
		temp = shs + angt * (angt * dhd - dhs - dhs);
		rednew = sth * (angt * dredg - sredg - half * sth * temp);
		if (rednew > redmax) {
			redmax = rednew;
			isav = i;
			rdprev = redsav;
		}
		else if (i == isav + 1) {
			rdnext = rednew;
		}
		redsav = rednew;
	}
	if (isav == 0) {
		goto L190;
	}
	if (isav < iu) {
		temp = (rdnext - rdprev) / (redmax + redmax - rdprev - rdnext);
		angt = angbd * ((double)isav + half * temp) / (double)iu;
	}
	cth = (one - angt * angt) / (one + angt * angt);
	sth = (angt + angt) / (one + angt * angt);
	temp = shs + angt * (angt * dhd - dhs - dhs);
	sdec = sth * (angt * dredg - sredg - half * sth * temp);
	if (sdec <= zero) {
		goto L190;
	}
	dredg = zero;
	gredsq = zero;
	LOOP(i, n) {
		gnew[i] = gnew[i] + (cth - one) * hred[i] + sth * hs[i];
		if (xbdi[i] == zero) {
			d[i] = cth * d[i] + sth * s[i];
			dredg += d[i] * gnew[i];
			gredsq += gnew[i] * gnew[i];
		}
		hred[i] = cth * hred[i] + sth * hs[i];
	}
	qred += sdec;
	if (iact > 0 && isav == iu) {
		++nact;
		xbdi[iact] = xsav;
		goto L100;
	}
	if (sdec > qred * 0.01) {
		goto L120;
	}
L190:
	*dsq = zero;
	LOOP(i, n) {
		temp = xopt[i] + d[i];
		temp = MIN(temp, su[i]);
		xnew[i] = MAX(temp, sl[i]);
		if (xbdi[i] == onemin) {
			xnew[i] = sl[i];
		}
		if (xbdi[i] == one) {
			xnew[i] = su[i];
		}
		d[i] = xnew[i] - xopt[i];
		*dsq += d[i] * d[i];
	}
	return;
L210:
	ih = 0;
	LOOP(j, n) {
		hs[j] = zero;
		LOOP(i, j) {
			++ih;
			if (i < j) {
				hs[j] += hq[ih] * s[i];
			}
			hs[i] += hq[ih] * s[j];
		}
	}
	LOOP(k, npt) {
		if (pq[k] != zero) {
			temp = zero;
			LOOP(j, n) {
				temp += XPT(k, j) * s[j];
			}
			temp *= pq[k];
			LOOP(i, n) {
				hs[i] += temp * XPT(k, i);
			}
		}
	}
	if (*crvmin != zero) {
		goto L50;
	}
	if (iterc > itcsav) {
		goto L150;
	}
	LOOP(i, n) {
		hred[i] = hs[i];
	}
	goto L120;
}

void update(long n, long npt, double* bmat,
	double* zmat, long ndim, double* vlag, double beta,
	double denom, long knew, double* w)
{

	double one = 1.0;
	double zero = 0.0;
	double alpha, tau, temp, tempa, tempb, ztest;
	long i, j, jp, k, nptm;
	zmat -= 1 + npt;
	bmat -= 1 + ndim;
	vlag -= 1;
	w -= 1;
	nptm = npt - n - 1;
	ztest = zero;
	LOOP(k, npt) {
		LOOP(j, nptm) {
			temp = abs(ZMAT(k, j));
			ztest = MAX(ztest, temp);
		}
	}
	ztest *= 1e-20;
	for (j = 2; j <= nptm; ++j) {
		if (abs(ZMAT(knew, j)) > ztest) {
			tempa = ZMAT(knew, 1);
			tempb = ZMAT(knew, j);
			temp = sqrt(tempa * tempa + tempb * tempb);
			tempa /= temp;
			tempb /= temp;
			LOOP(i, npt) {
				temp = tempa * ZMAT(i, 1) + tempb * ZMAT(i, j);
				ZMAT(i, j) = tempa * ZMAT(i, j) - tempb * ZMAT(i, 1);
				ZMAT(i, 1) = temp;
			}
		}
		ZMAT(knew, j) = zero;
	}
	LOOP(i, npt) {
		w[i] = ZMAT(knew, 1) * ZMAT(i, 1);
	}
	alpha = w[knew];
	tau = vlag[knew];
	vlag[knew] -= one;
	temp = sqrt(denom);
	tempb = ZMAT(knew, 1) / temp;
	tempa = tau / temp;
	LOOP(i, npt) {
		ZMAT(i, 1) = tempa * ZMAT(i, 1) - tempb * vlag[i];
	}
	LOOP(j, n) {
		jp = npt + j;
		w[jp] = BMAT(knew, j);
		tempa = (alpha * vlag[jp] - tau * w[jp]) / denom;
		tempb = (-beta * w[jp] - tau * vlag[jp]) / denom;
		LOOP(i, jp) {
			BMAT(i, j) = BMAT(i, j) + tempa * vlag[i] + tempb * w[i];
			if (i > npt) {
				BMAT(jp, i - npt) = BMAT(i, j);
			}
		}
	}
}

void rescue(long n, long npt,
	double (*objfun)(double*),
	double* xl, double* xu,
	double* xbase, double* xpt, double* fval, double* xopt,
	double* gopt, double* hq, double* pq, double* bmat, double* zmat,
	long ndim, double* sl, double* su, long* nf,
	double delta, long* kopt, double* vlag, double* ptsaux,
	double* ptsid, double* w)
{

	double half = 0.5;
	double one = 1.0;
	double zero = 0.0;
	double beta, bsum, den, denom, diff, distsq, dsqmin, f, fbase, hdiag,
		sfrac, sum, sumpq, temp, vlmxsq, vquad, winc, xp, xq;
	long i, ih, ihp, ihq, ip, iq, iw, j, jp, jpn, k, knew, kold, kpt,
		np, nptm, nrem;

	xp = 0;
	xq = 0;
	ihp = 0;
	ihq = 0;
	zmat -= 1 + npt;
	xpt -= 1 + npt;
	xl -= 1;
	xu -= 1;
	xbase -= 1;
	fval -= 1;
	xopt -= 1;
	gopt -= 1;
	hq -= 1;
	pq -= 1;
	bmat -= 1 + ndim;
	sl -= 1;
	su -= 1;
	vlag -= 1;
	ptsaux -= 3;
	ptsid -= 1;
	w -= 1;
	beta = zero;
	denom = zero;
	np = n + 1;
	sfrac = half / (double)np;
	nptm = npt - np;
	sumpq = zero;
	winc = zero;
	LOOP(k, npt) {
		distsq = zero;
		LOOP(j, n) {
			XPT(k, j) = XPT(k, j) - xopt[j];
			distsq += XPT(k, j) * XPT(k, j);
		}
		sumpq += pq[k];
		w[ndim + k] = distsq;
		winc = MAX(winc, distsq);
		LOOP(j, nptm) {
			ZMAT(k, j) = zero;
		}
	}
	ih = 0;
	LOOP(j, n) {
		w[j] = half * sumpq * xopt[j];
		LOOP(k, npt) {
			w[j] += pq[k] * XPT(k, j);
		}
		LOOP(i, j) {
			++ih;
			hq[ih] = hq[ih] + w[i] * xopt[j] + w[j] * xopt[i];
		}
	}
	LOOP(j, n) {
		xbase[j] += xopt[j];
		sl[j] -= xopt[j];
		su[j] -= xopt[j];
		xopt[j] = zero;
		PTSAUX(1, j) = MIN(delta, su[j]);
		PTSAUX(2, j) = MAX(-delta, sl[j]);
		if (PTSAUX(1, j) + PTSAUX(2, j) < zero) {
			temp = PTSAUX(1, j);
			PTSAUX(1, j) = PTSAUX(2, j);
			PTSAUX(2, j) = temp;
		}
		if (abs(PTSAUX(2, j)) < half * abs(PTSAUX(1, j))) {
			PTSAUX(2, j) = half * PTSAUX(1, j);
		}
		LOOP(i, ndim) {
			BMAT(i, j) = zero;
		}
	}
	fbase = fval[*kopt];
	ptsid[1] = sfrac;
	LOOP(j, n) {
		jp = j + 1;
		jpn = jp + n;
		ptsid[jp] = (double)j + sfrac;
		if (jpn <= npt) {
			ptsid[jpn] = (double)j / (double)np + sfrac;
			temp = one / (PTSAUX(1, j) - PTSAUX(2, j));
			BMAT(jp, j) = -temp + one / PTSAUX(1, j);
			BMAT(jpn, j) = temp + one / PTSAUX(2, j);
			BMAT(1, j) = -BMAT(jp, j) - BMAT(jpn, j);
			ZMAT(1, j) = sqrt(2.0) / abs(PTSAUX(1, j) * PTSAUX(2, j));
			ZMAT(jp, j) = ZMAT(1, j) * PTSAUX(2, j) * temp;
			ZMAT(jpn, j) = -ZMAT(1, j) * PTSAUX(1, j) * temp;
		}
		else {
			BMAT(1, j) = -one / PTSAUX(1, j);
			BMAT(jp, j) = one / PTSAUX(1, j);
			BMAT(j + npt, j) = -half * (PTSAUX(1, j) * PTSAUX(1, j));
		}
	}
	if (npt >= n + np) {
		for (k = 2 * np; k <= npt; ++k) {
			iw = (long)(((double)(k - np) - half) / (double)n);
			ip = k - np - iw * n;
			iq = ip + iw;
			if (iq > n) {
				iq -= n;
			}
			ptsid[k] = (double)ip + (double)iq / (double)np + sfrac;
			temp = one / (PTSAUX(1, ip) * PTSAUX(1, iq));
			ZMAT(1, k - np) = temp;
			ZMAT(ip + 1, k - np) = -temp;
			ZMAT(iq + 1, k - np) = -temp;
			ZMAT(k, k - np) = temp;
		}
	}
	nrem = npt;
	kold = 1;
	knew = *kopt;
L80:
	LOOP(j, n) {
		temp = BMAT(kold, j);
		BMAT(kold, j) = BMAT(knew, j);
		BMAT(knew, j) = temp;
	}
	LOOP(j, nptm) {
		temp = ZMAT(kold, j);
		ZMAT(kold, j) = ZMAT(knew, j);
		ZMAT(knew, j) = temp;
	}
	ptsid[kold] = ptsid[knew];
	ptsid[knew] = zero;
	w[ndim + knew] = zero;
	--nrem;
	if (knew != *kopt) {
		temp = vlag[kold];
		vlag[kold] = vlag[knew];
		vlag[knew] = temp;
		update(n, npt, &BMAT(1, 1), &ZMAT(1, 1), ndim, &vlag[1],
			beta, denom, knew, &w[1]);
		if (nrem == 0) {
			goto L350;
		}
		LOOP(k, npt) {
			w[ndim + k] = abs(w[ndim + k]);
		}
	}
L120:
	dsqmin = zero;
	LOOP(k, npt) {
		if (w[ndim + k] > zero) {
			if (dsqmin == zero || w[ndim + k] < dsqmin) {
				knew = k;
				dsqmin = w[ndim + k];
			}
		}
	}
	if (dsqmin == zero) {
		goto L260;
	}
	LOOP(j, n) {
		w[npt + j] = XPT(knew, j);
	}
	LOOP(k, npt) {
		sum = zero;
		if (k == *kopt) {
		}
		else if (ptsid[k] == zero) {
			LOOP(j, n) {
				sum += w[npt + j] * XPT(k, j);
			}
		}
		else {
			ip = (long)ptsid[k];
			if (ip > 0) {
				sum = w[npt + ip] * PTSAUX(1, ip);
			}
			iq = (long)((double)np * ptsid[k] - (double)(ip * np));
			if (iq > 0) {
				iw = 1;
				if (ip == 0) {
					iw = 2;
				}
				sum += w[npt + iq] * PTSAUX(iw, iq);
			}
		}
		w[k] = half * sum * sum;
	}
	LOOP(k, npt) {
		sum = zero;
		LOOP(j, n) {
			sum += BMAT(k, j) * w[npt + j];
		}
		vlag[k] = sum;
	}
	beta = zero;
	LOOP(j, nptm) {
		sum = zero;
		LOOP(k, npt) {
			sum += ZMAT(k, j) * w[k];
		}
		beta -= sum * sum;
		LOOP(k, npt) {
			vlag[k] += sum * ZMAT(k, j);
		}
	}
	bsum = zero;
	distsq = zero;
	LOOP(j, n) {
		sum = zero;
		LOOP(k, npt) {
			sum += BMAT(k, j) * w[k];
		}
		jp = j + npt;
		bsum += sum * w[jp];
		for (ip = npt + 1; ip <= ndim; ++ip) {
			sum += BMAT(ip, j) * w[ip];
		}
		bsum += sum * w[jp];
		vlag[jp] = sum;
		distsq += XPT(knew, j) * XPT(knew, j);
	}
	beta = half * distsq * distsq + beta - bsum;
	vlag[*kopt] += one;
	denom = zero;
	vlmxsq = zero;
	LOOP(k, npt) {
		if (ptsid[k] != zero) {
			hdiag = zero;
			LOOP(j, nptm) {
				hdiag += ZMAT(k, j) * ZMAT(k, j);
			}
			den = beta * hdiag + vlag[k] * vlag[k];
			if (den > denom) {
				kold = k;
				denom = den;
			}
		}
		temp = vlag[k] * vlag[k];
		vlmxsq = MAX(vlmxsq, temp);
	}
	if (denom <= vlmxsq * 0.01) {
		w[ndim + knew] = -w[ndim + knew] - winc;
		goto L120;
	}
	goto L80;
L260:
	LOOP(kpt, npt) {
		if (ptsid[kpt] == zero) {
			continue;
		}
		ih = 0;
		LOOP(j, n) {
			w[j] = XPT(kpt, j);
			XPT(kpt, j) = zero;
			temp = pq[kpt] * w[j];
			LOOP(i, j) {
				++ih;
				hq[ih] += temp * w[i];
			}
		}
		pq[kpt] = zero;
		ip = (long)ptsid[kpt];
		iq = (long)((double)np * ptsid[kpt] - (double)(ip * np))
			;
		if (ip > 0) {
			xp = PTSAUX(1, ip);
			XPT(kpt, ip) = xp;
		}
		if (iq > 0) {
			xq = PTSAUX(1, iq);
			if (ip == 0) {
				xq = PTSAUX(2, iq);
			}
			XPT(kpt, iq) = xq;
		}

		vquad = fbase;
		if (ip > 0) {
			ihp = (ip + ip * ip) / 2;
			vquad += xp * (gopt[ip] + half * xp * hq[ihp]);
		}
		if (iq > 0) {
			ihq = (iq + iq * iq) / 2;
			vquad += xq * (gopt[iq] + half * xq * hq[ihq]);
			if (ip > 0) {
				iw = MAX(ihp, ihq) - (ip >= iq ? ip - iq : iq - ip);
				vquad += xp * xq * hq[iw];
			}
		}
		LOOP(k, npt) {
			temp = zero;
			if (ip > 0) {
				temp += xp * XPT(k, ip);
			}
			if (iq > 0) {
				temp += xq * XPT(k, iq);
			}
			vquad += half * pq[k] * temp * temp;
		}

		LOOP(i, n) {
			temp = xbase[i] + XPT(kpt, i);
			temp = MAX(temp, xl[i]);
			w[i] = MIN(temp, xu[i]);
			if (XPT(kpt, i) == sl[i]) {
				w[i] = xl[i];
			}
			if (XPT(kpt, i) == su[i]) {
				w[i] = xu[i];
			}
		}
		++(*nf);
		f = objfun(w); // objfun(n, &w[1], data)
		fval[kpt] = f;
		if (f < fval[*kopt]) {
			*kopt = kpt;
		}
		diff = f - vquad;

		LOOP(i, n) {
			gopt[i] += diff * BMAT(kpt, i);
		}
		LOOP(k, npt) {
			sum = zero;
			LOOP(j, nptm) {
				sum += ZMAT(k, j) * ZMAT(kpt, j);
			}
			temp = diff * sum;
			if (ptsid[k] == zero) {
				pq[k] += temp;
			}
			else {
				ip = (long)ptsid[k];
				iq = (long)((double)np * ptsid[k] - (double)(ip * np));
				ihq = (iq * iq + iq) / 2;
				if (ip == 0) {
					hq[ihq] += temp * (PTSAUX(2, iq) * PTSAUX(2, iq));
				}
				else {
					ihp = (ip * ip + ip) / 2;
					hq[ihp] += temp * (PTSAUX(1, ip) * PTSAUX(1, ip));
					if (iq > 0) {
						hq[ihq] += temp * (PTSAUX(1, iq) * PTSAUX(1, iq));
						iw = MAX(ihp, ihq) - (ip >= iq ? ip - iq : iq - ip);
						hq[iw] += temp * PTSAUX(1, ip) * PTSAUX(1, iq);
					}
				}
			}
		}
		ptsid[kpt] = zero;
	}
L350:
	return;
}

int bobyqb(long n, long npt,
	double (*objfun)(double*),
	double* x, double* xl, double* xu,
	double rhobeg, double rhoend,
	double* xbase, double* xpt, double* fval,
	double* xopt, double* gopt, double* hq,
	double* pq, double* bmat, double* zmat,
	long ndim, double* sl, double* su, double* xnew,
	double* xalt, double* d, double* vlag, double* w)
{
	double half = 0.5;
	double one = 1.0;
	double ten = 10.0;
	double tenth = 0.1;
	double two = 2.0;
	double zero = 0.0;
	double adelt, alpha, bdtest, bdtol, beta, biglsq, bsum, cauchy, crvmin,
		curv, delsq, delta, den, denom, densav, diff, diffa, diffb, diffc,
		dist, distsq, dnorm, dsq, dx, errbig, f, fopt, fracsq, frhosq, fsave,
		gisq, gqsq, hdiag, pqold, ratio, rho, scaden, sum, suma, sumb, sumpq,
		sumw, sumz, temp, tempa, tempb, vquad, xoptsq;
	long i, ih, ip, itest, j, jj, jp, k, kbase, knew, kopt, ksav, nf,
		nfsav, nh, np, nptm, nresc, ntrits;
	int status = BOBYQA_SUCCESS;
	const char* reason = NULL;
	int iter = 0;
	beta = 0;
	x -= 1;
	xl -= 1;
	xu -= 1;
	xbase -= 1;
	xpt -= 1 + npt;
	fval -= 1;
	xopt -= 1;
	gopt -= 1;
	hq -= 1;
	pq -= 1;
	bmat -= 1 + ndim;
	zmat -= 1 + npt;
	sl -= 1;
	su -= 1;
	xnew -= 1;
	xalt -= 1;
	d -= 1;
	vlag -= 1;
	w -= 1;
	adelt = zero;
	alpha = zero;
	cauchy = zero;
	denom = zero;
	diff = zero;
	diffc = zero;
	f = zero;
	knew = 0;
	np = n + 1;
	nptm = npt - np;
	nh = n * np / 2;

	prelim(n, npt, objfun, &x[1], &xl[1], &xu[1], rhobeg, &xbase[1],
		&XPT(1, 1), &fval[1], &gopt[1], &hq[1], &pq[1], &BMAT(1, 1),
		&ZMAT(1, 1), ndim, &sl[1], &su[1], &nf, &kopt);
	xoptsq = zero;
	LOOP(i, n) {
		xopt[i] = XPT(kopt, i);
		xoptsq += xopt[i] * xopt[i];
	}
	fsave = fval[1];
	kbase = 1;

	rho = rhobeg;
	delta = rho;
	nresc = nf;
	ntrits = 0;
	diffa = zero;
	diffb = zero;
	itest = 0;
	nfsav = nf;

L20:
	if (kopt != kbase) {
		ih = 0;
		LOOP(j, n) {
			LOOP(i, j) {
				++ih;
				if (i < j) {
					gopt[j] += hq[ih] * xopt[i];
				}
				gopt[i] += hq[ih] * xopt[j];
			}
		}
		if (nf > npt) {
			LOOP(k, npt) {
				temp = zero;
				LOOP(j, n) {
					temp += XPT(k, j) * xopt[j];
				}
				temp = pq[k] * temp;
				LOOP(i, n) {
					gopt[i] += temp * XPT(k, i);
				}
			}
		}
	}

L60:
	trsbox(n, npt, &XPT(1, 1), &xopt[1], &gopt[1], &hq[1], &pq[1],
		&sl[1], &su[1], delta, &xnew[1], &d[1], &w[1], &w[np],
		&w[np + n], &w[np + 2 * n], &w[np + n * 3], &dsq, &crvmin);
	iter++;
	dnorm = sqrt(dsq);
	dnorm = MIN(dnorm, delta);
	if (dnorm < half * rho) {
		ntrits = -1;
		tempa = ten * rho;
		distsq = tempa * tempa;
		if (nf <= nfsav + 2) {
			goto L650;
		}
		errbig = MAX(diffa, diffb);
		errbig = MAX(errbig, diffc);
		frhosq = rho * 0.125 * rho;
		if (crvmin > zero && errbig > frhosq * crvmin) {
			goto L650;
		}
		bdtol = errbig / rho;
		LOOP(j, n) {
			bdtest = bdtol;
			if (xnew[j] == sl[j]) {
				bdtest = w[j];
			}
			if (xnew[j] == su[j]) {
				bdtest = -w[j];
			}
			if (bdtest < bdtol) {
				curv = hq[(j + j * j) / 2];
				LOOP(k, npt) {
					curv += pq[k] * (XPT(k, j) * XPT(k, j));
				}
				bdtest += half * curv * rho;
				if (bdtest < bdtol) {
					goto L650;
				}
			}
		}
		goto L680;
	}
	++ntrits;
L90:
	if (dsq <= xoptsq * 0.001) {
		fracsq = xoptsq * 0.25;
		sumpq = zero;
		LOOP(k, npt) {
			sumpq += pq[k];
			sum = -half * xoptsq;
			LOOP(i, n) {
				sum += XPT(k, i) * xopt[i];
			}
			w[npt + k] = sum;
			temp = fracsq - half * sum;
			LOOP(i, n) {
				w[i] = BMAT(k, i);
				vlag[i] = sum * XPT(k, i) + temp * xopt[i];
				ip = npt + i;
				LOOP(j, i) {
					BMAT(ip, j) = BMAT(ip, j) + w[i] * vlag[j] + vlag[i] * w[j];
				}
			}
		}

		LOOP(jj, nptm) {
			sumz = zero;
			sumw = zero;
			LOOP(k, npt) {
				sumz += ZMAT(k, jj);
				vlag[k] = w[npt + k] * ZMAT(k, jj);
				sumw += vlag[k];
			}
			LOOP(j, n) {
				sum = (fracsq * sumz - half * sumw) * xopt[j];
				LOOP(k, npt) {
					sum += vlag[k] * XPT(k, j);
				}
				w[j] = sum;
				LOOP(k, npt) {
					BMAT(k, j) = BMAT(k, j) + sum * ZMAT(k, jj);
				}
			}
			LOOP(i, n) {
				ip = i + npt;
				temp = w[i];
				LOOP(j, i) {
					BMAT(ip, j) = BMAT(ip, j) + temp * w[j];
				}
			}
		}

		ih = 0;
		LOOP(j, n) {
			w[j] = -half * sumpq * xopt[j];
			LOOP(k, npt) {
				w[j] += pq[k] * XPT(k, j);
				XPT(k, j) = XPT(k, j) - xopt[j];
			}
			LOOP(i, j) {
				++ih;
				hq[ih] = hq[ih] + w[i] * xopt[j] + xopt[i] * w[j];
				BMAT(npt + i, j) = BMAT(npt + j, i);
			}
		}
		LOOP(i, n) {
			xbase[i] += xopt[i];
			xnew[i] -= xopt[i];
			sl[i] -= xopt[i];
			su[i] -= xopt[i];
			xopt[i] = zero;
		}
		xoptsq = zero;
	}
	if (ntrits == 0) {
		goto L210;
	}
	goto L230;

L190:
	nfsav = nf;
	kbase = kopt;
	rescue(n, npt, objfun, &xl[1], &xu[1], &xbase[1],
		&XPT(1, 1), &fval[1], &xopt[1], &gopt[1], &hq[1],
		&pq[1], &BMAT(1, 1), &ZMAT(1, 1), ndim, &sl[1], &su[1],
		&nf, delta, &kopt, &vlag[1], &w[1], &w[n + np], &w[ndim + np]);

	xoptsq = zero;
	if (kopt != kbase) {
		LOOP(i, n) {
			xopt[i] = XPT(kopt, i);
			xoptsq += xopt[i] * xopt[i];
		}
	}
	nresc = nf;
	if (nfsav < nf) {
		nfsav = nf;
		goto L20;
	}
	if (ntrits > 0) {
		goto L60;
	}

L210:
	altmov(n, npt, &XPT(1, 1), &xopt[1], &BMAT(1, 1), &ZMAT(1, 1), ndim,
		&sl[1], &su[1], kopt, knew, adelt, &xnew[1], &xalt[1],
		&alpha, &cauchy, &w[1], &w[np], &w[ndim + 1]);
	LOOP(i, n) {
		d[i] = xnew[i] - xopt[i];
	}

L230:
	LOOP(k, npt) {
		suma = zero;
		sumb = zero;
		sum = zero;
		LOOP(j, n) {
			suma += XPT(k, j) * d[j];
			sumb += XPT(k, j) * xopt[j];
			sum += BMAT(k, j) * d[j];
		}
		w[k] = suma * (half * suma + sumb);
		vlag[k] = sum;
		w[npt + k] = suma;
	}
	beta = zero;
	LOOP(jj, nptm) {
		sum = zero;
		LOOP(k, npt) {
			sum += ZMAT(k, jj) * w[k];
		}
		beta -= sum * sum;
		LOOP(k, npt) {
			vlag[k] += sum * ZMAT(k, jj);
		}
	}
	dsq = zero;
	bsum = zero;
	dx = zero;
	LOOP(j, n) {
		dsq += d[j] * d[j];
		sum = zero;
		LOOP(k, npt) {
			sum += w[k] * BMAT(k, j);
		}
		bsum += sum * d[j];
		jp = npt + j;
		LOOP(i, n) {
			sum += BMAT(jp, i) * d[i];
		}
		vlag[jp] = sum;
		bsum += sum * d[j];
		dx += d[j] * xopt[j];
	}
	beta = dx * dx + dsq * (xoptsq + dx + dx + half * dsq) + beta - bsum;
	vlag[kopt] += one;
	if (ntrits == 0) {
		denom = vlag[knew] * vlag[knew] + alpha * beta;
		if (denom < cauchy && cauchy > zero) {
			LOOP(i, n) {
				xnew[i] = xalt[i];
				d[i] = xnew[i] - xopt[i];
			}
			cauchy = zero;
			goto L230;
		}
		if (denom <= half * (vlag[knew] * vlag[knew])) {
			if (nf > nresc) {
				goto L190;
			}
			goto cancellation_of_denominator;
		}

	}
	else {
		delsq = delta * delta;
		scaden = zero;
		biglsq = zero;
		knew = 0;
		LOOP(k, npt) {
			if (k == kopt) continue;
			hdiag = zero;
			LOOP(jj, nptm) {
				hdiag += ZMAT(k, jj) * ZMAT(k, jj);
			}
			den = beta * hdiag + vlag[k] * vlag[k];
			distsq = zero;
			LOOP(j, n) {
				tempa = XPT(k, j) - xopt[j];
				distsq += tempa * tempa;
			}
			temp = distsq / delsq;
			temp = temp * temp;
			temp = MAX(one, temp);
			if (temp * den > scaden) {
				scaden = temp * den;
				knew = k;
				denom = den;
			}
			temp *= vlag[k] * vlag[k];
			biglsq = MAX(biglsq, temp);
		}
		if (scaden <= half * biglsq) {
			if (nf > nresc) {
				goto L190;
			}
			goto cancellation_of_denominator;
		}
	}
L360:
	LOOP(i, n) {
		tempa = xbase[i] + xnew[i];
		tempa = MAX(tempa, xl[i]);
		x[i] = MIN(tempa, xu[i]);
		if (xnew[i] == sl[i]) {
			x[i] = xl[i];
		}
		if (xnew[i] == su[i]) {
			x[i] = xu[i];
		}
	}
	++nf;
	f = objfun(x); //objfun(n, &x[1], data);
	if (ntrits == -1) {
		fsave = f;
		goto done;
	}
	fopt = fval[kopt];
	vquad = zero;
	ih = 0;
	LOOP(j, n) {
		vquad += d[j] * gopt[j];
		LOOP(i, j) {
			++ih;
			temp = d[i] * d[j];
			if (i == j) {
				temp = half * temp;
			}
			vquad += hq[ih] * temp;
		}
	}
	LOOP(k, npt) {
		vquad += half * pq[k] * (w[npt + k] * w[npt + k]);
	}
	diff = f - fopt - vquad;
	diffc = diffb;
	diffb = diffa;
	diffa = abs(diff);
	if (dnorm > rho) {
		nfsav = nf;
	}
	if (ntrits > 0) {
		if (vquad >= zero) {
			goto step_failed;
		}
		ratio = (f - fopt) / vquad;
		if (ratio <= tenth) {
			delta *= half;
			delta = MIN(delta, dnorm);
		}
		else if (ratio <= 0.7) {
			delta *= half;
			delta = MAX(delta, dnorm);
		}
		else {
			tempa = dnorm + dnorm;
			delta *= half;
			delta = MAX(delta, tempa);
		}
		if (delta <= rho * 1.5) {
			delta = rho;
		}
		if (f < fopt) {
			ksav = knew;
			densav = denom;
			delsq = delta * delta;
			scaden = zero;
			biglsq = zero;
			knew = 0;
			LOOP(k, npt) {
				hdiag = zero;
				LOOP(jj, nptm) {
					hdiag += ZMAT(k, jj) * ZMAT(k, jj);
				}
				den = beta * hdiag + vlag[k] * vlag[k];
				distsq = zero;
				LOOP(j, n) {
					temp = XPT(k, j) - xnew[j];
					distsq += temp * temp;
				}
				temp = distsq / delsq;
				temp = temp * temp;
				temp = MAX(one, temp);
				if (temp * den > scaden) {
					scaden = temp * den;
					knew = k;
					denom = den;
				}
				temp *= (vlag[k] * vlag[k]);
				biglsq = MAX(biglsq, temp);
			}
			if (scaden <= half * biglsq) {
				knew = ksav;
				denom = densav;
			}
		}
	}
	update(n, npt, &BMAT(1, 1), &ZMAT(1, 1), ndim, &vlag[1],
		beta, denom, knew, &w[1]);
	ih = 0;
	pqold = pq[knew];
	pq[knew] = zero;
	LOOP(i, n) {
		temp = pqold * XPT(knew, i);
		LOOP(j, i) {
			++ih;
			hq[ih] += temp * XPT(knew, j);
		}
	}
	LOOP(jj, nptm) {
		temp = diff * ZMAT(knew, jj);
		LOOP(k, npt) {
			pq[k] += temp * ZMAT(k, jj);
		}
	}
	fval[knew] = f;
	LOOP(i, n) {
		XPT(knew, i) = xnew[i];
		w[i] = BMAT(knew, i);
	}
	LOOP(k, npt) {
		suma = zero;
		LOOP(jj, nptm) {
			suma += ZMAT(knew, jj) * ZMAT(k, jj);
		}
		sumb = zero;
		LOOP(j, n) {
			sumb += XPT(k, j) * xopt[j];
		}
		temp = suma * sumb;
		LOOP(i, n) {
			w[i] += temp * XPT(k, i);
		}
	}
	LOOP(i, n) {
		gopt[i] += diff * w[i];
	}
	if (f < fopt) {
		kopt = knew;
		xoptsq = zero;
		ih = 0;
		LOOP(j, n) {
			xopt[j] = xnew[j];
			xoptsq += xopt[j] * xopt[j];
			LOOP(i, j) {
				++ih;
				if (i < j) {
					gopt[j] += hq[ih] * d[i];
				}
				gopt[i] += hq[ih] * d[j];
			}
		}
		LOOP(k, npt) {
			temp = zero;
			LOOP(j, n) {
				temp += XPT(k, j) * d[j];
			}
			temp = pq[k] * temp;
			LOOP(i, n) {
				gopt[i] += temp * XPT(k, i);
			}
		}
	}
	if (ntrits > 0) {
		LOOP(k, npt) {
			vlag[k] = fval[k] - fval[kopt];
			w[k] = zero;
		}
		LOOP(j, nptm) {
			sum = zero;
			LOOP(k, npt) {
				sum += ZMAT(k, j) * vlag[k];
			}
			LOOP(k, npt) {
				w[k] += sum * ZMAT(k, j);
			}
		}
		LOOP(k, npt) {
			sum = zero;
			LOOP(j, n) {
				sum += XPT(k, j) * xopt[j];
			}
			w[k + npt] = w[k];
			w[k] = sum * w[k];
		}
		gqsq = zero;
		gisq = zero;
		LOOP(i, n) {
			sum = zero;
			LOOP(k, npt) {
				sum = sum + BMAT(k, i) * vlag[k] + XPT(k, i) * w[k];
			}
			if (xopt[i] == sl[i]) {
				tempa = MIN(zero, gopt[i]);
				gqsq += tempa * tempa;
				tempa = MIN(zero, sum);
				gisq += tempa * tempa;
			}
			else if (xopt[i] == su[i]) {
				tempa = MAX(zero, gopt[i]);
				gqsq += tempa * tempa;
				tempa = MAX(zero, sum);
				gisq += tempa * tempa;
			}
			else {
				gqsq += gopt[i] * gopt[i];
				gisq += sum * sum;
			}
			vlag[npt + i] = sum;
		}
		++itest;
		if (gqsq < ten * gisq) {
			itest = 0;
		}
		if (itest >= 3) {
			long i1 = MAX(npt, nh);
			LOOP(i, i1) {
				if (i <= n) {
					gopt[i] = vlag[npt + i];
				}
				if (i <= npt) {
					pq[i] = w[npt + i];
				}
				if (i <= nh) {
					hq[i] = zero;
				}
				itest = 0;
			}
		}
	}
	if (ntrits == 0) {
		goto L60;
	}
	if (f <= fopt + tenth * vquad) {
		goto L60;
	}
	tempa = two * delta;
	tempb = ten * rho;
	tempa = tempa * tempa;
	tempb = tempb * tempb;
	distsq = MAX(tempa, tempb);
L650:
	knew = 0;
	LOOP(k, npt) {
		sum = zero;
		LOOP(j, n) {
			tempa = XPT(k, j) - xopt[j];
			sum += tempa * tempa;
		}
		if (sum > distsq) {
			knew = k;
			distsq = sum;
		}
	}
	if (knew > 0) {
		dist = sqrt(distsq);
		if (ntrits == -1) {
			tempa = tenth * delta;
			tempb = half * dist;
			delta = MIN(tempa, tempb);
			if (delta <= rho * 1.5) {
				delta = rho;
			}
		}
		ntrits = 0;
		adelt = tenth * dist;
		adelt = MIN(adelt, delta);
		adelt = MAX(adelt, rho);
		dsq = adelt * adelt;
		goto L90;
	}
	if (ntrits == -1) {
		goto L680;
	}
	if (ratio > zero) {
		goto L60;
	}
	if (MAX(delta, dnorm) > rho) {
		goto L60;
	}
L680:
	if (rho > rhoend) {
		delta = half * rho;
		ratio = rho / rhoend;
		if (ratio <= 16.0) {
			rho = rhoend;
		}
		else if (ratio <= 250.0) {
			rho = sqrt(ratio) * rhoend;
		}
		else {
			rho = tenth * rho;
		}
		delta = MAX(delta, rho);

		cout << "Iteration " << iter << "\tRHO " << setprecision(1)
			<< scientific << rho << setprecision(defaultprecision)
			<< "\t" << (double)fval[kopt] << endl;

		ntrits = 0;
		nfsav = nf;
		goto L60;
	}
	if (ntrits == -1) {
		goto L360;
	}
done:
	if (fval[kopt] <= fsave) {
		LOOP(i, n) {
			tempa = xbase[i] + xopt[i];
			tempa = MAX(tempa, xl[i]);
			x[i] = MIN(tempa, xu[i]);
			if (xopt[i] == sl[i]) {
				x[i] = xl[i];
			}
			if (xopt[i] == su[i]) {
				x[i] = xu[i];
			}
		}
		f = fval[kopt];
	}

	cout << "Iteration " << iter << "\tRHO " << setprecision(1)
		<< scientific << rho << setprecision(defaultprecision)
		<< "\t" << (double)fval[kopt] << endl;// defaultfloat << endl;
	cout.unsetf(ios_base::floatfield);

	if (status == BOBYQA_SUCCESS) {
		xbase[1] = f;
	}
	return status;

cancellation_of_denominator:
	reason = "of much cancellation in a denominator";
	status = BOBYQA_ROUNDING_ERRORS;
	goto error;

step_failed:
	reason = "a trust region step has failed to reduce Q";
	status = BOBYQA_STEP_FAILED;
	goto error;

error:
	cout << "\n\tReturn from BOBYQA because " << reason << '.' << endl;
	goto done;

}

#undef ZMAT
#undef BMAT
#undef XPT
#undef PTSAUX

int bobyqa(long n, long npt,
	double (*objfun)(double*), double* x,
	double* xl, double* xu,
	double rhobeg, double rhoend,
	double* w)
{

	double zero = 0.0;
	double temp, tempa, tempb;
	long ibmat, id, ifv, igo, ihq, ipq, isl, isu, ivl, iw, ixa, ixb, ixn,
		ixo, ixp, izmat, j, jsl, jsu, ndim, np;

	w -= 1;
	//xu -= 1;
	//xl -= 1;
	//x -= 1;
	np = n + 1;
	if (npt < n + 2 || npt >(n + 2) * np / 2) {
		cout << "\n\tReturn from BOBYQA because NPT is not in the required interval." << endl;
		return BOBYQA_BAD_NPT;
	}

	ndim = npt + n;
	ixb = 1;
	ixp = ixb + n;
	ifv = ixp + n * npt;
	ixo = ifv + npt;
	igo = ixo + n;
	ihq = igo + n;
	ipq = ihq + n * np / 2;
	ibmat = ipq + npt;
	izmat = ibmat + ndim * n;
	isl = izmat + npt * (npt - np);
	isu = isl + n;
	ixn = isu + n;
	ixa = ixn + n;
	id = ixa + n;
	ivl = id + n;
	iw = ivl + ndim;

	LOOP(j, n) {
		temp = xu[j] - xl[j];
		if (temp < rhobeg + rhobeg) {
			cout << "\n\tReturn from BOBYQA because one of the differences XU(I)-XL(I) is less than 2*RHOBEG." << endl;
			return BOBYQA_TOO_CLOSE;
		}
		jsl = isl + j - 1;
		jsu = jsl + n;
		w[jsl] = xl[j] - x[j];
		w[jsu] = xu[j] - x[j];
		if (w[jsl] >= -rhobeg) {
			if (w[jsl] >= zero) {
				x[j] = xl[j];
				w[jsl] = zero;
				w[jsu] = temp;
			}
			else {
				x[j] = xl[j] + rhobeg;
				w[jsl] = -rhobeg;
				temp = xu[j] - x[j];
				w[jsu] = MAX(temp, rhobeg);
			}
		}
		else if (w[jsu] <= rhobeg) {
			if (w[jsu] <= zero) {
				x[j] = xu[j];
				w[jsl] = -temp;
				w[jsu] = zero;
			}
			else {
				x[j] = xu[j] - rhobeg;
				tempa = xl[j] - x[j];
				tempb = -rhobeg;
				w[jsl] = MIN(tempa, tempb);
				w[jsu] = rhobeg;
			}
		}
	}

	return bobyqb(n, npt, objfun, &x[1], &xl[1], &xu[1],
		rhobeg, rhoend,
		&w[ixb], &w[ixp], &w[ifv], &w[ixo], &w[igo],
		&w[ihq], &w[ipq], &w[ibmat], &w[izmat],
		ndim, &w[isl], &w[isu], &w[ixn], &w[ixa],
		&w[id], &w[ivl], &w[iw]);
}
