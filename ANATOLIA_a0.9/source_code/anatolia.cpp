/*****************************************************/
/***                ANATOLIA V1.0                  ***/
/***   Program for the total lineshape analysis    ***/
/***               of NMR spectra                  ***/
/*****************************************************/
/*** (C) 2017 Dmitry Cheshkov, Dmitry Sinitsyn,    ***/
/***             Kirill Sheberstov                 ***/
/*****************************************************/
/***             dcheshkov@gmail.com               ***/
/***          http://anatolia.nmrclub.ru           ***/
/***     https://github.com/dcheshkov/ANATOLIA     ***/
/*****************************************************/
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_multimin.h>

using namespace std;

#define IntensityThreshold 0.005

char Title[256];

#ifdef _WIN32
#define exit_ exit_with_pause()
bool ExitWithPause = true;
void exit_with_pause(void)
{

	if (!ExitWithPause) exit(0);
	char textline[8];
	cout << "Press any key..." << endl;
	cin.getline(textline, 8);
	exit(0);

};
#else
#define exit_ exit(0)
#endif

int getbitsum(unsigned int number)
{

	int bits = 0;
	for (int i = 0; i <= 8 * sizeof(unsigned int)-1; i++)
	if ((number & (1 << i)) > 0) { bits += 1; };
	return(bits);

};

int getbit(unsigned int number, int position)
{

	int result = 0;
	result = (number & (1 << (position - 1))) >> (position - 1);
	return(result);

};

bool isunsignint(const char* text)
{

	bool res = int(text[0]) != 0;
	for (int i = 0; int(text[i]) != 0; i++)
		res = res && isdigit(text[i]);
	return(res);

};

bool isunsignreal(const char* text)
{

	bool res = int(text[0]) != 0;
	int point = 0;
	for (int i = 0; int(text[i]) != 0; i++)
	{
		if (text[i] == '.') point++;
		res = res && (point < 2) && (isdigit(text[i]) || text[i] == '.');
	};
	return(res);

};

bool isreal(const char* text)
{

	int j = (text[0] == '-') || (text[0] == '+') ? 1 : 0;
	return isunsignreal(text + j);

};

class Hamiltonian
{
public:
	int nSpins;
	unsigned int nFunc;
	int nBlocks;
	double* Parameters;
	int nParams;
	int* Offs;
	int** Jcoup;
	unsigned int** bFunc;
	double*** Ham;
	int* BlockSize;
	int MaxBlockSize;
	double** EigenVal;
	bool*** Perturbation;
	int nFreqs;
	int nFreqsFiltered;
	double* Freqs;
	double* Intens;
	double* FreqsFiltered;
	double* IntensFiltered;

	Hamiltonian(ifstream& istr)
	{
		char textline[256];
		istr >> textline >> textline;
		if (!isunsignint(textline))  { cout << "Check the number of spins (Nspins)!" << endl; exit_; };
		nSpins = stoi(textline);
		istr.getline(textline, 256);
		istr.getline(textline, 256); // Shifts indices

		nParams = 0;
		Offs = new int[nSpins + 1]; Offs[0] = 0;
		Jcoup = new int*[nSpins + 1]; Jcoup[0] = NULL;
		for (int i = 1; i <= nSpins; i++)
		{
			Offs[i] = 0;
			Jcoup[i] = new int[nSpins + 1];
			for (int j = 0; j <= nSpins; j++)
				Jcoup[i][j] = 0;
		};

		// Spin offset indices reading
		for (int i = 1; i <= nSpins; i++)
		{
			istr >> textline;
			if (!isunsignint(textline)) { cout << "Wrong input of spin " << i << " chemical shift index." << endl; exit_; };
			Offs[i] = stoi(textline);
			if (nParams < Offs[i]) nParams = Offs[i];
		};

		istr.getline(textline, 256); // Rest of the line
		istr.getline(textline, 256); // Coupling Indices

		// Spin offset indices checking
		if (Offs[1] != 1) { cout << "Spin 1 should have chemical shift index 1." << endl; exit_; };
		int tmp = 1;
		for (int i = 2; i <= nSpins; i++)
		{
			if (Offs[i] < tmp) { cout << "Wrong input of spin " << i << " chemical shift index." << endl; exit_; };
			if (Offs[i] > tmp + 1) { cout << "Wrong input of spin " << i << " chemical shift index." << endl; exit_; };
			tmp = Offs[i];
		};
		tmp++;

		// Coupling constant indices reading
		for (int i = 1; i <= nSpins; i++)
		for (int j = i + 1; j <= nSpins; j++)
		{
			istr >> textline;
			if (!isunsignint(textline)) { cout << "Wrong index for J-coupling constant " << i << "," << j << "." << endl; exit_; };
			Jcoup[i][j] = stoi(textline);
			if (nParams < Jcoup[i][j]) nParams = Jcoup[i][j];
		};
		istr.getline(textline, 256); //  Rest of line

		// Coupling constant indices checking
		if (Jcoup[1][2] != tmp) { cout << "J-coupling constant 1,2 should have index " << tmp << " instead of " << Jcoup[1][2] << "." << endl; exit_; };
		for (int i = 1; i <= nSpins; i++)
		for (int j = i + 1; j <= nSpins; j++)
		{
			if (Jcoup[i][j] <= Offs[nSpins]) { cout << "Wrong index for J-coupling constant " << i << "," << j << "." << endl; exit_; };
			if (Jcoup[i][j] > tmp + 1) { cout << "Wrong index for J-coupling constant " << i << "," << j << "." << endl; exit_; };
			if (Jcoup[i][j] > tmp) tmp = Jcoup[i][j];
		};

		Parameters = new double[nParams + 1];
		for (int i = 0; i <= nParams; i++)
			Parameters[i] = 0;

		nFunc = 1 << nSpins;
		nBlocks = nSpins + 1;

		BlockSize = new int[nBlocks + 1];
		for (int i = 0; i <= nBlocks; i++)
			BlockSize[i] = 0;
		for (unsigned int i = 0; i <= nFunc - 1; i++)
			BlockSize[getbitsum(i) + 1]++;

		MaxBlockSize = 0;
		for (int i = 1; i <= nBlocks; i++)
		if (MaxBlockSize < BlockSize[i])
			MaxBlockSize = BlockSize[i];

		nFreqs = 0;
		nFreqsFiltered = 0;
		for (int i = 1; i <= nBlocks - 1; i++)
			nFreqs += BlockSize[i] * BlockSize[i + 1];
		Freqs = new double[nFreqs + 1];
		Intens = new double[nFreqs + 1];
		FreqsFiltered = new double[nFreqs + 1];
		IntensFiltered = new double[nFreqs + 1];

		for (int i = 0; i <= nFreqs; i++)
		{
			Freqs[i] = 0;
			Intens[i] = 0;
			FreqsFiltered[i] = 0;
			IntensFiltered[i] = 0;
		};

		bFunc = new unsigned int*[nBlocks + 1]; bFunc[0] = NULL;
		for (int i = 1; i <= nBlocks; i++)
		{
			int tmp = 1;
			bFunc[i] = new unsigned int[BlockSize[i] + 1];
			bFunc[i][0] = 0;
			for (unsigned int j = 0; j <= nFunc - 1; j++)
			if (getbitsum(j) + 1 == i) { bFunc[i][tmp] = j; tmp++; };
		};

		Ham = new double**[nBlocks + 1];
		Ham[0] = NULL;
		for (int i = 1; i <= nBlocks; i++)
		{
			Ham[i] = new double*[BlockSize[i] + 1];
			Ham[i][0] = NULL;
			for (int j = 1; j <= BlockSize[i]; j++)
			{
				Ham[i][j] = new double[BlockSize[i] + 1];
				for (int k = 0; k <= BlockSize[i]; k++)
					Ham[i][j][k] = 0;
			};
		};

		EigenVal = new double*[nBlocks + 1];
		EigenVal[0] = NULL;
		for (int i = 1; i <= nBlocks; i++)
		{
			EigenVal[i] = new double[BlockSize[i] + 1];
			for (int j = 0; j <= BlockSize[i]; j++)
				EigenVal[i][j] = 0;
		};

		Perturbation = new bool**[nBlocks];
		Perturbation[0] = NULL;
		for (int i = 1; i <= nBlocks - 1; i++)
		{
			Perturbation[i] = new bool*[BlockSize[i + 1] + 1];
			Perturbation[i][0] = NULL;
			for (int j = 1; j <= BlockSize[i + 1]; j++)
			{
				Perturbation[i][j] = new bool[BlockSize[i] + 1];
				Perturbation[i][j][0] = false;
			};
		};

		for (int i = 1; i <= nBlocks - 1; i++)
		for (int j = 1; j <= BlockSize[i + 1]; j++)
		{
			bool* Perturbij = Perturbation[i][j];
			unsigned int* bFi = bFunc[i];
			unsigned int bFippj = bFunc[i + 1][j];
			for (int k = 1; k <= BlockSize[i]; k++)
				Perturbij[k] = getbitsum(bFippj^bFi[k]) == 1;
		};

	};

	void Build(void)
	{
		int j_location[3] = { 0, 0, 0 };
		unsigned int offdiag = 0;
		int* sign = new int[nSpins + 1];
		int signk = 0;
		double element = 0;
		unsigned int bFij = 0;
		unsigned int* bFi = NULL;
		int bs = 0;
		double** Hami = NULL;

		for (int i = 1; i <= nBlocks; i++)
		{
			bs = BlockSize[i];
			Hami = Ham[i];
			bFi = bFunc[i];

			for (int j = 1; j <= bs; j++)
			{
				bFij = bFi[j];
				element = 0;

				for (int k = 1; k <= nSpins; k++)
					sign[k] = (getbit(bFij, nSpins - k + 1) == 1) ? -1 : 1;

				for (int k = 1; k <= nSpins; k++)
				{
					signk = sign[k];
					element += Parameters[Offs[k]] * signk;
					for (int kk = k + 1; kk <= nSpins; kk++) element += 0.5*signk*sign[kk] * Parameters[Jcoup[k][kk]];
				};

				Hami[j][j] = element;

				for (int k = 1; k <= j - 1; k++)
				{
					element = 0;
					offdiag = bFij^bFi[k];
					if (getbitsum(offdiag) == 2)
					{
						int tmp = 1;
						for (int kk = 1; kk <= nSpins; kk++)
						if (getbit(offdiag, nSpins - kk + 1) == 1)
							j_location[tmp++] = kk;
						element = Parameters[Jcoup[j_location[1]][j_location[2]]];
					};
					Hami[j][k] = element;
					//Hami[k][j] = element;
				};

				for (int k = j + 1; k <= bs; k++)
					Hami[j][k] = 0;
			};
		};

		delete[] sign;

	};

	void FindEigensystem(void)
	{

		for (int i = 1; i <= nBlocks; i++)
		{
			int bs = BlockSize[i];
			double** Hami = Ham[i];
			double* EigenVali = EigenVal[i];

			if (bs > 1)
			{

				gsl_matrix *ham  = gsl_matrix_alloc(bs, bs);
				gsl_matrix *evec = gsl_matrix_alloc(bs, bs);
				gsl_vector *eval = gsl_vector_alloc(bs);
				gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc(bs);

				int n = 0;
				for (int j = 1; j <= bs; j++)
					for (int k = 1; k <= bs; k++)
						ham->data[n++] = Hami[j][k];

				gsl_eigen_symmv(ham, eval, evec, w);

				n = 0;
				for (int j = 1; j <= bs; j++)
				{
					for (int k = 1; k <= bs; k++) Hami[j][k] = evec->data[n++];
					EigenVali[j] = eval->data[j-1];
				};

				gsl_eigen_symmv_free(w);
				gsl_matrix_free(ham);
				gsl_matrix_free(evec);
				gsl_vector_free(eval);

			}
			else
			{
				EigenVali[1] = Hami[1][1];
				Hami[1][1] = 1;
			};

		};

	};

	void HamTanspose(void)
	{
		double temp = 0;
		int tmp = 0;
		double** Hami = NULL;

		for (int i = 1; i <= nBlocks; i++)
		{
			tmp = BlockSize[i];
			if (tmp > 1)
			{
				Hami = Ham[i];
				for (int j = 1; j <= tmp; j++)
				for (int k = j + 1; k <= tmp; k++)
				{
					temp = Hami[j][k];
					Hami[j][k] = Hami[k][j];
					Hami[k][j] = temp;
				};
			};
		};
	};

	void CalcFreqIntens(void)
	{
		double** HamCurrBlock = NULL;
		double*  HamCurrBlockj = NULL;
		double** HamNextBlock = NULL;
		double*  HamNextBlockk = NULL;
		double*  EigenValCurr = NULL;
		double*	 EigenValNext = NULL;
		bool**   Perturbi = NULL;
		bool*    Perturbik = NULL;
		double   intensity = 0, maxintens = 0;
		double   threshold = 0;
		int	freqnum = 1,
			bs = 0,
			bs_next = 0;

		double* perturb = new double[MaxBlockSize + 1];
		perturb[0] = 0;

		HamTanspose();

		for (int i = 1; i <= nBlocks - 1; i++)
		{
			bs = BlockSize[i];
			bs_next = BlockSize[i + 1];
			HamCurrBlock = Ham[i];
			HamNextBlock = Ham[i + 1];
			EigenValCurr = EigenVal[i];
			EigenValNext = EigenVal[i + 1];
			Perturbi = Perturbation[i];

			for (int j = 1; j <= bs; j++)
			{
				HamCurrBlockj = HamCurrBlock[j];
				for (int k = 1; k <= bs_next; k++)
				{
					Perturbik = Perturbi[k];
					perturb[k] = 0;
					for (int ii = 1; ii <= bs; ii++)
					if (Perturbik[ii])
						perturb[k] += HamCurrBlockj[ii];
				};

				for (int k = 1; k <= bs_next; k++)
				{
					HamNextBlockk = HamNextBlock[k];
					intensity = 0;
					for (int ii = 1; ii <= bs_next; ii++)
						intensity += perturb[ii] * HamNextBlockk[ii];

					intensity = 0.25 * intensity * intensity;
					Intens[freqnum] = intensity;
					if (maxintens < intensity) maxintens = intensity;
					Freqs[freqnum] = 0.5 * (EigenValCurr[j] - EigenValNext[k]);
					freqnum++;
				};
			};
		};

		delete[] perturb;

		threshold = IntensityThreshold * maxintens;
		nFreqsFiltered = 0;
		for (int i = 1; i <= nFreqs; i++)
		{
			intensity = Intens[i];
			if (intensity > threshold)
			{
				nFreqsFiltered++;
				FreqsFiltered[nFreqsFiltered] = Freqs[i];
				IntensFiltered[nFreqsFiltered] = intensity;
			};
		};
	};

	void ComputeFreqIntens(void)
	{
		Build();
		FindEigensystem();
		CalcFreqIntens();
	};

	void PrintSpinSystem(ostream& ostr, double SF)
	{

		int oldprecision = int(ostr.precision());
		ostr.precision(3);

		ostr << "Chemical shifts (ppm):" << endl;
		for (int i = 1; i <= nSpins; i++)
			ostr << setw(10) << right << Parameters[Offs[i]] / SF;
		ostr << endl;
		ostr.precision(4);
		ostr << "J-coupling constants (Hz):" << endl;
		for (int i = 1; i < nSpins; i++)
		{
			for (int j = 1; j <= nSpins; j++)
				if (j <= i) ostr << setw(10) << ' ';
				else ostr << setw(10) << right << Parameters[Jcoup[i][j]];
			ostr << endl;
		};

		ostr.precision(oldprecision);

	};

};

class Data
{

public:
	int nPoints;
	double Magnitude;
	double MagnitudeMultiplier;
	double* Points;
	char* Filename;

	Data(void)
	{

		Filename = new char[256];
		strcpy(Filename, "");
		nPoints = 0;
		Magnitude = 0;
		MagnitudeMultiplier = 1;
		Points = NULL;

	};

	void Create(void)
	{

		if(nPoints > 0)
		{
			Points = new double[nPoints + 1];
			for(int i = 0; i <= nPoints; i++)
				Points[i] = 0;
		};

	};

	void Clean(void)
	{

		if (nPoints > 0)
		{
			delete[] Points;
			Points = NULL;
			nPoints = 0;
			Magnitude = 0;
			MagnitudeMultiplier = 1;
		};

	};

	void Rescale(int PivotPointIndex)
	{

		double rescaling = Magnitude * MagnitudeMultiplier / (PivotPointIndex == 0 ? CalcMagnitude() : Points[PivotPointIndex]);
		for (int i = 1; i <= nPoints; i++)
			Points[i] = rescaling*Points[i];

	};

	double CalcMagnitude(void)
	{

		double MaxIntens = 0;
		for(int i = 1; i <= nPoints; i++)
			if(MaxIntens < abs(Points[i])) MaxIntens = abs(Points[i]);
		return(MaxIntens);

	};

	void SaveSpecToFile(void)
	{

		if (Filename[0] == '-') return;
		int point = 0;
		ofstream ostr(Filename, ios::out | ios::binary);
		for(int i = 1; i <= nPoints; i++)
		{
			point = int(Points[i]);
			ostr.write((char*)&point,4);
		};
		ostr.close();

	};

	void LoadSpecFromFile(void)
	{

		int point = 0;

		Clean();
		
		if (int(Filename[0]) == 0) { cout << "Can't load spectrum, filename is empty!" << endl; exit_; };
		ifstream istr(Filename, ios::in | ios::binary | ios::ate);
		if (!istr) { cout << "File " << Filename << " does not exists!" << endl; exit_; };

		nPoints = (int)istr.tellg()/4;
		istr.seekg(0);
		Points = new double[nPoints+1];
		Points[0] = 0;

		for(int i = 1; i <= nPoints; i++)
		{
			istr.read((char*)&point,4);
			Points[i] = double(point);
			if(Magnitude < abs(Points[i])) Magnitude = abs(Points[i]);	
		};
		istr.close();

	};

	void ZeroData(void)
	{

		for(int i = 0; i <= nPoints; i++)
			Points[i] = 0;

	};

};

class PreOptimization
{

public:
	Hamiltonian* Hami;
	char* DatasetPath;
	int ExpProcNo;
	int BrExpProcNo;
	int TheorProcNo;
	Data ExperimentalSpec;
	Data ExperimentalSpecWithBroadening;
	Data TheoreticalSpec;
	double LineWidth;
	double LB;
	double SW_h;
	double O1_h;
	double BF;
	double SR;
	int nIntervals;
	int* StartPoint;
	int* EndPoint;
	int nPointsRated;
	double* ExpSpecPointsOnIntervals;
	double* TheorSpecPointsOnIntervals;
	double* FreqsOnIntervals;
	double SumOfExpSquaresOnIntervals;
	char* SpectraTextOutputFilename;

	PreOptimization(Hamiltonian* Ham, ifstream& istr)
	{

		Hami = Ham;
		DatasetPath = new char[256];
		SpectraTextOutputFilename = new char[256];
		ExpProcNo = 0;
		BrExpProcNo = 0;
		TheorProcNo = 0;
		char textline[256];
		BF = 0;
		O1_h = 0;
		SW_h = 0;
		SR = 0;
		LineWidth = 0;
		LB = 0;
		SumOfExpSquaresOnIntervals = 0;
		nIntervals = 0;
		nPointsRated = 0;

		StartPoint = NULL;
		EndPoint = NULL;
		ExpSpecPointsOnIntervals = NULL;
		TheorSpecPointsOnIntervals = NULL;
		FreqsOnIntervals = NULL;

		istr >> textline >> DatasetPath;
		for (int i = 0; int(DatasetPath[i]) != 0; i++) if (DatasetPath[i] == '\\') DatasetPath[i] = char('/');

		// Reading data from ACQUS
		sprintf(textline, "%s/acqus", DatasetPath);
		ifstream dataset(textline);
		if (!dataset) { cout << "File ACQUS does not exists!" << endl; exit_; };
		while (dataset)
		{
			dataset >> textline;
			if (strstr(textline, "$BF1=") != NULL)  dataset >> BF;
			if (strstr(textline, "$O1=") != NULL)   dataset >> O1_h;
			if (strstr(textline, "$SW_h=") != NULL) dataset >> SW_h;
			if (strstr(textline, "$SW_H=") != NULL) dataset >> SW_h;
		};
		dataset.close();

		istr >> textline >> textline;
		if (!isunsignint(textline)) { cout << "Wrong experimental spectrum processing number." << endl; exit_; };
		ExpProcNo = stoi(textline);
		// Reading data from PROCS of experimental spectrum
		sprintf(ExperimentalSpec.Filename, "%s/pdata/%i/procs", DatasetPath, ExpProcNo);
		dataset.open(ExperimentalSpec.Filename);
		if (!dataset) { cout << "File PROCS in processing number " << ExpProcNo << " does not exists!" << endl; exit_; };
		sprintf(ExperimentalSpec.Filename, "%s/pdata/%i/1r", DatasetPath, ExpProcNo);
		while (dataset)
		{
			dataset >> textline;
			if (strstr(textline, "$SF=") != NULL) dataset >> SR;
			if (strstr(textline, "$SI=") != NULL) dataset >> TheoreticalSpec.nPoints;
		};
		dataset.close();
		SR -= BF; SR *= 1e6;
		SR = double(int(SR * 1e4)) / 1e4;

		istr >> textline >> textline;
		if (textline[0] == '-')
			strcpy(ExperimentalSpecWithBroadening.Filename, "-");
		else
		{
			if (!isunsignint(textline)) { cout << "Check the processing number for experimental spectrum with broadening." << endl; exit_; };
			BrExpProcNo = stoi(textline);
			sprintf(ExperimentalSpecWithBroadening.Filename, "%s/pdata/%i/1r", DatasetPath, BrExpProcNo);
		};
		
		istr >> textline >> textline;
		if (textline[0] == '-')
			strcpy(TheoreticalSpec.Filename, "-");
		else
		{
			if (!isunsignint(textline)) { cout << "Check the processing number for theoretical spectrum." << endl; exit_; };
			TheorProcNo = stoi(textline);
			sprintf(TheoreticalSpec.Filename, "%s/pdata/%i/1r", DatasetPath, TheorProcNo);
		};
		istr.getline(textline, 256);   //  Rest of proc. no. line

		ExperimentalSpecWithBroadening.nPoints = TheoreticalSpec.nPoints;

	};

	void LoadIntervals(void)
	{
		stringstream parse;
		int position = 0;
		bool ok = false;
		char textline[256];
		sprintf(textline, "%s/pdata/%i/integrals.txt", DatasetPath, ExpProcNo);
		ifstream istr(textline);
		if (!istr) { cout << "File integrals.txt does not exist in the processing folder with experimental spectrum!" << endl; exit_; };
		
		nIntervals = 0;

		while (istr)
		{
			istr.getline(textline, 256);
			if (strstr(textline, "Integral") != NULL) { ok = true; break; };
		};
		if (!ok)  { cout << "File integrals.txt in the processing folder with experimental spectrum is corrupted!" << endl; exit_; };

		position = int(istr.tellg());

		int i = 0;
		istr.getline(textline, 256);
		while (istr)
		{
			i++;
			parse << textline;
			parse >> nIntervals;
			if (nIntervals != i)  { cout << "File integrals.txt in the processing folder with experimental spectrum is corrupted!" << endl; exit_; };
			parse.str(""); parse.clear();
			istr.getline(textline, 256);
		};

		if (nIntervals < 1)  { cout << "There is no defined intervals in the integrals.txt file for the experimental spectrum!" << endl; exit_; };

		istr.clear();
		istr.seekg(position);

		int j = 0;
		double StartPPM = 0;
		double EndPPM = 0;
		double IntVal = 0;

		double PPMOffset = (O1_h - SR + 0.5 * SW_h) / (BF + SR*1e-6);
		double PPMStep = SW_h / ((BF + SR*1e-6) * (ExperimentalSpec.nPoints - 1));

		StartPoint = new int[nIntervals + 1];
		EndPoint = new int[nIntervals + 1];
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
			if (StartPoint[i] > ExperimentalSpec.nPoints) StartPoint[i] = ExperimentalSpec.nPoints;
			EndPoint[i] = int((PPMOffset - EndPPM) / PPMStep) + 1;
			if (EndPoint[i] < 1) EndPoint[i] = 1;
			if (EndPoint[i] > ExperimentalSpec.nPoints) EndPoint[i] = ExperimentalSpec.nPoints;
			if (EndPoint[i] - StartPoint[i] <= 2) { cout << "Spectrum interval " << i << " is too small. Please increase the number of points." << endl; exit_; };
			nPointsRated += (EndPoint[i] - StartPoint[i] + 1);
		};

		istr.close();

		ExpSpecPointsOnIntervals = new double[nPointsRated + 1];
		TheorSpecPointsOnIntervals = new double[nPointsRated + 1];
		FreqsOnIntervals = new double[nPointsRated + 1];
		for (int i = 0; i <= nPointsRated; i++)
		{
			ExpSpecPointsOnIntervals[i] = 0;
			TheorSpecPointsOnIntervals[i] = 0;
			FreqsOnIntervals[i] = 0;
		};

		double Offset = O1_h + 0.5 * SW_h - SR;
		double FreqStep = SW_h / (TheoreticalSpec.nPoints - 1);
		int n = 1;
		for (int i = 1; i <= nIntervals; i++)
			for (int j = StartPoint[i]; j <= EndPoint[i]; j++)
				FreqsOnIntervals[n++] = Offset - FreqStep*(j - 1);

	};

	double GetIntervalStartFreq(int Inerval)
	{
		if (Inerval < 1 || Inerval > nIntervals) return(0.0);
		return(O1_h + 0.5 * SW_h - SR - (StartPoint[Inerval] - 1) * SW_h / (TheoreticalSpec.nPoints - 1));
	};

	double GetIntervalEndFreq(int Inerval)
	{
		if (Inerval < 1 || Inerval > nIntervals) return(0.0);
		return(O1_h + 0.5 * SW_h - SR - (EndPoint[Inerval] - 1) * SW_h / (TheoreticalSpec.nPoints - 1));
	};

	void CalcExpSpecMagnOnIntervals(void)
	{
		if (nIntervals > 0)
		{
			ExperimentalSpec.Magnitude = 0;
			ExperimentalSpec.MagnitudeMultiplier = 1;
			for (int i = 1; i <= nIntervals; i++)
				for (int j = StartPoint[i]; j <= EndPoint[i]; j++)
					if (ExperimentalSpec.Magnitude < abs(ExperimentalSpec.Points[j])) ExperimentalSpec.Magnitude = abs(ExperimentalSpec.Points[j]);
		};
	};

	void BroadOnIntervals(double LWBroad)
	{

		LB = LWBroad;
		ExperimentalSpecWithBroadening.ZeroData();

		double Betta = 2 * SW_h/((ExperimentalSpec.nPoints - 1) * LB); Betta *= Betta;
		int nPoints = int(ceill(sqrt(199 / Betta)));
		if(nPoints > 3)
		{
			double* Points = new double [nPoints + 1];
			for(int i = 0; i <= nPoints; i++)
				Points[i] = 1 / (1 + Betta * i * i);

			int startk = 0, stopk = 0;
			for(int i = 1; i <= nIntervals; i++)
				for(int j = StartPoint[i]; j <= EndPoint[i]; j++)
				{
					startk = j - nPoints; if(startk <= 0) startk = 1;
					stopk = j + nPoints;  if(stopk > ExperimentalSpec.nPoints) stopk = ExperimentalSpec.nPoints;
					for(int k = startk; k <= stopk; k++)
						ExperimentalSpecWithBroadening.Points[j] += ExperimentalSpec.Points[k] * Points[abs(j-k)];
				};
			delete[] Points;
			ExperimentalSpecWithBroadening.Rescale(0);
		}
		else
		{
			for(int i = 1; i <= nIntervals; i++)
				for(int j = StartPoint[i]; j <= EndPoint[i]; j++)
						ExperimentalSpecWithBroadening.Points[j] = ExperimentalSpec.Points[j];
		};

		SumOfExpSquaresOnIntervals = 0;
		int n = 1;
		double tmp = 0;
		for(int i = 1; i <= nIntervals; i++)
			for(int j = StartPoint[i]; j <= EndPoint[i]; j++)
			{
				tmp = ExperimentalSpecWithBroadening.Points[j];
				ExpSpecPointsOnIntervals[n++] = tmp; tmp *= tmp;
				SumOfExpSquaresOnIntervals += tmp;
			};

	};

	void CalcSpecOnIntervals(void)
	{

		double LW = (LineWidth + LB) / 2;
		double tmp = 0, tmp1 = 0;
		double CurrPointFreq = 0;
		double sqLW = LW * LW;
		for (int i = 1; i <= nPointsRated; i++)
		{
			CurrPointFreq = FreqsOnIntervals[i];
			tmp1 = 0;
			for (int j = 1; j <= Hami->nFreqsFiltered; j++)
			{
				tmp = CurrPointFreq - Hami->FreqsFiltered[j]; tmp *= tmp;
				tmp1 += Hami->IntensFiltered[j] * LW / (tmp + sqLW);
			};
			TheorSpecPointsOnIntervals[i] = tmp1;
		};

		//Rescale TheorSpec on intervals
		tmp = 0;
		for (int i = 1; i <= nPointsRated; i++)
		if (tmp < abs(TheorSpecPointsOnIntervals[i])) tmp = abs(TheorSpecPointsOnIntervals[i]);
		tmp1 = TheoreticalSpec.Magnitude * TheoreticalSpec.MagnitudeMultiplier / tmp;
		for (int i = 1; i <= nPointsRated; i++)
			TheorSpecPointsOnIntervals[i] *= tmp1;

	};

	void ComputeFullSpectrum(void)
	{

		double LW = (LineWidth + LB) / 2;
		double Offset = O1_h + 0.5 * SW_h - SR;
		double FreqStep = SW_h / (TheoreticalSpec.nPoints - 1);
		double tmp = 0, tmp1 = 0;
		double sqLW = LW * LW;
		double CurrFreq = 0;

		Hami->ComputeFreqIntens();

		for (int i = 1; i <= TheoreticalSpec.nPoints; i++)
		{
			CurrFreq = Offset - FreqStep*(i - 1);
			tmp1 = 0;
			for (int j = 1; j <= Hami->nFreqsFiltered; j++)
			{
				tmp = CurrFreq - Hami->FreqsFiltered[j]; tmp *= tmp;
				tmp1 += Hami->IntensFiltered[j] * LW / (tmp + sqLW);
			};
			TheoreticalSpec.Points[i] = tmp1;
		};

		// Rescaling of the full spectrum taking into account spectral intervals.
		int index = 0;
		tmp = 0;
		for (int i = 1; i <= nIntervals; i++)
			for (int j = StartPoint[i]; j <= EndPoint[i]; j++)
				if (tmp < TheoreticalSpec.Points[j]) { tmp = TheoreticalSpec.Points[j]; index = j; };
		TheoreticalSpec.Rescale(index);

	};

	void ComputeSpecOnIntervals(void)
	{

		Hami->ComputeFreqIntens();
		CalcSpecOnIntervals();

	};

	double SumOfSquareDeviationsOnIntervals(void)
	{
		double tmp = 0;
		double res = 0;
		for (int i = 1; i <= nPointsRated; i++)
		{
			tmp = ExpSpecPointsOnIntervals[i] - TheorSpecPointsOnIntervals[i]; tmp *= tmp;
			res += tmp;
		};
		return res;
	};

	double MeanSquareDeviation(Data Dat)
	{
		return (sqrt(SumOfSquareDeviationsOnIntervals() / nPointsRated));
	};

	double Badness(void)
	{

		ComputeSpecOnIntervals();
		return SumOfSquareDeviationsOnIntervals();

	};

	double BadnessScWd(void)
	{

		CalcSpecOnIntervals();
		return SumOfSquareDeviationsOnIntervals();

	};

	double CalcRFactor(void)
	{
		double rfactor = 1;
		if (SumOfExpSquaresOnIntervals > 0)
			rfactor = SumOfSquareDeviationsOnIntervals() / SumOfExpSquaresOnIntervals;
		return(100 * sqrt(rfactor));
	};

	void SaveSpecsOnIntervalsTXT(void)
	{

		if (SpectraTextOutputFilename[0] == '-') return;

		ofstream ostr(SpectraTextOutputFilename);
		ostr << setw(16) << left << "Freq(Hz)";
		ostr << setw(14) << right << "Exp.Intens.";
		ostr << setw(16) << right << "Theor.Intens.";
		ostr << setw(13) << right << "Diff." << endl;

		ostr.precision(6);
		ostr << fixed;
		for (int i = 1; i <= nPointsRated; i++)
		{
			ostr << setw(11) << left << FreqsOnIntervals[i];
			ostr << setw(16) << right << int(ExpSpecPointsOnIntervals[i]);
			ostr << setw(16) << right << int(TheorSpecPointsOnIntervals[i]);
			ostr << setw(16) << right << int(ExpSpecPointsOnIntervals[i]) - int(TheorSpecPointsOnIntervals[i]) << endl;
		};

		ostr.close();

	};

};

class OptHamiltonian
{
public:
	int nVarParams;
	int* VarParams;
	int nParams;
	int nBroadenings;
	char* InputParameters;
	char* OutputParameters;
	double* Parameters;
	double* ParameterErrors;
	double* LBs;
	char** ParNames;
	PreOptimization* PreOpt;
	bool ErrorsComputed;
	int nPreOptParams;
	int* PreOptParams;

	OptHamiltonian(PreOptimization* Pre)
	{

		PreOpt = Pre;
		nParams = Pre->Hami->nParams + 2;
		nVarParams = 0;
		nBroadenings = 0;
		LBs = NULL;
		ErrorsComputed = false;
		InputParameters = new char[256];
		OutputParameters = new char[256];
		Parameters = new double[nParams + 1];
		ParameterErrors = new double[nParams + 1];
		ParNames = new char*[nParams + 1];
		VarParams = new int[nParams + 1];
		ParNames[0] = NULL; Parameters[0] = 0;
		ParameterErrors[0] = 0; VarParams[0] = 0;
		for (int i = 1; i <= nParams; i++)
		{
			Parameters[i] = 0;
			ParameterErrors[i] = 0;
			VarParams[i] = 0;
			ParNames[i] = new char[256];
			strcpy(ParNames[i],"");
		};

		nPreOptParams = 0;
		PreOptParams = new int[3];
		for (int i = 0; i <= 2; i++)
			PreOptParams[i] = 0;

	};

	void LoadParameters(void)
	{

		char textline[256];
		int ParNum = 0;
		ifstream istr(InputParameters);
		if (!istr) { cout << "File " << InputParameters << " does not exists!" << endl; exit_; };

		for (int i = 1; i <= nParams - 2; i++)
		{
			istr >> textline;
			if (!isunsignint(textline)) { cout << "Incorrect number of parameter " << i << "." << endl; exit_; };
			ParNum = stoi(textline);
			if (ParNum != i) { cout << "Incorrect number of parameter " << i << "." << endl; exit_; };
			istr >> ParNames[i];
			istr >> textline;
			if (!isreal(textline)) { cout << "Incorrect value of parameter " << i << "." << endl; exit_; };
			Parameters[i] = stod(textline);
			istr.getline(textline, 256);
		};

		istr >> textline;
		if (!isunsignint(textline)) { cout << "Incorrect number of parameter " << nParams - 1 << " (linewidth)." << endl; exit_; };
		ParNum = stoi(textline);
		if (ParNum != nParams - 1) { cout << "Incorrect number of parameter " << nParams - 1 << " (linewidth)." << endl; exit_; };
		istr >> ParNames[nParams - 1];
		istr >> textline;
		if (!isunsignreal(textline)) { cout << "Incorrect value of parameter " << nParams - 1 << " (linewidth)." << endl; exit_; };
		Parameters[nParams - 1] = stod(textline);
		if (Parameters[nParams - 1] == 0) { cout << "Incorrect value of parameter " << nParams - 1 << " (linewidth)." << endl; exit_; };
		istr.getline(textline, 256);

		istr >> textline;
		if (!isunsignint(textline)) { cout << "Incorrect number of parameter " << nParams << " (spectrum magnitude)." << endl; exit_; };
		ParNum = stoi(textline);
		if (ParNum != nParams) { cout << "Incorrect number of parameter " << nParams << " (spectrum magnitude)." << endl; exit_; };
		istr >> ParNames[nParams];
		istr >> PreOpt->TheoreticalSpec.Magnitude;
		if (istr.fail() || PreOpt->TheoreticalSpec.Magnitude < 0) { cout << "Incorrect value of parameter " << nParams << " (spectrum magnitude)." << endl; exit_; };
		PreOpt->TheoreticalSpec.MagnitudeMultiplier = 1;
		Parameters[nParams] = 1;
		istr.close();

	};

	void CheckSpinOffsets(ostream& ostr)
	{

		bool check = false;
		int MaxOffs = PreOpt->Hami->Offs[PreOpt->Hami->nSpins];

		for (int i = 1; i <= MaxOffs; i++)
		{
			check = false;
			for (int j = 1; j <= PreOpt->nIntervals; j++)
				if (Parameters[i] <= PreOpt->GetIntervalStartFreq(j) && Parameters[i] >= PreOpt->GetIntervalEndFreq(j)) check = true;
			if (!check) ostr << "Warning! Chemical shift no. " << i << " (" << Parameters[i] << ") does not fall into any of defined spectral intervals." << endl;
		};

		for (int i = 1; i <= PreOpt->nIntervals; i++)
		{
			check = false;
			for (int j = 1; j <= MaxOffs; j++)
				if (Parameters[j] <= PreOpt->GetIntervalStartFreq(i) && Parameters[j] >= PreOpt->GetIntervalEndFreq(i)) check = true;
			if (!check) ostr << "Warning! Spectral interval no. " << i << " does not contain any chemical shift." << endl;
		};

		if (!check) ostr << endl;

	};

	void CheckBroadSequence(ostream& ostr)
	{

		bool printwarn = false;
		for (int i = 2; i <= nBroadenings; i++)
		if (LBs[i] >= LBs[i - 1]) printwarn = true;
		if (printwarn) ostr << "Warning! The sequence of additional broadenings should monotonically decrease!" << endl << endl;

	};

	void LoadBroadenings(ifstream& istr)
	{

		char textline[256];
		char textline1[256];
		stringstream parse;

		istr.getline(textline, 256);   //  LBs
		parse << textline;
		while (parse) { parse >> textline1;	nBroadenings++; };
		nBroadenings -= 2;
		if (nBroadenings < 1) { cout << "Wrong number of broadenings (" << nBroadenings << "), should be at least 1." << endl; exit_; };
		LBs = new double[nBroadenings + 1];
		LBs[0] = 0;
		parse.clear(); parse << textline;
		parse >> textline;
		for (int i = 1; i <= nBroadenings; i++)
		{
			parse >> textline;
			if (!isunsignreal(textline)) { cout << "Wrong value of broadening number " << i << "." << endl; exit_; };
			LBs[i] = stod(textline);
		};

	};

	void LoadVarParameters(ifstream& istr)
	{

		char textline[256];
		int ParNum = 0;

		int i = 1;
		istr >> textline;
		while (istr && i <= nParams)
		{
			if (!isunsignint(textline)) { cout << "Wrong index of varied parameter " << i << "." << endl; exit_; };
			VarParams[i] = stoi(textline);
			if (VarParams[i] > nParams) { cout << "Wrong index of varied parameter " << i << "." << endl; exit_; };
			i++;
			istr >> textline;
		};
		nVarParams = i - 1;

		if (nVarParams == 0) { cout << "There is no varied parameters!" << endl; exit_; };
		
		for (i = 2; i <= nVarParams; i++)
			if (VarParams[i] <= VarParams[i - 1]) { cout << "Wrong index of varied parameter " << i << "." << endl; exit_; };

		nPreOptParams = 0;
		bool preopt = false;
		for (i = 1; i <= nVarParams; i++)
			if ((VarParams[i] != nParams - 1) & (VarParams[i] != nParams)) preopt = true;
		if (preopt)
		{
			for (i = 1; i <= nVarParams; i++)
			{
				if (VarParams[i] == nParams - 1) { PreOptParams[++nPreOptParams] = nParams - 1; };
				if (VarParams[i] == nParams) { PreOptParams[++nPreOptParams] = nParams; };
			};
		};

	};

	void SetParametersToHamiltonian(void)
	{

		for(int i = 1; i <= nParams - 2; i++)
			PreOpt->Hami->Parameters[i] = Parameters[i];
		PreOpt->LineWidth = Parameters[nParams - 1];
		PreOpt->TheoreticalSpec.MagnitudeMultiplier = Parameters[nParams];

	};

	double Badness(void)
	{

		SetParametersToHamiltonian();
		return PreOpt->Badness();

	};

	double Badness(double* VarPars)
	{

		for(int i = 1; i <= nVarParams; i++)
			Parameters[VarParams[i]] = VarPars[i];
		return Badness();

	};

	void ComputeErrors()
	{

		ErrorsComputed = false;
		double** TheorSpecDerivativesOnIntervals = new double*[nParams + 1];
		TheorSpecDerivativesOnIntervals[0] = NULL;
		gsl_matrix *DTD = gsl_matrix_alloc(nParams, nParams);

		int nPoints = PreOpt->nPointsRated;

		for (int i = 1; i <= nParams; i++)
		{
			TheorSpecDerivativesOnIntervals[i] = new double[nPoints + 1];
			for (int j = 0; j <= nPoints; j++)
				TheorSpecDerivativesOnIntervals[i][j] = 0;
		};

		double step = 1.0e-7;
		double* Data = new double[nPoints + 1];
		Data[0] = 0;

		double Badn = Badness();

		for (int i = 1; i <= nPoints; i++)
			Data[i] = PreOpt->TheorSpecPointsOnIntervals[i];

		for (int i = 1; i <= nParams; i++)
		{
			double Par = Parameters[i];
			Parameters[i] += step;
			SetParametersToHamiltonian();
			PreOpt->ComputeSpecOnIntervals();
			Parameters[i] = Par;
			for (int j = 1; j <= nPoints; j++)
				TheorSpecDerivativesOnIntervals[i][j] = (PreOpt->TheorSpecPointsOnIntervals[j] - Data[j]) / step;
		};

		for (int i = 1; i <= nPoints; i++)
			PreOpt->TheorSpecPointsOnIntervals[i] = Data[i];

		for (int i = 1; i <= nParams; i++)
		for (int j = i; j <= nParams; j++)
		{
			DTD->data[nParams*(i - 1) + j - 1] = 0;
			for (int k = 1; k <= nPoints; k++)
				DTD->data[nParams*(i - 1) + j - 1] += TheorSpecDerivativesOnIntervals[i][k] * TheorSpecDerivativesOnIntervals[j][k];
		};

		for (int i = 0; i < nParams; i++)
			for (int j = i + 1; j < nParams; j++)
				DTD->data[nParams*j + i] = DTD->data[nParams*i + j];

		gsl_matrix *evec = gsl_matrix_alloc(nParams, nParams);
		gsl_vector *eval = gsl_vector_alloc(nParams);
		gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc(nParams);

		gsl_eigen_symmv(DTD, eval, evec, w);
		gsl_eigen_symmv_free(w);
		gsl_matrix_free(DTD);

		gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_ABS_ASC);

		gsl_matrix *temp = gsl_matrix_alloc(nParams, nParams);

		if (eval->data[0] > 1e-4)
		{
			for (int i = 0; i < nParams; i++)
				for (int j = 0; j < nParams; j++)
					temp->data[nParams*i + j] = evec->data[nParams*j + i] / eval->data[i];
			for (int i = 0; i < nParams; i++)
			{
				ParameterErrors[i+1] = 0.;
				for (int k = 0; k < nParams; k++)
					ParameterErrors[i+1] += evec->data[nParams*i + k] * temp->data[nParams*k + i];
			};
			for (int i = 1; i <= nParams; i++)
				ParameterErrors[i] = sqrt(ParameterErrors[i] * Badn / (nPoints - nParams));
			ParameterErrors[nParams] *= PreOpt->TheoreticalSpec.Magnitude;

			ErrorsComputed = true;
		};
		
		gsl_matrix_free(evec);
		gsl_matrix_free(temp);
		gsl_vector_free(eval);

	};

	void SaveParameters(void)
	{

		int parnamelen = 0;
		ofstream ostr(OutputParameters);
		
		for (int i = 1; i <= nParams; i++)
		{
			int tmp = strlen(ParNames[i]);
			if (parnamelen < tmp) parnamelen = tmp;
		};
		parnamelen += 4;

		ostr << fixed;
		for (int i = 1; i <= nParams; i++)
		{
			double Param = Parameters[i];
			if (i == nParams) {	ostr << scientific; Param *= PreOpt->TheoreticalSpec.Magnitude; };
			ostr << setw(4) << left << i << setw(parnamelen) << left << ParNames[i];
			ostr << setw(13) << right << Param;
			if (ErrorsComputed) ostr << setw(4) << ' ' << "+-" << ParameterErrors[i];
			ostr << endl;
		};
		ostr << fixed;
		if (!ErrorsComputed) ostr << "Errors are not calculated due to singularity of normal equations matrix!" << endl;
		ostr << endl;

		ostr << Title << endl << endl;
		CheckSpinOffsets(ostr);
		CheckBroadSequence(ostr);

		ostr << "Line Broadening: " << fixed << PreOpt->LB << endl;
		ostr << "Theoretical spectrum linewidth: " << PreOpt->LB + Parameters[nParams - 1] << endl;
		ostr << "RSS Value: " << scientific << Badness() << endl;
		int oldprecision = int(ostr.precision());
		ostr.precision(2);
		ostr << "R-Factor: " << fixed << PreOpt->CalcRFactor() << " %" << endl;
		ostr.precision(oldprecision);
		ostr << endl;

		PreOpt->Hami->PrintSpinSystem(ostr, (PreOpt->BF + PreOpt->SR * 1e-6));

		ostr.close();

	};

};

OptHamiltonian* HamOpt = NULL;

double GlobalBadnessScWdGSL(const gsl_vector *v, void *params)
{

	for (int i = 1; i <= HamOpt->nPreOptParams; i++)
		HamOpt->Parameters[HamOpt->PreOptParams[i]] = v->data[i - 1];

	HamOpt->PreOpt->LineWidth = HamOpt->Parameters[HamOpt->nParams - 1];
	HamOpt->PreOpt->TheoreticalSpec.MagnitudeMultiplier = HamOpt->Parameters[HamOpt->nParams];
	return HamOpt->PreOpt->BadnessScWd();

};

double GlobalBadnessGSL(const gsl_vector *v, void *params)
{

	for (int i = 1; i <= HamOpt->nVarParams; i++)
		HamOpt->Parameters[HamOpt->VarParams[i]] = v->data[i - 1];
	return HamOpt->Badness();

};

int main (void)
{

#ifdef _WIN32                                                     // Running under Windows Platform
	char *env = NULL;
	if (getenv("XWINNMRHOME") != NULL) ExitWithPause = false;     // Running under TopSpin
	if ((env = getenv("ExitWithoutPause")) != NULL)               // Checking ExitWithoutPause enviroment variable
		if (strcmp(env, "1") == 0) ExitWithPause = false;
#endif

	bool SimMode = true;
	bool MagnitudeFromExpSpec = true;
	char textline[256];

	// Parsing input control file and initialize relevant data structures.
	ifstream input("Input_Data.txt");
	if (!input) { cout << "File Input_Data.txt does not exists!" << endl; exit_; };
	input.getline(Title, 256);      //  Title
	input.getline(textline, 256);   //  Empty line
	input >> textline >> SimMode;   //  Sim mode
	if (input.fail()) { cout << "Wrong simulation mode value, sould be 0 or 1." << endl; exit_; };
	input.getline(textline, 256);   //  Rest of sim mode line
	input.getline(textline, 256);   //  Empty line
	input.getline(textline, 256);   //  Spin System
	Hamiltonian Hami(input);
	input.getline(textline, 256);   //  Empty line
	input.getline(textline, 256);   //  Spectra parameters
	PreOptimization Pre(&Hami, input);
	HamOpt = new OptHamiltonian(&Pre);
	input.getline(textline, 256);   //  Empty line
	input.getline(textline, 256);   //  Optimization parameters
	input >> textline >> HamOpt->InputParameters; // Input parameters filename
	HamOpt->LoadParameters();       // Reading input parameters from file
	Pre.TheoreticalSpec.Create();
	input >> textline >> HamOpt->OutputParameters; // Output parameters filename
	input.getline(textline, 256);   //  Rest of output parameters filename line
	input >> textline >> Pre.SpectraTextOutputFilename; // Filename for spectra in ASCII text format
	input.getline(textline, 256);   //  Rest of SpectraTextOutputFilename line
	HamOpt->LoadBroadenings(input); //  Parsing of broadenings list
	HamOpt->CheckBroadSequence(cout);
	input >> textline >> MagnitudeFromExpSpec;   //  MagnitudeFromExpSpec
	if (input.fail()) { cout << "Wrong MagnitudeFromExpSpec flag value, sould be 0 or 1." << endl; exit_; };
	input.getline(textline, 256);   //  Rest of MagnitudeFromExpSpec line
	Pre.ExperimentalSpec.LoadSpecFromFile();
	if (Pre.ExperimentalSpec.nPoints != Pre.TheoreticalSpec.nPoints) { cout << "File with experimental spectrum is corrupted." << endl; exit_; };
	Pre.LoadIntervals();            // Loading and computing of points intervals on the basis of intergal regions from integrals.txt.
	HamOpt->CheckSpinOffsets(cout);

	Pre.CalcExpSpecMagnOnIntervals();
	Pre.ExperimentalSpecWithBroadening.Magnitude = Pre.ExperimentalSpec.Magnitude;
	Pre.ExperimentalSpecWithBroadening.MagnitudeMultiplier = 1;

	if (MagnitudeFromExpSpec)
	{
		Pre.TheoreticalSpec.Magnitude = Pre.ExperimentalSpec.Magnitude;
		Pre.TheoreticalSpec.MagnitudeMultiplier = 1;
		HamOpt->Parameters[HamOpt->nParams] = 1;
	};

	Pre.ExperimentalSpecWithBroadening.Create();

	if (SimMode)
	{
		HamOpt->SetParametersToHamiltonian();
		Pre.BroadOnIntervals(HamOpt->LBs[1]);
		Pre.ExperimentalSpecWithBroadening.SaveSpecToFile();
		Pre.ComputeSpecOnIntervals();
		Pre.SaveSpecsOnIntervalsTXT();
		Pre.LB = 0;
		Pre.ComputeFullSpectrum();
		Pre.TheoreticalSpec.SaveSpecToFile();
		exit_;
	};

	input.getline(textline, 256);   //  Empty line
	input.getline(textline, 256);   //  List of optimized parameters
	HamOpt->LoadVarParameters(input);
	input.close();
	// End of INPUT file parsing

	int iter = 0, n = HamOpt->nVarParams, npre = HamOpt->nPreOptParams;

	const gsl_multimin_fminimizer_type *T =	gsl_multimin_fminimizer_nmsimplex2;
	gsl_multimin_function minex_func;

	gsl_multimin_fminimizer *s = gsl_multimin_fminimizer_alloc(T, n);
	gsl_vector *ss = gsl_vector_alloc(n), *x = gsl_vector_alloc(n);

	gsl_multimin_fminimizer *sscwd = NULL;
	gsl_vector *ssscwd = NULL, *xscwd = NULL;

	if (npre > 0)
	{
		sscwd = gsl_multimin_fminimizer_alloc(T, npre);
		ssscwd = gsl_vector_alloc(npre);
		xscwd = gsl_vector_alloc(npre);
	};

	int status;
	double size;

	for (int k = 1; k <= HamOpt->nBroadenings; k++)
	{

		Pre.BroadOnIntervals(HamOpt->LBs[k]);

		cout << "Broadening:  " << HamOpt->LBs[k] << endl;
		iter = 0;

		// Scale and/or Linewidth preoptimization
		if (npre > 0)
		{
			if(npre == 2 ) cout << "Scale & Linewidth preoptimization:" << endl;
			else if(HamOpt->PreOptParams[1] == HamOpt->nParams - 1) cout << "Linewidth preoptimization:" << endl;
			else cout << "Scale preoptimization:" << endl;

			HamOpt->SetParametersToHamiltonian();
			Hami.ComputeFreqIntens();

			gsl_vector_set_all(ssscwd, 1.0);

			for (int i = 0; i < npre; i++)
				gsl_vector_set(xscwd, i, HamOpt->Parameters[HamOpt->PreOptParams[i + 1]]);

			minex_func.n = npre;
			minex_func.f = GlobalBadnessScWdGSL;
			minex_func.params = NULL;
			gsl_multimin_fminimizer_set(sscwd, &minex_func, xscwd, ssscwd);

			do
			{
				iter++;
				status = gsl_multimin_fminimizer_iterate(sscwd);

				if (status)
					break;

				size = gsl_multimin_fminimizer_size(sscwd);
				status = gsl_multimin_test_size(size, 1e-10);

				if (status == GSL_SUCCESS)
					cout << "converged to minimum at" << endl;

				cout << "iteration " << iter << '\t' << sscwd->fval << '\t' << size << endl;

			} while (status == GSL_CONTINUE);

			for (int i = 0; i < npre; i++)
				HamOpt->Parameters[HamOpt->PreOptParams[i + 1]] = gsl_vector_get(sscwd->x, i);

			if (HamOpt->Parameters[HamOpt->nParams - 1] < 0 && HamOpt->Parameters[HamOpt->nParams] < 0)
			{
				HamOpt->Parameters[HamOpt->nParams - 1] *= -1; HamOpt->Parameters[HamOpt->nParams] *= -1;
			};

			Pre.TheoreticalSpec.Magnitude *= HamOpt->Parameters[HamOpt->nParams];
			HamOpt->Parameters[HamOpt->nParams] = 1;

			cout << "Main optimization:" << endl;
			iter = 0;

		};

		// Main optimization
		gsl_vector_set_all(ss, 1.0);
		for (int i = 0; i < n; i++)
			gsl_vector_set(x, i, HamOpt->Parameters[HamOpt->VarParams[i + 1]]);

		minex_func.n = n;
		minex_func.f = GlobalBadnessGSL;
		minex_func.params = NULL;
		gsl_multimin_fminimizer_set(s, &minex_func, x, ss);

		do
		{
			iter++;
			status = 
				gsl_multimin_fminimizer_iterate(s);

			if (status)
				break;

			size = gsl_multimin_fminimizer_size(s);
			status = gsl_multimin_test_size(size, 1e-10);

			if (status == GSL_SUCCESS)
				cout << "converged to minimum at" << endl;

			cout << "iteration " << iter << '\t' << s->fval << '\t' << size << endl;
			
		} while (status == GSL_CONTINUE);

		for (int i = 0; i < n; i++)
			HamOpt->Parameters[HamOpt->VarParams[i + 1]] = gsl_vector_get(s->x, i);

	};


	if (npre > 0)
	{
		gsl_vector_free(xscwd);
		gsl_vector_free(ssscwd);
		gsl_multimin_fminimizer_free(sscwd);
	};

	gsl_vector_free(x);
	gsl_vector_free(ss);
	gsl_multimin_fminimizer_free(s);

	cout.precision(2); cout << fixed;
	cout << "R-Factor: " << Pre.CalcRFactor() << " %" << endl;
	cout.precision(6);

	Pre.SaveSpecsOnIntervalsTXT();
	HamOpt->SetParametersToHamiltonian();
	Pre.ComputeFullSpectrum();
	Pre.TheoreticalSpec.SaveSpecToFile();
	Pre.ExperimentalSpecWithBroadening.SaveSpecToFile();
	HamOpt->ComputeErrors();
	HamOpt->SetParametersToHamiltonian();
	HamOpt->SaveParameters();
	exit_;

};
