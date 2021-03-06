#include<iostream>
#include<fstream>
#include<cmath>
#include"vector3.h"
#include<cstdlib>
#include"Functions.h"
#include <inttypes.h>
using namespace std;

//walker
//----------------------fast_drand40()----------------------------------?
static double d40d = (double)(1LL << 40);
static uint64_t d40m1 = (1LL << 40) - 1LL;
uint64_t prev1 = 12LL;

double randoms()
{
	prev1 = (prev1 * 762939453125LL) & d40m1;
	return (double)prev1 / d40d;
}
//??????? ?????????? ???????? ??????????? q ? ????? alias ??? Walker
int walker_prep(vector<double>& p, double* q, int* alias)
{
	int p_cnt = p.size();
	int* over = new int[p_cnt];
	int* under = new int[2 * p_cnt];
	memset(over, 0, sizeof(int) * p_cnt);
	memset(under, 0, sizeof(int) * 2 * p_cnt);
	int i, j, l, k, k_o = 0, k_u = 0;
	for (i = 0; i < p_cnt; i++)
	{
		alias[i] = 0;
		q[i] = 0.;
	}
	for (i = 0; i < p_cnt; i++)
	{
		q[i] = p[i] * (double)p_cnt;
		// assign a label to q[i]: overfull q>1, underfull g<=1
		if (q[i] > 1.) { // overfull
			over[k_o] = i;
			k_o++;
		}
		else { // underfull
			under[k_u] = i;
			k_u++;
		}

	}
	for (i = 0, j = 0; i < k_u && j < k_o; i++)
	{
		l = under[i];
		k = over[j];
		alias[l] = k; // number of q which >1
		q[k] = (q[k] + q[l]) - 1.; // new value of q which was >1
		if (q[k] < 1.) {
			j++; // take next over element
			under[k_u] = k; // put new under element
			k_u++; // add the number of under elements
		}
		else if ((q[k] - 1.) <= 1e-15)
			alias[k] = k;
	}
	return 0;
}
//??????? Walker ?????????? ????????? ??????, ????????? ? ???????????? ? ???????????? ?????????? :
int walker(int p_cnt, double* q, int* alias)
{
	int i, num;

	double rnd = randoms() * p_cnt;
	i = (int)rnd;
	rnd = randoms();
	if (rnd <= q[i]) { num = i; }
	else { num = alias[i]; }

	return num;
}
extern const double Pi;
class Calculation
{
	//variables
public:
	//add drift velocity
	double V;
	//diffusion coefficient
	double D;
	//average particle lifetime
	double tau;
	//Peclet number
	double Pe;
	//constant mu
	double mu;
	//epsilon boundaru size
	double epsilon;
	//layer thickness 
	double H;
	//photon flight length factor 
	double l;
	//reverse recombination probability
	double rrp;
	//recombination probability
	double rp;
	//add drift velocity vector
	Vector3 v;
	//particles number
	int N;

private:
	double survive_probability(double R)
	{
		if (tau != INFINITY)
		{
			double r = R / sqrt(tau * D);
			double k = norm(v) * R / (2 * D);
			if (k == 0)
				return (mu * r / sinh(mu * r));
			else
				return (mu * r / sinh(mu * r)) * (sinh(k) / k);
		}
		else
		{
			return 1;
		}
	}
	Vector3 exit_point_on_sphere(double R)
	{
		double k = norm(v) * R / (2 * D);
		double cos_tetta;
		if (k != 0)
			cos_tetta = 1 + log(1 - (1 - exp(-2 * k)) * uniDistrib()) / k;
		else
			cos_tetta = 1 - 2 * uniDistrib();
		double sin_tetta = sqrt(1 - pow(cos_tetta, 2));
		double phi = 2 * Pi * uniDistrib();
		return Vector3(R * sin_tetta * cos(phi), R * sin_tetta * sin(phi), R * cos_tetta);
	}
	double first_passage_time(double R)
	{
		double L = sqrt(tau * D);
		double r = R / L;
		while (1)
		{			
			int n = 0;
			double S = 0;
			double U = uniDistrib();
			double c = (Pi * Pi) / (r * r) + 1 + (V * V * tau) / (4 * D);
			double x = expDistrib(1 / c);
			while (1)
			{
				n += 1;
				S += (n + 1) * (n + 1) * exp(-Pi * Pi * x * ((n + 1) * (n + 1) - 1) / (r * r));
				if (U >= S)
					return x;
				n += 1;
				S -= (n + 1) * (n + 1) * exp(-Pi * Pi * x * ((n + 1) * (n + 1) - 1) / (r * r));
				if (U < S)
					break;
			}
		}
	}
	double time_if_absorbed(double R)
	{
		double k = norm(v) * R / (2 * D);
		while (1)
		{
			int n = 0;
			double S = 0;
			double U = uniDistrib();
			double c = Pi * Pi * D / (R * R) + 1 / tau + V * V / (4 * D);
			double x = expDistrib(1 / c);
			while (1)
			{
				n += 1;
				S += (n + 1) * (n + 1) / (k * k + (n + 1) * (n + 1) * Pi * Pi) * (k * k + Pi * Pi) * exp(-Pi * Pi * D * x * ((n + 1) * (n + 1) - 1) / (R * R));
				if (U >= S)
					return x;
				n += 1;
				S -= (n + 1) * (n + 1) / (k * k + (n + 1) * (n + 1) * Pi * Pi) * (k * k + Pi * Pi) * exp(-Pi * Pi * D * x * ((n + 1) * (n + 1) - 1) / (R * R));
				if (U < S)
					break;
			}
		}
	}
	double time_if_absorbed_density(double R, double t)
	{
		double k = norm(v) * R / (2 * D);
		double F;
		if (k == 0)
			F = 1 - (mu * R / sqrt(D)) / sinh(mu * R / sqrt(D));
		else
			F = 1 - (mu * R / sqrt(D)) / sinh(mu * R / sqrt(D)) * sinh(k) / k;
		double sum = 0;
		for (int n = 1; n <= 1000; ++n)
		{
			double c = -(n * n * Pi * Pi * D / (R * R) + 1 / tau + norm(v) * norm(v) / (4 * D));
			double p = 2 * n * n * Pi * Pi / (k * k + n * n * Pi * Pi);
			sum += pow(-1, n + 1) * p * exp(c * t);
		}
		if (k == 0)
			return sum / F / tau;
		else
			return sum / F / tau * (sinh(k) / k);
	}
	void arrayInDensity(vector<double>& array, ofstream& file, double dx)
	{
		file.setf(ios::fixed);
		int size = array.size();
		sort(array.begin(), array.end());
		double delta = *array.begin() + dx;
		auto next = array.begin()++;
		auto cur = array.begin();
		while (cur != array.end())
		{
			int count = 0;
			while (*next <= delta && next != array.end())
			{
				++count;
				next++;
				cur++;
			}
			file << delta << "\t" << (double)count / size / dx << endl;
			delta += dx;
		}
	}
	Vector3 vectorSphereIntersection(Vector3& start, Vector3& direction, double R)
	{
		double r = norm(start);
		double cosTetta;
		if (r != 0)
			cosTetta = -scalProd(start, direction) / (r * norm(direction));
		else
			cosTetta = 0;
		double vectorLength = r * cosTetta + sqrt(r * r * (cosTetta * cosTetta - 1) + R * R);
		//cout << norm(start + direction * vectorLength) << endl;
		return start + direction* vectorLength;
	}
	double doubleLayerDiffusionEquationSolution(Calculation &calc1, Calculation &calc2, double layerHeight, double sourseHeight)
	{
		double D1 = calc1.D, D2 = calc2.D, L1 = sqrt(calc1.D * calc1.tau), L2 = sqrt(calc2.D * calc2.tau), H = calc1.H, h = layerHeight, z = sourseHeight;
		double k1 = D2 / L2 * (exp(h / L2) + exp((2 * H - h) / L2)) + D1 / L1 * (exp(h / L2) - exp((2 * H - h) / L2)) * (exp(-h / L1) + exp(h / L1)) / (exp(-h / L1) - exp(h / L1));
		double k2 = D1 / L1 * ((exp((H - h) / L2) - exp(h / L1)) * (exp(-h / L1) + exp(h / L1)) / (exp(-h / L1) - exp(h / L1)) - exp(h / L1)) - D2 / L2 * exp((H - h) / L2);
		double c = k2 / k1;
		double b = c * (exp(h / L2) - exp((2 * H - h) / L2)) / (exp(-h / L1) - exp(h / L1)) - (exp((H - h) / L2) - exp(h / L1)) / (exp(-h / L1) - exp(h / L1));
		double a = -b - 1;
		double d = -c * exp(2 * H / L2) - exp(H / L2);
		if (z < h)
			return a * exp(z / L1) + b * exp(-z / L1) + 1;
		else
			return c * exp(z / L2) + d * exp(-z / L2) + 1;
	}
	void guessFromDistributionGrid(vector<vector<double>>& grid, int& time, int& coord)
	{
		double guess = uniDistrib();
		double sum = 0;
		for (int i = 0; i < grid.size();++i)
			for (int j = 0; j < grid[0].size(); ++j)
			{
				sum += grid[i][j];
				if (sum >= guess)
				{
					time = i;
					coord = j;
					return;
				}
			}
	}
	void createDislocationArray(vector<Cylinder>& dislocations, double R, int disNum,double squareLen)
	{
		for (int i = 1; i <= disNum; ++i)
		{
			int maxTryCount = 0;
			while (1)
			{
				double x = uniDistrib() * squareLen - squareLen / 2;
				double y = uniDistrib() * squareLen - squareLen / 2;
				Cylinder cyl(x, y, R);
				if (!isCylindersIntersec(cyl, dislocations))
				{
					dislocations.push_back(cyl);
					break;
				}
				maxTryCount++;
				if (maxTryCount >= 1000)
				{
					cout << "I build only " << i - 1 << " cylinders\n";
					return;
				}
			}
		}
	}
public:
	void print()
	{
		cout << N << "\n" << V << "\n" << D << "\n" << tau << "\n" << epsilon << "\n" << H << "\n" << Pe << "\n" << mu << "\n";
	}
	void excitonBoundaryDensityFromNonCentralPoint(double R, double r, double tetta)
	{
		ofstream graph;
		//vector<double> times;
		vector<double> cosTettas;
		vector<double> phis;
		for (int i = 1; i <= N; ++i)
		{
			Vector3 u(r * sin(tetta), 0, r * cos(tetta));
			double t = 0;
			bool isAbsorbed = false;
			while (true)
			{
				double sphereRadius = R - norm(u);
				if (sphereRadius < epsilon)
					break;
				if (uniDistrib() >= survive_probability(sphereRadius))
				{
					isAbsorbed = true;
					break;
				}			
				u += exit_point_on_sphere(sphereRadius);
				//t += first_passage_time(sphereRadius);
			}
			if (!isAbsorbed)
			{
				//times.push_back(t);
				cosTettas.push_back(1 - u.getZ() / norm(u));
				phis.push_back(vectorsAzimut(u));
			}
			if (i % 10000 == 0)
				cout << i << endl;
		}
		//graph.open("times density from non-central point.txt");
		//arrayInDensity(times, graph, 0.01);
		//graph.close();
		graph.open("cos tetta density from non-central point.txt");
		arrayInDensity(cosTettas, graph, 0.01);
		graph.close();
		graph.open("phi density from non-central point.txt");
		arrayInDensity(phis, graph, 0.01);
		graph.close();
	}
	void readParametres()
	{
		ifstream file;
		file.open("parametres.txt");
		file >> N >> V >> D >> tau >> epsilon >> rp >> rrp >> H >> l;
		v.setZ(V);
		if (tau <= 0)
			tau = INFINITY;
		else
		{
			Pe = norm(v) * sqrt(tau) / sqrt(D);
			mu = sqrt(1 + Pe * Pe / 4);
		}
	}
	void laplasEquationSimaulationFromNonCentralPoint(double R, double r, double tetta)
	{
		ofstream graph;
		vector<double> cosTettas;
		vector<double> phis;
		for (int i = 1; i <= N; ++i)
		{
			Vector3 u(r * sin(tetta), 0, r * cos(tetta));
			Vector3 v = randUnitVector();
			Vector3 variant1 = vectorSphereIntersection(u, v, R);
			Vector3 v_ = -v;
			Vector3 variant2 = vectorSphereIntersection(u, v_, R);
			double a = norm(variant1 - u);
			double b = norm(variant2 - u);
			if (uniDistrib() < b / (a + b))
			{
				cosTettas.push_back(variant1.getZ() / norm(variant1));
				phis.push_back(vectorsAzimut(variant1));
			}
			else
			{
				cosTettas.push_back(variant2.getZ() / norm(variant2));
				phis.push_back(vectorsAzimut(variant2));
			}
			if (i % 10000 == 0)
				cout << i << endl;
		}
		graph.open("cos tetta density from non-central point(simplified).txt");
		arrayInDensity(cosTettas, graph, 0.01);
		graph.close();
		graph.open("phi density from non-central point(simplified).txt");
		arrayInDensity(phis, graph, 0.01);
		graph.close();
	}
	void twoLayersDiffusionCheck(double D1, double D2, double layerHeight)
	{
		ofstream calcGraph;
		ofstream analitGraph;
		calcGraph.open("calcGraph.txt");
		analitGraph.open("analitGraph.txt");
		double dz = H / 50;
		Calculation calc1 = *this, calc2 = *this;
		calc1.D = D1;
		calc2.D = D2;
		if (calc1.tau <= 0)
			calc1.tau = INFINITY;
		else
		{
			calc1.Pe = norm(calc1.v) * sqrt(calc1.tau) / sqrt(calc1.D);
			calc1.mu = sqrt(1 + calc1.Pe * calc1.Pe / 4);
		}
		if (calc2.tau <= 0)
			calc2.tau = INFINITY;
		else
		{
			calc2.Pe = norm(calc2.v) * sqrt(calc2.tau) / sqrt(calc2.D);
			calc2.mu = sqrt(1 + calc2.Pe * calc2.Pe / 4);
		}
		for (double sourseHeight = dz; sourseHeight <= H - dz; sourseHeight += dz)
		{
			int absorbed = 0;
			for (int i = 1; i <= N; ++i)
			{
				Vector3 u(0, 0, sourseHeight);
				bool counted = false;
				while (!counted)
				{
					if (abs(layerHeight - u.getZ()) <= epsilon)
					{
						if (uniDistrib() < D1 / (D1 + D2))
							u += Vector3(0, 0, -sqrt(epsilon));
						else
							u += Vector3(0, 0, sqrt(epsilon));
					}
					if (u.getZ() >= H - epsilon || u.getZ() <= epsilon)
						break;
					double R;
					//D1
					if (u.getZ() < layerHeight)
					{
						R = min(abs(u.getZ()), abs(layerHeight - u.getZ()));
						if (uniDistrib() > calc1.survive_probability(R))
						{
							absorbed++;
							break;
						}
						else
						{
							u += calc1.exit_point_on_sphere(R);
						}
					}
					//D2
					else
					{
						R = min(abs(H - u.getZ()), abs(layerHeight - u.getZ()));
						if (uniDistrib() > calc2.survive_probability(R))
						{
							absorbed++;
							break;
						}
						else
						{
							u += calc2.exit_point_on_sphere(R);
						}
					}
				}
				if (i %10000 == 0)
				cout << sourseHeight << "\t" << i << endl;
			}
			calcGraph << sourseHeight << "\t" << double(absorbed) / N << endl;
			analitGraph << sourseHeight << "\t" << doubleLayerDiffusionEquationSolution(calc1, calc2, layerHeight, sourseHeight) << endl;
		}
	}
	void timeIfAbsorbedDensityCheck(double R)
	{
		double dt = 0.00001;
		double tmax = 0.001;
		vector<int> times((int)(tmax / dt), 0);
		for (int i = 1; i <= N; ++i)
		{
			double guessTime = time_if_absorbed(R);
			if (guessTime < tmax)
				times[(int)(guessTime / dt)]++;
		}
		ofstream checkSim;
		checkSim.open("checkSim.txt");
		for (int i = 0; i < times.size(); ++i)
			checkSim << i * dt << "\t" << (double)times[i] / N / dt << endl;
		checkSim.close();


		ofstream checkAn;
		checkAn.open("checkAn.txt");
		for (int i = 0; i < (int)(tmax / dt); ++i)
			checkAn << i * dt << "\t" << time_if_absorbed_density(R, i * dt) << endl;
		checkAn.close();
	}
	void equationSystemWithRecombinationsCheck()
	{
		struct ParticleParametres
		{
			double z;
			double t;
			ParticleParametres(double z, double t)
			{
				this->z = z;
				this->t = t;
			}
		};
		vector<ParticleParametres> excitonsPos;
		vector<ParticleParametres> photonsPos;
		double h;
		int recombinationsNumber = 0;
		int stepNumber = 0;
		for (int i = 1; i <= N; ++i)
		{
			//???????? ?????????
			h =  uniDistrib()* (H - 2 * epsilon) + epsilon;
			Vector3 u(0, 0, h);
			bool counted = false;
			double t = expDistrib(0.05);
			while (!counted)
			{
			reverseRecombination:		//????? ???????? ????????????
				double R = min(abs(H - u.getZ()), abs(u.getZ()),H/20);

				//???????? ?? ?????????
				if (uniDistrib() > survive_probability(R))	//??????? ?? ??????
				{
					t += time_if_absorbed(R);
					excitonsPos.push_back(ParticleParametres(u.getZ(), t));
					//???????? ?? ???????????? ? ?????
					if (uniDistrib() <= rp)			//???????????? ? ?????
					{
						++recombinationsNumber;
						while (!counted)
						{
							u += randUnitVector() * expDistrib(l);					
							++stepNumber;
							if (u.getZ() >= H || u.getZ() <= 0)
							{
								counted = true;
							}
							if (!counted)
								photonsPos.push_back(ParticleParametres(u.getZ(), t));
							if (!counted && uniDistrib() < rrp)
							{
								++recombinationsNumber;
								goto reverseRecombination;
							}
						}
					}
					//?????????? ???????
					else
					{
						counted = true;
					}
				}
				else
				{
					//???????? ???????? ? ????? ?????
					u += exit_point_on_sphere(R);
					++stepNumber;
					t += first_passage_time(R);
					//???????????? ????? ?? ???????
					if (abs(H - u.getZ()) <= epsilon || abs(u.getZ()) <= epsilon)
					{
						counted = true;
					}
					if (!counted)
						excitonsPos.push_back(ParticleParametres(u.getZ(), t));
				}
			}
			if (i % 10000 == 0)
				cout << i << endl;
		}
		cout << (double)recombinationsNumber / N << endl;
		cout << (double)stepNumber / N << endl;
		cout << "start sorting\n";
		sort(excitonsPos.begin(), excitonsPos.end(), [](ParticleParametres &left, ParticleParametres &right)
			{
				if (left.t < right.t)
					return true;
				else if (left.t > right.t)
					return false;
				else return left.z < right.z;
			});
		sort(photonsPos.begin(), photonsPos.end(), [](ParticleParametres &left, ParticleParametres &right)
			{
				if (left.t < right.t)
					return true;
				else if (left.t > right.t)
					return false;
				else return left.z < right.z;
			});
		ofstream concentrationGrid;
		ofstream intensityGrid;
		concentrationGrid.open("concentrationGrid.txt");
		intensityGrid.open("intensityGrid.txt");
		double tmax = min(excitonsPos[(int)excitonsPos.size() * 0.95].t, photonsPos[(int)photonsPos.size() * 0.95].t);
		cout << tmax << endl;
		double dz = 1;
		double dt = 0.005;
		vector<vector<double>> excitonNet((int)(tmax / dt) + 1, vector<double>((int)(H / dz),0));
		auto cur = excitonsPos.begin();
		while (cur != excitonsPos.end())
		{
			if (cur->t > tmax)
				break;
			excitonNet[(int)(cur->t / dt)][(int)(cur->z / dz)]++;
			++cur;
		}
		for (int i = 0; i < excitonNet.size(); ++i)
			for (int j = 0; j < excitonNet[0].size(); ++j)
				excitonNet[i][j] = excitonNet[i][j] / N;
				//concentrationGrid << (j + 1) * dz << "\t" << (i + 1) * dt << "\t" << (double)excitonNet[i][j] / N / (dz * dt) << endl;
		vector<vector<double>> photonNet((int)(tmax / dt) + 1, vector<double>((int)(H / dz), 0));
		cur = photonsPos.begin();
		while (cur != photonsPos.end())
		{
			
			if (cur->t > tmax)
				break;
			photonNet[(int)(cur->t / dt)][(int)(cur->z / dz)]++;
			++cur;
		}
		
		for (int i = 0; i < photonNet.size(); ++i)
			for (int j = 0; j < photonNet[0].size(); ++j)
				photonNet[i][j] = photonNet[i][j] / N;
				//intensityGrid << (j + 1) * dz << "\t" << (i + 1) * dt << "\t" << (double)excitonNet[i][j] / N / (dz * dt) << endl;
		
#define n excitonNet
#define I photonNet
		vector<vector<double>> result((int)(tmax / dt) + 1, vector<double>((int)(H / dz), 0));
		for (int i = 1; i < result.size() - 1; ++i)
			for (int j = 1; j < result[0].size() - 1; ++j)
			{
				result[i][j] = (n[i + 1][j] - n[i - 1][j]) / (2 * dt) - D * (n[i][j + 1] - 2 * n[i][j] + n[i][j - 1]) / (dz * dz) + n[i][j] / tau - V * (n[i][j + 1] - n[i][j - 1]) / (2 * dz) -I[i][j] * rrp / tau;
			}
		ofstream resultFile;
		resultFile.open("equationRightPart.txt");
		for (int i = 0; i < result.size() - 1; ++i)
			for (int j = 3; j < result[0].size() - 3; ++j)
				resultFile << (j + 1) * dz << "\t" << (i + 1) * dt << "\t" << result[i][j] << endl;
		cout << "Done\n";
		cin.get();
	}
	void equationSystemWithRecombinationsMonteKarloCheck()
	{

		double h;
		int recombinationsNumber = 0;
		int stepNumber = 0;
		double sourceIntensity = 1;
		double dz = 2;
		double dt = 0.01;
		double tmax = 100 * dt;		
		vector<vector<double>> excitonNet((int)(tmax / dt) + 1, vector<double>((int)(H / dz), 0));
		vector<vector<double>> photonNet((int)(tmax / dt) + 1, vector<double>((int)(H / dz), 0));
		for (int i = 1; i <= N; ++i)
		{
			int lastTimeCell = 0, curTimeCell, lastCoordCell;
			//???????? ?????????
			h = uniDistrib() * (H - 2 * epsilon) + epsilon;
			lastCoordCell = (int)(h / dz);
			Vector3 u(0, 0, h);
			bool counted = false;
			double t = 0;
			while (!counted)
			{
			reverseRecombination:		//????? ???????? ????????????
				double R = min(abs(H - u.getZ()), abs(u.getZ()), H / 20);
				//???????? ?? ?????????
				if (uniDistrib() > survive_probability(R))	//??????? ?? ??????
				{
					t += time_if_absorbed(R);
					if (t < tmax)
						excitonNet[(int)(t / dt)][(int)(u.getZ() / dz)]++;
					//???????? ?? ???????????? ? ?????
					if (uniDistrib() <= rp)			//???????????? ? ?????
					{
						++recombinationsNumber;
						while (!counted)
						{
							u += randUnitVector() * expDistrib(l);
							++stepNumber;
							if (u.getZ() >= H || u.getZ() <= 0)
							{
								counted = true;
							}
							if (!counted && t < tmax)
								photonNet[(int)(t / dt)][(int)(u.getZ() / dz)] += 1.0;
							if (!counted && uniDistrib() < rrp)
							{
								++recombinationsNumber;
								goto reverseRecombination;
							}
						}
					}
					//?????????? ???????
					else
					{
						counted = true;
					}
				}
				else
				{
					//???????? ???????? ? ????? ?????
					u += exit_point_on_sphere(R);
					++stepNumber;
					t += first_passage_time(R);				
					//???????????? ????? ?? ???????
					if (abs(H - u.getZ()) <= epsilon || abs(u.getZ()) <= epsilon)
					{
						counted = true;
					}
					if (!counted && t < tmax)
					{
						excitonNet[(int)(t / dt)][(int)(u.getZ() / dz)]++;
						curTimeCell = (int)(t / dt);
						if (curTimeCell - lastTimeCell > 1)
						{
							for (int timeLost = lastTimeCell + 1; timeLost < curTimeCell; ++timeLost)
								excitonNet[timeLost][lastCoordCell]++;
						}
						lastTimeCell = (int)(t / dt);
						lastCoordCell = (int)(u.getZ() / dz);
					}
				}
			}
			if (i % 10000 == 0)
				cout << i << endl;
		}
		for (int i = 0; i < photonNet.size(); ++i)
		{
			int count = 0;
			for (int j = 0; j < photonNet[0].size(); ++j)
			{
				photonNet[i][j] = photonNet[i][j] / N / (dz * dt) * sourceIntensity;
				excitonNet[i][j] = excitonNet[i][j] / N / (dz * dt) * sourceIntensity;
				count += excitonNet[i][j];
			}
			cout << count << endl;
		}
		//?????????? ??????? ????????? ? ????????? ???????
		for (int i = 0; i < photonNet.size(); ++i)
			for (int j = 0; j < photonNet[0].size(); ++j)
				photonNet[i][j] = photonNet[i][j] * rrp;
		for (int j = 0; j < photonNet[0].size(); ++j)
			photonNet[0][j] += sourceIntensity;
		double gridSum = 0;
		for (int i = 0; i < photonNet.size(); ++i)
			for (int j = 0; j < photonNet[0].size(); ++j)
				gridSum += photonNet[i][j];
		for (int i = 0; i < photonNet.size(); ++i)
			for (int j = 0; j < photonNet[0].size(); ++j)
				photonNet[i][j] = photonNet[i][j] / gridSum;
		cout << gridSum << "\n\n";
		vector<double> walkerVector;
		for (int i = 0; i < photonNet.size(); ++i)
			for (int j = 0; j < photonNet[0].size(); ++j)
				walkerVector.push_back(photonNet[i][j]);
		int* alias_w = new int[walkerVector.size()]; // alias for Walker
		double* q_w = new double[walkerVector.size()]; // q = prob*p_cnt for Walker
		walker_prep(walkerVector, q_w, alias_w);
		cout << "Preparation ends\n";
		//?????? ??????? ????????? ?????? ???????? ??????????? ?????????? ???????
		vector<vector<double>> excitonNet2((int)(tmax / dt) + 1, vector<double>((int)(H / dz), 0));
		for (int i = 1; i <= N; ++i)
		{
			//???????? ?????????
			int time, coord, guess;
			int lastTimeCell = 0, curTimeCell, lastCoordCell;
			guess = walker(walkerVector.size(), q_w, alias_w);
			time = guess / (int)photonNet[0].size();
			coord = guess % (int)photonNet[0].size();
			h = coord * dz;
			if (h < epsilon)
				h = epsilon * 2;
			if (h > H - epsilon)
				h = H - epsilon * 2;
			Vector3 u(0, 0, h);
			lastCoordCell = (int)(h / dz);
			bool counted = false;
			double t = time * dt;
			while (!counted)
			{
				double R = min(abs(H - u.getZ()), abs(u.getZ()),H / 20);
				//???????? ?? ?????????
				if (uniDistrib() > survive_probability(R))	//??????? ?? ??????
				{
					t += time_if_absorbed(R);
					if (t < tmax)
						excitonNet2[(int)(t / dt)][(int)(u.getZ() / dz)]++;
					counted = true;
				}
				else
				{
					//???????? ???????? ? ????? ?????
					u += exit_point_on_sphere(R);
					++stepNumber;
					t += first_passage_time(R);
					//???????????? ????? ?? ???????
					if (abs(H - u.getZ()) <= epsilon || abs(u.getZ()) <= epsilon)
					{
						counted = true;
					}
					if (!counted && t < tmax)
					{
						excitonNet2[(int)(t / dt)][(int)(u.getZ() / dz)]++;
						curTimeCell = (int)(t / dt);
						if (curTimeCell - lastTimeCell > 1)
						{
							for (int timeLost = lastTimeCell + 1; timeLost < curTimeCell; ++timeLost)
								excitonNet2[timeLost][lastCoordCell]++;
						}
						lastTimeCell = (int)(t / dt);
						lastCoordCell = (int)(u.getZ() / dz);
					}
				}
			}
			if (i % 10000 == 0)
				cout << i << endl;
		}
		for (int i = 0; i < excitonNet2.size(); ++i)
			for (int j = 0; j < excitonNet2[0].size(); ++j)
				excitonNet2[i][j] = excitonNet2[i][j] / N / (dz * dt) * sourceIntensity;// *photonNet[i][j] * gridSum;
		ofstream excitonGrid, excitonGrid2;
		excitonGrid.open("excitonGrid.txt");
		excitonGrid2.open("excitonGrid2.txt");
		for (int i = 1; i < excitonNet.size() - 1; ++i)
			for (int j = 1; j < excitonNet[0].size() - 1; ++j)
			{
				excitonGrid << (j + 1) * dz << "\t" << (i + 1) * dt << "\t" << (double)excitonNet[i][j] << endl;
				excitonGrid2 << (j + 1) * dz << "\t" << (i + 1) * dt << "\t" << (double)excitonNet2[i][j] << endl;
			}
	}
	void equationSystemWithRecombinationsMonteKarloCheck2()
	{
		double dt = 0.01, h;
		double maxTime = 2;
		double sourceIntencity = 1;
		vector<double> excitonTimes((int)(maxTime / dt), 0);
		vector<double> photonTimes((int)(maxTime / dt), 0);
		vector<double> excitonsOnBoundaryTimes((int)(maxTime / dt), 0);
		for (int j = 1; j <= N; ++j)
		{
			//???????? ?????????
			h = uniDistrib() * H;
			Vector3 u(0, 0, h);
			bool counted = false;
			double t = uniDistrib() * 2;

			vector<pair<double, bool>> particlesTimes;//?????? ????? - ????? ???????????? ???????, ?????? - ?? ???(1 - ???????, 0 - ?????)

			while (!counted)
			{
			reverseRecombination1:		//????? ???????? ????????????
				particlesTimes.push_back(make_pair(t, 1));
				double R = min(abs(H - u.getZ()), abs(u.getZ()));

				//???????? ?? ?????????
				if (uniDistrib() > survive_probability(R))	//??????? ?? ??????
				{
					t += time_if_absorbed(R);
					while (1)
					{
						double t_abs = time_if_absorbed(R) / first_passage_time(R);
						if (t_abs < 1)
						{
							u += exit_point_on_sphere(t_abs * R);
							break;
						}
					}
					//???????? ?? ???????????? ? ?????
					if (uniDistrib() <= rp)			//???????????? ? ?????
					{
						while (!counted)
						{
							u += randUnitVector() * expDistrib(l);
							particlesTimes.push_back(make_pair(t, 0));
							if (u.getZ() >= H || u.getZ() <= 0)
							{
								counted = true;
							}
							if (!counted && uniDistrib() < rrp)
							{
								goto reverseRecombination1;
							}
						}
					}
					//?????????? ???????
					else
					{
						counted = true;
						particlesTimes.push_back(make_pair(t, 1));
					}
				}
				else
				{
					//???????? ???????? ? ????? ?????
					u += exit_point_on_sphere(R);
					particlesTimes.push_back(make_pair(t, 1));
					t += first_passage_time(R);
					//???????????? ????? ?? ???????
					if (abs(H - u.getZ()) <= epsilon || abs(u.getZ()) <= epsilon)
					{
						counted = true;
						if (t < maxTime)
							excitonsOnBoundaryTimes[(int)(t / dt)]++;
					}
				}
			}
			if (j % 10000 == 0)
				cout << j << endl;

			//????????? ??????? ??????
			if (!particlesTimes.empty())
			{
				vector<int> uniquePhotonTimes;
				for (int i = 0; i < particlesTimes.size(); ++i)
				{
					if (particlesTimes[i].second == 0)
					{
						int num = (int)(particlesTimes[i].first / dt);
						if (uniquePhotonTimes.size() == 0 || uniquePhotonTimes.back() != num)
							uniquePhotonTimes.push_back(num);
					}
				}
				int lastTime = min((int)(particlesTimes.back().first / dt), (int)(maxTime / dt) - 1);
				int firstTime = (int)(particlesTimes[0].first / dt);
				for (int i = firstTime; i <= lastTime; ++i)
					excitonTimes[i]++;
				for (int num : uniquePhotonTimes)
				{
					if (num < lastTime)
					{
						photonTimes[num]++;
						excitonTimes[num]--;
					}
				}
			}
		}
		for (int i = 1; i < excitonsOnBoundaryTimes.size(); ++i)
		{
			excitonsOnBoundaryTimes[i] += excitonsOnBoundaryTimes[i - 1];
			photonTimes[i] += photonTimes[i - 1];
		}
		cout << excitonTimes[0] << endl;
		cout << photonTimes[0] << endl;
		for (int num = 0; num < photonTimes.size(); ++num)
		{
			photonTimes[num] = photonTimes[num] / N * sourceIntencity;
			excitonTimes[num] = excitonTimes[num] / N * sourceIntencity;
			excitonsOnBoundaryTimes[num] = excitonsOnBoundaryTimes[num] / N * sourceIntencity;
		}
		vector<double> concDeriv(excitonTimes.size());
		for (int i = 1; i < excitonTimes.size() - 1; ++i)
			concDeriv[i] = (excitonTimes[i + 1] - excitonTimes[i - 1]) / (2 * dt);
		ofstream file;
		file.open("equationSystemCheckWithIntegration.txt");
		for (int i = 2; i < concDeriv.size() - 2; ++i)
			file << dt * i << "\t" << concDeriv[i]
			+ (excitonsOnBoundaryTimes[i + 1] - excitonsOnBoundaryTimes[i - 1]) / (2 * dt)
			+ excitonTimes[i] / tau
			- (photonTimes[i + 1] - photonTimes[i - 1]) / (2 * dt) << endl;
		file.close();
	}
	void doubleRecombinationSimulation(int recombinationNumberParam) //0 - ??? ????????????,1 - ??????????? ???? ??? ??????????????? ? ??????????????? ???????, 
		//2 - ????????????? ????????????
	{
		double h;
		vector<double> excitonTimes;
		vector<double> photonTimes;
		vector<double> absorbedTimes;
		vector<double> trajectoryEndTimes;
		int recombinationsCount = 0;
		for (int i = 1; i <= N; ++i)
		{
			//???????? ?????????
			h = H / 2;
			Vector3 u(0, 0, h);
			bool counted = false;
			bool reverseRecombinated = false;
			double t = 0;
			while (!counted)
			{
			reverseRecombination:		//????? ???????? ????????????
				double R = min(abs(H - u.getZ()), abs(u.getZ()));

				//???????? ?? ?????????
				if (uniDistrib() > survive_probability(R))	//??????? ?? ??????
				{
					t += time_if_absorbed(R);
					while (1)
					{
						double t_abs = time_if_absorbed(R) / first_passage_time(R);
						if (t_abs < 1)
						{
							u += exit_point_on_sphere(t_abs * R);
							break;
						}
					}
					//???????? ?? ???????????? ? ?????
					if (recombinationNumberParam != 0 && uniDistrib() <= rp && !(recombinationNumberParam == 1 && reverseRecombinated))			//???????????? ? ?????
					{
						recombinationsCount++;
						while (!counted)
						{
							u += randUnitVector() * expDistrib(l);
							if (u.getZ() >= H || u.getZ() <= 0)
							{
								counted = true;
								photonTimes.push_back(t);
								trajectoryEndTimes.push_back(t);
							}
							if (!counted && uniDistrib() < rrp)
							{
								recombinationsCount++;
								reverseRecombinated = true;
								goto reverseRecombination;
							}
						}
					}
					//?????????? ???????
					else
					{
						counted = true;
						absorbedTimes.push_back(t);
						trajectoryEndTimes.push_back(t);
					}
				}
				else
				{
					//???????? ???????? ? ????? ?????
					u += exit_point_on_sphere(R);
					t += first_passage_time(R);
					//???????????? ????? ?? ???????
					if (abs(H - u.getZ()) <= epsilon || abs(u.getZ()) <= epsilon)
					{
						counted = true;
						excitonTimes.push_back(t);
						trajectoryEndTimes.push_back(t);
					}
				}
			}
			if (i % 10000 == 0)
				cout << i << endl;
		}
		sort(excitonTimes.begin(), excitonTimes.end());
		sort(photonTimes.begin(), photonTimes.end());
		sort(absorbedTimes.begin(), absorbedTimes.end());
		sort(trajectoryEndTimes.begin(), trajectoryEndTimes.end());

		ofstream excitonIntensityOverTime;
		ofstream photonIntensityOverTime;

		string paramString;
		if (recombinationNumberParam == 0)
			paramString = "NoRecombinations";
		else if (recombinationNumberParam == 1)
			paramString = "OneDoubleRecombination";
		else if (recombinationNumberParam == 2)
			paramString = "MultyDoubleRecombinations";

		excitonIntensityOverTime.open("excitonsIntensityAtBoundaryOverTime" + paramString + ".txt");
		photonIntensityOverTime.open("photonsIntensityAtBoundaryOverTime" + paramString + ".txt");

		excitonTimes.erase(excitonTimes.begin() + (int)excitonTimes.size() * 0.98, excitonTimes.end());
		if (!photonTimes.empty())
			photonTimes.erase(photonTimes.begin() + (int)photonTimes.size() * 0.98, photonTimes.end());
		absorbedTimes.erase(absorbedTimes.begin() + (int)absorbedTimes.size() * 0.98, absorbedTimes.end());
		trajectoryEndTimes.erase(trajectoryEndTimes.begin() + (int)trajectoryEndTimes.size() * 0.98, trajectoryEndTimes.end());

		auto cur = excitonTimes.begin();
		double dt = excitonTimes.back() / 100;
		double t = dt;
		while (cur != excitonTimes.end())
		{
			int count = 0;
			while (*cur < t && cur != excitonTimes.end())
			{
				count++;
				cur++;
			}
			excitonIntensityOverTime << t << "\t" << (double)count / N << endl;
			t += dt;
		}
		cur = photonTimes.begin();
		t = dt;
		while (cur != photonTimes.end())
		{
			int count = 0;
			while (*cur < t && cur != photonTimes.end())
			{
				count++;
				cur++;
			}
			photonIntensityOverTime << t << "\t" << (double)count / N << endl;
			t += dt;
		}

		ofstream absorbedExcitonsOverTime;
		absorbedExcitonsOverTime.open("absorbedExcitonsOverTime" + paramString + ".txt");
		cur = absorbedTimes.begin();
		t = dt;
		int count = 0;
		while (cur != absorbedTimes.end())
		{	
			while (*cur < t && cur != absorbedTimes.end())
			{
				count++;
				cur++;
			}
			absorbedExcitonsOverTime << t << "\t" << (double)count / N << endl;
			t += dt;
		}

		ofstream aliveExcitonsOverTime;
		aliveExcitonsOverTime.open("aliveExcitonsOverTime" + paramString + ".txt");
		cur = trajectoryEndTimes.begin();
		t = dt;
		count = N;
		while (cur != trajectoryEndTimes.end())
		{
			
			while (*cur < t && cur != trajectoryEndTimes.end())
			{
				count--;
				cur++;
			}
			aliveExcitonsOverTime << t << "\t" << (double)count / N << endl;
			t += dt;
		}
		excitonIntensityOverTime.close();
		photonIntensityOverTime.close();
		absorbedExcitonsOverTime.close();
		aliveExcitonsOverTime.close();
		cout << (double)recombinationsCount / N << endl;
	}
	void doubleRecombinationWithCylindricalDislocationsSimulation()
	{
		//??????? ?????????????? ??????????(??????? ??????? ?? ???? ???? ? ???????????? ?????? (x,y) ? ???????? 
		double cylX = 0, cylY = 0.1, cylR = 0.01;
		Cylinder cyl(cylX, cylY, cylR);
		Cylinder epscyl(cylX, cylY, cylR + epsilon);
		double h;
		vector<double> excitonTimes;
		vector<double> photonTimes;
		for (int i = 1; i <= N; ++i)
		{
			//???????? ?????????
			h = H / 2;
			Vector3 u(0, 0, h);
			bool counted = false;
			double t = 0;
			while (!counted)
			{
			reverseRecombination1:		//????? ???????? ????????????
				double R = min(abs(H - u.getZ()), abs(u.getZ()), abs(cylinderDistance(u,cyl)),0.05);

				//???????? ?? ?????????
				if (uniDistrib() > survive_probability(R))	//??????? ?? ??????
				{
					t += time_if_absorbed(R);
					//???????? ?? ???????????? ? ?????
					if (uniDistrib() <= rp)			//???????????? ? ?????
					{
						while (!counted)
						{
							u += randUnitVector() * expDistrib(l);
							if (u.getZ() >= H || u.getZ() <= 0)
							{
								counted = true;
								photonTimes.push_back(t);
							}
							if (isInsideCylinder(u, epscyl))
							{
								counted = true;
							}
							if (!counted && uniDistrib() < rrp)
							{
								goto reverseRecombination1;
							}
						}
					}
					//?????????? ???????
					else
					{
						counted = true;
					}
				}
				else
				{
					//???????? ???????? ? ????? ?????
					u += exit_point_on_sphere(R);
					t += first_passage_time(R);
					//???????????? ????? ?? ???????
					if (abs(H - u.getZ()) <= epsilon || abs(u.getZ()) <= epsilon)
					{
						counted = true;
						excitonTimes.push_back(t);
					}
					if (isInsideCylinder(u,epscyl))
						counted = true;
				}
			}
			if (i % 10000 == 0)
				cout << i << endl;
		}
		sort(excitonTimes.begin(), excitonTimes.end());
		sort(photonTimes.begin(), photonTimes.end());
		ofstream excitonIntensityOverTime;
		ofstream photonIntensityOverTime;
		excitonIntensityOverTime.open("excitonsIntensityAtBoundaryOverTimeWithCylinderDislocation.txt");
		photonIntensityOverTime.open("photonsIntensityAtBoundaryOverTimeWithCylinderDislocation.txt");

		excitonTimes.erase(excitonTimes.begin() + (int)excitonTimes.size() * 0.98, excitonTimes.end());
		photonTimes.erase(photonTimes.begin() + (int)photonTimes.size() * 0.98, photonTimes.end());

		auto cur = excitonTimes.begin();
		double dt = min(excitonTimes.back(), photonTimes.back()) / 100;
		double t = dt;
		while (cur != excitonTimes.end())
		{
			int count = 0;
			while (*cur < t && cur != excitonTimes.end())
			{
				count++;
				cur++;
			}
			excitonIntensityOverTime << t << "\t" << (double)count / N << endl;
			t += dt;
		}
		cur = photonTimes.begin();
		t = dt;
		while (cur != photonTimes.end())
		{
			int count = 0;
			while (*cur < t && cur != photonTimes.end())
			{
				count++;
				cur++;
			}
			photonIntensityOverTime << t << "\t" << (double)count / N << endl;
			t += dt;
		}
		excitonIntensityOverTime.close();
		photonIntensityOverTime.close();
	}
	void doubleRecombinationWithCylindricalDislocationHitMap()
	{
		//??????? ?????????????? ??????????(??????? ??????? ?? ???? ???? ? ???????????? ?????? (x,y) ? ???????? 
		double cylX = 3, cylY = 4, cylR = 2;
		Cylinder cyl(cylX, cylY, cylR);
		Cylinder epscyl(cylX, cylY, cylR + epsilon);
		double h;
		double dx = 0.1, dy = 0.1;
		int netSize = 80;
		vector<vector<double>> excitonNet(netSize*2, vector<double>(netSize * 2, 0));
		vector<vector<double>> photonNet(netSize * 2, vector<double>(netSize * 2, 0));
		int WTFCOUNT = 0;
		for (int i = 1; i <= N; ++i)
		{
			//???????? ?????????
			h = 7 * H / 8;
			Vector3 u(0, 0, h);
			bool counted = false;
			Vector3 lastPhotonPos;
			while (!counted)
			{
			reverseRecombination1:		//????? ???????? ????????????
				double R = min(abs(H - u.getZ()), abs(u.getZ()), abs(cylinderDistance(u, cyl)), H / 20);

				//???????? ?? ?????????
				if (uniDistrib() > survive_probability(R))	//??????? ?? ??????
				{
					//???????? ?? ???????????? ? ?????
					if (uniDistrib() <= rp)			//???????????? ? ?????
					{
						while (!counted)
						{
							lastPhotonPos = u;
							u += randUnitVector() * expDistrib(l);
							if (u.getZ() >= H)
							{
								counted = true;
								double t = (H - lastPhotonPos.getZ()) / (u.getZ() - lastPhotonPos.getZ());
								u.setX(lastPhotonPos.getX() + t * (u.getX() - lastPhotonPos.getX()));
								u.setY(lastPhotonPos.getY() + t * (u.getY() - lastPhotonPos.getY()));
								if (u.getX() > -netSize*dx && u.getX() < netSize * dx && u.getY() > -netSize * dy && u.getY() < netSize * dy)
									photonNet[(int)(u.getX() / dx) + netSize][(int)(u.getY() / dy) + netSize]++;
							}
							if (u.getZ() <= 0)
								counted = true;
							if (isInsideCylinder(u, epscyl))
							{
								counted = true;
							}
							if (!counted && uniDistrib() < rrp)
							{
								goto reverseRecombination1;
							}
						}
					}
					//?????????? ???????
					else
					{
						counted = true;
					}
				}
				else
				{
					//???????? ???????? ? ????? ?????
					u += exit_point_on_sphere(R);
					//???????????? ????? ?? ???????
					if (abs(H - u.getZ()) <= epsilon)
					{
						counted = true;
						if (u.getX() > -netSize * dx && u.getX() < netSize * dx && u.getY() > -netSize * dy && u.getY() < netSize * dy)
						{
							excitonNet[(int)(u.getX() / dx) + netSize][(int)(u.getY() / dy) + netSize]++;
							WTFCOUNT++;
						}
					}
					if (abs(u.getZ()) <= epsilon)
						counted = true;
					if (isInsideCylinder(u, epscyl))
						counted = true;
				}
			}
			if (i % 10000 == 0)
				cout << i << endl;
		}
		ofstream excitonIntensity;
		ofstream photonIntensity;
		excitonIntensity.open("excitonsIntensityAtBoundaryHitMap.txt");
		photonIntensity.open("photonsIntensityAtBoundaryHitMap.txt");

		for (int i = 0; i < excitonNet.size(); ++i)
			for (int j = 0; j < excitonNet[0].size(); ++j)
				excitonIntensity << (i - netSize) * dx << "\t" << (j - netSize) * dy << "\t" << excitonNet[i][j] << endl;
		for (int i = 0; i < photonNet.size(); ++i)
			for (int j = 0; j < photonNet[0].size(); ++j)
				photonIntensity << (i - netSize) * dx << "\t" << (j - netSize) * dy << "\t" << photonNet[i][j] << endl;
		excitonIntensity.close();
		photonIntensity.close();
		cout << WTFCOUNT << endl;
	}
	void doubleRecombinationWithMultipleDislocations(double cylR,int dislocationNum)
	{
		//???????? ??????? ??????????
		vector<Cylinder> dislocations;
		vector<Cylinder> epsDislocations;

		int sqrtNum = (int)sqrt(dislocationNum);
		double squareLen = 3 * sqrt(D * tau);
		double d = (squareLen - 2 * cylR * sqrtNum) / (sqrtNum + 1);

		for (int i = 1; i <= sqrtNum; ++i)
			for (int j = 1; j <= sqrtNum; ++j)
			{
				dislocations.push_back(Cylinder(i * d + (2 * i - 1) * cylR - squareLen / 2, j * d + (2 * j - 1) * cylR - squareLen / 2, cylR));
			}
		//createDislocationArray(dislocations, cylR, dislocationNum, squareLen);
		for (int i = 0; i < dislocations.size(); ++i)
		{
			epsDislocations.push_back(Cylinder(dislocations[i].centerX, dislocations[i].centerY, cylR + epsilon));
		}

		double dt = 0.01;
		double maxTime = 2;
		vector<double> excitonTimes((int)(maxTime / dt), 0);
		vector<double> photonTimes((int)(maxTime / dt),0);
		for (int i = 1; i <= N; ++i)
		{
			//???????? ?????????
			Vector3 u(uniDistrib()*squareLen/2 - squareLen/4, uniDistrib() * squareLen / 2 - squareLen / 4, uniDistrib()*H);
			bool counted = false;
			if (isInsideCylinders(u, epsDislocations))
				counted = true;
			double t = 0;
			
			vector<pair<double, bool>> particlesTimes;//?????? ????? - ????? ???????????? ???????, ?????? - ?? ???(1 - ???????, 0 - ?????)
			
			while (!counted)
			{
			reverseRecombination1:		//????? ???????? ????????????
				particlesTimes.push_back(make_pair(t, 1));
				double R = min(abs(H - u.getZ()), abs(u.getZ()), minDistance(u, dislocations), H / 20);

				//???????? ?? ?????????
				if (uniDistrib() > survive_probability(R))	//??????? ?? ??????
				{
					t += time_if_absorbed(R);
					//???????? ?? ???????????? ? ?????
					if (uniDistrib() <= rp)			//???????????? ? ?????
					{
						while (!counted)
						{
							u += randUnitVector() * expDistrib(l);
							particlesTimes.push_back(make_pair(t, 0));
							if (u.getZ() >= H || u.getZ() <= 0)
							{
								counted = true;
							}
							if (!counted && uniDistrib() < rrp)
							{
								goto reverseRecombination1;
							}
						}
					}
					//?????????? ???????
					else
					{
						counted = true;
					}
				}
				else
				{
					//???????? ???????? ? ????? ?????
					u += exit_point_on_sphere(R);
					particlesTimes.push_back(make_pair(t, 1));
					t += first_passage_time(R);
					//???????????? ????? ?? ???????
					if (abs(H - u.getZ()) <= epsilon || abs(u.getZ()) <= epsilon)
					{
						counted = true;					
					}
					if (isInsideCylinders(u, epsDislocations))
						counted = true;
				}
			}
			if (i % 10000 == 0)
				cout << i << endl;

			//????????? ??????? ??????
			if (!particlesTimes.empty())
			{
				vector<int> uniquePhotonTimes;
				for (int i = 0; i < particlesTimes.size(); ++i)
				{
					if (particlesTimes[i].second == 0)
					{
						int num = (int)(particlesTimes[i].first / dt);
						if (uniquePhotonTimes.size() == 0 || uniquePhotonTimes.back() != num)
							uniquePhotonTimes.push_back(num);
					}
				}
				int lastTime = min((int)(particlesTimes.back().first / dt), (int)(maxTime / dt) - 1);
				for (int i = 0; i < lastTime; ++i)
					excitonTimes[i]++;
				for (int num : uniquePhotonTimes)
				{
					if (num < lastTime)
					{
						photonTimes[num]++;
						excitonTimes[num]--;
					}
				}
			}
		}
		ofstream excitonNumberOverTime;
		ofstream photonNumberOverTime;
		excitonNumberOverTime.open("excitonsNumberAtBoundaryOverTimeWithMultipleDislocations.txt");
		photonNumberOverTime.open("photonsNumberAtBoundaryOverTimeWithMultipleDislocations.txt");
		for (int i = 0; i < excitonTimes.size(); ++i)
		{
			excitonNumberOverTime << i * dt << "\t" << excitonTimes[i] / N << endl;
		}
		for (int i = 0; i < photonTimes.size(); ++i)
		{
			photonNumberOverTime << i * dt << "\t" << photonTimes[i] / N << endl;
		}
		
		excitonNumberOverTime.close();
		photonNumberOverTime.close();
	}
};
void main()
{
	Calculation calc;
	calc.readParametres();
	double R, r, tetta;
	ifstream conditions;
	conditions.open("initial conditions.txt");
	conditions >> R >> r >> tetta;
	conditions.close();
	//calc.timeIfAbsorbedDensityCheck(100);
	//calc.equationSystemWithRecombinationsMonteKarloCheck2();
	//calc.doubleRecombinationWithMultipleDislocations(5, 50);
	//calc.doubleRecombinationWithCylindricalDislocationHitMap();
	//calc.equationSystemWithRecombinationsMonteKarloCheck();
	calc.doubleRecombinationSimulation(0);
	calc.doubleRecombinationSimulation(1);
	calc.doubleRecombinationSimulation(2);
	//calc.doubleRecombinationWithCylindricalDislocationsSimaulation();
	//calc.equationSystemWithRecombinationsCheck();
	//calc.twoLayersDiffusionCheck(10, 10000, 20);
	//calc.excitonBoundaryDensityFromNonCentralPoint(R , r, tetta);
	//calc.laplasEquationSimaulationFromNonCentralPoint(R, r, tetta);
	cout << "ready\n";
	cin.get();
}