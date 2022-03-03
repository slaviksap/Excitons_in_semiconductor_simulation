#include<iostream>
#include<fstream>
#include<cmath>
#include"vector3.h"
#include<cstdlib>
#include"Functions.h"
#include <inttypes.h>
using namespace std;

//walker
//----------------------fast_drand40()----------------------------------�
static double d40d = (double)(1LL << 40);
static uint64_t d40m1 = (1LL << 40) - 1LL;
uint64_t prev1 = 12LL;

double randoms()
{
	prev1 = (prev1 * 762939453125LL) & d40m1;
	return (double)prev1 / d40d;
}
//������� ���������� �������� ����������� q � ����� alias ��� Walker
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
//������� Walker ���������� ��������� ������, ��������� � ������������ � ������������ ���������� :
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
			if (norm(v) == 0)
				return r / sinh(r);
			double k = norm(v) * R / (2 * D);
			return (mu * r / sinh(mu * r)) * sinh(k) / k;
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
		while (1)
		{
			int n = 0;
			double S = 0;
			double U = uniDistrib();
			double c = (Pi * Pi) / (R * R) + 1 + (V * V * tau) / (4 * D);
			double x = expDistrib(1 / c);
			while (1)
			{
				n += 1;
				S += (n + 1) * (n + 1) * exp(-Pi * Pi * x * ((n + 1) * (n + 1) - 1) / (R * R));
				if (U >= S)
					return x;
				n += 1;
				S -= (n + 1) * (n + 1) * exp(-Pi * Pi * x * ((n + 1) * (n + 1) - 1) / (R * R));
				if (U < S)
					break;
			}
		}
	}
	double time_if_absorbed(double R)
	{
		double k = R / 2 * Pe;
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
		double dz = 0.05;
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
			//�������� ���������
			h =  uniDistrib()* (H - 2 * epsilon) + epsilon;
			Vector3 u(0, 0, h);
			bool counted = false;
			double t = expDistrib(0.05);
			while (!counted)
			{
			reverseRecombination:		//����� �������� ������������
				double R = min(abs(H - u.getZ()), abs(u.getZ()),0.05);

				//�������� �� ���������
				if (uniDistrib() > survive_probability(R))	//������� �� ������
				{
					t += time_if_absorbed(R);
					excitonsPos.push_back(ParticleParametres(u.getZ(), t));
					//�������� �� ������������ � �����
					if (uniDistrib() <= rp)			//������������ � �����
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
					//���������� �������
					else
					{
						counted = true;
					}
				}
				else
				{
					//�������� �������� � ����� �����
					u += exit_point_on_sphere(R);
					++stepNumber;
					t += first_passage_time(R);
					//������������ ����� �� �������
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
		double dz = 0.05;
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
		double dz = 0.05;
		double dt = 0.005;
		double tmax = 50 * dt;		
		vector<vector<double>> excitonNet((int)(tmax / dt) + 1, vector<double>((int)(H / dz), 0));
		vector<vector<double>> photonNet((int)(tmax / dt) + 1, vector<double>((int)(H / dz), 0));
		for (int i = 1; i <= N; ++i)
		{
			int lastTimeCell = 0, curTimeCell, lastCoordCell;
			//�������� ���������
			h = uniDistrib() * (H - 2 * epsilon) + epsilon;
			lastCoordCell = (int)(h / dz);
			Vector3 u(0, 0, h);
			bool counted = false;
			double t = 0;
			while (!counted)
			{
			reverseRecombination:		//����� �������� ������������
				double R = min(abs(H - u.getZ()), abs(u.getZ()),0.05);
				//�������� �� ���������
				if (uniDistrib() > survive_probability(R))	//������� �� ������
				{
					t += time_if_absorbed(R);
					if (t < tmax)
						excitonNet[(int)(t / dt)][(int)(u.getZ() / dz)]++;
					//�������� �� ������������ � �����
					if (uniDistrib() <= rp)			//������������ � �����
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
								photonNet[(int)(t / dt)][(int)(u.getZ() / dz)] += 1.0 / 10;
							if (!counted && uniDistrib() < rrp)
							{
								++recombinationsNumber;
								goto reverseRecombination;
							}
						}
					}
					//���������� �������
					else
					{
						counted = true;
					}
				}
				else
				{
					//�������� �������� � ����� �����
					u += exit_point_on_sphere(R);
					++stepNumber;
					t += first_passage_time(R);				
					//������������ ����� �� �������
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
		//���������� ������� ��������� � ��������� �������
		for (int i = 0; i < photonNet.size(); ++i)
			for (int j = 0; j < photonNet[0].size(); ++j)
				photonNet[i][j] = photonNet[i][j] * rrp;
		for (int j = 0; j < photonNet[0].size(); ++j)
			photonNet[0][j] += H/sourceIntensity;
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
		int* alias_w = new int[N]; // alias for Walker
		double* q_w = new double[N]; // q = prob*p_cnt for Walker
		walker_prep(walkerVector, q_w, alias_w);
		cout << "Preparation ends\n";
		//������ ������� ��������� ������ �������� ����������� ���������� �������
		vector<vector<double>> excitonNet2((int)(tmax / dt) + 1, vector<double>((int)(H / dz), 0));
		for (int i = 1; i <= N; ++i)
		{
			//�������� ���������
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
				double R = min(abs(H - u.getZ()), abs(u.getZ()),0.05);
				//�������� �� ���������
				if (uniDistrib() > survive_probability(R))	//������� �� ������
				{
					t += time_if_absorbed(R);
					if (t < tmax)
						excitonNet2[(int)(t / dt)][(int)(u.getZ() / dz)]++;
					counted = true;
				}
				else
				{
					//�������� �������� � ����� �����
					u += exit_point_on_sphere(R);
					++stepNumber;
					t += first_passage_time(R);
					//������������ ����� �� �������
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
				excitonNet2[i][j] = excitonNet2[i][j] / N / (dz * dt) * photonNet[i][j] * gridSum;
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
	void doubleRecombinationSimulation()
	{
		double h;
		vector<double> excitonTimes;
		vector<double> photonTimes;
		for (int i = 1; i <= N; ++i)
		{
			//�������� ���������
			h = H / 2;
			Vector3 u(0, 0, h);
			bool counted = false;
			double t = 0;
			while (!counted)
			{
			reverseRecombination:		//����� �������� ������������
				double R = min(abs(H - u.getZ()), abs(u.getZ()),0.05);

				//�������� �� ���������
				if (uniDistrib() > survive_probability(R))	//������� �� ������
				{
					t += time_if_absorbed(R);
					//�������� �� ������������ � �����
					if (uniDistrib() <= rp)			//������������ � �����
					{
						while (!counted)
						{
							u += randUnitVector() * expDistrib(l);
							if (u.getZ() >= H || u.getZ() <= 0)
							{
								counted = true;
								photonTimes.push_back(t);
							}
							if (!counted && uniDistrib() < rrp)
							{
								goto reverseRecombination;
							}
						}
					}
					//���������� �������
					else
					{
						counted = true;
					}
				}
				else
				{
					//�������� �������� � ����� �����
					u += exit_point_on_sphere(R);
					t += first_passage_time(R);
					//������������ ����� �� �������
					if (abs(H - u.getZ()) <= epsilon || abs(u.getZ()) <= epsilon)
					{
						counted = true;
						excitonTimes.push_back(t);
					}
				}
			}
			if (i % 10000 == 0)
				cout << i << endl;
		}
		sort(excitonTimes.begin(), excitonTimes.end());
		sort(photonTimes.begin(), photonTimes.end());
		ofstream excitonIntensityOverTime;
		ofstream photonIntensityOverTime;
		excitonIntensityOverTime.open("excitonsIntensityAtBoundaryOverTime.txt");
		photonIntensityOverTime.open("photonsIntensityAtBoundaryOverTime.txt");

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
	void doubleRecombinationWithCylindricalDislocationsSimaulation()
	{
		//������� �������������� ����������(������� ������� �� ���� ���� � ������������ ������ (x,y) � �������� 
		double cylX = 0, cylY = 0.1, cylR = 0.01;
		double h;
		vector<double> excitonTimes;
		vector<double> photonTimes;
		for (int i = 1; i <= N; ++i)
		{
			//�������� ���������
			h = H / 2;
			Vector3 u(0, 0, h);
			bool counted = false;
			double t = 0;
			while (!counted)
			{
			reverseRecombination1:		//����� �������� ������������
				double R = min(abs(H - u.getZ()), abs(u.getZ()), abs(cylinderDistance(u, cylX, cylY, cylR)),0.05);

				//�������� �� ���������
				if (uniDistrib() > survive_probability(R))	//������� �� ������
				{
					t += time_if_absorbed(R);
					//�������� �� ������������ � �����
					if (uniDistrib() <= rp)			//������������ � �����
					{
						while (!counted)
						{
							u += randUnitVector() * expDistrib(l);
							if (u.getZ() >= H || u.getZ() <= 0)
							{
								counted = true;
								photonTimes.push_back(t);
							}
							if (isInsideCylinder(u, cylX, cylY, cylR + epsilon))
							{
								counted = true;
							}
							if (!counted && uniDistrib() < rrp)
							{
								goto reverseRecombination1;
							}
						}
					}
					//���������� �������
					else
					{
						counted = true;
					}
				}
				else
				{
					//�������� �������� � ����� �����
					u += exit_point_on_sphere(R);
					t += first_passage_time(R);
					//������������ ����� �� �������
					if (abs(H - u.getZ()) <= epsilon || abs(u.getZ()) <= epsilon)
					{
						counted = true;
						excitonTimes.push_back(t);
					}
					if (isInsideCylinder(u, cylX, cylY, cylR + epsilon))
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
		//������� �������������� ����������(������� ������� �� ���� ���� � ������������ ������ (x,y) � �������� 
		double cylX = 0.05, cylY = 0.1, cylR = 0.07;
		double h;
		double dx = 0.01, dy = 0.01;
		int netSize = 30;
		vector<vector<double>> excitonNet(netSize*2, vector<double>(netSize * 2, 0));
		vector<vector<double>> photonNet(netSize * 2, vector<double>(netSize * 2, 0));
		for (int i = 1; i <= N; ++i)
		{
			//�������� ���������
			h = 3 * H / 4;
			Vector3 u(0, 0, h);
			bool counted = false;
			while (!counted)
			{
			reverseRecombination1:		//����� �������� ������������
				double R = min(abs(H - u.getZ()), abs(u.getZ()), abs(cylinderDistance(u, cylX, cylY, cylR)), 0.05);

				//�������� �� ���������
				if (uniDistrib() > survive_probability(R))	//������� �� ������
				{
					//�������� �� ������������ � �����
					if (uniDistrib() <= rp)			//������������ � �����
					{
						while (!counted)
						{
							u += randUnitVector() * expDistrib(l);
							if (u.getZ() >= H)
							{
								counted = true;
								if (u.getX() > -netSize*dx && u.getX() < netSize * dx && u.getY() > -netSize * dy && u.getY() < netSize * dy)
									photonNet[(int)(u.getX() / dx) + netSize][(int)(u.getY() / dy) + netSize]++;
							}
							if (u.getZ() <= 0)
								counted = true;
							if (isInsideCylinder(u, cylX, cylY, cylR + epsilon))
							{
								counted = true;
							}
							if (!counted && uniDistrib() < rrp)
							{
								goto reverseRecombination1;
							}
						}
					}
					//���������� �������
					else
					{
						counted = true;
					}
				}
				else
				{
					//�������� �������� � ����� �����
					u += exit_point_on_sphere(R);
					//������������ ����� �� �������
					if (abs(H - u.getZ()) <= epsilon)
					{
						counted = true;
						if (u.getX() > -netSize * dx && u.getX() < netSize * dx && u.getY() > -netSize * dy && u.getY() < netSize * dy)
							excitonNet[(int)(u.getX() / dx) + netSize][(int)(u.getY() / dy) + netSize]++;
					}
					if (abs(u.getZ()) <= epsilon)
						counted = true;
					if (isInsideCylinder(u, cylX, cylY, cylR + epsilon))
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
				excitonIntensity << (i - 30) * dx << "\t" << (j - 30) * dy << "\t" << excitonNet[i][j] / N << endl;
		for (int i = 0; i < photonNet.size(); ++i)
			for (int j = 0; j < photonNet[0].size(); ++j)
				photonIntensity << (i - 30) * dx << "\t" << (j - 30) * dy << "\t" << photonNet[i][j] / N << endl;
		excitonIntensity.close();
		photonIntensity.close();
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
	calc.doubleRecombinationWithCylindricalDislocationHitMap();
	//calc.equationSystemWithRecombinationsMonteKarloCheck();
	//calc.doubleRecombinationSimulation();
	//calc.doubleRecombinationWithCylindricalDislocationsSimaulation();
	//calc.equationSystemWithRecombinationsCheck();
	//calc.twoLayersDiffusionCheck(1, 0.5, 0.5);
	//calc.excitonBoundaryDensityFromNonCentralPoint(R , r, tetta);
	//calc.laplasEquationSimaulationFromNonCentralPoint(R, r, tetta);
}