// Baseline algorithm for haplotype inference
// By Zeeshan Pirzada


// inputs -> genotype vectors (0, 1, 2)
// outputs -> haplotype vectors (0,1)

#include <vector>
#include <iostream>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <algorithm>
#include <set>

using namespace std;

/* 
   helper function for getHap()
   converts number to binary vector
   todo: this is inefficient use of space,
   make better
*/

vector<int> convertToBin(int a, int width)
{
	vector<int> ret;

	for (int i = 0; i < width; i++)
	{
		if (a % 2 == 0)
			ret.push_back(0);
		else
			ret.push_back(1);

		a = a/2;
	}

	reverse(ret.begin(), ret.end());
	return ret;
}


/* 
   returns all possible haplotypes of width m
   ie, if m = 3 returns:
   000, 001, 010, 011, 100, 101, 110, 111
*/

vector<vector<int> > getHaplotypes(int m)
{
	int num = pow(2, m);

	vector<vector<int> > ret;

	for (int i = 0; i < num; i++)
	{
		vector<int> toAdd = convertToBin(i, m);
		ret.push_back(toAdd);
	}

	return ret;
}

/*
	Given a set of haplotypes, determine set of
	semantically valid solutions.
	Basically, just calculate powerset.
*/

vector<vector<vector<int> > > getSolutions(vector<vector<int> > haplotypes)
{
	vector<vector<vector<int> > > ret;

	// insert empty set;
	vector<vector<int> > empty;
	ret.push_back(empty);

	// iterate through haplotypes;
	
	int s = haplotypes.size();
	for (int i = 0; i < s; i++)
	{
		int num_so_far = ret.size();

		for (int j = 0; j < num_so_far; j++)
		{
			vector<vector<int> > temp = ret[j];
			temp.push_back(haplotypes[i]);
			ret.push_back(temp);
		}

	}

	return ret;
}

/*
	determines if h1 and h2 
	explain g
*/

bool explains(vector<int> g, vector<int> h1, vector<int> h2)
{
	for (int i = 0; i < g.size(); i++)
	{
		if (g[i] == 0)
		{
			if (h1[i] != 0 || h2[i] != 0)
				return false;
		}

		if (g[i] == 1)
		{
			if (h1[i] != 1 || h2[i] != 1)
				return false;
		}

		if (g[i] == 2) // ambiguous
		{
			// g_i = 2 so the haplotype has to be 01 or 10

			if (!(h1[i] ^ h2[i]))
				return false;
		}
	}

	return true;

}

/*
	checks if haplotypes explains
	genotypes
	todo: make this not O(N^2)
*/

bool isValidSolution(vector<vector<int> > genotypes, vector<vector<int> > haplotypes)
{
	

	for (int i = 0; i < genotypes.size(); i++)
	{

		bool explained = false;

		// find haplotypes that match current genotype
		// enumerate all pairs
	
		for (int j = 0; j < haplotypes.size(); j++)
		{
			for (int k = j; k < haplotypes.size(); k++)
			{
				if (explains(genotypes[i], haplotypes[j], haplotypes[k]))
				{
					explained = true;
				}
			}
		}

		if (!explained)
			return false;

	}

	return true;
}

/*
	Baseline method: iterate all possible solutions
	and pick optimal. 
	INCREDIBLY INEFFICIENT: O(2^(2MN))
*/

vector<vector<int> > baseline(vector<vector<int> > genotypes)
{
	int min_so_far = INT_MAX;
	
	vector<vector<int> > ret;

	vector<vector<int> > haps = getHaplotypes(genotypes[0].size());

	vector<vector<vector<int> > > solutions = getSolutions(haps);

	for (int i = 0; i < solutions.size(); i++)
	{
		if (isValidSolution(genotypes, solutions[i]) && solutions[i].size() < min_so_far)
		{
			min_so_far = solutions[i].size();
			ret = solutions[i];
		}

	}

	return ret;
}

/*
	See if genotype is ambigious
	ie, does it have a 2

*/

bool ambigious(vector<int> g)
{
	for (int i = 0; i < g.size(); i++)
	{
		if (g[i] == 2)
			return true;
	}
	return false;
}

/*	
	Can haplotype infer a pair from
	genotype. store inferred pairing in
	result and return true

*/

bool infer(vector<int> genotype, vector<int> haplotype, vector<int>& result)
{
	for (int i = 0; i < genotype.size(); i++)
	{

		if ((genotype[i] == 0 && haplotype[i] != 0)
		    	|| (genotype[i] == 1 && haplotype[i] != 1))
			return false;
		else
		{

			if (genotype[i] == 2)
			{
				if (haplotype[i] == 1)
					result.push_back(0);
				else
					result.push_back(1);
			}	
			else
			{
				result.push_back(genotype[i]);
			}
		}
		
	}
	return true;
}

void printVector(vector<int> a)
{
	for (int i = 0; i < a.size(); i++)
		cout << a[i];
	cout << endl;
}

/*
	Uses Clarks Method, this is better but not
	great, O(m!) worst case
*/

vector<vector<int> > clarks(vector<vector<int> > genotypes)
{
	set<vector<int> > unexplained;
	set<vector<int> > explained;

	vector<vector<int> > known_haplotypes;

	// filter out all the unambigious.
	for (int i = 0; i < genotypes.size(); i++)
	{
		if (ambigious(genotypes[i]))
		{
			unexplained.insert(genotypes[i]);
		}
		else
		{
			explained.insert(genotypes[i]);
			known_haplotypes.push_back(genotypes[i]);
		}
	}



	int size = unexplained.size();
	while (size > 0)
	{
		// cycle through known haplotypes

		for (int i = 0; i < known_haplotypes.size(); i++)
		{


			// apply inference rule on unexplained
			set<vector<int> >::iterator it;
			for (it = unexplained.begin(); it != unexplained.end();)
			{
				vector<int> g = *it;
				vector<int> h;

				if (infer(g, known_haplotypes[i], h))
				{

					known_haplotypes.push_back(h);
					explained.insert(g);
					unexplained.erase(it++); 
				}
				else
				{
					it++;
				}

			}
		}


		if (unexplained.size() == size)
			break;
		else
			size = unexplained.size();

	}

	return known_haplotypes;
}

vector<vector<int> > generateGenotypes(int numGenotypes, int numSNPS)
{
	vector<vector<int> > ret;

	for (int i = 0; i < numGenotypes; i++)
	{
		vector<int> a;
		for (int j = 0; j < numSNPS; j++)
		{
			int pos = rand() % 3;
			a.push_back(pos);
		}
		ret.push_back(a);
	}
	return ret;
}

int main(int argc, char *argv[])
{

	srand(time(NULL));

	int numSNPs = atoi(argv[1]);
	
	vector<vector<int> > genotypes = generateGenotypes(10, numSNPs);
	
	vector<vector<int> > opt = baseline(genotypes);

	cout << numSNPs << " SNPs, optimal solution: " << opt.size() << " haplotypes." << endl;
	
	

}


